// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/search_scales/split.hpp"

#include "quicksynth/search_scales/scale_bound.hpp"
#include "quicksynth/support/clock.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <utility>

namespace quicksynth {

namespace {

// Returns each walk's share of `points` leaves, at least 1. Returns 0 (no cap)
// for an infinite budget or a share >= 1e17 (node cap is 64 leaf caps, int64).
int64_t leaf_cap(double points, double walks) {
    const double share = points / walks;
    if (!(share < 1e17)) return 0;
    return std::max<int64_t>(1, static_cast<int64_t>(share));
}

// Copies a set into memory allocated with the split, for locality.
std::shared_ptr<const PrefixSet> own_copy(const std::shared_ptr<const PrefixSet>& set) {
    return std::make_shared<const PrefixSet>(*set);
}

// Copies of the fixed sets, then the extra sets uncopied. run_model copies
// those lazily, since a certified scan never splits.
PrefixSets split_sets(const PrefixSets& fixed, const PrefixSets& extra) {
    PrefixSets out;
    for (const auto& set : fixed) out.push_back(own_copy(set));
    out.insert(out.end(), extra.begin(), extra.end());
    return out;
}

// Transform-cache size: 4 set-ups per prefix of each distinct scale.
int64_t transform_cache_for(const PrefixSets& sets) {
    std::map<int, int64_t> sizes;  // per distinct scale
    for (const auto& s : sets) sizes[s->scale] = static_cast<int64_t>(s->prefixes.size());
    int64_t total = 0;
    for (const auto& [scale, n] : sizes) total += n;
    return std::max<int64_t>(kSplitTransformCache, 4 * total + 1024);
}

// One remainder problem per prefix and determinant class, targeting P^-1 V.
std::vector<Enumerator::Problem> remainders(const PrefixSet& set, const Unitary& V,
                                            int det_classes) {
    std::vector<Enumerator::Problem> problems;
    problems.reserve(set.prefixes.size() * static_cast<size_t>(det_classes));
    for (size_t p = 0; p < set.prefixes.size(); ++p) {
        const Unitary rotated = matmul(set.inverses[p], V);
        for (int det_class = 0; det_class < det_classes; ++det_class)
            problems.push_back({rotated, det_class});
    }
    return problems;
}

SearchLog composite_log(int scale, int prefix_scale, double seconds, double estimated,
                        int64_t points, std::optional<int64_t> best_cost, int sub = 0) {
    SearchLog e;
    e.scale = scale;
    e.sub = sub;
    e.prefix_scale = prefix_scale;
    e.seconds = seconds;
    e.estimated = estimated;
    e.points = points;
    e.best_cost = best_cost;
    return e;
}

}  // namespace

PrefixSplit::PrefixSplit(double epsilon, GateSet gateset, Price price, int threads,
                         const PrefixSets& fixed, const PrefixSets& extra, Settings settings)
    : epsilon_(epsilon),
      price_(price),
      pricer_(Pricer::for_gateset(gateset, price)),
      settings_(settings),
      sets_(split_sets(fixed, extra)),
      fixed_(fixed.size()),
      copied_(fixed.size()),
      enumerator_(epsilon, std::move(gateset), price, threads,
                  Enumerator::Settings{transform_cache_for(sets_)}) {}

PrefixSplit::Result PrefixSplit::run(const Unitary& V, int scales, std::optional<int> first,
                                     std::optional<int64_t> best_cost) {
    Result out;
    for (size_t s = 0; s < fixed_; ++s) {
        if (s > 0 && out.best) break;
        Result more = run_one(V, *sets_[s], scales, first, best_cost);
        more.log.insert(more.log.begin(), out.log.begin(), out.log.end());
        out = std::move(more);
    }
    return out;
}

PrefixSplit::Result PrefixSplit::run_one(const Unitary& V, const PrefixSet& set, int scales,
                                         std::optional<int> first,
                                         std::optional<int64_t> best_cost) {
    Result out;
    const std::vector<CyclotomicMatrix>& prefixes = set.prefixes;
    // with no prefixes the budget test below never stops the loop
    if (prefixes.empty()) return out;
    const int k = set.scale;
    const GateSet& gateset = enumerator_.gateset();
    const int det_classes = GateSet::kDeterminantClasses;
    const int64_t walks = static_cast<int64_t>(det_classes) * static_cast<int64_t>(prefixes.size());
    // Per-walk leaf cap bounds work where a rotated region's volume estimate is
    // far off.
    const int64_t walk_budget =
        leaf_cap(settings_.point_budget, static_cast<double>(std::max<int64_t>(walks, 1)));
    // composite level l = L K' - sub' of scale K'; its remainder is level l - L k
    const int L = settings_.levels;
    int l =
        first ? *first : (L == 1 ? enumerator_.first_scale() : enumerator_.region().first_level(L));
    const int l_max = l + L * kSplitScaleGuard;
    const Vec4 target = su2_vector(V);
    const std::vector<Enumerator::Problem> problems = remainders(set, V, det_classes);
    int after_hit = 0;
    bool hit = false;
    for (;; ++l) {
        if (l_max < l) break;
        const LatticeLevel remainder = level_at(std::max(l - L * k, 0), L);
        const double est = enumerator_.estimate_points(remainder.K, remainder.sub);
        if (static_cast<double>(walks) * est > settings_.point_budget) break;
        if (l - L * k < 0) continue;
        // Levels ascend, so stop once the bound at the scale K reached the best cost.
        const int K = level_at(l, L).K;
        const std::optional<double> bound = split_scale_bound(gateset, price_, K);
        if (out.best && bound && *bound >= static_cast<double>(out.best->cost.value)) break;
        if (settings_.incumbent && bound && *bound >= static_cast<double>(*settings_.incumbent))
            break;
        if (out.best && (!best_cost || out.best->cost.value < *best_cost))
            best_cost = out.best->cost.value;
        size_t walked = prefixes.size();
        if (settings_.rate && best_cost)
            walked = affordable_prefixes(est, prefixes.size(), *best_cost);
        if (walked == 0) break;
        const size_t count = walked * static_cast<size_t>(det_classes);
        const double composite_start = steady_seconds();
        const Composite c = walk_composite(set, problems, target, remainder.K + k, walk_budget,
                                           out.best, count, remainder.sub);
        const double spent = steady_seconds() - composite_start;
        const std::optional<int64_t> best =
            out.best ? std::optional<int64_t>(out.best->cost.value) : std::nullopt;
        out.log.push_back(composite_log(remainder.K + k, k, spent, static_cast<double>(count) * est,
                                        c.visited, best, remainder.sub));
        if (walked < prefixes.size()) break;
        hit = hit || c.found;
        if (hit && ++after_hit >= scales) break;
    }
    return out;
}

PrefixSplit::Prefixed PrefixSplit::prefixed(const Unitary& V, size_t set) const {
    const GateSet& gateset = enumerator_.gateset();
    return {sets_[set].get(), remainders(*sets_[set], V, GateSet::kDeterminantClasses),
            search_direction(V, 0, gateset)};
}

std::optional<PrefixSplit::Level> PrefixSplit::level(const Prefixed& p, int l) const {
    const int L = settings_.levels;
    const int r = l - L * p.set->scale;
    if (r < 0) return std::nullopt;
    const LatticeLevel remainder = level_at(r, L);
    return Level{remainder, static_cast<double>(p.problems.size()),
                 enumerator_.estimate_points(remainder.K, remainder.sub)};
}

SearchLog PrefixSplit::walk_level(const Prefixed& p, int l, std::optional<Enumerator::Hit>& best) {
    const Level at = level(p, l).value();
    const int k = p.set->scale;
    const double start = steady_seconds();
    const Composite c = walk_composite(*p.set, p.problems, p.target, at.remainder.K + k,
                                       leaf_cap(settings_.point_budget, std::max(at.walks, 1.0)),
                                       best, p.problems.size(), at.remainder.sub);
    const double spent = steady_seconds() - start;
    const std::optional<int64_t> cost =
        best ? std::optional<int64_t>(best->cost.value) : std::nullopt;
    return composite_log(at.remainder.K + k, k, spent, at.walks * at.est, c.visited, cost,
                         at.remainder.sub);
}

size_t PrefixSplit::affordable_prefixes(double est, size_t prefixes, int64_t best_cost) const {
    const TimeModel& time = time_model(enumerator_.gateset());
    // time allowance kScaleSaving x best / rate over the time model's per-prefix time
    const double per_prefix =
        GateSet::kDeterminantClasses * (time.split_overhead + est / time.split_throughput);
    const double allowance = kScaleSaving * static_cast<double>(best_cost) / *settings_.rate;
    if (allowance > static_cast<double>(prefixes) * per_prefix) return prefixes;
    if (!settings_.partial) return 0;
    return std::min(prefixes, static_cast<size_t>(allowance / per_prefix));
}

PrefixSplit::Composite PrefixSplit::walk_composite(const PrefixSet& set,
                                                   const std::vector<Enumerator::Problem>& problems,
                                                   const Vec4& target, int K, int64_t leaf_cap,
                                                   std::optional<Enumerator::Hit>& best,
                                                   size_t count, int sub) {
    Composite out;
    const int r = K - set.scale;
    const std::vector<WalkResult> walks =
        count < problems.size()
            ? enumerator_.survivors_batch({problems.begin(), problems.begin() + count}, r, leaf_cap,
                                          sub)
            : enumerator_.survivors_batch(problems, r, leaf_cap, sub);
    const size_t n = static_cast<size_t>(enumerator_.gateset().dim());
    for (size_t w = 0; w < walks.size(); ++w) {
        out.visited += walks[w].visited;
        const std::vector<int64_t>& rows = walks[w].coefficients;
        const size_t p = w / static_cast<size_t>(GateSet::kDeterminantClasses);
        for (const CyclotomicMatrix& e :
             enumerator_.matrices(rows.data(), rows.size() / n, r, problems[w].det_class)) {
            const CyclotomicMatrix composite = compose(set.prefixes[p], e);
            const double distance = chord_distance(matrix_direction(composite), target);
            if (distance > epsilon_) continue;
            out.found = true;
            const Cost cost = pricer_.price(composite);
            if (!best || cost.value < best->cost.value) best = Enumerator::Hit{cost, composite};
        }
    }
    return out;
}

PrefixSplit::Result PrefixSplit::run_model(const Unitary& V, const StoppingModel& model,
                                           const ScaleScan::Result& scan, double rate,
                                           bool walk_first) {
    for (; copied_ < sets_.size(); ++copied_) sets_[copied_] = own_copy(sets_[copied_]);
    Result out;
    const GateSet& gateset = enumerator_.gateset();
    const int det_classes = GateSet::kDeterminantClasses;
    // Candidates: the first set of each prefix scale k >= 2, sorted by k.
    struct Candidate {
        const PrefixSet* set;
        std::array<double, 3> phi;
        std::vector<Enumerator::Problem> problems;
    };
    std::vector<Candidate> candidates;
    std::set<int> scales_seen;
    for (size_t i = 0; i < sets_.size(); ++i) {
        const PrefixSet& set = *sets_[i];
        if (set.scale < 2 || set.prefixes.empty() || !scales_seen.insert(set.scale).second)
            continue;
        candidates.push_back({&set,
                              split_coverage(gateset, model.price(), set.scale,
                                             static_cast<int64_t>(set.prefixes.size())),
                              remainders(set, V, det_classes)});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) { return a.set->scale < b.set->scale; });
    const int L = scan.last_full_scale ? *scan.last_full_scale : model.start() - 1;
    std::optional<int64_t> best;
    if (scan.best) best = scan.best->cost.value;
    const Vec4 target = su2_vector(V);
    const TimeModel& time = time_model(gateset);
    const double a = time.split_overhead;
    const double theta = time.split_throughput;
    const auto seconds = [&](int Kc, const Candidate& c) {
        return static_cast<double>(c.problems.size()) *
               (a + enumerator_.estimate_points(Kc - c.set->scale) / theta);
    };
    // Returns the best-value candidate at Kc and its value (-1: none within
    // budget; `shallow`: none deep enough yet). With no answer yet, picks the
    // highest hit chance per second instead.
    const auto pick = [&](int Kc, bool& shallow) {
        std::pair<int, double> choice{-1, -std::numeric_limits<double>::infinity()};
        shallow = true;
        for (size_t i = 0; i < candidates.size(); ++i) {
            const Candidate& c = candidates[i];
            if (Kc - c.set->scale < 0) continue;
            shallow = false;
            const double work = static_cast<double>(c.problems.size()) *
                                enumerator_.estimate_points(Kc - c.set->scale);
            if (work > settings_.point_budget) continue;
            const double v =
                best ? model.split_saving(Kc, c.phi, *best, L) - rate * seconds(Kc, c)
                     : model.split_p_cheaper(Kc, c.phi, model.cost_cap() + 1, L) / seconds(Kc, c);
            if (choice.first < 0 || v > choice.second) choice = {static_cast<int>(i), v};
        }
        return choice;
    };
    bool first = true;
    for (int Kc = L + 1; Kc <= L + 1 + kSplitScaleGuard; ++Kc) {
        // proved: nothing past the scan's scales is cheaper
        if (best && model.c_int(L + 1) >= *best) break;
        bool shallow = false;
        const auto [chosen, value] = pick(Kc, shallow);
        if (chosen < 0) {
            if (shallow) continue;
            break;
        }
        if (best && !(first && walk_first) && !(value > 0.0)) break;
        first = false;
        const Candidate& c = candidates[static_cast<size_t>(chosen)];
        const int k = c.set->scale;
        const double composite_start = steady_seconds();
        const double walks = static_cast<double>(c.problems.size());
        const Composite composite =
            walk_composite(*c.set, c.problems, target, Kc, leaf_cap(settings_.point_budget, walks),
                           out.best, c.problems.size());
        if (out.best && (!best || out.best->cost.value < *best)) best = out.best->cost.value;
        const double spent = steady_seconds() - composite_start;
        out.log.push_back(composite_log(Kc, k, spent, walks * enumerator_.estimate_points(Kc - k),
                                        composite.visited, best));
    }
    return out;
}

}  // namespace quicksynth
