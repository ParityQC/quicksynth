// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/search_scales/level_search.hpp"

#include "quicksynth/search_scales/scale_bound.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/enumerate_lattice/search_region.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace quicksynth {

namespace {

// Predicted LLL set-up seconds per walk for one search's directions. The
// enumerator reduces sub 0 afresh, derives each finer sub from the next coarser
// and keeps every sub on the way, so subs 0..finest are built.
class Setups {
  public:
    double seconds(int sub) const {
        if (sub <= finest_) return 0.0;
        const double fresh = finest_ < 0 ? kFreshSetupSeconds : 0.0;
        return fresh + kDerivedSetupSeconds * (sub - std::max(finest_, 0));
    }
    void built(int sub) { finest_ = std::max(finest_, sub); }
    int finest() const { return finest_; }

  private:
    int finest_ = -1;
};

std::vector<size_t> every_prefix(const PrefixSplit::Prefixed& p) {
    std::vector<size_t> out(p.set->prefixes.size());
    std::iota(out.begin(), out.end(), size_t{0});
    return out;
}

}  // namespace

LevelSearch::LevelSearch(ScaleScan& scan, PrefixSplit* split, Price price, Settings settings)
    : scan_(scan), split_(split), price_(price), settings_(settings) {}

LevelSearch::Result LevelSearch::run(const Unitary& V, const ScaleScan::Progress& progress) {
    const int det_classes = GateSet::kDeterminantClasses;
    Enumerator& enumerator = scan_.enumerator();
    const GateSet& gateset = enumerator.gateset();
    const int L = scan_.levels();
    const double scan_budget = scan_.settings().point_budget;
    const double split_budget = split_ ? split_->settings().point_budget : 0.0;
    Result out;
    ScaleScan::Result scan;
    std::optional<Enumerator::Hit> split_best;
    std::optional<PrefixSplit::Prefixed> prefixed;  // the split's prefix set
    Setups scan_setups, split_setups;
    bool split_ran = false;
    const auto walk_scan = [&](int l, double est) {
        const size_t logged = scan.log.size();
        scan_.walk_level(V, l, est, scan, progress);
        scan_setups.built(level_at(l, L).sub);
        for (size_t i = logged; i < scan.log.size(); ++i) {
            out.log.push_back(scan.log[i]);
            count_other(out.log.back(), split_best);
        }
    };
    // Walks composite level l on `prefixes`; returns the remainder's sub.
    const auto walk_split = [&](int l, const std::vector<size_t>& prefixes) {
        out.log.push_back(split_->walk_level(*prefixed, l, prefixes, split_best));
        count_other(out.log.back(), scan.best);
        split_ran = true;
        return out.log.back().sub;
    };
    // The cheaper answer (the scan's on a tie) and the summary of what was walked.
    const auto finish = [&] {
        scan_.conclude(scan);
        out.summary = scan_summary(scan);
        out.summary.split_ran = split_ran;
        out.best = scan.best;
        if (split_best && (!out.best || split_best->cost.value < out.best->cost.value)) {
            out.best = split_best;
            out.summary.source = Source::split;
        }
        return std::move(out);
    };

    // Whether the split leaves out level y before the first answer: its set-ups
    // take longer than its walks and the next level, which holds its elements,
    // stays a factor kSkipHeadroom under the split budget.
    const auto setup_bound = [&](int y, const PrefixSplit::Level& at) {
        if (!(split_setups.seconds(at.remainder.sub) > walk_seconds(at.est))) return false;
        const std::optional<PrefixSplit::Level> next = split_->level(*prefixed, y + 1);
        return next && !(next->walks * next->est * kSkipHeadroom > split_budget);
    };

    // Before the first answer: the scan's levels in order.
    const int first = enumerator.region().first_level(L);
    for (int l = first;; ++l) {
        const auto [K, sub] = level_at(l, L);
        const double est = enumerator.estimate_points(K, sub);
        if (est > scan_budget) break;
        walk_scan(l, est);
        if (scan.best) break;
    }
    int cov = scan.last_level.value_or(first - 1);  // the scan's last level
    int scov = cov;                                 // the split's; none while <= cov
    std::optional<int> answer_level = scan.first_hit_level;
    // Without a scan answer: the split's levels past the scan's, each fixed set in turn.
    for (size_t set = 0; !scan.best && split_ && !split_best && set < split_->fixed_sets(); ++set) {
        prefixed = split_->prefixed(V, set);
        split_setups = Setups{};
        scov = cov;
        for (int y = cov + 1; !prefixed->problems.empty() && y <= cov + L * kSplitScaleGuard; ++y) {
            const std::optional<PrefixSplit::Level> at = split_->level(*prefixed, y);
            if (!at) continue;
            if (at->walks * at->est > split_budget) break;
            if (settings_.skip_setup_levels && setup_bound(y, *at)) continue;
            split_setups.built(walk_split(y, every_prefix(*prefixed)));
            scov = y;
            if (split_best) {
                answer_level = y;
                break;
            }
        }
    }
    if (!answer_level) return finish();

    // After it: the scan's most valuable level or the split's most valuable batch
    // of prefix walks, while one is worth its time.
    const bool split_after = split_ && (!scan.best || enumerator.epsilon() < settings_.split_below);
    for (size_t set = 0; split_after && !prefixed && set < split_->fixed_sets(); ++set) {
        PrefixSplit::Prefixed p = split_->prefixed(V, set);
        if (!p.problems.empty()) prefixed = std::move(p);
    }
    const PrefixSet* prefix_set = split_after && prefixed ? prefixed->set : nullptr;
    const std::unique_ptr<LevelGain> gain =
        level_gain(gateset, price_, enumerator.epsilon(), V, enumerator.first_scale(),
                   *answer_level, prefix_set);
    const LevelRates rates = gain->rates(settings_.rate);
    const size_t prefixes = prefix_set ? prefix_set->prefixes.size() : 0;
    std::vector<int> depth(prefixes, scov);  // each prefix's last composite level
    std::vector<Setups> prefix_setups(prefixes, split_setups);
    const auto best_cost = [&] {
        if (!split_best) return scan.best->cost.value;
        if (!scan.best) return split_best->cost.value;
        return std::min(scan.best->cost.value, split_best->cost.value);
    };
    struct Walk {
        size_t prefix;
        int level;
        double seconds, value;
    };
    // Prefix p's composite level of largest positive value past both searches.
    const auto most_valuable_walk = [&](size_t p, int64_t best) {
        std::optional<Walk> pick;
        const int lo = std::max(cov, depth[p]);
        for (int y = lo + 1; y <= lo + L * kSplitScaleGuard; ++y) {
            const std::optional<PrefixSplit::Level> at = split_->level(*prefixed, y);
            if (!at) continue;
            if (at->walks * at->est > split_budget) break;
            const std::optional<double> bound =
                split_scale_bound(gateset, price_, level_at(y, L).K);
            if (bound && *bound >= static_cast<double>(best)) break;
            const double seconds = det_classes * (walk_seconds(at->est) +
                                                  prefix_setups[p].seconds(at->remainder.sub));
            const double v = gain->split(cov, {{p, depth[p], y}}, best) - rates.split * seconds;
            if (v > (pick ? pick->value : 0.0)) pick = Walk{p, y, seconds, v};
        }
        return pick;
    };
    for (;;) {
        const int64_t best = best_cost();
        std::optional<int> scan_level;
        double scan_est = 0.0, scan_value = 0.0;
        for (int a = cov + 1;; ++a) {
            const auto [K, sub] = level_at(a, L);
            if (scan_.scale_lower_bound(K) >= static_cast<double>(best)) break;
            const double est = enumerator.estimate_points(K, sub);
            if (est > scan_budget) break;
            const double seconds = det_classes * (walk_seconds(est) + scan_setups.seconds(sub));
            const double v = gain->scan(cov, a, depth, best) - rates.scan * seconds;
            if (v > scan_value) {
                scan_value = v;
                scan_level = a;
                scan_est = est;
            }
        }
        // Each prefix's most valuable walk, priced once per class and state; the
        // batch of the most valuable, or the first alone where together they are not
        // worth their time.
        std::map<std::tuple<int, int, int>, std::optional<Walk>> priced;
        std::vector<Walk> batch;
        for (size_t p = 0; p < prefixes; ++p) {
            const auto key =
                std::make_tuple(gain->prefix_class(p), depth[p], prefix_setups[p].finest());
            auto it = priced.find(key);
            if (it == priced.end()) it = priced.emplace(key, most_valuable_walk(p, best)).first;
            if (it->second) batch.push_back({p, it->second->level, it->second->seconds,
                                             it->second->value});
        }
        std::stable_sort(batch.begin(), batch.end(),
                         [](const Walk& a, const Walk& b) { return a.value > b.value; });
        batch.resize(std::min(batch.size(), gain->batch(prefixes)));
        double split_value = 0.0;
        if (!batch.empty()) {
            std::vector<PrefixWalk> walks;
            double seconds = 0.0;
            for (const Walk& w : batch) {
                walks.push_back({w.prefix, depth[w.prefix], w.level});
                seconds += w.seconds;
            }
            split_value = gain->split(cov, walks, best) - rates.split * seconds;
            if (!(split_value > 0.0)) {
                batch.resize(1);
                split_value = batch.front().value;
            }
        }
        // ties go to the scan
        if (scan_level && scan_value >= split_value) {
            walk_scan(*scan_level, scan_est);
            cov = *scan_level;
            continue;
        }
        if (batch.empty()) break;
        std::map<int, std::vector<size_t>> by_level;
        for (const Walk& w : batch) by_level[w.level].push_back(w.prefix);
        for (auto& [level, walked] : by_level) {
            std::sort(walked.begin(), walked.end());
            const int sub = walk_split(level, walked);
            for (const size_t p : walked) {
                depth[p] = level;
                prefix_setups[p].built(sub);
            }
        }
    }
    return finish();
}

}  // namespace quicksynth
