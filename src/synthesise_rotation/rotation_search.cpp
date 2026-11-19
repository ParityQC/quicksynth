// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/synthesise_rotation/rotation_search.hpp"

#include "quicksynth/synthesise_rotation/rotation_body.hpp"
#include "quicksynth/solve_norm_equation/relative_norm.hpp"
#include "quicksynth/search_scales/scale_bound.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/number_field.hpp"
#include "quicksynth/support/clock.hpp"
#include "quicksynth/support/host.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>

namespace quicksynth {

namespace {

// New body points expected per class at level l: the segment's area
// (2/3) eps^3 R^2 times the conjugate disks over the covolume, halved for
// "new" (test_rotation_body.cpp checks both against counts).
double expected_points(const GateSet& gateset, double epsilon, int level) {
    const double per_level = gateset.step == 1 ? std::pow(M_PI, 3) / 768 : M_PI / 12;
    return per_level * epsilon * epsilon * epsilon * std::ldexp(1.0, level);
}

int64_t checked_int64(__int128 x) {
    if (x > INT64_MAX || x < INT64_MIN)
        throw std::overflow_error("RotationSearch: a coefficient passes 64 bits");
    return static_cast<int64_t>(x);
}

// The least cost of an element whose u is new at level l, in class j.
//
// Clifford+T: u = delta^sub u' with delta not dividing u' has
// v_delta(u u*) = 2 sub, and v_delta(2) = 4, so the Bloch entry
// (u u* - t t*) / 2^K = 2 u u* / 2^K - 1 has sqrt2-denominator exponent l - 2.
// Clifford Bloch matrices are integral and T's has denominator sqrt2, so a word
// for the element holds at least l - 2 T gates [cf. Kliuchnikov, Maslov, Mosca,
// arXiv:1206.5236]. The element is M(u) T^j up to a power of omega, so
// its determinant is omega^j times an even power of omega, and a word's is
// omega^(T-count) times one: the T-count has the parity of j.
//
// Clifford+sqrt(T): the scale bound of K(l).
int64_t level_bound(const GateSet& gateset, const Price& price, int l, int j) {
    const int L = gateset.degree() / 2;
    const int64_t scale = integer_scale_bound(gateset, price, (l + L - 1) / L);
    if (gateset.step == 1) return scale;
    int t_count = std::max(l - 2, 0);
    if ((t_count - j) % 2 != 0) ++t_count;
    return std::max(scale, price.t * t_count);
}

void add(RotationWalk& total, const RotationWalk& w) {
    total.candidates += w.candidates;
    total.solvable += w.solvable;
    total.undecided += w.undecided;
    total.capped += w.capped;
    total.pricings += w.pricings;
}

template <size_t N>
struct Search {
    static constexpr int L = static_cast<int>(N / 2);
    using Ring = RingElement<N, __int128>;

    const GateSet& gateset;
    const Price& price;
    const RotationSearch::Settings& settings;
    double epsilon;
    Vec4 target;  // the hit filter's direction, as Enumerator::enumerate
    std::array<Vec4, GateSet::kDeterminantClasses> directions;
    Enumerator completer;
    Pricer pricer;
    // One LLL per (class, sub), reused at every scale.
    std::array<std::array<std::unique_ptr<RotationBody<N>>, L>, GateSet::kDeterminantClasses>
        bodies;

    Search(const GateSet& g, const Price& p, const RotationSearch::Settings& s, const Unitary& D,
           double eps)
        : gateset(g),
          price(p),
          settings(s),
          epsilon(eps),
          target(search_direction(D, 0, g)),
          completer(eps, g, p, 1, {}),
          pricer(Pricer::for_gateset(g, p)) {
        for (int j = 0; j < GateSet::kDeterminantClasses; ++j)
            directions[j] = search_direction(D, j, g);
    }

    // Prime decompositions recur from one norm equation to the next and
    // depend on p alone, so one cache per thread serves every search.
    static PrimeIdeals<N>& primes() {
        thread_local PrimeIdeals<N> cache;
        return cache;
    }

    const RotationBody<N>& body(int j, int sub) {
        auto& slot = bodies[j][sub];
        if (!slot) slot = std::make_unique<RotationBody<N>>(directions[j], sub, epsilon);
        return *slot;
    }

    // Least weight of P^i, P = diag(1, zeta) on Clifford+sqrt(T): a phase-i
    // completion costs at least phase 0's cost less twice this.
    int64_t phase_weight(int i) const {
        return i % 2 != 0 ? price.q : std::min(price.t, 2 * price.q);
    }

    // Whether phase i of a choice whose phase 0 costs `first` can beat the best.
    bool worth_pricing(int i, int64_t first, int64_t best_cost) const {
        return first - 2 * phase_weight(i) < best_cost;
    }

    // Prices the completion [u | t] at (K, class j) into `best` if it is within
    // epsilon; its cost.
    std::optional<int64_t> price_completion(const Ring& u, const Ring& t, int K, int j,
                                            RotationWalk& w, std::optional<Enumerator::Hit>& best) {
        int64_t row[2 * N];
        for (size_t k = 0; k < N; ++k) {
            row[k] = checked_int64(u[k]);
            row[N + k] = checked_int64(t[k]);
        }
        for (const CyclotomicMatrix& e : completer.matrices(row, 1, K, j)) {
            if (chord_distance(matrix_direction(e), target) > epsilon) continue;
            const Cost cost = pricer.price(e);
            ++w.pricings;
            if (!best || cost.value < best->cost.value) best = Enumerator::Hit{cost, e};
            return cost.value;
        }
        return std::nullopt;
    }

    // Prices the completions of u at (K, class j) within epsilon into `best`,
    // choice by choice, phases in order; on Clifford+sqrt(T) with lazy_phases
    // only the phases worth_pricing() keeps.
    void complete(const Ring& u, int K, int j, RotationWalk& w,
                  std::optional<Enumerator::Hit>& best) {
        RealElement<L, __int128> power{};
        power[0] = static_cast<__int128>(1) << K;
        const NormSolutions<N> solutions = solve_relative_norm<N>(
            real_sub(power, relative_norm(u)), primes(), kChoiceCap, settings.transpose_halving);
        if (solutions.outcome == NormOutcome::undecided) {
            ++w.undecided;
            return;
        }
        if (solutions.outcome == NormOutcome::none) return;
        ++w.solvable;
        w.capped += solutions.capped;
        const bool lazy = N == 8 && settings.lazy_phases;
        for (const Ring& t : solutions.choices) {
            const std::optional<int64_t> first = price_completion(u, t, K, j, w, best);
            for (int i = 1; i < L; ++i)
                if (!lazy || !first || worth_pricing(i, *first, best->cost.value))
                    price_completion(u, ring_times_zeta(t, i), K, j, w, best);
        }
    }

    RotationSearch::Result run() {
        RotationSearch::Result out;
        SearchSummary& summary = out.summary;
        summary.source = Source::rotation;
        const int useful = RotationSearch::first_useful_level(gateset, epsilon);
        const int last = useful + settings.levels_past_useful;
        std::optional<int> anchor;     // the hit cap's: first solvable level, at least `useful`
        std::optional<int> first_cut;  // first level with a walk past the point limit
        for (int l = 0;; ++l) {
            const int K = (l + L - 1) / L, sub = L * K - l;
            if (K > kRotationMaxScale) {
                out.stop = RotationStop::scale_cap;
                return finish(std::move(out), l - 1, false, false, first_cut);
            }
            if (summary.scales.empty() || summary.scales.back() != K) summary.scales.push_back(K);
            for (int j = 0; j < GateSet::kDeterminantClasses; ++j) {
                // the walk's candidates cost at least `floor`, and only a lower
                // cost replaces the best
                const int64_t floor = level_bound(gateset, price, l, j);
                const auto settled = [&] { return out.best && out.best->cost.value <= floor; };
                if (settled()) continue;
                poll_interrupt();
                if (settings.deadline && wall_seconds() > *settings.deadline) {
                    out.stop = RotationStop::time_limit;
                    return finish(std::move(out), l, false, true, first_cut);
                }
                const double start = steady_seconds();
                const RotationBody<N>& b = body(j, sub);
                RotationWalk w;
                // Clifford+sqrt(T): a walk with kRotationWalkSolvable times the
                // count rule's u is left, and the rule then stops the search
                // after this level; a slab (rotation_body.hpp) holds millions.
                const auto counted = [&] {
                    return N == 8 &&
                           w.solvable >= int64_t{kRotationWalkSolvable} * settings.solvable;
                };
                const int64_t points = b.walk(K, settings.point_limit, false,
                                              [&](const typename RotationBody<N>::Point& p) {
                                                  if (b.contains(p.u, K)) {
                                                      ++w.candidates;
                                                      complete(p.u, K, j, w, out.best);
                                                  }
                                                  return !settled() && !counted();
                                              });
                const bool cut = points > settings.point_limit || counted();
                SearchLog entry;
                entry.scale = K;
                entry.det_class = j;
                entry.sub = sub;
                entry.estimated = expected_points(gateset, epsilon, l);
                entry.points = points;
                if (out.best) entry.best_cost = out.best->cost.value;
                entry.rotation = w;
                entry.seconds = steady_seconds() - start;
                summary.estimated_points += entry.estimated;
                add(out.counts, w);
                out.log.push_back(std::move(entry));
                // The rest of a cut walk is lost, so the search goes on
                // uncertified: its later levels still hold candidates.
                if (cut && !first_cut) first_cut = l;
            }
            if (out.best && !summary.first_hit_scale) summary.first_hit_scale = K;
            int64_t bound = level_bound(gateset, price, l + 1, 0);
            for (int j = 1; j < GateSet::kDeterminantClasses; ++j)
                bound = std::min(bound, level_bound(gateset, price, l + 1, j));
            if (out.best && bound >= out.best->cost.value) {
                out.stop = RotationStop::bound;
                return finish(std::move(out), l, true, false, first_cut);
            }
            if (N == 8) {
                if (!anchor && out.counts.solvable > 0) anchor = std::max(l, useful);
                if (out.best && out.counts.solvable >= settings.solvable) {
                    out.stop = RotationStop::count;
                    return finish(std::move(out), l, false, false, first_cut);
                }
                if (anchor && l >= *anchor + settings.levels_past_hit) {
                    out.stop = RotationStop::hit_cap;
                    return finish(std::move(out), l, false, false, first_cut);
                }
            }
            if (l >= last) {
                out.stop = RotationStop::level_cap;
                return finish(std::move(out), l, false, false, first_cut);
            }
        }
    }

    // Summary fields once level l is the last one walked, `partial` if not to
    // its end; `cut` is the first level with a cut walk.
    RotationSearch::Result finish(RotationSearch::Result out, int l, bool bounded, bool partial,
                                  std::optional<int> cut) const {
        SearchSummary& summary = out.summary;
        const auto scale = [](int level) { return (level + L - 1) / L; };
        // a scale is complete once its sub = 0 level is walked, every walk whole
        int full = l == L * scale(l) && !partial ? scale(l) : scale(l) - 1;
        if (cut) full = std::min(full, scale(*cut) - 1);
        if (full >= 0) summary.last_full_scale = full;
        out.certified = bounded && !cut && out.counts.undecided == 0 && out.counts.capped == 0;
        if (out.best) {
            summary.scan_cost = out.best->cost.value;
            if (out.certified) summary.lower_bound = static_cast<double>(out.best->cost.value);
        }
        summary.scan_certified = out.certified;
        return out;
    }
};

}  // namespace

RotationSearch::RotationSearch(GateSet gateset, Price price, Settings settings)
    : gateset_(std::move(gateset)), price_(price), settings_(settings) {}

int RotationSearch::first_useful_level(const GateSet& gateset, double epsilon) {
    return static_cast<int>(std::ceil(-std::log2(expected_points(gateset, epsilon, 0))));
}

RotationSearch::Result RotationSearch::run(const Unitary& D, double epsilon) const {
    if (D[2] != 0 || D[3] != 0 || D[4] != 0 || D[5] != 0)
        throw std::invalid_argument("RotationSearch: the target is not diagonal");
    if (!(epsilon > 0 && epsilon < 1))
        throw std::invalid_argument("RotationSearch: epsilon outside (0, 1)");
    if (gateset_.step == 1) return Search<8>(gateset_, price_, settings_, D, epsilon).run();
    return Search<4>(gateset_, price_, settings_, D, epsilon).run();
}

}  // namespace quicksynth
