// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The cheapest candidate near a diagonal target D = diag(z, z*), level by
// level [Ross, Selinger, arXiv:1403.2975; Kliuchnikov,
// Bocharov, Roetteler, Yard, arXiv:1510.03888].
//
// A candidate [[u, -t* w], [t, u* w]] / sqrt2^K, u, t in O_K = Z[zeta], w the
// determinant phase of class j, has its distance to D from u alone. Levels
// l = L K - sub, L = N / 2, run from 0 up, both classes per level. At each, the
// u of the rotation body (rotation_body.hpp) that are new at the level and pass
// its exact test give xi = 2^K - u u*; the relative norm equation t t* = xi
// (relative_norm.hpp) gives every t up to roots of unity, and each t zeta^i,
// i < L, completes a candidate through Enumerator::matrices (exact unit-norm
// check, class mapping, reduction). Completions within epsilon of D (the joint
// walk's hit filter) are priced. The phases zeta^i, i >= L, are the S- and
// Z-conjugates t -> i t, -t, and the transpose halving of the solutions keeps
// one of each transposed pair; neither changes the distance or the cost.
//
// A u new at (K, sub) has t of the same delta-valuation, since u u* + t t* =
// 2^K, so every candidate first appears at exactly one level and the element's
// sqrt2-exponent is K. Each level and class has a cost bound, nondecreasing in
// the level: on Clifford+T, T-count >= l - 2 with the parity of the class; on
// Clifford+sqrt(T), the bound of scale K (integer_scale_bound, scale_bound.hpp).
// A walk is skipped, or left early, once the best cost is at most its bound,
// and the search stops after level l once level l + 1's bound reaches the best
// cost.
//
// On Clifford+sqrt(T) that bound is about half the costs found, so a count
// rule stops the search instead: whole levels are walked until Settings::solvable
// u were solvable. Within a level the costs spread by about 10, and the least
// lies up to about 10 levels past the first useful one; the rule is
// deterministic.
//
// Phase sweep on Clifford+sqrt(T). The phase t -> t zeta^i is conjugation by
// P^i, P = diag(1, zeta) = sqrt(T), so the costs of the phases of one choice
// differ by at most twice the cost of P^i: q_cost for odd i, t_cost for i = 2
// (the cost is the least weight of a word [paper: cor:cost], which is
// subadditive). Phases 1..3 of a choice are priced only if phase 0's cost less
// that difference is below the best cost; a skipped phase cannot be cheaper
// than the best, so the sweep finds the element the full sweep finds. Since
// every phase that could win is priced, the transpose halving stays exact.
#pragma once

#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/data_structures/cost.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"
#include "quicksynth/support/diagnostics.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace quicksynth {

#pragma GCC visibility push(hidden)

// Ellipsoid points one walk may list. Directions within about epsilon of a field
// direction hold a slab of points that grows with the scale (rotation_body.hpp);
// a walk past this many is cut, which clears the certificate, and the search
// goes on with the next walk.
inline constexpr int64_t kRotationPointLimit = int64_t{1} << 20;
// Levels searched past the first useful one (one new body point expected per
// class) before the search gives up; near-exact angles answer up to about 11
// levels past it.
inline constexpr int kRotationLevelsPastUseful = 32;
// Clifford+sqrt(T) count rule: solvable u before the search stops, by default
// and under `fast`, and levels past the first solvable one, or the first useful
// one if later, before it gives up. A u can be solvable up to about 12 levels
// below the first useful level, and from there on the solvable u per level
// about double. The counts were chosen on 100 random angles at epsilon 1e-3 to
// 1e-10, price (1, 3), against the reference costs of golden_rotation.json's
// sqrt_t table: from 16 to 4096 solvable u the mean cost falls by 0.4 to 1.7
// and the time grows in proportion. At 256 the mean is 1.1 to 1.6 below the
// reference and every angle but one at most 2 above; at 40 the mean is 0.6 to
// 1.1 below, in a fifth of the time or less.
inline constexpr int kRotationSolvable = 256;
inline constexpr int kRotationSolvableFast = 40;
inline constexpr int kRotationLevelsPastHit = 16;
// Solvable u of one walk, in multiples of the count rule's, before the walk is
// cut. A level's walks hold about the count rule's u when it stops, except in
// a slab (rotation_body.hpp).
inline constexpr int kRotationWalkSolvable = 16;
// Largest scale walked: u and t have coefficients at most 2^(K/2), below the
// 2^60 that Enumerator::matrices takes.
inline constexpr int kRotationMaxScale = 118;

// Why the search stopped. Only bound and count are answers on
// Clifford+sqrt(T).
enum class RotationStop {
    bound,       // the next level's cost bound reached the best cost
    count,       // Clifford+sqrt(T): Settings::solvable u were solvable
    level_cap,   // kRotationLevelsPastUseful (Settings::levels_past_useful) passed
    hit_cap,     // Clifford+sqrt(T): Settings::levels_past_hit passed before the count
    scale_cap,   // the next scale passes kRotationMaxScale
    time_limit,  // Settings::deadline passed before a walk
};

class RotationSearch {
  public:
    struct Settings {
        int64_t point_limit = kRotationPointLimit;
        int levels_past_useful = kRotationLevelsPastUseful;
        // Clifford+sqrt(T): the count rule
        int solvable = kRotationSolvable;
        int levels_past_hit = kRotationLevelsPastHit;
        // Clifford+sqrt(T); false: every phase of every choice is priced
        bool lazy_phases = true;
        bool transpose_halving = true;   // solve_relative_norm
        std::optional<double> deadline;  // wall_seconds(); checked before each walk
    };
    struct Result {
        std::optional<Enumerator::Hit> best;  // first of the least cost, in walk order
        RotationStop stop = RotationStop::bound;
        // Stopped by the bound with every candidate decided, no choice list
        // cut and no walk cut: no element of the gate set within epsilon of D
        // costs less.
        bool certified = false;
        std::vector<SearchLog> log;  // one per walked (K, sub, class), in walk order
        RotationWalk counts;         // totals over the log
        SearchSummary summary;       // source rotation; scan_certified = certified
    };

    // `price` is the searched price: Price::clifford_t(t) on Clifford+T.
    RotationSearch(GateSet gateset, Price price, Settings settings);

    // The search around D, a diagonal unitary (std::invalid_argument
    // otherwise), at epsilon in (0, 1). Polls the host's interrupt check per walk.
    Result run(const Unitary& D, double epsilon) const;

    // The first level with one new body point expected per class:
    // pi eps^3 2^l / 12 on Clifford+T, pi^3 eps^3 2^l / 768 on Clifford+sqrt(T).
    static int first_useful_level(const GateSet& gateset, double epsilon);

  private:
    GateSet gateset_;
    Price price_;
    Settings settings_;
};

#pragma GCC visibility pop

}  // namespace quicksynth
