// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Prefix split: searches a deeper denominator scale for the same work.
//
// An element U of scale K factors as P R with P of scale k from a cost ball and
// R of scale K - k. Since d(R, P^-1 V) = d(U, V), R is enumerated around the
// rotated target in a region 2^(deg k) times smaller. Heuristic: U need not
// factor with a prefix of the chosen scale. Each composite is checked against
// epsilon.
#pragma once

#include "quicksynth/search_scales/scan.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/cost_ball.hpp"
#include "quicksynth/data_structures/prefix_sets.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"
#include "quicksynth/support/diagnostics.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace quicksynth {

constexpr double kSplitPointBudget = 3e7;
// Minimum transform-cache size; each prefix is a separate target and set-up.
constexpr int64_t kSplitTransformCache = 8192;
// Termination guard, not a setting: no budget reaches this many scales past
// the first.
constexpr int kSplitScaleGuard = 64;

class PrefixSplit {
  public:
    struct Settings {
        // per composite scale, walks x est(K - k); split into per-walk leaf caps
        double point_budget = kSplitPointBudget;
        // run(): rate rule, magic states per second; with an answer in hand, a composite
        // scale is walked iff kScaleSaving x best cost > rate x its time_model time.
        // None: off.
        std::optional<double> rate;
        // rate rule: a composite scale over its allowance is walked on the cheapest
        // prefixes the allowance covers, then the split stops; false: not walked
        bool partial = false;
        // run(): levels per scale (SearchRegion::first_level); 1: whole scales
        int levels = 1;
        // run(): cost to beat; stop where the scale bound reaches it
        std::optional<int64_t> incumbent;
    };
    struct Result {
        // first of the cheapest composites P R
        std::optional<Enumerator::Hit> best;
        std::vector<SearchLog> log;  // one entry per composite scale
    };
    // Holds copies of CostBall::prefix_sets. run() tries `fixed` (prefix scale,
    // then fallbacks); run_model() also uses `extra`, copied on first run.
    PrefixSplit(double epsilon, GateSet gateset, Price price, int threads, const PrefixSets& fixed,
                const PrefixSets& extra, Settings settings);

    // Tries the fixed sets in order until one finds an element; logs accumulate.
    // Walks `scales` composite levels from the first hit, fewer under the rate
    // rule. `first` is the first composite level walked (none: first feasible);
    // `best_cost` the scan's answer, for the rate rule. Levels are scales when
    // Settings::levels is 1.
    Result run(const Unitary& V, int scales, std::optional<int> first,
               std::optional<int64_t> best_cost = std::nullopt);

    // Model-rule split (stopping_rule.hpp), from composite scale L + 1 (L: the scan's
    // last full scale), starting from the scan's answer.
    // - Per scale K', picks the prefix scale k >= 2 maximising
    //   S_split - rate sec_split among those with det_classes x M_k x
    //   est(K' - k) within budget; sec_split = det_classes x M_k
    //   (a + est(K' - k) / theta_split).
    // - Walks while the value is positive and c_int(L + 1) < best cost. A
    //   one-scale look-ahead (next scale 16-256x dearer) gave the same scales
    //   on the calibration set.
    // - With no answer yet, walks every affordable scale, with the prefix scale
    //   of highest hit chance per second.
    // - `walk_first` walks the first composite scale regardless of value.
    Result run_model(const Unitary& V, const StoppingModel& model, const ScaleScan::Result& scan,
                     double rate, bool walk_first);

    // Per-level access for LevelSearch, on composite levels (Settings::levels).
    // One fixed set's remainder problems around a target.
    struct Prefixed {
        const PrefixSet* set = nullptr;
        std::vector<Enumerator::Problem> problems;
        Vec4 target;
    };
    // A composite level: its remainder level, walks and estimated points per walk.
    struct Level {
        LatticeLevel remainder;
        double walks = 0.0;
        double est = 0.0;
    };
    const Settings& settings() const { return settings_; }
    size_t fixed_sets() const { return fixed_; }
    Prefixed prefixed(const Unitary& V, size_t set) const;
    // none below the prefix scale
    std::optional<Level> level(const Prefixed& p, int l) const;
    // Walks composite level l on every prefix, each walk capped at its share of
    // the point budget; keeps the first cheapest composite within epsilon in
    // `best`. Returns the level's log entry.
    SearchLog walk_level(const Prefixed& p, int l, std::optional<Enumerator::Hit>& best);

  private:
    Result run_one(const Unitary& V, const PrefixSet& set, int scales, std::optional<int> first,
                   std::optional<int64_t> best_cost);

    // Prefixes of `prefixes` the rate rule affords per walk of `est` points: all, a
    // part (partial), or none.
    size_t affordable_prefixes(double est, size_t prefixes, int64_t best_cost) const;

    struct Composite {
        bool found = false;
        int64_t visited = 0;
    };
    // Walks composite scale K: one problem per prefix and determinant class at
    // scale K - set.scale (the first `count` problems), each capped at `leaf_cap`
    // leaves. Keeps the first
    // cheapest P R within epsilon in `best`; returns whether any was found and
    // the leaves visited.
    Composite walk_composite(const PrefixSet& set, const std::vector<Enumerator::Problem>& problems,
                             const Vec4& target, int K, int64_t leaf_cap,
                             std::optional<Enumerator::Hit>& best, size_t count, int sub = 0);

    double epsilon_;
    Price price_;
    Pricer pricer_;
    Settings settings_;
    PrefixSets sets_;
    size_t fixed_;   // sets_ tried by run()
    size_t copied_;  // sets_ already copied
    Enumerator enumerator_;
};

}  // namespace quicksynth
