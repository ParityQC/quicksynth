// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// QuickSynth's search: the scale scan, then, unless it
// certified, the prefix split from the first scale not ruled out. The cheaper
// word wins.
//
// Fixed stopping rule: the split runs below split_below or when the scan finds
// nothing. Prefixes are the max_prefixes cheapest of the prefix scale, then
// every shallower scale as a fallback. With finer levels on Clifford+sqrt(T)
// the split starts past the scan's last level, and with level_search and a rate
// LevelSearch (level_search.hpp) runs both searches instead.
//
// Model stopping rule (stopping_rule.hpp): the split runs whenever the scan did not certify,
// from the first scale not walked in full, choosing prefix scales from 2 to
// max(3, prefix scale) by value. Below split_below its first composite scale is
// always walked. A target on stopping_fallback gets the fixed stopping rule.
#pragma once

#include "quicksynth/search_scales/level_search.hpp"
#include "quicksynth/search_scales/scan.hpp"
#include "quicksynth/search_scales/split.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/cost_ball.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"
#include "quicksynth/support/diagnostics.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace quicksynth {

constexpr double kSplitBelow = 1e-5;
constexpr int64_t kMaxPrefixes = 600;

// Below this epsilon the split uses prefix scale 3 or 4, a 3e8 budget and a
// prefix ball of kMaxPrefixes prefixes (split_radius).
constexpr double kDeepEpsilon = 1e-7;

// Below this epsilon Clifford+T splits at prefix scale 6 with every coset of
// that scale (split_radius, no max_prefixes cap). A composite scale K costs
// 3 est(K) / 4^k points with all 1.5 4^k cosets of scale k. The first answer
// lies ~5 scales (10 T) deeper per decade; a fixed budget and prefix scale
// reach ~2.5 scales deeper, so scale 4 at 3e8 stops short below 1e-8.
constexpr double kDeeperEpsilon = 1e-8;

// Returns the prefix scale: 2 for epsilon >= kDeepEpsilon, else 3 on
// Clifford+sqrt(T) and 4 on Clifford+T. Scale k shrinks the remainder region by
// 2^(deg k) but needs ~16^k (sqrt(T)) or ~4^k (T) prefixes. Empirical
// thresholds. Clifford+T below kDeeperEpsilon: 6.
int prefix_scale_for(double epsilon, const GateSet& gateset);

// Returns the split's point budget: 3e7 for epsilon >= kDeepEpsilon, else 3e8.
// Below kDeepEpsilon the box is held at its minimum half-thickness, so 3e7 stops
// short of the first hit on most targets. Empirical.
double split_budget_for(double epsilon);

class ScanThenSplit {
  public:
    struct Settings {
        ScaleScan::Settings scan;
        std::optional<int> prefix_scale;     // none: prefix_scale_for(epsilon)
        std::optional<double> split_budget;  // split point budget; none: split_budget_for
        std::optional<int64_t> max_prefixes = kMaxPrefixes;  // per scale; none: no cap
        // false: no shallower fallback sets. They serve targets near a Clifford,
        // and with no cap each would walk up to split_budget per composite scale.
        bool fallback_scales = true;
        double split_below = kSplitBelow;  // fixed stopping rule: split below this epsilon
        // fixed stopping rule: the split's rate rule (split.hpp); none: kFixedSplitScales
        std::optional<double> split_rate;
        bool split_partial = false;  // PrefixSplit::Settings::partial
        // finer levels (scan.fine_levels): composite levels walked from the first
        // hit; none: L (kFixedSplitScales - 1) + 1, the same depth
        std::optional<int> split_levels;
        // finer levels on Clifford+sqrt(T): LevelSearch at level_search_rate (none:
        // scan.rate; neither: the fixed depth); false: the fixed depth
        bool level_search = true;
        std::optional<double> level_search_rate;
        bool skip_setup_levels = true;  // LevelSearch::Settings
    };
    struct Result {
        std::optional<Enumerator::Hit> best;
        SearchSummary summary;
        std::vector<SearchLog> log;  // scan walks, then split composite scales
    };

    // `ball` supplies the split's prefixes (null: scan only). The model stopping rule
    // also takes stopping_extra_scales sets from `stopping_ball` (null: `ball`).
    ScanThenSplit(double epsilon, GateSet gateset, Price price, int threads,
                  std::shared_ptr<const CostBall> ball, Settings settings,
                  std::shared_ptr<const CostBall> stopping_ball);

    // Runs the scan, then the split as described above. `progress` goes to the
    // scan.
    Result run(const Unitary& V, const ScaleScan::Progress& progress = {});

  private:
    Price price_;
    double split_below_;
    std::optional<int> split_levels_;
    std::optional<double> level_search_rate_;  // none: no LevelSearch
    bool skip_setup_levels_;
    ScaleScan scan_;
    std::unique_ptr<PrefixSplit> split_;
};

// The split's prefix scales: prefix_scales(scale), or `scale` alone without
// fallback_scales.
std::vector<int> split_scales(const ScanThenSplit::Settings& settings, int scale);

}  // namespace quicksynth
