// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// QuickSynth's search on finer levels (SearchRegion::first_level): the scan and the
// prefix split, interleaved level by level.
//
// Before the first answer the scan walks its levels in order, as
// ScaleScan::run does. Without a scan answer the split walks the levels past
// the scan's last, trying the fixed prefix sets in order as PrefixSplit::run
// does. With skip_setup_levels it leaves out a level whose LLL set-ups would
// take longer than its walks, since the next level holds its elements, unless
// the next level is within a factor kSkipHeadroom of the point budget.
//
// After the first answer each round prices the scan's levels past its last (cov)
// and, for each of the split's prefixes, the composite levels past its own last
// and cov. A candidate's value is its expected saving (LevelGain) minus a rate
// times its predicted seconds: walk_seconds per walk plus the LLL set-ups of subs
// not yet built. The split's candidate is a batch: the LevelGain::batch prefixes
// of largest positive value, each to its most valuable level, priced together.
// The round walks the scan's level or the batch, whichever is worth more (ties
// to the scan); the search stops when neither is worth its time. Candidates stop
// at the point budgets and where a scale bound reaches the best cost. Split
// candidates after a scan answer only below split_below.
//
// Every answer is a walked element, checked against epsilon.
#pragma once

#include "quicksynth/search_scales/scan.hpp"
#include "quicksynth/search_scales/split.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/su2.hpp"
#include "quicksynth/support/diagnostics.hpp"

#include <optional>
#include <vector>

namespace quicksynth {

// A skipped level's next level stays this far under the split budget, so its
// per-walk leaf cap leaves head-room over a low rotated-region estimate.
constexpr double kSkipHeadroom = 4.0;

class LevelSearch {
  public:
    struct Settings {
        // magic states per predicted second; LevelGain::rates splits it
        double rate = 0.0;
        double split_below = 0.0;
        bool skip_setup_levels = true;
    };
    struct Result {
        std::optional<Enumerator::Hit> best;
        SearchSummary summary;
        std::vector<SearchLog> log;  // scan and split walks in the order walked
    };

    // `split` may be null (scan only). Both must use finer levels.
    LevelSearch(ScaleScan& scan, PrefixSplit* split, Price price, Settings settings);

    Result run(const Unitary& V, const ScaleScan::Progress& progress = {});

  private:
    ScaleScan& scan_;
    PrefixSplit* split_;
    Price price_;
    Settings settings_;
};

}  // namespace quicksynth
