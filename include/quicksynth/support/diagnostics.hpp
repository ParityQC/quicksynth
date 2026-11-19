// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Diagnostics reported with an answer, never read by the search: walk log,
// answering-stage summary, verbose walk line.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace quicksynth {

// Search that found the answer; rotation: the z-rotation path (stage 0).
enum class Source { scan, split, rotation };

// Counts of one walk of the z-rotation path (synthesise_rotation/).
struct RotationWalk {
    int64_t candidates = 0;  // new body points, after the exact body test
    int64_t solvable = 0;    // candidates whose norm equation has a solution
    int64_t undecided = 0;   // candidates whose norm equation the solver left open
    int64_t capped = 0;      // solvable candidates with more choices than the cap
    int64_t pricings = 0;    // completions within epsilon, priced
};

// One scan walk, one split composite scale or one walk of the z-rotation path;
// calibrates the time model.
struct SearchLog {
    int scale = 0;  // K, or the composite scale
    std::optional<int> det_class;
    std::optional<int> prefix_scale;
    int sub = 0;  // sublattice delta^sub of the (remainder) scale; finer levels only
    double seconds = 0.0;
    double estimated = 0.0;  // estimated points: est(K), or walks x est(K - k)
    int64_t points = 0;      // leaves walked
    // best cost after this step; includes the scan's answer except in
    // PrefixSplit::run
    std::optional<int64_t> best_cost;
    std::optional<RotationWalk> rotation;  // the z-rotation path's walks only
};

// What QuickSynth's answering stage searched and which search answered.
struct SearchSummary {
    std::vector<int> scales;  // walked scales, increasing, each once
    std::optional<int> first_hit_scale, last_full_scale;
    double lower_bound = 0.0;       // lower bound on the optimum
    double estimated_points = 0.0;  // scan work in stopping-rule units
    std::optional<int64_t> scan_cost;
    bool scan_certified = false;  // every scale below the answer's bound completed
    bool split_ran = false;
    Source source = Source::scan;
};

// Progress line for one scan walk; lists the eight lowest hit costs.
std::string walk_line(int K, int det_class, int64_t points, std::vector<int64_t> costs,
                      double estimated, double seconds);

}  // namespace quicksynth
