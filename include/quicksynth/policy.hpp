// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Constants and escalation stages of the one-call synthesis.
#pragma once

#include "quicksynth/search_scales/scan_then_split.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/data_structures/gate_set.hpp"

#include <string>
#include <vector>

namespace quicksynth {

// Double-precision distances err by about 1e-16, so the searches run
// at epsilon minus this, and the final check of an answer at epsilon minus
// kDistanceError (5e-16, distance.hpp). It must exceed kDistanceError, or
// the check rejects search hits (static_assert in synthesise.cpp).
constexpr double kEpsilonMargin = 1e-15;
constexpr int kPrefixRadius = 6;  // split prefix ball radius down to 1e-7
// Last stage's split budget on Clifford+T below kDeeperEpsilon: composite scale
// 45 takes 4.1e9 points at 1e-9 and 1.3e10 at 1.47e-9.
constexpr double kDeeperLastBudget = 1.5e10;

// One stage's settings; `changes` is its diff from stage 1 as a Python dict,
// printed under verbose.
struct Stage {
    ScanThenSplit::Settings settings;
    std::string changes;
};

// Returns the escalation stages; `first` is stage 1.
// 1. QuickSynth defaults.
// 2. Each deeper prefix scale up to 4 (on Clifford+T below 1e-7 one past the
//    split's: 5, or 7 below kDeeperEpsilon), split budget 3e8. Near a Clifford
//    the low scales hold no answer.
// 3. Default prefix scale, split budget 3e9, scan budget 1e9. Clifford+T below
//    kDeeperEpsilon: prefix scale 8 (98304 cosets), split budget kDeeperLastBudget.
// Stages 2+ use the fixed stopping rule, the default point budget and no prefix cap.
std::vector<Stage> stages(double epsilon, const GateSet& gateset,
                          const ScanThenSplit::Settings& first);

}  // namespace quicksynth
