// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Prefix sets for the prefix split: per denominator scale
// k, the cost-ball cosets of scale k, cheapest first, with their inverses.
#pragma once

#include "quicksynth/data_structures/cost_ball.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace quicksynth {

// Prefixes of one denominator scale and their inverses P^dagger, cheapest first.
struct PrefixSet {
    int scale = 0;
    std::vector<CyclotomicMatrix> prefixes;
    std::vector<Unitary> inverses;
};

// Returns k, k - 1, ..., 1. Shallower scales are fallbacks for targets with no
// prefix of scale k.
std::vector<int> prefix_scales(int scale);

// Returns one prefix set per entry of `scales`, in order. Levels are read from
// cost 0 until every set holds `cap` prefixes; sets are then cut to `cap`
// (none: all levels, no cut). Only kept prefixes are multiplied out.
std::vector<PrefixSet> collect_prefixes(const CostBall& ball, const std::vector<int>& scales,
                                        std::optional<int64_t> cap);

}  // namespace quicksynth
