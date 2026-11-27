// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Prefix sets for the prefix split: per denominator scale
// k, the cost-ball cosets of scale k, cheapest first, with their inverses and
// normal forms.
#pragma once

#include "quicksynth/data_structures/cost_ball.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

namespace quicksynth {

// A prefix P's normal form: its cost, its syllables by type (0: T, 1: Q), the
// type of its first syllable when that has no connector, and per type b whether
// the set also holds a child P C b, which reaches every element P C b R sooner.
struct PrefixForm {
    int64_t cost = 0;
    int n_t = 0, n_q = 0;
    std::optional<int> bare_first;
    std::array<bool, 2> extended{};
};

// Prefixes of one denominator scale, their inverses P^dagger and their normal
// forms, cheapest first.
struct PrefixSet {
    int scale = 0;
    std::vector<CyclotomicMatrix> prefixes;
    std::vector<Unitary> inverses;
    std::vector<PrefixForm> forms;
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
