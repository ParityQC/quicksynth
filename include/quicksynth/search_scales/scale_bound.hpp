// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Lower bounds on the cost of any element of denominator scale K or deeper, at
// the price (t, q): t (2K - 3) on Clifford+T and min(t/2, q/3) (4K - 7) on
// Clifford+sqrt(T). A scan stops, with an optimal answer, once the next scale's
// bound reaches its best cost.
#pragma once

#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/gate_set.hpp"

#include <algorithm>
#include <cstdint>
#include <optional>

namespace quicksynth {

// Cost of a T syllable: min(t, 2q) on Clifford+sqrt(T) (Q Q = T), else t.
inline int64_t t_syllable_cost(const GateSet& gateset, const Price& price) {
    return gateset.step == 1 ? std::min(price.t, 2 * price.q) : price.t;
}

// Fixed-rule bound: t max(0, 2K - 3) on Clifford+T; min(t/2, q/3) (4K - 7)
// on Clifford+sqrt(T), since t n_T + q n_Q >= min(t/2, q/3) (2 n_T + 3 n_Q).
double scan_scale_bound(const GateSet& gateset, const Price& price, int K);

// Prefix-split bound: the Clifford+T bound; none on Clifford+sqrt(T), where the
// split stops on scale count and budget.
std::optional<double> split_scale_bound(const GateSet& gateset, const Price& price, int K);

// Model-rule bound c_int(K): t max(0, 2K - 3) on Clifford+T; on
// Clifford+sqrt(T) the least t n_T + q n_Q with 2 n_T + 3 n_Q >= 4K - 7 (t from
// t_syllable_cost). At least scan_scale_bound.
int64_t integer_scale_bound(const GateSet& gateset, const Price& price, int K);

}  // namespace quicksynth
