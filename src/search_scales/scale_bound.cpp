// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/search_scales/scale_bound.hpp"

#include <algorithm>
#include <limits>

namespace quicksynth {

double scan_scale_bound(const GateSet& gateset, const Price& price, int K) {
    const double t = static_cast<double>(price.t);
    const double q = static_cast<double>(price.q);
    if (gateset.step != 1) return t * std::max(0, 2 * K - 3);
    return std::min(t / 2, q / 3) * (4 * K - 7);
}

std::optional<double> split_scale_bound(const GateSet& gateset, const Price& price, int K) {
    if (gateset.step == 1) return std::nullopt;
    return static_cast<double>(price.t) * std::max(0, 2 * K - 3);
}

int64_t integer_scale_bound(const GateSet& gateset, const Price& price, int K) {
    if (gateset.step != 1) return price.t * std::max(0, 2 * K - 3);
    const int64_t t = t_syllable_cost(gateset, price);
    const int64_t need = 4 * static_cast<int64_t>(K) - 7;
    if (need <= 0) return 0;
    int64_t best = std::numeric_limits<int64_t>::max();
    for (int64_t n_q = 0; n_q <= need / 3 + 1; ++n_q)
        best = std::min(best, price.q * n_q + t * ((std::max<int64_t>(0, need - 3 * n_q) + 1) / 2));
    return best;
}

}  // namespace quicksynth
