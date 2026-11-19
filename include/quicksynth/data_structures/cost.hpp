// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Gate prices and word costs, used by the descent
// (descent.hpp) and the Pricer (price.hpp).
#pragma once

#include <cstdint>

namespace quicksynth {

// Cost of a word with its T and sqrt(T) counts.
struct Cost {
    int64_t value = 0, n_t = 0, n_q = 0;
};

// Price per T (t) and per sqrt(T) (q) gate; defaults to the catalysis price.
struct Price {
    int64_t t = 1, q = 3;

    // Clifford+T price; q = 3 is not searched (the T-only ball has no sqrt(T)).
    static Price clifford_t(int64_t t) { return {t, 3}; }

    // The cost of n_t T and n_q sqrt(T) gates.
    Cost cost(int64_t n_t, int64_t n_q) const { return {t * n_t + q * n_q, n_t, n_q}; }
};

}  // namespace quicksynth
