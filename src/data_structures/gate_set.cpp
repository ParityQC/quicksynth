// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/gate_set.hpp"

#include <stdexcept>

namespace quicksynth {

GateSet GateSet::sqrt_t() { return GateSet{{1, 3, 5, 7}, 1}; }
GateSet GateSet::clifford_t() { return GateSet{{1, 3}, 2}; }

GateSet GateSet::from_name(const std::string& name) {
    if (name == "sqrt_t") return sqrt_t();
    if (name == "clifford_t") return clifford_t();
    throw std::invalid_argument("unknown gate set " + name);
}

}  // namespace quicksynth
