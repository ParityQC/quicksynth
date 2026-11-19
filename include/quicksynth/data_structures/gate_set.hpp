// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Supported gate sets: the Cliffords plus diag(1, zeta_n), Q for
// Clifford+sqrt(T) (n = 16) or T for Clifford+T (n = 8). Z[zeta_n] sits in
// Z[zeta16] as the powers zeta16^(step i), step = 16 / n.
#pragma once

#include <string>
#include <vector>

namespace quicksynth {

struct GateSet {
    // Embedding exponents a, zeta_n -> zeta_n^a; one per conjugate pair, identity first.
    std::vector<int> places;
    int step = 1;

    static GateSet sqrt_t();
    static GateSet clifford_t();
    // "sqrt_t" or "clifford_t"; std::invalid_argument for any other name.
    static GateSet from_name(const std::string& name);
    const char* name() const { return step == 1 ? "sqrt_t" : "clifford_t"; }

    int degree() const { return 8 / step; }
    int dim() const { return 2 * degree(); }  // two ring elements
    int conj_dim() const { return dim() - 4; }
    int num_places() const { return static_cast<int>(places.size()); }
    // Number of determinant classes p (det-exponent parity) in
    // zeta^(step l) M(u) Q^p / sqrt2^K, det M(u) = 2^K.
    static constexpr int kDeterminantClasses = 2;
};

}  // namespace quicksynth
