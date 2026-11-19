// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Targets that are Clifford-equivalent to a z-rotation: U = C1 D C2 up to
// phase, C1 and C2 Cliffords and D = diag(z, z*). The Bloch image of such a U
// maps some coordinate axis e_i to a signed axis s e_j. With R(C2) e_i = e_z
// and R(C1) e_z = s e_j, R(D) = R(C1)^T R(U) R(C2)^T fixes e_z, so D is
// diagonal; no separate X case is needed.
//
// The word of U is left + w + right for a word w of D: word_matrix multiplies
// the letters left to right.
#pragma once

#include "quicksynth/data_structures/su2.hpp"

#include <array>
#include <optional>
#include <string>

namespace quicksynth {

struct DiagonalForm {
    std::string left;   // C1, its shortest word (spelling_tables().words)
    std::string right;  // C2
    Unitary d{};        // diag(z, z*), |z| = 1 to rounding
    // Chord distance from U to C1 D C2, both in doubles: the off-axis part of
    // U plus rounding (about 1e-16 for an exact Clifford-equivalent target).
    double residual = 0.0;
};

// R(U)_jk = tr(sigma_j U sigma_k U^dagger) / 2 in doubles, rows and columns
// x, y, z: the convention of clifford_rotation.
using BlochRotation = std::array<std::array<double, 3>, 3>;
BlochRotation bloch_rotation(const Unitary& U);

// The diagonal form of U with the least residual over the axis pairs (i, j)
// with |R(U)_ji| >= 1/2, ties to the first pair in the order z, x, y of i,
// then of j; nullopt if that residual exceeds `tolerance`.
std::optional<DiagonalForm> diagonal_form(const Unitary& U, double tolerance);

}  // namespace quicksynth
