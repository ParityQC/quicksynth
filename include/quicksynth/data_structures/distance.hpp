// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Exact elements and words evaluated in floating point.
//
// precise_distance() is the diamond distance between the word's unitary W and
// the target U, which arrives unitary (the caller's nearest_unitary): W from
// its exact element in double precision (matrix_unitary), and the distance
// from the shorter chord between the two SU(2) vectors (chord_distance), as
// the searches measure it. Its error is absolute, typically about 1e-16 and
// measured at most 4.4e-16 at every distance (kDistanceError);
// sqrt(1 - |tr(U^dagger W)|^2 / 4) would lose every digit below 1e-8 to
// cancellation. A distance within about 5e-16 of epsilon therefore does not
// decide which side of epsilon the word lies on.
//
// The SU(2) vector of the target is that of its first column, so the target
// must be unitary: one farther than kTargetUnitarity from it is refused.
#pragma once

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <array>
#include <string>

namespace quicksynth {

// The error of precise_distance() the final check allows for (accept() in
// synthesise.cpp takes an answer within epsilon - kDistanceError). It is
// measured, not bounded: against a 40-digit evaluation the error was at most
// 4.4e-16 over 33,000 word and target pairs, typically about 1e-16, and flat
// in the distance and the length of the word. The searches run at epsilon -
// kEpsilonMargin (1e-15), which must exceed it.
constexpr double kDistanceError = 5e-16;

// The largest ||U^dagger U - I|| (Frobenius) of a target precise_distance()
// takes. nearest_unitary's polar factor is within about 2e-15 of unitary; a
// target farther off than this would give a distance off by about as much.
constexpr double kTargetUnitarity = 1e-12;

// The distance from the word's unitary to the unitary target (see the top);
// std::invalid_argument on a letter outside the alphabet, a target entry that
// is not finite, or a target farther than kTargetUnitarity from unitary.
double precise_distance(const std::string& word, const Unitary& target);

}  // namespace quicksynth
