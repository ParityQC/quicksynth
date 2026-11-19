// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Gate word rewriting. A word is a matrix product, read left to right, over
// H S X Y Z T Q and the inverses s t q; T, Q, t, q are magic letters.
//
// normal_form() writes a word as A_1 M_1 ... A_n M_n C: M_i magic, A_i in
// {"", "H", "SH"} (left coset representatives of <S, X> in the Clifford group),
// C the shortest word of the closing Clifford. A Clifford before M factors as
// A D with D in <S, X>, and D M = M^(+-1) D, so only A is written. Cliffords are
// tracked as CliffordRotation (clifford.hpp, modulo phase), in a fixed
// search order. The result is verified to match the word up to a phase zeta^j.
//
// merge_letters() and t_as_qq() apply exact identities, phase included.
#pragma once

#include <string>

namespace quicksynth {

// Applies S S S S = H H = X X = Y Y = Z Z = 1, T T = S, Q Q = T, S S = Z and
// orders diagonal letters T, Q, S, Z, until a fixed point.
std::string merge_letters(std::string word);

// Replaces every T by Q Q.
std::string t_as_qq(const std::string& word);

// Normal-form spelling; std::invalid_argument on an unknown letter.
std::string normal_form(const std::string& word);

}  // namespace quicksynth
