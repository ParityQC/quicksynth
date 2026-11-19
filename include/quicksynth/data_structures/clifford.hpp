// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The Clifford group modulo phase as Bloch images.
#pragma once

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace quicksynth {

// Bloch-sphere rotation of a Clifford (modulo phase): a 3 x 3 signed permutation
// matrix.
using CliffordRotation = std::array<std::array<int, 3>, 3>;

constexpr CliffordRotation kIdentityRotation{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};

CliffordRotation rotation_product(const CliffordRotation& a, const CliffordRotation& b);
CliffordRotation rotation_transpose(const CliffordRotation& a);

// The Bloch image of a Clifford letter H S s X Y Z; none for any other letter.
std::optional<CliffordRotation> clifford_rotation(char letter);

// Bloch image of a word, skipping non-Clifford letters.
CliffordRotation rotation_of(const std::string& word);

// Lookup tables for normal_form (normal_form.hpp), which writes a word as
// A_1 M_1 ... A_n M_n C. It accumulates the Clifford P between magic letters.
// At a magic letter M it splits P = A D, with A a coset representative and D in
// K = <S, X>, writes A M and moves D past M. S commutes with M and
// X M X = M^-1 up to phase, so M is inverted when D contains an X. The Clifford
// left at the end is written as the word C.
struct SpellingTables {
    // The left coset representatives "", "H", "SH" of K in the 24 Cliffords
    // (24 = 3 x 8), each with its Bloch image.
    std::vector<std::pair<std::string, CliffordRotation>> cosets;
    // The 8 elements X^b S^i of K, each with b: whether it holds an X.
    std::vector<std::pair<CliffordRotation, bool>> sx;
    // The 24 Cliffords, each with the first word to reach it breadth first over
    // H S s X Z Y: a shortest word, the ties broken by that letter order, which
    // the recorded spellings depend on.
    std::vector<std::pair<CliffordRotation, std::string>> words;
};
// Built on the first call and kept; thread-safe.
const SpellingTables& spelling_tables();

// A Clifford P split as P = A D, A a coset representative and D in K.
struct CosetSplit {
    std::string representative;  // A: "", "H" or "SH"
    CliffordRotation remainder;  // D
    bool holds_x = false;        // whether D holds an X
};
// Every Clifford has exactly one split.
CosetSplit split_coset(const CliffordRotation& p);

// The shortest word of a Clifford (SpellingTables::words).
const std::string& clifford_word(const CliffordRotation& p);

// The 24 Cliffords, each with its first spelling among S^l, then S^l H S^k,
// S^l H S^k H S^m and S^l H S^k H S^m H S^n, in that order: the fewest H, and
// the S between them in runs that merge_letters folds (SS = Z, SSSS = 1). The
// descent spells its leftover Clifford with it.
const std::vector<std::pair<CliffordRotation, std::string>>& s_h_spellings();

// Connectors moving column pairs {0, 1}, {0, 2}, {1, 2} to the front; the descent's
// kColumnOrders and syllable_generators (cost_ball.hpp) follow this order.
inline constexpr const char* kConnectors[3] = {"", "SH", "HSSS"};

}  // namespace quicksynth
