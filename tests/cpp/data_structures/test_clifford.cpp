// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Clifford group modulo phase and its spelling tables.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/data_structures/clifford.hpp"

#include <doctest/doctest.h>

#include <set>
#include <string>
#include <vector>

namespace quicksynth::test {

TEST_CASE("each table spells each of the 24 Cliffords modulo phase once") {
    // A Clifford's Bloch image determines it modulo phase.
    const auto distinct = [](const std::vector<std::string>& words) {
        std::set<CliffordRotation> images;
        for (const std::string& w : words) images.insert(rotation_of(w));
        return images.size();
    };
    std::vector<std::string> s_h, tables;
    for (const auto& [rotation, word] : s_h_spellings()) s_h.push_back(word);
    for (const auto& [rotation, word] : spelling_tables().words) tables.push_back(word);
    CHECK(s_h.size() == 24);
    CHECK(distinct(s_h) == 24);
    CHECK(distinct(tables) == 24);
    CHECK(spelling_tables().cosets.size() == 3);
    CHECK(spelling_tables().sx.size() == 8);
}

TEST_CASE("a word's signed permutation is the product of its letters'") {
    Rng rng(7);
    for (int i = 0; i < 100; ++i) {
        const std::string a = random_word(rng, "HSsXYZ", 12), b = random_word(rng, "HSsXYZ", 12);
        CHECK(rotation_of(a + b) == rotation_product(rotation_of(a), rotation_of(b)));
        CHECK(rotation_product(rotation_of(a), rotation_transpose(rotation_of(a))) ==
              kIdentityRotation);
    }
    CHECK_FALSE(clifford_rotation('T'));
    CHECK(rotation_of("TQHt") == rotation_of("H"));
}

}  // namespace quicksynth::test
