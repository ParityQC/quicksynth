// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Word elements against golden_words.json (earlier Python reference).

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"

#include <doctest/doctest.h>

#include <stdexcept>
#include <string>

namespace quicksynth::test {

namespace {

// zeta^j I, 0 <= j < 8.
CyclotomicMatrix scalar(int j) {
    CyclotomicMatrix e;
    e.matrix[0][0][j] = 1;
    e.matrix[1][1][j] = 1;
    return e;
}

bool same(const CyclotomicMatrix& a, const CyclotomicMatrix& b) { return row_of(a) == row_of(b); }

}  // namespace

TEST_CASE("phase words are the scalars they spell") {
    CHECK(same(word_matrix("HSHSHS"), scalar(2)));  // (HS)^3 = omega
    CHECK(same(word_matrix("XQXQ"), scalar(1)));    // X Q X Q = zeta
    CHECK(same(word_matrix("Qq"), scalar(0)));
    CHECK(same(word_matrix(""), scalar(0)));
}

TEST_CASE("the elements of every letter and a few hundred words are the reference's") {
    for (const json& e : golden("words")["elements"]) {
        const std::string word = e["word"];
        INFO(word);
        CHECK(row_of(word_matrix(word)) == e["element"]);
        CHECK(decimal_row(word_matrix_digits(word)) == decimal_row(digits_of_row(e["element"])));
    }
}

TEST_CASE("an element past 64 bits is an overflow in int64 and exact in digits") {
    for (const json& r : golden("words")["past_64_bits"]["rows"]) {
        const std::string word = r["word"];
        INFO(word);
        CHECK_THROWS_AS(word_matrix(word), std::overflow_error);
        CHECK(decimal_row(word_matrix_digits(word)) == decimal_row(digits_of_row(r["element"])));
    }
}

TEST_CASE("a letter outside the alphabet is named") {
    CHECK_THROWS_WITH_AS(word_matrix("HTR"), doctest::Contains("unknown letter 'R'"),
                         std::invalid_argument);
    CHECK_THROWS_AS(word_matrix_digits("HTR"), std::invalid_argument);
}

TEST_CASE("compose is the reduced product of words") {
    Rng rng(41);
    for (int i = 0; i < 200; ++i) {
        const std::string a = random_word(rng, "HSXYZTQstq", 40);
        const std::string b = random_word(rng, "HSXYZTQstq", 40);
        INFO(a, " ", b);
        CHECK(same(compose(word_matrix(a), word_matrix(b)), word_matrix(a + b)));
    }
}

TEST_CASE("the first column of a word's element has unit norm at its own scale only") {
    Rng rng(43);
    for (int i = 0; i < 200; ++i) {
        const std::string word = random_word(rng, "HSXYZTQstq", 40);
        INFO(word);
        const CyclotomicMatrix e = word_matrix(word);
        CHECK(unit_norm(e.matrix[0][0], e.matrix[1][0], e.denominator_exponent));
        CHECK_FALSE(unit_norm(e.matrix[0][0], e.matrix[1][0], e.denominator_exponent + 1));
        // u1 + zeta^t changes the norm's constant coefficient by an odd number
        Z16 shifted = e.matrix[0][0];
        ++shifted[i % kZetaDegree];
        CHECK_FALSE(unit_norm(shifted, e.matrix[1][0], e.denominator_exponent));
    }
}

}  // namespace quicksynth::test
