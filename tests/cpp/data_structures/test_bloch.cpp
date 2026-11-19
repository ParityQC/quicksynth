// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Bloch images over Z[lambda]: normal-form residues and exponents
// [paper: prop:blochexp], and reference images in int64 and arbitrary digits.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/data_structures/bloch.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"

#include <doctest/doctest.h>

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace quicksynth::test {

namespace {

using Residue = std::array<std::array<int, 3>, 3>;

// Exponent and image mod lambda; Z[lambda]/lambda = F_2 is the parity of each
// entry's constant digit.
std::pair<int64_t, Residue> residue(const std::string& word) {
    const CyclotomicMatrix e = word_matrix(word);
    const BlochDigits b = bloch_digits(e.matrix, e.denominator_exponent);
    Residue r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            const std::string& d = at(b, i, j, 0);
            r[i][j] = (d.back() - '0') % 2;
        }
    return {b.exponent, r};
}

std::vector<std::string> image_digits(const json& list) {
    std::vector<std::string> out;
    for (const json& d : list) out.push_back(decimal(d));
    return out;
}

std::vector<std::string> image_digits(const BlochDigits& b) {
    return std::vector<std::string>(b.m.begin(), b.m.end());
}

}  // namespace

TEST_CASE("every syllable reduces to the same rank-one residue") {
    // lambda^e B(b) = J = (e1 + e2)(e1 + e2)^T mod lambda
    const Residue J{{{1, 1, 0}, {1, 1, 0}, {0, 0, 0}}};
    // clang-format off
    for (const auto& [syllable, exponent] :
         std::vector<std::pair<std::string, int>>{{"Q", 3}, {"QX", 3}, {"QS", 3}, {"QQQ", 3},
                                                  {"QSX", 3}, {"T", 2}, {"TX", 2}, {"TS", 2}}) {
        const std::string& label = syllable;  // C++17: INFO cannot capture a binding
        INFO(label);
        CHECK(residue(syllable) == std::make_pair(int64_t{exponent}, J));
    }
    // clang-format on
}

TEST_CASE("the connectors send z off its axis and K = <S X> does not") {
    for (const auto& [connector, image] :
         std::vector<std::pair<std::string, int>>{{"H", 0}, {"SH", 1}, {"S", 2}, {"X", 2}}) {
        const std::string& label = connector;  // C++17: INFO cannot capture a binding
        INFO(label);
        const auto [k, r] = residue(connector);
        CHECK(k == 0);
        for (int i = 0; i < 3; ++i) CHECK(r[i][2] == (i == image ? 1 : 0));
    }
}

TEST_CASE("the exponent of a normal form is three per sqrt(T) and two per T") {
    for (uint64_t seed = 0; seed < 20; ++seed) {
        Rng rng(seed);
        const std::vector<std::string> first{"", "H", "SH"}, connectors{"H", "SH"},
            syllables{"Q", "T", "QQQ"}, last{"", "H", "S", "HS"};
        std::string word = first[static_cast<size_t>(rng.below(3))];
        int64_t expected = 0;
        const int64_t n = 1 + rng.below(29);
        for (int64_t i = 0; i < n; ++i) {
            if (i) word += connectors[static_cast<size_t>(rng.below(2))];
            const std::string& s = syllables[static_cast<size_t>(rng.below(3))];
            word += s;
            expected += s == "T" ? 2 : 3;
        }
        word += last[static_cast<size_t>(rng.below(4))];
        INFO(word);
        const auto [k, r] = residue(word);
        CHECK(k == expected);
        bool nonzero = false;
        for (const auto& row : r)
            for (int x : row) nonzero = nonzero || x != 0;
        CHECK(nonzero);
    }
}

TEST_CASE("the Bloch images of the reference's elements") {
    for (const json& e : golden("words")["elements"]) {
        INFO(e["word"].get<std::string>());
        const CyclotomicMatrix element = matrix_of_row(e["element"]);
        const BlochDigits b = bloch_digits(element.matrix, element.denominator_exponent);
        CHECK(image_digits(b) == image_digits(e["bloch"]));
        CHECK(b.exponent == e["bloch_exponent"]);
        CHECK(bloch_exponent(element.matrix, element.denominator_exponent) == b.exponent);
    }
}

TEST_CASE("the Bloch images of elements past 64 bits") {
    for (const json& r : golden("words")["past_64_bits"]["rows"]) {
        INFO(r["word"].get<std::string>());
        const CyclotomicMatrixDigits element = digits_of_row(r["element"]);
        const BlochDigits b = bloch_digits(element);
        CHECK(image_digits(b) == image_digits(r["bloch"]));
        CHECK(b.exponent == r["bloch_exponent"]);
        CHECK(bloch_exponent(element) == b.exponent);
    }
}

TEST_CASE("an unreduced element past 64 bits has the image of the element") {
    const json& elements = golden("words")["elements"];
    for (size_t i = 0; i < elements.size(); i += 5) {
        INFO(elements[i]["word"].get<std::string>());
        const CyclotomicMatrix element = matrix_of_row(elements[i]["element"]);
        const CyclotomicMatrixDigits wide = widened(element, 70);
        CHECK(image_digits(bloch_digits(wide)) ==
              image_digits(bloch_digits(element.matrix, element.denominator_exponent)));
        CHECK(bloch_exponent(wide) == bloch_exponent(element.matrix, element.denominator_exponent));
    }
}

}  // namespace quicksynth::test
