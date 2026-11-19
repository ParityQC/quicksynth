// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Exact elements and words in floating point, against the frozen reference of
// golden_words.json: its complex matrices (CPython's complex arithmetic) and the
// distances of its spellings (mpmath at 40 digits). The core evaluates both in
// double precision, the distance from the shorter chord, so each is held to
// kTolerance rather than to the bit. Elements past 64 bits, which the reference
// holds without matrices, against their words multiplied out in long double.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/distance.hpp"
#include "quicksynth/data_structures/su2.hpp"
#include "quicksynth/exact_synthesis/normal_form.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <string>

namespace quicksynth::test {

namespace {

// Rounding in double precision: typically about 1e-16, measured at most
// 4.4e-16 against 40 digits (kDistanceError), and at most 3.5e-16 on these
// matrices and 2.2e-16 on these distances. The targets are decoded, not taken
// through an SVD, so the values are the same on every host.
constexpr double kTolerance = 4e-16;

}  // namespace

TEST_CASE("the matrices of the reference's elements") {
    for (const json& e : golden("words")["elements"]) {
        INFO(e["word"].get<std::string>());
        const CyclotomicMatrix element = matrix_of_row(e["element"]);
        const Unitary m = matrix_unitary(element);
        for (size_t i = 0; i < 8; ++i) CHECK(std::abs(m[i] - number(e["matrix"][i])) <= kTolerance);
    }
}

TEST_CASE("the distances of the reference's words past 64 bits") {
    for (const json& r : golden("words")["past_64_bits"]["rows"]) {
        const std::string word = r["word"];
        INFO(word);
        CHECK(precise_distance(word, unitary_of(letters_multiplied(word))) <= kTolerance);
    }
}

TEST_CASE("the spelling and the distance of the reference's words") {
    // random words against Haar targets, targets within 1e-16..1e-6 of the word
    // (where the trace formula would cancel every digit) and diagonal targets
    for (const json& s : golden("words")["spelling"]) {
        const std::string word = s["word"];
        INFO(word);
        const Unitary V = decode_target(s["target"]);
        const std::string simplified = s["simplified"];
        const double distance = number(s["distance"]);
        CHECK(normal_form(word) == simplified);
        CHECK(std::abs(precise_distance(simplified, V) - distance) <= kTolerance);
    }
}

}  // namespace quicksynth::test
