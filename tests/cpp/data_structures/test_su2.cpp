// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Double-precision SU(2) helpers: unitary vectors, chord distance, element
// matrices and epsilon-ball measure.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <string>

namespace quicksynth::test {

TEST_CASE("every letter is unitary numerically") {
    for (char letter : std::string("HSXYZTQstq")) {
        INFO(letter);
        CHECK(unitarity_error(matrix_of(matrix_unitary(word_matrix(std::string(1, letter))))) <
              1e-15);
    }
}

TEST_CASE("an element's matrix and direction are its exact value's") {
    for (const json& e : golden("words")["elements"]) {
        const CyclotomicMatrix element = word_matrix(e["word"]);
        const Matrix exact = matrix_value(element);
        const Matrix got = matrix_of(matrix_unitary(element));
        for (int i = 0; i < 4; ++i) CHECK(std::abs(got[i] - exact[i]) < 1e-14L);
        const Vec4 d = matrix_direction(element);
        CHECK(gap_up_to_sign(d, su2_of(exact)) < 1e-14L);
    }
}

TEST_CASE("the SU(2) vector does not see the global phase") {
    Rng rng(0);
    for (int i = 0; i < 20; ++i) {
        const Matrix U = haar_unitary(rng);
        const Vec4 a = su2_vector(unitary_of(U));
        const Vec4 b = su2_vector(unitary_of(scaled(U, std::polar(1.0L, 0.7L))));
        CHECK(gap_up_to_sign(a, b) < 1e-12);
        CHECK(search_direction(unitary_of(U), 0, GateSet::sqrt_t()) == a);
    }
}

TEST_CASE("the chord distance is accurate where the trace formula is not") {
    const Unitary I = unitary_of(kI);
    for (double theta : {0.3, 1e-3, 1e-6, 1e-9}) {
        INFO(theta);
        const Matrix rz = rotation(theta, 'z');
        const double exact = std::abs(std::sin(theta / 2));
        const double chord = chord_distance(su2_vector(unitary_of(rz)), su2_vector(I));
        CHECK(std::abs(chord - exact) <= 1e-12 * exact);
        CHECK(std::abs(chord_diamond_distance(rz, kI) - exact) <= 1e-12 * exact);
        if (theta >= 1e-3) CHECK(std::abs(trace_diamond_distance(rz, kI) - exact) <= 1e-6 * exact);
    }
    CHECK(chord_distance(su2_vector(I), su2_vector(I)) == 0.0);
    CHECK(chord_distance(su2_vector(I), su2_vector(unitary_of(scaled(kI, -1)))) == 0.0);
}

TEST_CASE("the measure of an epsilon-ball is the window study's") {
    for (const json& pair : golden("window")["ball_measure"]) {
        const double eps = number(pair[0]), mu = number(pair[1]);
        INFO(eps);
        CHECK(close(ball_measure(eps), mu, 1e-14, 0.0));
    }
}

}  // namespace quicksynth::test
