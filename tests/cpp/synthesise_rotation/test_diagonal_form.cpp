// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Diagonal forms against targets built in long double: every Clifford pair
// around a z-rotation, the x and y axes, anti-diagonal targets, a measured
// off-axis tilt, Haar targets, and the word order left + w + right.

#include "quicksynth/synthesise_rotation/diagonal_form.hpp"

#include <doctest/doctest.h>

#include "fixtures.hpp"
#include "oracles.hpp"
#include "quicksynth/data_structures/clifford.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/distance.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

namespace quicksynth::test {

namespace {

const double kPi = std::acos(-1.0);

// phi with D = Rz(phi) = diag(exp(-i phi / 2), exp(i phi / 2)) up to sign.
double z_angle(const Unitary& d) { return -2 * std::atan2(d[1], d[0]); }

// Distance of x to the nearest multiple of `period`.
double off_grid(double x, double period) { return std::fabs(std::remainder(x, period)); }

// R e_z and R^T e_z of a Clifford, as signed axes.
std::array<int, 3> column_z(const CliffordRotation& p) { return {p[0][2], p[1][2], p[2][2]}; }
std::array<int, 3> row_z(const CliffordRotation& p) { return p[2]; }

std::array<int, 3> negated(std::array<int, 3> a) {
    for (int& x : a) x = -x;
    return a;
}

}  // namespace

TEST_CASE("the Bloch rotation of each Clifford letter is its Bloch image") {
    for (char letter : std::string("HSsXYZ")) {
        INFO(letter);
        const BlochRotation r =
            bloch_rotation(unitary_of(letters_multiplied(std::string(1, letter))));
        const CliffordRotation p = *clifford_rotation(letter);
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) CHECK(std::fabs(r[j][k] - p[j][k]) < 1e-15);
    }
}

TEST_CASE("every Clifford pair around a z-rotation is recovered") {
    const auto& words = spelling_tables().words;
    std::mt19937_64 rng(20261116);
    std::uniform_real_distribution<double> angle(-kPi, kPi);
    // no structured bindings: INFO's lambda cannot capture them in C++17 (Clang)
    for (const auto& entry1 : words)
        for (const auto& entry2 : words)
            for (int trial = 0; trial < 3; ++trial) {
                const CliffordRotation& p1 = entry1.first;
                const CliffordRotation& p2 = entry2.first;
                const std::string& w1 = entry1.second;
                const std::string& w2 = entry2.second;
                const double theta = angle(rng);
                INFO(w1, " Rz(", theta, ") ", w2);
                const Unitary U = unitary_of(
                    mul(mul(letters_multiplied(w1), rotation(theta, 'z')), letters_multiplied(w2)));
                const std::optional<DiagonalForm> f = diagonal_form(U, 1e-14);
                REQUIRE(f.has_value());
                CHECK(f->residual <= 1e-15);
                CHECK(f->d[2] == 0);
                CHECK(f->d[4] == 0);
                // Unique up to the z-axis stabiliser S^k and an X on both sides,
                // which turns theta into -theta.
                const CliffordRotation q1 = rotation_of(f->left), q2 = rotation_of(f->right);
                const bool same = column_z(q1) == column_z(p1);
                CHECK(column_z(q1) == (same ? column_z(p1) : negated(column_z(p1))));
                CHECK(row_z(q2) == (same ? row_z(p2) : negated(row_z(p2))));
                const double phi = z_angle(f->d);
                CHECK(off_grid(phi - (same ? theta : -theta), kPi / 2) < 1e-12);
            }
}

TEST_CASE("x and y rotations and anti-diagonal targets have a diagonal form") {
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> angle(-kPi, kPi);
    for (int trial = 0; trial < 50; ++trial) {
        const double theta = angle(rng), alpha = angle(rng), beta = angle(rng);
        const Matrix anti{Complex(0), std::polar(1.0L, static_cast<long double>(alpha)),
                          std::polar(1.0L, static_cast<long double>(beta)), Complex(0)};
        for (const Matrix& m :
             {rotation(theta, 'x'), rotation(theta, 'y'), anti, mul(kX, rotation(theta, 'z'))}) {
            const std::optional<DiagonalForm> f = diagonal_form(unitary_of(m), 1e-14);
            REQUIRE(f.has_value());
            CHECK(f->residual <= 1e-15);
        }
    }
}

TEST_CASE("an off-axis tilt is measured as the residual") {
    // Rz(theta) Rx(delta) has SU(2) vector (a, b), |b| = sin(delta / 2); the
    // nearest diagonal a / |a| is at chord distance exactly sin(delta / 2).
    std::mt19937_64 rng(11);
    std::uniform_real_distribution<double> angle(-kPi, kPi);
    const auto& words = spelling_tables().words;
    for (double delta : {1e-12, 1e-9, 1e-6, 1e-3}) {
        for (int trial = 0; trial < 20; ++trial) {
            const std::string& w1 = words[rng() % words.size()].second;
            const std::string& w2 = words[rng() % words.size()].second;
            const double theta = angle(rng);
            INFO(delta, " ", w1, " ", theta, " ", w2);
            const Unitary U = unitary_of(
                mul(mul(mul(letters_multiplied(w1), rotation(theta, 'z')), rotation(delta, 'x')),
                    letters_multiplied(w2)));
            const double expected = std::sin(delta / 2);
            const std::optional<DiagonalForm> f = diagonal_form(U, 2 * expected);
            REQUIRE(f.has_value());
            CHECK(close(f->residual, expected, 1e-6, 1e-15));
            CHECK_FALSE(diagonal_form(U, 0.5 * expected).has_value());
        }
    }
}

TEST_CASE("Haar targets have no diagonal form") {
    for (size_t i = 0; i < 20; ++i) {
        INFO("calibration target ", i);
        CHECK_FALSE(diagonal_form(calibration_targets()[i], 1e-3).has_value());
    }
}

TEST_CASE("the word of a diagonal form is left then the diagonal then right") {
    // Exact targets C_a M C_b, M a magic letter; D is M or M^-1 (an X on both
    // sides inverts it) times the z-stabiliser S^b, so a word x of D from
    // M^(0, +-1) S^b closes the form exactly.
    const auto& words = spelling_tables().words;
    int reversed_differs = 0;
    for (const std::string magic : {"T", "Q"})
        for (const auto& entry : words)
            for (size_t b = 0; b < words.size(); b += 5) {
                const std::string& wa = entry.second;
                const std::string& wb = words[b].second;
                INFO(wa, magic, wb);
                const Unitary U = matrix_unitary(word_matrix(wa + magic + wb));
                const std::optional<DiagonalForm> f = diagonal_form(U, 1e-14);
                REQUIRE(f.has_value());
                std::string best;
                double gap = 1;
                const std::string inverse(1, static_cast<char>(magic[0] - 'A' + 'a'));
                for (const std::string& m : {std::string(), magic, inverse})
                    for (int s = 0; s < 4; ++s) {
                        const std::string x = m + std::string(s, 'S');
                        const double g = chord_distance(su2_vector(matrix_unitary(word_matrix(x))),
                                                        su2_vector(f->d));
                        if (g < gap) {
                            gap = g;
                            best = x;
                        }
                    }
                REQUIRE(gap < 1e-15);
                CHECK(precise_distance(f->left + best + f->right, U) < 1e-15);
                if (precise_distance(f->right + best + f->left, U) > 1e-3) ++reversed_differs;
            }
    CHECK(reversed_differs > 0);
}

}  // namespace quicksynth::test
