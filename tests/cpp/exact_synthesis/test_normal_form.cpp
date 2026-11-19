// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// normal_form, the merge_letters and t_as_qq identities, and words
// past 128 bits.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/exact_synthesis/normal_form.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/distance.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace quicksynth::test {

TEST_CASE("the spelling of small words") {
    CHECK(normal_form("") == "");
    CHECK(normal_form("XTX") == "t");
    CHECK(normal_form("XQX") == "q");
    CHECK(normal_form("SSSS") == "");
}

TEST_CASE("the spelling keeps the gate and its magic letters") {
    // and puts at most two Clifford letters before each magic letter
    Rng rng(4);
    for (int i = 0; i < 200; ++i) {
        const std::string w = random_word(rng, "HSXYZTQsqt", 51);
        INFO(w);
        const std::string s = normal_form(w);
        CHECK(precise_distance(s, matrix_unitary(word_matrix(w))) < 1e-14);
        CHECK(magic_counts(s) == magic_counts(w));
        size_t run = 0;
        for (char c : s) {
            if (std::string("TtQq").find(c) != std::string::npos) {
                CHECK(run <= 2);
                run = 0;
            } else {
                ++run;
            }
        }
    }
}

TEST_CASE("merging letters and spelling T as Q Q keep the element with its phase") {
    Rng rng(3);
    for (int i = 0; i < 300; ++i) {
        const std::string w = random_word(rng, "HSXYZTQ", 60);
        INFO(w);
        const CyclotomicMatrix e = word_matrix(w);
        CHECK(row_of(word_matrix(merge_letters(w))) == row_of(e));
        CHECK(row_of(word_matrix(t_as_qq(w))) == row_of(e));
        CHECK(t_as_qq(w).find('T') == std::string::npos);
    }
    CHECK(merge_letters("SSSS") == "");
    CHECK(merge_letters("TT") == "S");
    CHECK(merge_letters("QQ") == "T");
}

TEST_CASE("a word past 128 bits: exact spelling and distance") {
    // its numerator passes 2^120 and its Bloch image 2^128; the word times its
    // inverse is exactly the identity element (its spelling up to a phase), at
    // distance 0 up to rounding in double precision
    Rng rng(4);
    const std::vector<std::string> syllables{"HT", "HQ", "SHT", "HQX"};
    std::string word;
    for (int i = 0; i < 400; ++i) word += syllables[static_cast<size_t>(rng.below(4))];
    std::string inverse;
    for (auto c = word.rbegin(); c != word.rend(); ++c)
        inverse += *c == 'T' ? 't' : *c == 'Q' ? 'q' : *c == 'S' ? 's' : *c;
    CHECK_THROWS_AS(word_matrix(word), std::overflow_error);
    size_t widest = 0;
    for (const std::string& h : word_matrix_digits(word).hex)
        widest = std::max(widest, h.size() - (h[0] == '-' ? 1 : 0));
    CHECK(widest > 30);  // hex digits; past 2^120
    CHECK(magic_counts(normal_form(word)) == magic_counts(word));
    const CyclotomicMatrix identity = word_matrix("");
    CHECK(decimal_row(word_matrix_digits(word + inverse)) == decimal_row(digits_of(identity)));
    const std::vector<std::string> spelt_identity =
        decimal_row(word_matrix_digits(normal_form(word + inverse)));
    bool phase = false;  // the spelling keeps the gate up to a phase zeta^j
    for (int j = 0; j < 16; ++j)
        phase = phase || spelt_identity == decimal_row(digits_of(times_zeta(identity, j)));
    CHECK(phase);
    const Unitary I = unitary_of(kI);
    CHECK(precise_distance(word + inverse, I) <= 1e-15);
    CHECK(precise_distance(normal_form(word + inverse), I) <= 1e-15);
    const Unitary V = unitary_of(mul(rotation(0.3, 'x'), rotation(1.1, 'z')));
    const double spelt = precise_distance(normal_form(word), V);
    const double raw = precise_distance(word, V);
    CHECK(std::abs(spelt - raw) <= 1e-12 * raw);
}

TEST_CASE("a word past 2^1024: a finite distance") {
    // its coefficients pass the range of a double; the distance is that of its
    // letters in long double all the same (test_exact.py checks its matrix at 50
    // digits)
    Rng rng(5);
    const std::vector<std::string> syllables{"HT", "HQ", "SHT", "HQX"};
    std::string word;
    for (int i = 0; i < 3400; ++i) word += syllables[static_cast<size_t>(rng.below(4))];
    const CyclotomicMatrixDigits digits = word_matrix_digits(word);
    CHECK(digits.denominator_exponent > 2100);
    size_t widest = 0;
    for (const std::string& h : digits.hex)
        widest = std::max(widest, h.size() - (h[0] == '-' ? 1 : 0));
    CHECK(widest > 256);  // hexadecimal digits: past 2^1024
    CHECK(precise_distance(word, unitary_of(letters_multiplied(word))) <=
          1e-15 + letters_multiplied_error(word.size()));
    const Unitary V = unitary_of(mul(rotation(0.3, 'x'), rotation(1.1, 'z')));
    const double d = precise_distance(normal_form(word), V);
    CHECK(std::isfinite(d));
    CHECK(std::abs(d - precise_distance(word, V)) <= 1e-15);
}

TEST_CASE("a letter outside the alphabet and a target off unitary are refused") {
    CHECK_THROWS_WITH_AS(normal_form("HTR"), doctest::Contains("unknown letter 'R'"),
                         std::invalid_argument);
    CHECK_THROWS_WITH_AS(precise_distance("HTR", unitary_of(kI)),
                         doctest::Contains("unknown letter"), std::invalid_argument);
    // the distance reads the target's first column, so a shear of 1e-9 would
    // otherwise move it by about 1e-9 unnoticed
    const Matrix W = matrix_value(word_matrix("HTSHTHQHTSHtHQ"));
    for (long double shear : {1e-9L, 1e-3L}) {
        const Unitary U = unitary_of(mul(W, Matrix{1, shear, 0, 1}));
        CHECK_THROWS_WITH_AS(precise_distance("HT", U), doctest::Contains("not unitary"),
                             std::invalid_argument);
    }
    const Unitary near = unitary_of(mul(W, rotation(1e-9, 'x')));
    CHECK(std::abs(precise_distance("HTSHTHQHTSHtHQ", near) - 5e-10) <= 1e-15);
}

}  // namespace quicksynth::test
