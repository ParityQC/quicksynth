// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Descent words and prices against golden_words.json and golden_balls.json.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/data_structures/bloch.hpp"
#include "quicksynth/data_structures/cost_ball.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/exact_synthesis/descent.hpp"
#include "quicksynth/exact_synthesis/normal_form.hpp"
#include "quicksynth/price_candidates/price.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

namespace quicksynth::test {

namespace {

json cost_triple(const Cost& c) { return json::array({c.value, c.n_t, c.n_q}); }

std::string word_of(const CyclotomicMatrix& e, int64_t t, int64_t q) {
    return Descent(t, q).word(e.matrix, e.denominator_exponent);
}

json price_of(const CyclotomicMatrix& e, int64_t t, int64_t q) {
    return cost_triple(Pricer::descent(t, q).price(e));
}

// T-count price; nullopt off Clifford+T.
std::optional<json> t_count_of(const CyclotomicMatrixDigits& e) {
    if (bloch_exponent(e) % 2 != 0) return std::nullopt;
    try {
        return cost_triple(Pricer::t_count(1).price(e));
    } catch (const std::invalid_argument&) {
        return std::nullopt;
    }
}

}  // namespace

TEST_CASE("the descent spells and prices every coset of the reference's balls") {
    // (1, 3) to radius 8, (1, 1) to 4, (2, 3) to 9, (2, 5) to 11
    for (const json& b : golden("balls")["balls"]) {
        const int64_t t = b["price"][0], q = b["price"][1];
        INFO(t, ",", q);
        const json& words = b["words"];
        for (size_t n = 0; n < words.size(); ++n) {
            const std::string w = words[n];
            INFO(w);
            const CyclotomicMatrix e = word_matrix(w);
            CHECK(word_of(e, t, q) == b["descent"][n]);
            CHECK(price_of(e, t, q) == b["cost"][n]);
            const std::string key = std::to_string(n);
            if (b["simplified"].contains(key))
                CHECK(normal_form(b["descent"][n]) == b["simplified"][key]);
        }
    }
}

TEST_CASE("the descent spells every element of the small cosets") {
    // 12 radius-6 cosets times all 24 Cliffords: every leftover Clifford occurs
    for (const json& section : golden("words")["small_cosets"]) {
        const int64_t t = section["price"][0], q = section["price"][1];
        for (size_t n = 0; n < section["words"].size(); ++n) {
            INFO(section["words"][n].get<std::string>());
            CHECK(word_of(word_matrix(section["words"][n]), t, q) == section["descent"][n]);
        }
    }
}

TEST_CASE("Clifford+T words past 128 bits in their Bloch digits") {
    // T-only ball and T-counts 40, 60, 80 (Bloch digits past 128 bits)
    const json& ct = golden("words")["clifford_t"];
    const int64_t t = ct["price"][0], q = ct["price"][1];
    for (size_t n = 0; n < ct["words"].size(); ++n) {
        const std::string w = ct["words"][n];
        INFO(w);
        const CyclotomicMatrix e = word_matrix(w);
        const std::string descent = ct["descent"][n];
        CHECK(word_of(e, t, q) == descent);
        CHECK(descent.find('Q') == std::string::npos);
        const auto [n_t, n_q] = magic_counts(descent);
        CHECK(price_of(e, t, q) == json::array({n_t + 3 * n_q, n_t, n_q}));
        CHECK(cost_triple(Pricer::t_count(1).price(e)) == json::array({n_t, n_t, 0}));
    }
}

TEST_CASE("random words at every price") {
    // includes two prices where T is spelt Q Q
    for (const json& section : golden("words")["random_words"]) {
        const int64_t t = section["price"][0], q = section["price"][1];
        for (size_t n = 0; n < section["words"].size(); ++n) {
            INFO(section["words"][n].get<std::string>(), " at ", t, ",", q);
            const CyclotomicMatrix e = word_matrix(section["words"][n]);
            CHECK(word_of(e, t, q) == section["descent"][n]);
            CHECK(price_of(e, t, q) == section["cost"][n]);
        }
    }
    for (const json& section : golden("words")["priced_words"]) {
        const int64_t t = section["price"][0], q = section["price"][1];
        for (size_t n = 0; n < section["words"].size(); ++n) {
            INFO(section["words"][n].get<std::string>(), " at ", t, ",", q);
            CHECK(price_of(word_matrix(section["words"][n]), t, q) == section["cost"][n]);
        }
    }
}

TEST_CASE("the words prices T-counts and cosets of elements past 64 bits") {
    // coefficients of 74 bits (Clifford+sqrt(T)), 71 and 125 bits (Clifford+T)
    const json& section = golden("words")["past_64_bits"];
    for (const json& r : section["rows"]) {
        const std::string w = r["word"];
        INFO(w);
        const CyclotomicMatrixDigits e = digits_of_row(r["element"]);
        for (size_t p = 0; p < section["prices"].size(); ++p) {
            const int64_t t = section["prices"][p][0], q = section["prices"][p][1];
            CHECK(Descent(t, q).word(e) == r["descent"][p]);
            CHECK(cost_triple(Pricer::descent(t, q).price(e)) == r["cost"][p]);
        }
        const std::optional<json> t_count = t_count_of(e);
        CHECK((t_count ? *t_count : json(nullptr)) == r["t_count"]);
        for (const char* clifford : {"H", "SH", "XS"})
            CHECK(coset_key(word_matrix_digits(w + clifford)) == coset_key(e));
        CHECK(coset_key(word_matrix_digits(w + "T")) != coset_key(e));
    }
}

TEST_CASE("an unreduced element past 64 bits prices and spells as the element") {
    // every fifth reference element as 2^70 M / sqrt2^(k + 140)
    const json& elements = golden("words")["elements"];
    for (size_t i = 0; i < elements.size(); i += 5) {
        INFO(elements[i]["word"].get<std::string>());
        const CyclotomicMatrix e = matrix_of_row(elements[i]["element"]);
        const CyclotomicMatrixDigits wide = widened(e, 70);
        CHECK(t_count_of(wide) == t_count_of(digits_of(e)));
        for (const auto& [t, q] : {std::pair<int64_t, int64_t>{1, 3}, {6, 3}}) {
            CHECK(cost_triple(Pricer::descent(t, q).price(wide)) == price_of(e, t, q));
            CHECK(Descent(t, q).word(wide) == word_of(e, t, q));
        }
    }
}

TEST_CASE("the descent's words have the oracle's cost and reproduce the coset") {
    auto ball = std::make_shared<CostBall>(Price{1, 3});
    ball->build(7);
    int checked = 0;
    for (int c = 0; c <= 7; ++c)
        for (int i = 0; i < std::min<int>(static_cast<int>(ball->level(c).size()), 40); ++i) {
            const CyclotomicMatrix e = word_matrix(ball->word(c, i));
            const std::string word = word_of(e, 1, 3);
            const auto [n_t, n_q] = magic_counts(word);
            CHECK(Pricer::descent(1, 3).price(e).value == c);
            CHECK(n_t + 3 * n_q == c);
            CHECK(coset_key(word_matrix(word)) == coset_key(e));
            ++checked;
        }
    CHECK(checked > 100);
}

TEST_CASE("an element that is not unitary is refused by the descent") {
    // [[a, a], [c, d]] has a Bloch image but is not unitary
    const CyclotomicMatrix good = word_matrix("HTSHQ");
    CyclotomicMatrix bad = good;
    bad.matrix[0][1] = good.matrix[0][0];
    bloch_digits(bad.matrix, bad.denominator_exponent);
    const CyclotomicMatrixDigits wide = widened(bad, 70);
    const Descent descent(1, 3);
    CHECK_THROWS_WITH_AS(descent.word(bad.matrix, bad.denominator_exponent),
                         doctest::Contains("descent"), std::invalid_argument);
    CHECK_THROWS_WITH_AS(descent.word(wide), doctest::Contains("descent"), std::invalid_argument);
    CHECK_THROWS_WITH_AS(Pricer::descent(1, 3).price(bad), doctest::Contains("descent"),
                         std::invalid_argument);
    CHECK_THROWS_WITH_AS(Pricer::descent(1, 3).price(wide), doctest::Contains("descent"),
                         std::invalid_argument);
}

}  // namespace quicksynth::test
