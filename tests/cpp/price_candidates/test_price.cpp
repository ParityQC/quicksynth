// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Prices of cyclotomic matrices: normal-form weight on Clifford+sqrt(T),
// t * T-count on Clifford+T; int64 and arbitrary-size digits.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"

#include <doctest/doctest.h>

#include <stdexcept>
#include <string>

namespace quicksynth::test {

namespace {

bool same(const Cost& a, const Cost& b) {
    return a.value == b.value && a.n_t == b.n_t && a.n_q == b.n_q;
}

}  // namespace

TEST_CASE("each gate set chooses its pricer") {
    const Pricer sqrt_t = Pricer::for_gateset(GateSet::sqrt_t(), Price{2, 3});
    const Pricer clifford_t = Pricer::for_gateset(GateSet::clifford_t(), Price{2, 3});
    const Price p = Price::clifford_t(5);
    CHECK((p.t == 5 && p.q == 3));
    for (const json& e : golden("words")["elements"]) {
        const std::string word = e["word"];
        INFO(word);
        const CyclotomicMatrix element = word_matrix(word);
        CHECK(same(sqrt_t.price(element), Pricer::descent(2, 3).price(element)));
        CHECK(same(sqrt_t.price(digits_of(element)), sqrt_t.price(element)));
        if (word.find_first_of("Qq") != std::string::npos) continue;
        const auto [n_t, n_q] = magic_counts(word);
        const Cost c = clifford_t.price(element);
        CHECK(c.value == 2 * c.n_t);
        CHECK(c.n_t <= n_t);
        CHECK(c.n_q == 0);
        CHECK(same(clifford_t.price(digits_of(element)), c));
    }
}

TEST_CASE("the T-count refuses an element off Clifford+T") {
    // H Q: its Bloch exponent is odd
    const CyclotomicMatrix e = word_matrix("HQ");
    CHECK_THROWS_WITH_AS(Pricer::t_count(1).price(e), doctest::Contains("not Clifford+T"),
                         std::invalid_argument);
    CHECK_THROWS_WITH_AS(Pricer::t_count(1).price(digits_of(e)),
                         doctest::Contains("not Clifford+T"), std::invalid_argument);
}

}  // namespace quicksynth::test
