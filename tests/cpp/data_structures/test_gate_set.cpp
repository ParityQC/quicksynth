// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Gate set rings inside Z[zeta16] and gate set names.

#include "quicksynth/data_structures/gate_set.hpp"

#include <doctest/doctest.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace quicksynth::test {

TEST_CASE("the generator determines the gate set") {
    const GateSet sqrt_t = GateSet::sqrt_t(), clifford_t = GateSet::clifford_t();
    CHECK(sqrt_t.places == std::vector<int>{1, 3, 5, 7});
    CHECK(clifford_t.places == std::vector<int>{1, 3});
    CHECK(sqrt_t.step == 1);
    CHECK(clifford_t.step == 2);
    CHECK(GateSet::kDeterminantClasses == 2);
    for (const GateSet& g : {sqrt_t, clifford_t}) {
        CHECK(2 * g.num_places() == g.degree());
        CHECK(g.dim() == 2 * g.degree());
        CHECK(GateSet::from_name(g.name()).places == g.places);
        CHECK(GateSet::from_name(g.name()).step == g.step);
    }
    CHECK(std::string(sqrt_t.name()) == "sqrt_t");
    CHECK(std::string(clifford_t.name()) == "clifford_t");
}

TEST_CASE("an unknown gate set is refused") {
    CHECK_THROWS_AS(GateSet::from_name("sqrtT"), std::invalid_argument);
    CHECK_THROWS_AS(GateSet::from_name(""), std::invalid_argument);
}

}  // namespace quicksynth::test
