// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Prefix-split prefix sets against golden_prefixes.json.

#include "fixtures.hpp"
#include "search_scales/replay.hpp"

#include "quicksynth/data_structures/cost_ball.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/prefix_sets.hpp"
#include "quicksynth/price_candidates/price.hpp"

#include <doctest/doctest.h>

#include <optional>
#include <string>
#include <vector>

namespace quicksynth::test {

namespace {

// The set with each prefix's cost (its ball level).
json as_json(const Pricer& pricer, const PrefixSet& s) {
    json costs = json::array(), elements = json::array();
    for (const CyclotomicMatrix& p : s.prefixes) {
        costs.push_back(pricer.price(p).value);
        elements.push_back(row_of(p));
    }
    return {{"scale", s.scale}, {"costs", costs}, {"elements", elements}};
}

}  // namespace

TEST_CASE("the prefix sets are the reference's capped and uncapped") {
    // radius-6 balls, with and without fallback scales
    for (const json& p : golden("prefixes")["prefixes"]) {
        CHECK(p["radius"] == 6);  // matches recorded_prefix_ball
        const auto ball = recorded_prefix_ball(p["use_q"].get<bool>() ? "sqrtT" : "T");
        const std::vector<int> scales = p["scales"];
        const std::optional<int64_t> cap = optional_of(p["cap"]);
        const std::vector<PrefixSet> sets = collect_prefixes(*ball, scales, cap);
        const PrefixSets& kept = ball->prefix_sets(scales, cap);
        const Pricer pricer =
            Pricer::for_gateset(p["use_q"] ? GateSet::sqrt_t() : GateSet::clifford_t(), {1, 3});
        REQUIRE(sets.size() == p["sets"].size());
        REQUIRE(kept.size() == sets.size());
        for (size_t i = 0; i < sets.size(); ++i) {
            json want = p["sets"][i];
            want.erase("words");
            CHECK(as_json(pricer, sets[i]) == want);
            CHECK(as_json(pricer, *kept[i]) == want);
        }
    }
}

TEST_CASE("a prefix set holds the elements of its denominator scale") {
    const auto ball = recorded_prefix_ball("sqrtT");
    for (const auto& [scale, count] : std::vector<std::pair<int, size_t>>{{0, 4}, {1, 24}}) {
        const std::vector<PrefixSet> sets = collect_prefixes(*ball, {scale}, std::nullopt);
        REQUIRE(sets.size() == 1);
        CHECK(sets[0].prefixes.size() == count);
        for (const CyclotomicMatrix& p : sets[0].prefixes) CHECK(p.denominator_exponent == scale);
    }
}

TEST_CASE("denominator zero means monomial and not Clifford") {
    // T itself has no root of two
    const auto ball = recorded_prefix_ball("sqrtT");
    const std::vector<PrefixSet> sets = collect_prefixes(*ball, {0}, std::nullopt);
    REQUIRE(sets.size() == 1);
    for (const CyclotomicMatrix& p : sets[0].prefixes) {
        int zeros = 0;
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b) {
                bool zero = true;
                for (int64_t c : p.matrix[a][b]) zero = zero && c == 0;
                zeros += zero;
            }
        CHECK(p.denominator_exponent == 0);
        CHECK(zeros == 2);
    }
}

TEST_CASE("the scales of a split are its own and then every shallower one") {
    CHECK(prefix_scales(3) == std::vector<int>{3, 2, 1});
    CHECK(prefix_scales(1) == std::vector<int>{1});
}

}  // namespace quicksynth::test
