// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The set-up cache: a batch leaves it as its keys taken one at a time would.

#include "fixtures.hpp"

#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/enumerate_lattice/setup.hpp"

#include <doctest/doctest.h>

#include <functional>
#include <vector>

namespace quicksynth::test {

namespace {

using IntMap = ClearingMap<int, int, std::hash<int>>;

int ten_times(int key) { return 10 * key; }

}  // namespace

TEST_CASE("a batch leaves the clearing map as inserting one key at a time would") {
    IntMap map(2);
    {
        auto batch = map.batch();
        batch.insert(1);
        batch.insert(2);
        batch.commit(ten_times);
    }
    // One at a time: {1, 2}, + 3 -> {1, 2, 3}, 4 finds 3 > 2 entries and
    // clears -> {4}, + 1 + 5 -> {4, 1, 5}.
    auto batch = map.batch();
    REQUIRE(batch.stored(1) != nullptr);
    CHECK(*batch.stored(1) == 10);
    batch.insert(3);
    batch.insert(4);
    CHECK(batch.stored(1) == nullptr);
    CHECK_FALSE(batch.inserted(3));
    CHECK(batch.inserted(4));
    batch.insert(1);
    batch.insert(5);
    batch.commit(ten_times);
    const auto after = map.batch();
    for (int key : {4, 1, 5}) {
        REQUIRE(after.stored(key) != nullptr);
        CHECK(*after.stored(key) == 10 * key);
    }
    for (int key : {2, 3}) CHECK(after.stored(key) == nullptr);
}

TEST_CASE("a batch never committed leaves the clearing map untouched") {
    IntMap map(1);
    {
        auto batch = map.batch();
        batch.insert(1);
        batch.insert(2);
        batch.commit(ten_times);
    }
    {
        auto batch = map.batch();
        batch.insert(3);  // clears
        batch.insert(4);
    }
    const auto after = map.batch();
    CHECK(after.stored(1) != nullptr);
    CHECK(after.stored(2) != nullptr);
    CHECK(after.stored(3) == nullptr);
}

TEST_CASE("a set-up lost to a clear and met again in the same batch is cached again") {
    // Capacity 1: the third problem clears the first two set-ups; the second
    // problem, met again, is inserted again, as taken one at a time.
    const SearchRegion region(1e-3, GateSet::sqrt_t());
    const std::vector<Unitary>& targets = calibration_targets();
    REQUIRE(targets.size() >= 3);
    const int K = region.first_scale();
    SetupCache setups(1);
    const auto first = setups.for_problems(
        region, 1, {{targets[0], 0}, {targets[1], 0}, {targets[2], 0}, {targets[1], 0}}, K, 0);
    CHECK(first[3].get() == first[1].get());
    const auto again = setups.for_problems(region, 1, {{targets[1], 0}, {targets[2], 0}}, K, 0);
    CHECK(again[0].get() == first[1].get());
    CHECK(again[1].get() == first[2].get());
}

}  // namespace quicksynth::test
