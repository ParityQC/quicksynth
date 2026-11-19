// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// LevelSearch through ScanThenSplit: answers, stopping, the split without a scan
// answer, and Clifford+T left on the fixed depth.

#include "fixtures.hpp"
#include "oracles.hpp"
#include "search_scales/replay.hpp"

#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/search_scales/level_search.hpp"
#include "quicksynth/search_scales/scan_then_split.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace quicksynth::test {

namespace {

ScanThenSplit::Settings level_settings(double rate) {
    ScanThenSplit::Settings s;
    s.scan.fine_levels = true;
    s.scan.rate = rate;
    s.split_below = std::numeric_limits<double>::infinity();
    return s;
}

ScanThenSplit::Result search_at(double eps, const ScanThenSplit::Settings& s, const Unitary& V,
                                const GateSet& gateset = GateSet::sqrt_t()) {
    const bool sqrt_t = gateset.step == 1;
    ScanThenSplit search(eps, gateset, sqrt_t ? Price{1, 3} : Price{1, 1}, 1,
                         recorded_prefix_ball(sqrt_t ? "sqrtT" : "T"), s, nullptr);
    return search.run(V);
}

// Best cost after the walk that first found an element.
std::optional<int64_t> first_answer(const ScanThenSplit::Result& r) {
    for (const SearchLog& e : r.log)
        if (e.best_cost) return e.best_cost;
    return std::nullopt;
}

int level_of(const SearchLog& e) { return levels_per_scale(GateSet::sqrt_t()) * e.scale - e.sub; }

}  // namespace

TEST_CASE("the level search answers within epsilon and never above its first answer") {
    int split_walks = 0;
    for (const double eps : {1e-3, 1e-6}) {
        for (size_t i = 0; i < 4; ++i) {
            INFO("epsilon ", eps, " calibration target ", i);
            const Unitary& V = calibration_targets()[i];
            const ScanThenSplit::Result r = search_at(eps, level_settings(6.0), V);
            REQUIRE(r.best.has_value());
            const Matrix U = matrix_value(r.best->element);
            CHECK(chord_diamond_distance(U, matrix_of(V)) <= eps);
            CHECK(Pricer::descent(1, 3).price(r.best->element).value == r.best->cost.value);
            REQUIRE(first_answer(r).has_value());
            CHECK(r.best->cost.value <= *first_answer(r));
            // best costs only fall, and the last is the answer
            std::optional<int64_t> before;
            for (const SearchLog& e : r.log) {
                if (before) CHECK((e.best_cost && *e.best_cost <= *before));
                before = e.best_cost;
            }
            CHECK(r.log.back().best_cost == best_cost(r.best));
            if (r.summary.source == Source::split)
                CHECK((!r.summary.scan_cost || r.best->cost.value < *r.summary.scan_cost));
            else
                CHECK(r.best->cost.value == r.summary.scan_cost);
            // each search walks its levels in increasing order
            int scan_last = -1, split_last = -1;
            for (const SearchLog& e : r.log) {
                int& last = e.prefix_scale ? split_last : scan_last;
                if (!e.prefix_scale && level_of(e) == last) continue;  // second class
                CHECK(level_of(e) > last);
                last = level_of(e);
                split_walks += e.prefix_scale ? 1 : 0;
            }
        }
    }
    CHECK(split_walks > 0);
}

TEST_CASE(
    "the level search stops sooner at a higher rate and at its first answer at an "
    "unbounded one") {
    const Unitary& V = calibration_targets()[5];
    const ScanThenSplit::Result slow = search_at(1e-4, level_settings(0.6), V);
    const ScanThenSplit::Result fast = search_at(1e-4, level_settings(6.0), V);
    const ScanThenSplit::Result none = search_at(1e-4, level_settings(1e12), V);
    REQUIRE(none.best.has_value());
    CHECK(none.log.size() <= fast.log.size());
    CHECK(fast.log.size() <= slow.log.size());
    CHECK(best_cost(slow.best) <= best_cost(fast.best));
    CHECK(best_cost(fast.best) <= best_cost(none.best));
    // no walk past the level of the first answer
    CHECK(none.summary.source == Source::scan);
    const auto hit = std::find_if(none.log.begin(), none.log.end(),
                                  [](const SearchLog& e) { return e.best_cost.has_value(); });
    REQUIRE(hit != none.log.end());
    CHECK(level_of(none.log.back()) == level_of(*hit));
    CHECK(none.log.back().best_cost == first_answer(none));
}

TEST_CASE("without a scan answer the level search's split answers past the scan's levels") {
    ScanThenSplit::Settings s = level_settings(6.0);
    // the scan stops short of any element
    s.scan.point_budget = 2.0;
    int answered = 0;
    for (size_t i = 0; i < 3; ++i) {
        INFO("calibration target ", i);
        const Unitary& V = calibration_targets()[i];
        const ScanThenSplit::Result r = search_at(1e-6, s, V);
        CHECK_FALSE(r.summary.scan_cost.has_value());
        if (!r.best) continue;
        ++answered;
        CHECK(r.summary.source == Source::split);
        CHECK(chord_diamond_distance(matrix_value(r.best->element), matrix_of(V)) <= 1e-6);
        int scan_last = -1;
        for (const SearchLog& e : r.log)
            if (!e.prefix_scale) scan_last = std::max(scan_last, level_of(e));
        for (const SearchLog& e : r.log)
            if (e.prefix_scale) CHECK(level_of(e) > scan_last);
    }
    CHECK(answered > 0);
}

TEST_CASE("before the first answer the split leaves out levels its set-ups would dominate") {
    ScanThenSplit::Settings walk_all = level_settings(6.0);
    walk_all.scan.point_budget = 2.0;
    walk_all.skip_setup_levels = false;
    ScanThenSplit::Settings skip = walk_all;
    skip.skip_setup_levels = true;
    const GateSet g = GateSet::sqrt_t();
    const SearchRegion region(1e-6, g);
    const auto split_log = [](const ScanThenSplit::Result& r) {
        std::vector<SearchLog> out;
        for (const SearchLog& e : r.log)
            if (e.prefix_scale) out.push_back(e);
        return out;
    };
    for (size_t i = 0; i < 3; ++i) {
        INFO("calibration target ", i);
        const Unitary& V = calibration_targets()[i];
        const std::vector<SearchLog> all = split_log(search_at(1e-6, walk_all, V));
        const std::vector<SearchLog> some = split_log(search_at(1e-6, skip, V));
        REQUIRE_FALSE(all.empty());
        REQUIRE_FALSE(some.empty());
        CHECK(level_of(some.front()) > level_of(all.front()));
        // the first walk builds every set-up afresh, so its walks take longer
        const SearchLog& first = some.front();
        CHECK(walk_seconds(region.estimate_points(first.scale - *first.prefix_scale, first.sub)) >=
              kFreshSetupSeconds);
        // with the next level over budget, the last one within it is walked
        ScanThenSplit::Settings tight = skip;
        tight.split_budget = 2 * all.front().estimated;
        const ScanThenSplit::Result last = search_at(1e-6, tight, V);
        const std::vector<SearchLog> one = split_log(last);
        REQUIRE(one.size() == 1);
        CHECK(level_of(one.front()) == level_of(all.front()));
    }
}

TEST_CASE("Clifford+T keeps the fixed depth under level_search") {
    ScanThenSplit::Settings fixed = level_settings(3.0);
    fixed.level_search = false;
    for (size_t i = 0; i < 3; ++i) {
        INFO("calibration target ", i);
        const Unitary& V = calibration_targets()[i];
        const ScanThenSplit::Result a =
            search_at(1e-4, level_settings(3.0), V, GateSet::clifford_t());
        const ScanThenSplit::Result b = search_at(1e-4, fixed, V, GateSet::clifford_t());
        CHECK(best_cost(a.best) == best_cost(b.best));
        REQUIRE(a.log.size() == b.log.size());
        for (size_t j = 0; j < a.log.size(); ++j) {
            CHECK(a.log[j].scale == b.log[j].scale);
            CHECK(a.log[j].sub == b.log[j].sub);
            CHECK(a.log[j].points == b.log[j].points);
        }
    }
}

}  // namespace quicksynth::test
