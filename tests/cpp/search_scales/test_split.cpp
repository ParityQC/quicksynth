// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Prefix split: golden_split.json replay and comparison with a complete scan.

#include "fixtures.hpp"
#include "oracles.hpp"
#include "search_scales/replay.hpp"

#include "quicksynth/data_structures/prefix_sets.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/search_scales/scan.hpp"
#include "quicksynth/search_scales/split.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace quicksynth::test {

namespace {

// Leaves over all composite scales walked.
int64_t points_of(const PrefixSplit::Result& r) {
    int64_t points = 0;
    for (const SearchLog& e : r.log) points += e.points;
    return points;
}

// Composite scales walked by the last prefix set tried.
std::vector<int> last_set_scales(const PrefixSplit::Result& r) {
    std::vector<int> scales;
    for (auto e = r.log.rbegin(); e != r.log.rend() && e->prefix_scale == r.log.back().prefix_scale;
         ++e)
        scales.insert(scales.begin(), e->scale);
    return scales;
}

ScaleScan::Settings fixed_scales(int scales) {
    ScaleScan::Settings s;
    s.extra_scales = scales - 1;
    s.rate = std::nullopt;
    return s;
}

// Split with shallower fallbacks, <= 600 prefixes each, 3 scales from the first hit.
void replay(const std::string& ring, int64_t t, int64_t q, int prefix_scale) {
    const std::shared_ptr<const CostBall> ball = recorded_prefix_ball(ring);
    const std::vector<Recorded> runs =
        recorded_runs("split", ring, {t, q}, "prefix_scale", prefix_scale);
    CHECK(runs.size() == 24);
    for (const Recorded& record : runs) {
        const json& run = record.run;
        const json& want = record.want;
        const Unitary& V = record.V;
        INFO(run.dump());
        const PrefixSets& sets = ball->prefix_sets(prefix_scales(prefix_scale), 600);
        PrefixSplit split(number(run["epsilon"]), run_gateset(run), run_price(run), 1, sets, {},
                          {});
        const PrefixSplit::Result r = split.run(V, 3, std::nullopt);
        CHECK(best_cost(r.best) == optional_of(want["best_cost"]));
        CHECK(json(last_set_scales(r)) == want["scales"]);
        CHECK(points_of(r) == want["points"]);
        // one remainder per prefix, class and composite scale
        int64_t prefixes_tried = 0;
        for (const SearchLog& e : r.log)
            for (const auto& set : sets)
                if (set->scale == e.prefix_scale)
                    prefixes_tried +=
                        GateSet::kDeterminantClasses * static_cast<int64_t>(set->prefixes.size());
        CHECK(prefixes_tried == want["prefixes_tried"]);
        if (want["best_cost"].is_null()) {
            CHECK_FALSE(r.best.has_value());
            continue;
        }
        REQUIRE(r.best.has_value());
        check_same_circuit(run, want, r.best->element);
        // recompute the split's distance test
        CHECK(chord_distance(matrix_direction(r.best->element), su2_vector(V)) ==
              number(want["distance"]));
    }
}

}  // namespace

TEST_SUITE("golden") {
    TEST_CASE("split sqrtT t=1 q=3 prefix 1") { replay("sqrtT", 1, 3, 1); }
    TEST_CASE("split sqrtT t=1 q=3 prefix 2") { replay("sqrtT", 1, 3, 2); }
    TEST_CASE("split T t=1 q=1 prefix 2") { replay("T", 1, 1, 2); }
}

TEST_CASE("the split never beats the complete scan and walks fewer points") {
    const auto ball = recorded_prefix_ball("sqrtT");
    for (size_t i = 38; i < 40; ++i) {
        INFO("calibration target ", i);
        const Unitary& V = calibration_targets()[i];
        const ScaleScan::Result plain =
            ScaleScan(0.02, GateSet::sqrt_t(), {1, 3}, 1, fixed_scales(2)).run(V);
        PrefixSplit split(0.02, GateSet::sqrt_t(), {1, 3}, 1, ball->prefix_sets({1}, std::nullopt),
                          {}, {});
        const PrefixSplit::Result r = split.run(V, 2, std::nullopt);
        if (!plain.best || !r.best) continue;
        CHECK(r.best->cost.value >= plain.best->cost.value);  // heuristic
        int64_t walked = 0;
        for (const SearchLog& e : plain.log) walked += e.points;
        if (walked > 5000) CHECK(points_of(r) < walked);  // small regions: no saving
        const Matrix U = matrix_value(r.best->element);
        CHECK(unitarity_error(U) < 1e-9);
        CHECK(chord_diamond_distance(U, matrix_of(V)) <= 0.02);
    }
}

TEST_CASE("prefix scale zero is the ordinary search") {
    const auto ball = recorded_prefix_ball("sqrtT");
    const Unitary& V = calibration_targets()[40];
    PrefixSplit split(0.05, GateSet::sqrt_t(), {1, 3}, 1, ball->prefix_sets({0}, std::nullopt), {},
                      {});
    const PrefixSplit::Result r = split.run(V, 1, std::nullopt);
    const ScaleScan::Result plain =
        ScaleScan(0.05, GateSet::sqrt_t(), {1, 3}, 1, fixed_scales(1)).run(V);
    if (r.best && plain.best) CHECK(r.best->cost.value <= plain.best->cost.value);
}

TEST_CASE("a split budget past the leaf range is no cap") {
    const auto ball = recorded_prefix_ball("sqrtT");
    const Unitary& V = calibration_targets()[41];
    std::vector<json> results;
    for (double budget : {1e15, 1e30, std::numeric_limits<double>::infinity()}) {
        PrefixSplit::Settings settings;
        settings.point_budget = budget;
        PrefixSplit split(3e-2, GateSet::sqrt_t(), {1, 3}, 1, ball->prefix_sets({2}, std::nullopt),
                          {}, settings);
        const PrefixSplit::Result r = split.run(V, 3, std::nullopt);
        REQUIRE(r.best.has_value());
        std::vector<int> scales;
        for (const SearchLog& e : r.log) scales.push_back(e.scale);
        results.push_back({r.best->cost.value, scales, points_of(r)});
    }
    CHECK(results[0] == results[1]);
    CHECK(results[1] == results[2]);
}

TEST_CASE("a model split finds nothing cheaper than complete enumeration of its scales") {
    // 240 scale-2 prefixes, budget for L+1..L+3; complete enumeration of
    // L..L+3 is at least as cheap
    CostBall ball({1, 1});
    ball.build(6);
    int walked = 0;
    for (double eps : {3e-2, 1e-2}) {
        Enumerator enumerator(eps, GateSet::sqrt_t(), {1, 1}, 1, {});
        const int L = enumerator.first_scale();
        for (size_t i = 42; i < 46; ++i) {
            INFO(eps, " calibration target ", i);
            const Unitary& V = calibration_targets()[i];
            ScaleScan::Result scan;  // scale L walked in full
            scan.last_full_scale = L;
            std::optional<Enumerator::Hit> complete;  // cheapest over L..L+3
            for (int K = L; K < L + 4; ++K) {
                for (int det_class : {0, 1})
                    for (const Enumerator::Hit& h : enumerator.enumerate(V, K, det_class).hits)
                        if (!complete || h.cost.value < complete->cost.value) complete = h;
                if (K == L) scan.best = complete;
            }
            const PrefixSets sets = ball.prefix_sets({2}, std::nullopt);
            REQUIRE(sets[0]->prefixes.size() == 240);
            const double walks = 2.0 * 240;
            PrefixSplit::Settings settings;
            settings.point_budget = walks * enumerator.estimate_points(L + 1) * 1.001;
            PrefixSplit split(eps, GateSet::sqrt_t(), {1, 1}, 1, sets, {}, settings);
            const auto model = stopping_model(GateSet::sqrt_t(), {1, 1}, eps, L);
            const PrefixSplit::Result r = split.run_model(V, *model, scan, 0.0, false);
            for (const SearchLog& e : r.log)
                CHECK(std::set<int>{L + 1, L + 2, L + 3}.count(e.scale) == 1);
            if (!r.log.empty()) ++walked;
            std::optional<int64_t> final_cost = best_cost(scan.best);
            const std::optional<int64_t> split_cost = best_cost(r.best);
            if (split_cost && (!final_cost || *split_cost < *final_cost)) final_cost = split_cost;
            if (!final_cost) continue;
            REQUIRE(complete.has_value());
            CHECK(complete->cost.value <= *final_cost);
        }
    }
    CHECK(walked > 0);
}

}  // namespace quicksynth::test
