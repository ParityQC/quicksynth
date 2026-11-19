// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Scan-then-split search: golden_combined.json replay, result fields, defaults
// and log.

#include "fixtures.hpp"
#include "oracles.hpp"
#include "search_scales/replay.hpp"

#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/search_scales/scan.hpp"
#include "quicksynth/search_scales/scan_then_split.hpp"
#include "quicksynth/search_scales/split.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace quicksynth::test {

namespace {

// Rate rule at the recorded throughput. Split when uncertified, or below the
// default threshold for "split": "default".
void replay(const std::string& ring, int64_t t, int64_t q, const std::string& split, int count) {
    // "always" runs have no split field
    const std::vector<Recorded> runs =
        recorded_runs("combined", ring, {t, q}, "split", split == "always" ? json() : json(split));
    CHECK(runs.size() == count);
    for (const Recorded& record : runs) {
        const json& run = record.run;
        const json& want = record.want;
        INFO(run.dump());
        ScanThenSplit::Settings s;
        s.scan.throughput = recorded_throughput(ring);
        if (split == "always") s.split_below = std::numeric_limits<double>::infinity();
        ScanThenSplit search(number(run["epsilon"]), run_gateset(run), run_price(run), 1,
                             recorded_prefix_ball(ring), s, nullptr);
        const ScanThenSplit::Result r = search.run(record.V);
        const std::optional<CyclotomicMatrix> element =
            r.best ? std::optional<CyclotomicMatrix>(r.best->element) : std::nullopt;
        CHECK(best_cost(r.best) == optional_of(want["best_cost"]));
        const json source =
            r.best ? json(r.summary.source == Source::scan ? "scan" : "split") : json(nullptr);
        CHECK(source == want["source"]);
        CHECK(r.summary.scan_cost == optional_of(want["scan_cost"]));
        CHECK(r.summary.scan_certified == want["certified"]);
        if (!want["best_cost"].is_null()) check_same_circuit(run, want, element);
    }
}

}  // namespace

TEST_SUITE("golden") {
    TEST_CASE("combined sqrtT t=1 q=3 throughput") { replay("sqrtT", 1, 3, "always", 64); }
    TEST_CASE("combined T t=1 q=1 throughput") { replay("T", 1, 1, "always", 64); }
    TEST_CASE("combined sqrtT t=1 q=3 throughput split default") {
        replay("sqrtT", 1, 3, "default", 64);
    }
    TEST_CASE("combined T t=1 q=1 throughput split default") { replay("T", 1, 1, "default", 64); }
}

TEST_CASE("QuickSynth keeps the cheaper of the two searches") {
    for (const GateSet& gateset : {GateSet::sqrt_t(), GateSet::clifford_t()}) {
        INFO(gateset.name());
        const bool sqrt_t = gateset.step == 1;
        const Price price = sqrt_t ? Price{1, 3} : Price{1, 1};
        ScanThenSplit search(1e-2, gateset, price, 1, recorded_prefix_ball(sqrt_t ? "sqrtT" : "T"),
                             {}, nullptr);
        for (size_t i = 46; i < 49; ++i) {
            INFO("calibration target ", i);
            const Unitary& V = calibration_targets()[i];
            const ScanThenSplit::Result r = search.run(V);
            REQUIRE(r.best.has_value());
            if (r.summary.source == Source::split)
                CHECK((!r.summary.scan_cost || r.best->cost.value < *r.summary.scan_cost));
            else
                CHECK(r.best->cost.value == r.summary.scan_cost);
            const Matrix U = matrix_value(r.best->element);
            CHECK(unitarity_error(U) < 1e-9);
            CHECK(chord_diamond_distance(U, matrix_of(V)) <= 1e-2);
            if (sqrt_t)
                CHECK(Pricer::descent(1, 3).price(r.best->element).value == r.best->cost.value);
            else
                CHECK(Pricer::t_count(1).price(r.best->element).n_t == r.best->cost.value);
        }
    }
}

TEST_CASE("the split's rate rule at 0 walks the fixed scales and at 1e12 keeps the scan's answer") {
    const Unitary& V = calibration_targets()[46];
    const auto run = [&](std::optional<double> rate) {
        ScanThenSplit::Settings s;
        s.split_below = std::numeric_limits<double>::infinity();
        s.split_rate = rate;
        ScanThenSplit search(1e-3, GateSet::sqrt_t(), {1, 3}, 1, recorded_prefix_ball("sqrtT"), s,
                             nullptr);
        return search.run(V);
    };
    const ScanThenSplit::Result none = run(std::nullopt);
    const ScanThenSplit::Result zero = run(0.0);
    const ScanThenSplit::Result huge = run(1e12);
    REQUIRE(none.best.has_value());
    CHECK(best_cost(zero.best) == best_cost(none.best));
    CHECK(zero.log.size() == none.log.size());
    // the scan answered, so no composite scale is worth its time
    REQUIRE(none.summary.scan_cost.has_value());
    const auto split_walks = [](const ScanThenSplit::Result& r) {
        return std::count_if(r.log.begin(), r.log.end(),
                             [](const SearchLog& e) { return e.prefix_scale.has_value(); });
    };
    CHECK(split_walks(none) > 0);
    CHECK(split_walks(huge) == 0);
    CHECK(best_cost(huge.best) == none.summary.scan_cost);
}

TEST_CASE("on Clifford+T a certified scan leaves out the split") {
    // certified by the T-count bound
    ScanThenSplit search(1e-2, GateSet::clifford_t(), {1, 1}, 1, recorded_prefix_ball("T"), {},
                         nullptr);
    int certified = 0;
    for (size_t i = 49; i < 53; ++i) {
        const ScanThenSplit::Result r = search.run(calibration_targets()[i]);
        if (!r.summary.scan_certified) continue;
        ++certified;
        CHECK_FALSE(r.summary.split_ran);
        REQUIRE(r.best.has_value());
        CHECK(r.summary.source == Source::scan);
    }
    CHECK(certified > 0);
}

TEST_CASE("without a ball the scan runs alone") {
    ScanThenSplit search(1e-2, GateSet::sqrt_t(), {1, 3}, 1, nullptr, {}, nullptr);
    const ScanThenSplit::Result r = search.run(calibration_targets()[53]);
    CHECK_FALSE(r.summary.split_ran);
    REQUIRE(r.best.has_value());
    CHECK(r.summary.source == Source::scan);
    CHECK(r.best->cost.value == r.summary.scan_cost);
}

TEST_CASE("the result carries what the scan walked and whether the split ran") {
    // at fixed throughput a standalone scan walks the same scales
    for (const double eps : {1e-3, 1e-6}) {
        INFO("epsilon ", eps);
        ScanThenSplit::Settings settings;
        settings.scan.throughput = 7.9e6;
        ScanThenSplit search(eps, GateSet::sqrt_t(), {1, 3}, 1, recorded_prefix_ball("sqrtT"),
                             settings, nullptr);
        ScaleScan alone(eps, GateSet::sqrt_t(), {1, 3}, 1, settings.scan);
        for (size_t i = 0; i < 4; ++i) {
            INFO("calibration target ", i);
            const Unitary& V = calibration_targets()[i];
            const ScanThenSplit::Result r = search.run(V);
            const ScaleScan::Result s = alone.run(V);
            std::vector<int> scales;
            for (const SearchLog& e : s.log) scales.push_back(e.scale);
            std::sort(scales.begin(), scales.end());
            scales.erase(std::unique(scales.begin(), scales.end()), scales.end());
            CHECK(r.summary.scales == scales);
            CHECK(r.summary.first_hit_scale == s.first_hit_scale);
            CHECK(r.summary.last_full_scale == s.last_full_scale);
            CHECK(r.summary.lower_bound == s.lower_bound);
            CHECK(r.summary.estimated_points == s.estimated_points);
            CHECK(r.summary.split_ran == (!s.certified && (!s.best || eps < settings.split_below)));
        }
    }
}

TEST_CASE("the defaults are the values the fixtures were recorded at") {
    // former Python wrapper defaults
    const ScaleScan::Settings scan;
    CHECK(scan.rate == 0.3);
    CHECK(scan.extra_scales == 2);
    CHECK(scan.point_budget == 3e6);
    CHECK(kScaleSaving == 0.06);
    CHECK(PrefixSplit::Settings{}.point_budget == 3e7);
    const ScanThenSplit::Settings both;
    CHECK(both.max_prefixes == 600);
    CHECK(both.split_below == 1e-5);
    CHECK(Enumerator::Settings{}.transform_cache == 4096);
}

TEST_CASE("a split entry of the log counts the scan's answer") {
    // a split entry's best cost is the search's, never above the scan's
    const auto ball = recorded_prefix_ball("sqrtT");
    int checked = 0;
    for (StoppingRule stopping_rule : {StoppingRule::fixed, StoppingRule::model}) {
        ScanThenSplit::Settings s;
        s.split_below = std::numeric_limits<double>::infinity();
        s.scan.stopping_rule = stopping_rule;
        ScanThenSplit search(1e-3, GateSet::sqrt_t(), {1, 3}, 1, ball, s, nullptr);
        for (size_t i = 54; i < 57; ++i) {
            INFO("calibration target ", i);
            const ScanThenSplit::Result r = search.run(calibration_targets()[i]);
            for (const SearchLog& e : r.log)
                if (e.prefix_scale && r.summary.scan_cost) {
                    REQUIRE(e.best_cost.has_value());
                    CHECK(*e.best_cost <= *r.summary.scan_cost);
                    ++checked;
                }
            REQUIRE_FALSE(r.log.empty());
            CHECK(r.log.back().best_cost == best_cost(r.best));
        }
    }
    CHECK(checked > 0);
}

TEST_CASE("with finer levels on Clifford+sqrt(T) the split starts past the scan's last level") {
    ScanThenSplit::Settings s;
    s.split_below = std::numeric_limits<double>::infinity();
    s.scan.fine_levels = true;
    s.split_levels = 2;
    s.level_search = false;
    ScanThenSplit search(1e-3, GateSet::sqrt_t(), {1, 3}, 1, recorded_prefix_ball("sqrtT"), s,
                         nullptr);
    const int L = levels_per_scale(GateSet::sqrt_t());
    int checked = 0;
    for (size_t i = 0; i < 4; ++i) {
        INFO("calibration target ", i);
        const ScanThenSplit::Result r = search.run(calibration_targets()[i]);
        REQUIRE(r.summary.scan_cost.has_value());
        int scan_last = 0;
        for (const SearchLog& e : r.log) {
            const int level = L * e.scale - e.sub;
            if (!e.prefix_scale) {
                scan_last = std::max(scan_last, level);
                continue;
            }
            CHECK(level > scan_last);
            ++checked;
        }
    }
    CHECK(checked > 0);
}

}  // namespace quicksynth::test
