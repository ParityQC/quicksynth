// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Scale scan: golden_scan.json replay, scale bound, rate rule, power price.

#include "fixtures.hpp"
#include "search_scales/replay.hpp"

#include "quicksynth/search_scales/scan.hpp"

#include <doctest/doctest.h>

#include <limits>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace quicksynth::test {

namespace {

// Settings of a recorded variant; rate rule at the recorded throughput or off.
ScaleScan::Settings variant_settings(const json& run) {
    ScaleScan::Settings s;
    const std::string variant = run["variant"];
    if (variant == "throughput") s.throughput = recorded_throughput(run["ring"]);
    // certify: no rate rule, no scale cap
    if (variant == "certify") s.extra_scales = std::numeric_limits<int>::max();
    if (variant == "scales2") s.extra_scales = 1;
    if (variant != "throughput") s.rate = std::nullopt;
    return s;
}

// (K, determinant class, points, hit costs) per walk.
using Scale = std::tuple<int, int, int64_t, std::map<int64_t, int64_t>>;

std::vector<Scale> scales_of(const std::vector<Enumerator::Result>& walks) {
    std::vector<Scale> out;
    for (const Enumerator::Result& s : walks) {
        std::map<int64_t, int64_t> histogram;
        for (const Enumerator::Hit& h : s.hits) ++histogram[h.cost.value];
        out.emplace_back(s.K, s.det_class, s.points, histogram);
    }
    return out;
}

std::vector<Scale> scales_of(const json& want) {
    std::vector<Scale> out;
    for (const json& s : want["scales"]) {
        std::map<int64_t, int64_t> histogram;
        for (const auto& [cost, count] : s[3].items()) histogram[std::stoll(cost)] = count;
        out.emplace_back(s[0], s[1], s[2], histogram);
    }
    return out;
}

std::optional<CyclotomicMatrix> best_matrix(const ScaleScan::Result& r) {
    if (!r.best) return std::nullopt;
    return r.best->element;
}

void check_scan(const json& run, const json& want, const ScaleScan::Result& got,
                const std::vector<Enumerator::Result>& walks) {
    CHECK(scales_of(walks) == scales_of(want));
    CHECK(best_cost(got.best) == optional_of(want["best_cost"]));
    CHECK(got.first_hit_scale == optional_of(want["first_hit_scale"]));
    CHECK(got.certified == want["certified"]);
    CHECK(got.lower_bound == number(want["lower_bound"]));
    CHECK(got.last_full_scale == optional_of(want["last_full_scale"]));
    CHECK(got.estimated_points == number(want["estimated_points"]));
    if (want["best_cost"].is_null()) {
        CHECK_FALSE(got.best.has_value());
        return;
    }
    // same scales, same first cheapest hit; ties not accepted
    check_same_circuit(run, want, best_matrix(got));
}

// Runs of one gate set, price and variant: 12 targets 1e-2..1e-6, 4 at 1e-7.
void replay(const std::string& ring, int64_t t, int64_t q, const std::string& variant) {
    const std::vector<Recorded> runs = recorded_runs("scan", ring, {t, q}, "variant", variant);
    CHECK(runs.size() == 64);
    for (const Recorded& r : runs) {
        INFO(r.run.dump());
        ScaleScan scan(number(r.run["epsilon"]), run_gateset(r.run), run_price(r.run), 1,
                       variant_settings(r.run));
        std::vector<Enumerator::Result> walks;
        const ScaleScan::Result got =
            scan.run(r.V, [&](const Enumerator::Result& w, double) { walks.push_back(w); });
        check_scan(r.run, r.want, got, walks);
    }
}

}  // namespace

TEST_SUITE("golden") {
    TEST_CASE("scan sqrtT t=1 q=1 throughput") { replay("sqrtT", 1, 1, "throughput"); }
    TEST_CASE("scan sqrtT t=1 q=1 certify") { replay("sqrtT", 1, 1, "certify"); }
    TEST_CASE("scan sqrtT t=1 q=3 throughput") { replay("sqrtT", 1, 3, "throughput"); }
    TEST_CASE("scan sqrtT t=1 q=3 certify") { replay("sqrtT", 1, 3, "certify"); }
    TEST_CASE("scan sqrtT t=1 q=3 scales2") { replay("sqrtT", 1, 3, "scales2"); }
    TEST_CASE("scan sqrtT t=1 q=3 rate_none") { replay("sqrtT", 1, 3, "rate_none"); }
    TEST_CASE("scan sqrtT t=2 q=3 throughput") { replay("sqrtT", 2, 3, "throughput"); }
    TEST_CASE("scan sqrtT t=2 q=3 certify") { replay("sqrtT", 2, 3, "certify"); }
    TEST_CASE("scan T t=1 q=1 throughput") { replay("T", 1, 1, "throughput"); }
    TEST_CASE("scan T t=1 q=1 certify") { replay("T", 1, 1, "certify"); }
    TEST_CASE("scan T t=1 q=1 scales2") { replay("T", 1, 1, "scales2"); }
    TEST_CASE("scan T t=1 q=1 rate_none") { replay("T", 1, 1, "rate_none"); }
}

TEST_CASE("the scale bound follows the price") {
    for (const auto& [price, bound] :
         std::vector<std::pair<Price, double>>{{{2, 3}, 33}, {{1, 3}, 16.5}, {{1, 1}, 11}}) {
        ScaleScan::Settings s;
        s.rate = std::nullopt;
        CHECK(ScaleScan(1e-3, GateSet::sqrt_t(), price, 1, s).scale_lower_bound(10) == bound);
    }
}

TEST_CASE("a larger rate stops sooner and never answers cheaper") {
    // rate in magic states/s; never stops before a hit
    for (size_t i = 34; i < 37; ++i) {
        INFO("calibration target ", i);
        const Unitary& V = calibration_targets()[i];
        ScaleScan::Settings patient, hasty;
        patient.rate = 0.01;
        hasty.rate = 1e6;
        const ScaleScan::Result p = ScaleScan(0.02, GateSet::sqrt_t(), {}, 1, patient).run(V);
        const ScaleScan::Result h = ScaleScan(0.02, GateSet::sqrt_t(), {}, 1, hasty).run(V);
        int64_t walked_p = 0, walked_h = 0;
        for (const SearchLog& e : p.log) walked_p += e.points;
        for (const SearchLog& e : h.log) walked_h += e.points;
        CHECK(walked_h <= walked_p);
        REQUIRE(h.best.has_value());
        if (p.best) CHECK(p.best->cost.value <= h.best->cost.value);
    }
}

TEST_CASE("an incumbent stops the scan where the scale bound reaches it") {
    const Unitary& V = calibration_targets()[0];
    const Price price = Price::clifford_t(1);
    const ScaleScan::Result free = ScaleScan(1e-3, GateSet::clifford_t(), price, 1, {}).run(V);
    REQUIRE(free.best.has_value());
    REQUIRE(!free.log.empty());
    const int first = free.log.front().scale;
    const int64_t best = free.best->cost.value;
    for (const int64_t incumbent : {int64_t{0}, int64_t{2 * first - 3}, best, best + 100}) {
        INFO("incumbent ", incumbent);
        ScaleScan::Settings s;
        s.incumbent = incumbent;
        ScaleScan scan(1e-3, GateSet::clifford_t(), price, 1, s);
        const ScaleScan::Result r = scan.run(V);
        for (const SearchLog& e : r.log) CHECK(scan.scale_lower_bound(e.scale) < incumbent);
        if (incumbent <= scan.scale_lower_bound(first)) CHECK(r.log.empty());
        if (incumbent > best) CHECK(r.log.size() == free.log.size());
    }
}

TEST_CASE("the certified scan matches the certified optima at the power price" *
          doctest::test_suite("slow")) {
    // CertSynth optima at (2, 3), 1e-3
    const std::vector<int64_t> certified{28, 28, 27, 26, 27};
    ScaleScan::Settings s;
    s.extra_scales = std::numeric_limits<int>::max();
    s.rate = std::nullopt;
    s.point_budget = 2e8;
    ScaleScan scan(1e-3, GateSet::sqrt_t(), {2, 3}, 4, s);
    for (size_t i = 0; i < certified.size(); ++i) {
        INFO("calibration target ", i);
        const ScaleScan::Result r = scan.run(calibration_targets()[i]);
        CHECK(r.certified);
        REQUIRE(r.best.has_value());
        CHECK(r.best->cost.value == certified[i]);
    }
}

}  // namespace quicksynth::test
