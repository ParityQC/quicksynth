// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// One-call policy: constants, stages and stage settings. The margin bound is a
// static_assert in synthesise.cpp.

#include "fixtures.hpp"

#include "quicksynth/policy.hpp"
#include "quicksynth/search_scales/scan.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/synthesise.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace quicksynth::test {

namespace {

std::vector<int> prefix_scales_of(double eps, const GateSet& g) {
    QuickSynth::Settings s;
    s.gateset = g;
    std::vector<int> out;
    for (const Stage& stage : stages(eps, g, QuickSynth(s).first_stage(eps, false)))
        out.push_back(stage.settings.prefix_scale.value());
    return out;
}

}  // namespace

TEST_CASE("the policy's constants") {
    CHECK(kEpsilonMargin == 1e-15);
    CHECK(kSplitBelow == 1e-5);
    CHECK(kDefaultRate == 0.3);
    CHECK(kPrefixRadius == 6);
    CHECK(kMaxPrefixes == 600);
    CHECK(kDeeperEpsilon == 1e-8);
    CHECK(kDeeperLastBudget == 1.5e10);
    CHECK(reference_throughput(GateSet::sqrt_t()) == 8.3e6);
    CHECK(reference_throughput(GateSet::clifford_t()) == 2.1e7);
}

TEST_CASE("Clifford+T splits at scale five below 1e-7") {
    // stage 1 splits at scale 4 (384 prefixes, T-count <= 8), stage 2 at 5
    using Scales = std::vector<int>;
    CHECK(prefix_scales_of(1e-8, GateSet::clifford_t()) == Scales{4, 5, 4});
    CHECK(prefix_scales_of(1e-7, GateSet::clifford_t()) == Scales{2, 3, 4, 2});
    CHECK(prefix_scales_of(1e-8, GateSet::sqrt_t()) == Scales{3, 4, 3});
}

TEST_CASE("Clifford+T splits at scale six below 1e-8") {
    // every coset of the prefix scale: 6144 of scale 6 at radius 12, 24576 of 7 at 14
    using Scales = std::vector<int>;
    CHECK(prefix_scales_of(9.9e-9, GateSet::clifford_t()) == Scales{6, 7, 8});
    CHECK(prefix_scales_of(1e-9, GateSet::sqrt_t()) == Scales{3, 4, 3});
    QuickSynth::Settings s;
    s.gateset = GateSet::clifford_t();
    const ScanThenSplit::Settings deep = QuickSynth(s).first_stage(1e-9, false);
    CHECK(!deep.max_prefixes);
    CHECK(split_scales(deep, 6) == std::vector<int>{6});
    const ScanThenSplit::Settings first = QuickSynth(s).first_stage(1e-8, false);
    CHECK(first.max_prefixes == kMaxPrefixes);
    CHECK(split_scales(first, 4) == std::vector<int>{4, 3, 2, 1});
    // later stages keep stage 1's uncapped sets without fallbacks
    const std::vector<Stage> all = stages(1e-9, s.gateset, deep);
    REQUIRE(all.size() == 3);
    for (const Stage& stage : all) CHECK(!stage.settings.fallback_scales);
    CHECK(all[1].changes == "{'prefix_scale': 7, 'split_budget': 300000000.0}");
    CHECK(all[2].settings.split_budget == kDeeperLastBudget);
    CHECK(all[2].changes ==
          "{'prefix_scale': 8, 'split_budget': 15000000000.0, 'point_budget': 1000000000.0}");
    CHECK(stages(1e-8, s.gateset, first).back().settings.split_budget == 3e9);
    s.gateset = GateSet::sqrt_t();
    const ScanThenSplit::Settings sqrt_first = QuickSynth(s).first_stage(1e-9, false);
    CHECK(sqrt_first.max_prefixes == kMaxPrefixes);
    CHECK(sqrt_first.fallback_scales);
}

TEST_CASE("the policy reaches the stage settings") {
    // fixed stopping rule by default; the model rule in stage 1 only, fixed
    // afterwards and on fallback targets
    for (const GateSet& g : {GateSet::sqrt_t(), GateSet::clifford_t()}) {
        INFO(g.name());
        QuickSynth::Settings s;
        s.gateset = g;
        CHECK(QuickSynth(s).first_stage(1e-6, false).scan.stopping_rule == StoppingRule::fixed);
        s.stopping_rule = StoppingRule::model;
        const ScanThenSplit::Settings first = QuickSynth(s).first_stage(1e-6, false);
        CHECK(!first.scan.throughput);
        CHECK(first.scan.rate == kDefaultRate);
        CHECK(first.split_below == kSplitBelow - kEpsilonMargin);
        CHECK(first.scan.stopping_rule == StoppingRule::model);
        const std::vector<Stage> all = stages(1e-6, g, first);
        REQUIRE(all.size() > 1);
        for (size_t i = 1; i < all.size(); ++i)
            CHECK(all[i].settings.scan.stopping_rule == StoppingRule::fixed);
        CHECK(QuickSynth(s).first_stage(1e-6, true).scan.stopping_rule == StoppingRule::fixed);
    }
}

TEST_CASE("fast runs LevelSearch at its own rate in every stage") {
    const auto first = [](const GateSet& g, std::optional<double> rate) {
        QuickSynth::Settings s;
        s.gateset = g;
        s.fast = true;
        s.rate = rate;
        return QuickSynth(s).first_stage(1e-6, false);
    };
    CHECK(kFineFastRate == 24.0);
    for (const GateSet& g : {GateSet::sqrt_t(), GateSet::clifford_t()}) {
        INFO(g.name());
        CHECK(first(g, std::nullopt).scan.rate == kFastRate);
        CHECK(first(g, std::nullopt).split_rate == kFastRate);
        CHECK(first(g, std::nullopt).level_search_rate == kFineFastRate);
        // an explicit rate wins
        CHECK(first(g, 2.0).scan.rate == 2.0);
        CHECK(first(g, 2.0).level_search_rate == 2.0);
    }
    // the model stopping rule walks whole scales in stage 1; the fixed later
    // stages run LevelSearch at kFineFastRate
    QuickSynth::Settings s;
    s.fast = true;
    s.stopping_rule = StoppingRule::model;
    const ScanThenSplit::Settings model_first = QuickSynth(s).first_stage(1e-6, false);
    REQUIRE(model_first.scan.stopping_rule == StoppingRule::model);
    const std::vector<Stage> all = stages(1e-6, s.gateset, model_first);
    REQUIRE(all.size() > 1);
    CHECK(all[1].settings.scan.stopping_rule == StoppingRule::fixed);
    CHECK(all[1].settings.scan.fine_levels);
    CHECK(all[1].settings.level_search_rate == kFineFastRate);
}

TEST_CASE("finer levels are on unless fine_levels is false") {
    const auto fine = [](const GateSet& g, bool fast, std::optional<bool> fine_levels) {
        QuickSynth::Settings s;
        s.gateset = g;
        s.fast = fast;
        s.fine_levels = fine_levels;
        return QuickSynth(s).first_stage(1e-6, false);
    };
    const GateSet sqrt_t = GateSet::sqrt_t();
    CHECK(kFineDefaultRate == 0.3);
    const ScanThenSplit::Settings fine_default = fine(sqrt_t, false, std::nullopt);
    CHECK(fine_default.scan.fine_levels);
    CHECK(fine_default.scan.rate == kDefaultRate);
    CHECK(fine_default.level_search_rate == kFineDefaultRate);
    CHECK(fine(sqrt_t, false, true).level_search_rate == kFineDefaultRate);
    CHECK(fine(sqrt_t, true, std::nullopt).scan.fine_levels);
    // whole scales at kDefaultRate, the default before the finer levels
    const ScanThenSplit::Settings whole_default = fine(sqrt_t, false, false);
    CHECK(!whole_default.scan.fine_levels);
    CHECK(whole_default.scan.rate == kDefaultRate);
    CHECK(!whole_default.split_rate);
    // the fast setting on whole scales, for ablation
    const ScanThenSplit::Settings whole = fine(sqrt_t, true, false);
    CHECK(!whole.scan.fine_levels);
    CHECK(whole.scan.rate == kFastRate);
    CHECK(whole.split_rate == kFastRate);
    // Clifford+T never walks finer levels
    for (const bool fast : {false, true})
        for (const std::optional<bool> f : {std::optional<bool>{}, std::optional<bool>{true}})
            CHECK(!fine(GateSet::clifford_t(), fast, f).scan.fine_levels);
    // the model stopping rule walks whole scales in stage 1; the fixed later
    // stages run LevelSearch at kFineDefaultRate
    QuickSynth::Settings s;
    s.stopping_rule = StoppingRule::model;
    const std::vector<Stage> all = stages(1e-6, s.gateset, QuickSynth(s).first_stage(1e-6, false));
    REQUIRE(all.size() > 1);
    CHECK(all[0].settings.scan.stopping_rule == StoppingRule::model);
    CHECK(all[0].settings.scan.rate == kDefaultRate);
    CHECK(all[1].settings.scan.fine_levels);
    CHECK(all[1].settings.level_search_rate == kFineDefaultRate);
}

TEST_CASE("the stages follow the requested epsilon not the margin") {
    // epsilon - epsilon_margin is below the 1e-7 thresholds
    QuickSynth::Settings s;
    const ScanThenSplit::Settings first = QuickSynth(s).first_stage(1e-7, false);
    CHECK(first.prefix_scale == 2);
    CHECK(first.split_budget == 3e7);
    CHECK(first.split_below == 1e-5 - 1e-15);
    CHECK(first.max_prefixes == 600);
    s.gateset = GateSet::clifford_t();
    const std::vector<Stage> all = stages(1e-8, s.gateset, QuickSynth(s).first_stage(1e-8, false));
    CHECK(all[1].settings.prefix_scale == 5);
    CHECK(!all[1].settings.max_prefixes);
}

TEST_CASE("each stage line prints what the stage changes") {
    QuickSynth::Settings s;
    std::vector<std::string> changes;
    for (const Stage& stage : stages(1e-6, s.gateset, QuickSynth(s).first_stage(1e-6, false)))
        changes.push_back(stage.changes);
    const std::vector<std::string> expected{
        "{}",
        "{'max_prefixes': None, 'prefix_scale': 3, 'split_budget': 300000000.0}",
        "{'max_prefixes': None, 'prefix_scale': 4, 'split_budget': 300000000.0}",
        "{'max_prefixes': None, 'split_budget': 3000000000.0, 'point_budget': 1000000000.0}",
    };
    CHECK(changes == expected);
}

TEST_CASE("a slower throughput walks fewer scales") {
    const Unitary& V = calibration_targets()[0];
    std::vector<int> deepest;
    for (double throughput : {reference_throughput(GateSet::sqrt_t()), 1.0}) {
        ScaleScan::Settings s;
        s.throughput = throughput;
        const ScaleScan::Result r = ScaleScan(1e-6, GateSet::sqrt_t(), {}, 1, s).run(V);
        int K = 0;
        for (const SearchLog& e : r.log) K = std::max(K, e.scale);
        deepest.push_back(K);
    }
    CHECK(deepest[1] < deepest[0]);
}

}  // namespace quicksynth::test
