// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// One-call synthesis result: the hit in normal form at the hit's cost, priced
// from the word's magic counts (no sqrt(T) on Clifford+T). The distance is
// within kDistanceError of the long-double oracle and below
// epsilon - kDistanceError. Prefix balls: prefix radius against
// golden_words.json, the Clifford+T prefix ball below 1e-7, the model
// split's prefix sets, and prepare(). Stage 0 against golden_rotation.json:
// the Clifford+T optimum of the corrected Ross-Selinger search, Clifford+sqrt(T)
// against the reference costs and its own recorded costs, the diagonal runs of
// golden_synthesise.json against the joint walk, and the benchmark ZZ angles;
// near-diagonal targets up to the rotation tolerance.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/synthesise.hpp"
#include "quicksynth/policy.hpp"
#include "quicksynth/synthesise_rotation/rotation_search.hpp"
#include "quicksynth/search_scales/scan_then_split.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/distance.hpp"
#include "quicksynth/data_structures/prefix_sets.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/support/host.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace quicksynth::test {

namespace {

// Whether the prefix ball of `radius` levels is cached (Clifford+sqrt(T) with
// use_q). A build polls the interrupt check before its first level and an
// interrupted build is not cached, so the probe never adds a ball.
bool ball_is_cached(const Price& price, int radius, bool use_q = true) {
    set_interrupt_check([] { throw std::runtime_error("ball not cached"); });
    bool cached = true;
    try {
        prefix_ball(price, use_q, radius);
    } catch (const std::runtime_error&) {
        cached = false;
    }
    set_interrupt_check(nullptr);
    return cached;
}

// Rz(theta) = diag(exp(-i theta/2), exp(i theta/2)).
Unitary z_rotation(double theta) {
    const double c = std::cos(theta / 2), s = std::sin(theta / 2);
    return {c, -s, 0, 0, 0, 0, c, s};
}

// Rows of a {"columns", "rows"} fixture table, keyed by column.
std::vector<json> table_rows(const json& table) {
    std::vector<json> out;
    for (const json& row : table["rows"]) {
        json r;
        for (size_t c = 0; c < table["columns"].size(); ++c)
            r[table["columns"][c].get<std::string>()] = row[c];
        out.push_back(std::move(r));
    }
    return out;
}

QuickSynth synthesiser(const GateSet& gateset, bool fast = false, bool rotation_path = true) {
    QuickSynth::Settings settings;
    settings.gateset = gateset;
    settings.fast = fast;
    settings.rotation_path = rotation_path;
    return QuickSynth(settings);
}

// The answer's cost re-priced from its word at the default price, and its
// distance to U from the long-double oracle.
void check_answer(const Synthesised& a, const Unitary& U, double eps, const GateSet& gateset) {
    const CyclotomicMatrix e = word_matrix(a.word);
    const Cost check = Pricer::for_gateset(gateset, Price{}).price(e);
    CHECK(check.value == a.cost.value);
    CHECK(check.n_t == a.cost.n_t);
    CHECK(check.n_q == a.cost.n_q);
    const auto [n_t, n_q] = magic_counts(a.word);
    CHECK(n_t == a.cost.n_t);
    CHECK(n_q == a.cost.n_q);
    if (gateset.step != 1) CHECK(n_q == 0);
    const double exact = chord_diamond_distance(matrix_value(e), matrix_of(U));
    CHECK(std::abs(static_cast<long double>(a.distance) - exact) <= kDistanceError);
    CHECK(exact <= eps);
    CHECK(a.distance <= eps - kDistanceError);
}

// A verified stage-0 answer, certified exactly on Clifford+T.
void check_stage_0(const Synthesised& a, const Unitary& U, double eps, const GateSet& gateset) {
    CHECK(a.stages == 0);
    CHECK(a.certified == (gateset.step != 1));
    REQUIRE(a.search.has_value());
    CHECK(a.search->source == Source::rotation);
    check_answer(a, U, eps, gateset);
}

// A Clifford+sqrt(T) stage-0 answer at its recorded cost, or a pure-T word of
// the Clifford+T search below it (source clifford_t), which the recording
// predates.
void check_recorded_sqrt_t(const Synthesised& a, const Unitary& U, double eps, int64_t recorded) {
    REQUIRE(a.search.has_value());
    if (a.search->source != Source::clifford_t) {
        check_stage_0(a, U, eps, GateSet::sqrt_t());
        CHECK(a.cost.value == recorded);
        return;
    }
    CHECK(a.stages == 0);
    check_answer(a, U, eps, GateSet::sqrt_t());
    CHECK(a.cost.value < recorded);
}

// Golden Clifford+sqrt(T) rows of ids below `ids` at `epsilons` (all if empty):
// default and fast at their recorded costs, or below them by a pure-T word,
// the default within two of the reference on each angle and below it in the
// mean of each epsilon, fast within half a gate of it.
void check_sqrt_t_rotations(int ids, const std::vector<double>& epsilons) {
    const QuickSynth synth = synthesiser(GateSet::sqrt_t());
    const QuickSynth fast = synthesiser(GateSet::sqrt_t(), true);
    std::map<double, std::array<int64_t, 4>> sums;  // n, default, fast, reference
    for (const json& row : table_rows(golden("rotation")["sqrt_t"])) {
        const double eps = row["epsilon"];
        if (row["id"].get<int>() >= ids) continue;
        if (!epsilons.empty() && std::find(epsilons.begin(), epsilons.end(), eps) == epsilons.end())
            continue;
        INFO("id ", row["id"].get<int>(), " eps ", eps);
        const Unitary U = z_rotation(row["theta"]);
        const Synthesised a = synth.synthesise(U, eps);
        const Synthesised f = fast.synthesise(U, eps);
        check_recorded_sqrt_t(a, U, eps, row["cost"]);
        check_recorded_sqrt_t(f, U, eps, row["fast_cost"]);
        const int64_t reference = row["reference"];
        CHECK(a.cost.value <= reference + 2);
        auto& s = sums[eps];
        s[0] += 1;
        s[1] += a.cost.value;
        s[2] += f.cost.value;
        s[3] += reference;
    }
    CHECK(sums.size() == (epsilons.empty() ? 6 : epsilons.size()));
    // no structured bindings: INFO's lambda cannot capture them in C++17 (Clang)
    for (const auto& entry : sums) {
        const double eps = entry.first;
        const auto& s = entry.second;
        INFO("eps ", eps, " rows ", s[0]);
        CHECK(s[1] <= s[3]);
        CHECK(2 * s[2] <= 2 * s[3] + s[0]);
    }
}

// Near-exact rows of golden_rotation.json, whose bodies hold slabs of points:
// Clifford+T at the reference T-count, Clifford+sqrt(T) default and fast at
// their recorded costs.
void check_near_exact(bool sqrt_t) {
    const QuickSynth synth = synthesiser(sqrt_t ? GateSet::sqrt_t() : GateSet::clifford_t());
    const QuickSynth fast = synthesiser(GateSet::sqrt_t(), true);
    int64_t rows = 0, certified = 0;
    for (const json& row : table_rows(golden("rotation")["near_exact"])) {
        const double eps = row["epsilon"];
        INFO("theta ", row["theta"].get<double>(), " eps ", eps);
        const Unitary U = z_rotation(row["theta"]);
        const Synthesised a = synth.synthesise(U, eps);
        CHECK(a.stages == 0);
        check_answer(a, U, eps, sqrt_t ? GateSet::sqrt_t() : GateSet::clifford_t());
        ++rows;
        if (!sqrt_t) {
            CHECK(a.cost.value == row["t_count"].get<int64_t>());
            certified += a.certified;
            continue;
        }
        CHECK_FALSE(a.certified);
        CHECK(a.cost.value == row["cost"].get<int64_t>());
        const Synthesised f = fast.synthesise(U, eps);
        check_stage_0(f, U, eps, GateSet::sqrt_t());
        CHECK(f.cost.value == row["fast_cost"].get<int64_t>());
    }
    CHECK(rows == 24);
    // a walk cut in a slab leaves the optimum uncertified
    if (!sqrt_t) CHECK(certified == 14);
}

}  // namespace

TEST_CASE("an answer is its hit spelt at the hit's cost" * doctest::test_suite("slow")) {
    for (const GateSet& gateset : {GateSet::sqrt_t(), GateSet::clifford_t()}) {
        QuickSynth::Settings settings;
        settings.gateset = gateset;
        const QuickSynth synth(settings);
        // search price (T-count on Clifford+T)
        const Price price =
            gateset.step == 1 ? settings.price : Price::clifford_t(settings.price.t);
        for (const double eps : {1e-3, 1e-5}) {
            for (int i = 0; i < 8; ++i) {
                const Unitary& U = calibration_targets()[i];
                const Synthesised r = synth.synthesise(U, eps);
                INFO(gateset.name(), " ", eps, " ", i, " ", r.word);
                const auto [n_t, n_q] = magic_counts(r.word);
                CHECK(r.cost.n_t == n_t);
                CHECK(r.cost.n_q == n_q);
                CHECK(r.cost.value == price.cost(n_t, n_q).value);
                if (gateset.step != 1) {
                    CHECK(n_q == 0);
                    CHECK(r.word.find_first_of("Qq") == std::string::npos);
                }
                const double exact =
                    chord_diamond_distance(matrix_value(word_matrix(r.word)), matrix_of(U));
                CHECK(std::abs(static_cast<long double>(r.distance) - exact) <= kDistanceError);
                CHECK(r.distance <= eps - kDistanceError);
            }
        }
    }
}

TEST_CASE("the split's prefix radius is the reference's at each price") {
    for (const json& r : golden("words")["prefix_radius"]) {
        INFO(r.dump());
        const Price price{r["price"][0].get<int64_t>(), r["price"][1].get<int64_t>()};
        CHECK(prefix_radius(price, r["use_q"], r["scale"]) == r["radius"]);
    }
}

TEST_CASE("the Clifford+T prefix ball holds every prefix of scale 4 by T-count 8") {
    const auto ball = prefix_ball({1, 3}, false, 10);
    std::vector<int64_t> found;
    for (int c = 0; c <= 10; ++c) {
        int64_t n = 0;
        for (size_t i = 0; i < ball->level(c).size(); ++i)
            n += word_matrix(ball->word(c, static_cast<int>(i))).denominator_exponent == 4;
        found.push_back(n);
    }
    int64_t total = 0;
    for (int64_t n : found) total += n;
    CHECK(total == 384);
    CHECK(found[9] == 0);
    CHECK(found[10] == 0);
    CHECK(split_radius({1, 3}, GateSet::clifford_t(), 1e-8, 4) == 8);
    CHECK(split_radius({1, 3}, GateSet::clifford_t(), 1e-7, 4) == kPrefixRadius);
    CHECK(split_radius({1, 3}, GateSet::clifford_t(), 1e-9, 6) == 12);
    CHECK(prefix_ball({1, 3}, false, 12)->scale_count(6) == 6144);
    CHECK_THROWS_AS(split_radius({int64_t{1} << 30, 3}, GateSet::clifford_t(), 1e-9, 6),
                    std::invalid_argument);
}

TEST_CASE("every set the model split is given has a measured coverage") {
    // stage-1 sets (fixed ball, then stopping ball of stopping_extra_scales,
    // <= 600 each) must have a measured split_coverage, not the default
    struct Case {
        GateSet gateset;
        Price price;
    };
    std::map<int, std::array<double, 3>> fallback;
    for (int k : {2, 3, 4}) fallback[k] = split_coverage(GateSet::sqrt_t(), {5, 7}, k, 1);
    for (const Case& c : {Case{GateSet::sqrt_t(), {1, 1}}, Case{GateSet::sqrt_t(), {1, 3}},
                          Case{GateSet::sqrt_t(), {2, 3}}, Case{GateSet::clifford_t(), {1, 3}}}) {
        const bool use_q = c.gateset.step == 1;
        for (double eps : {1e-4, 1e-8}) {
            INFO(c.gateset.name(), " ", c.price.t, ",", c.price.q, " ", eps);
            const int scale = prefix_scale_for(eps, c.gateset);
            const int radius = split_radius(c.price, c.gateset, eps, scale);
            PrefixSets sets =
                prefix_ball(c.price, use_q, radius)->prefix_sets(prefix_scales(scale), 600);
            const std::vector<int> more = stopping_extra_scales(scale);
            if (!more.empty()) {
                const int r = prefix_radius(c.price, use_q, more.back());
                const PrefixSets extra = prefix_ball(c.price, use_q, r)->prefix_sets(more, 600);
                sets.insert(sets.end(), extra.begin(), extra.end());
            }
            std::set<int> seen;
            for (const auto& s : sets) {
                const int k = s->scale;
                const int64_t size = static_cast<int64_t>(s->prefixes.size());
                if (k < 2 || size == 0 || seen.count(k)) continue;
                seen.insert(k);
                CHECK(split_coverage(c.gateset, c.price, k, size) != fallback[k]);
            }
            CHECK(seen.count(2) == 1);
            CHECK(seen.count(3) == 1);
        }
    }
}

TEST_CASE("prepare builds the prefix balls of a call") {
    // No other test uses price (1, 4), so its balls start uncached. At 1e-7 the
    // split takes radius kPrefixRadius and the model stopping rule's scale-3 prefixes radius 9.
    QuickSynth::Settings settings;
    settings.price = {1, 4};
    settings.stopping_rule = StoppingRule::fixed;
    REQUIRE(!ball_is_cached(settings.price, kPrefixRadius));
    QuickSynth(settings).prepare(1e-7);
    CHECK(ball_is_cached(settings.price, kPrefixRadius));
    CHECK(!ball_is_cached(settings.price, 9));
    settings.stopping_rule = StoppingRule::model;
    const QuickSynth synth(settings);
    synth.prepare(1e-7);
    CHECK(ball_is_cached(settings.price, 9));
    // a call after prepare builds no larger ball
    CHECK(synth.synthesise(calibration_targets()[0], 1e-7).distance <= 1e-7);
    CHECK(!ball_is_cached(settings.price, 10));
}

TEST_CASE("stage 0 answers Clifford+T z-rotations at the joint walk's certified cost") {
    // Price t = 5 is used nowhere else, so its Clifford+T balls start uncached.
    QuickSynth::Settings settings;
    settings.gateset = GateSet::clifford_t();
    settings.price = {5, 3};
    const QuickSynth synth(settings);
    settings.rotation_path = false;
    const QuickSynth joint(settings);
    REQUIRE(!ball_is_cached(settings.price, kPrefixRadius, false));
    std::vector<std::pair<Unitary, Synthesised>> answers;
    for (const double theta : {0.7, 2.5, -1.1}) {
        const Unitary U{std::cos(theta / 2), -std::sin(theta / 2), 0, 0, 0, 0,
                        std::cos(theta / 2), std::sin(theta / 2)};
        for (const double eps : {1e-3, 1e-7}) {
            INFO("theta ", theta, " eps ", eps);
            const Synthesised r = synth.synthesise(U, eps);
            CHECK(r.stages == 0);
            CHECK(r.certified);
            CHECK(r.search->source == Source::rotation);
            const auto [n_t, n_q] = magic_counts(r.word);
            CHECK(r.cost.value == 5 * n_t);
            CHECK(r.cost.n_t == n_t);
            CHECK(n_q == 0);
            const double exact =
                chord_diamond_distance(matrix_value(word_matrix(r.word)), matrix_of(U));
            CHECK(std::abs(static_cast<long double>(r.distance) - exact) <= kDistanceError);
            CHECK(r.distance <= eps - kDistanceError);
            if (eps == 1e-3) answers.emplace_back(U, r);
        }
    }
    // stage 0 builds no prefix ball
    CHECK(!ball_is_cached(settings.price, kPrefixRadius, false));
    for (const auto& [U, r] : answers) {
        const Synthesised w = joint.synthesise(U, 1e-3);
        CHECK(w.stages >= 1);
        CHECK(w.certified);
        CHECK(w.cost.value == r.cost.value);
    }
}

TEST_CASE("stage 0 answers Clifford+sqrt(T) z-rotations uncertified by the count rule") {
    QuickSynth::Settings settings;
    settings.gateset = GateSet::sqrt_t();
    settings.price = {1, 3};
    CHECK(QuickSynth(settings).rotation_solvable() == kRotationSolvable);
    settings.fast = true;
    CHECK(QuickSynth(settings).rotation_solvable() == kRotationSolvableFast);
    settings.rotation_solvable = 7;
    CHECK(QuickSynth(settings).rotation_solvable() == 7);
    settings.rotation_solvable = std::nullopt;
    const QuickSynth fast(settings);
    settings.fast = false;
    const QuickSynth synth(settings);
    const Pricer pricer = Pricer::for_gateset(settings.gateset, settings.price);
    for (const double theta : {0.7, 2.5, -1.1}) {
        const Unitary U{std::cos(theta / 2), -std::sin(theta / 2), 0, 0, 0, 0,
                        std::cos(theta / 2), std::sin(theta / 2)};
        for (const double eps : {1e-3, 1e-7}) {
            INFO("theta ", theta, " eps ", eps);
            const Synthesised r = synth.synthesise(U, eps);
            const Synthesised f = fast.synthesise(U, eps);
            for (const Synthesised* a : {&r, &f}) {
                CHECK(a->stages == 0);
                CHECK_FALSE(a->certified);
                CHECK(a->search->source == Source::rotation);
                const Cost check = pricer.price(word_matrix(a->word));
                CHECK(check.value == a->cost.value);
                CHECK(check.n_t == a->cost.n_t);
                CHECK(check.n_q == a->cost.n_q);
                const double exact =
                    chord_diamond_distance(matrix_value(word_matrix(a->word)), matrix_of(U));
                CHECK(std::abs(static_cast<long double>(a->distance) - exact) <= kDistanceError);
                CHECK(a->distance <= eps - kDistanceError);
            }
            int64_t solvable = 0, fast_solvable = 0;
            for (const auto& [stage, e] : r.log) solvable += e.rotation->solvable;
            for (const auto& [stage, e] : f.log) fast_solvable += e.rotation->solvable;
            CHECK(solvable >= kRotationSolvable);
            CHECK(fast_solvable >= kRotationSolvableFast);
            CHECK(fast_solvable < solvable);
            CHECK(r.cost.value <= f.cost.value);
        }
    }
}

TEST_CASE("a Clifford+sqrt(T) stage 0 word stopped by a cap bounds the stages' answer" *
          doctest::test_suite("slow")) {
    // without the count rule stage 0 stops at its hit cap at 1e-5 (by the level
    // bound at 1e-4 and above); the stages then run as without stage 0, and
    // its word answers only if cheaper
    QuickSynth::Settings settings;
    settings.gateset = GateSet::sqrt_t();
    settings.rotation_solvable = INT32_MAX;
    const QuickSynth capped(settings);
    const QuickSynth joint = synthesiser(GateSet::sqrt_t(), false, false);
    const double eps = 1e-5;
    const Unitary U = z_rotation(0.7);
    const Synthesised a = capped.synthesise(U, eps);
    const Synthesised w = joint.synthesise(U, eps);
    check_answer(a, U, eps, GateSet::sqrt_t());
    CHECK_FALSE(a.certified);
    std::set<int> stages;
    for (const auto& entry : a.log) stages.insert(entry.first);
    CHECK(stages.count(0) == 1);
    CHECK(stages.count(w.stages) == 1);
    if (a.stages == 0) {
        CHECK(a.cost.value < w.cost.value);
        CHECK(a.search->source == Source::rotation);
    } else {
        CHECK(a.word == w.word);
        CHECK(a.stages == w.stages);
    }
}

TEST_CASE("stage 0 answers every golden Clifford+T z-rotation certified at the reference T-count") {
    // the corrected Ross-Selinger optimum (golden_rotation.json)
    const QuickSynth synth = synthesiser(GateSet::clifford_t());
    int64_t rows = 0;
    for (const json& row : table_rows(golden("rotation")["clifford_t"])) {
        const double eps = row["epsilon"];
        INFO("id ", row["id"].get<int>(), " eps ", eps);
        const Unitary U = z_rotation(row["theta"]);
        const Synthesised a = synth.synthesise(U, eps);
        check_stage_0(a, U, eps, GateSet::clifford_t());
        CHECK(a.cost.value == row["t_count"].get<int64_t>());
        ++rows;
    }
    CHECK(rows == 720);
}

TEST_CASE("stage 0 answers near-exact Clifford+T z-rotations at the reference T-count") {
    check_near_exact(false);
}

TEST_CASE("stage 0 answers near-exact Clifford+sqrt(T) z-rotations at their recorded costs" *
          doctest::test_suite("slow")) {
    check_near_exact(true);
}

TEST_CASE("stage 0 answers golden Clifford+sqrt(T) z-rotations within two of the reference") {
    check_sqrt_t_rotations(20, {1e-3, 1e-6, 1e-10});
}

TEST_CASE("stage 0 answers every golden Clifford+sqrt(T) z-rotation within two of the reference" *
          doctest::test_suite("slow")) {
    check_sqrt_t_rotations(50, {});
}

TEST_CASE("stage 0 answers near-diagonal targets up to the rotation tolerance") {
    // C1 Rz C2 with an off-axis tilt of f epsilon between Rz and C2: within
    // the tolerance (1e-2) stage 0 answers uncertified, beyond it the call
    // falls through to the stages unchanged
    const std::array<std::pair<Matrix, Matrix>, 2> cliffords{std::pair{mul(kH, kS), kH},
                                                             std::pair{kS, mul(mul(kH, kS), kH)}};
    for (const GateSet& gateset : {GateSet::clifford_t(), GateSet::sqrt_t()}) {
        const QuickSynth synth = synthesiser(gateset);
        const QuickSynth joint = synthesiser(gateset, false, false);
        for (const auto& [c1, c2] : cliffords) {
            for (const double eps : {1e-3, 1e-4}) {
                for (const double f : {1e-3, 1e-1}) {
                    INFO(gateset.name(), " eps ", eps, " f ", f);
                    const Matrix tilt = rotation(2 * std::asin(f * eps), 1, 1, 0);
                    const Unitary U = unitary_of(mul(mul(mul(c1, rotation(0.7, 'z')), tilt), c2));
                    const Synthesised a = synth.synthesise(U, eps);
                    check_answer(a, U, eps, gateset);
                    if (f < 1e-2) {
                        CHECK(a.stages == 0);
                        CHECK_FALSE(a.certified);
                        continue;
                    }
                    const Synthesised w = joint.synthesise(U, eps);
                    CHECK(a.stages >= 1);
                    CHECK(a.word == w.word);
                    CHECK(a.cost.value == w.cost.value);
                    CHECK(a.stages == w.stages);
                    CHECK(a.certified == w.certified);
                }
            }
        }
    }
}

TEST_CASE("stage 0 costs no more than the joint walk on the diagonal golden runs") {
    // the runs of golden_synthesise.json that stage 0 answers, at their
    // recorded rotation cost; on Clifford+T the joint walk's cost where it
    // certified it
    const json& runs = golden("synthesise")["runs"];
    const std::vector<json> diagonal = table_rows(golden("rotation")["diagonal"]);
    CHECK(diagonal.size() == 35);
    for (const json& d : diagonal) {
        const json& run = runs[d["run"].get<size_t>()];
        INFO("run ", d["run"].get<int>());
        for (const auto& [key, value] : run["config"].items())
            CHECK((key == "gateset" || key == "method" || key == "workers" ||
                   key == "deterministic"));
        const GateSet gateset =
            GateSet::from_name(run["config"].value("gateset", std::string("sqrt_t")));
        const double eps = run["epsilon"];
        const Unitary U = decode_target(run["target"]);
        const Synthesised a = synthesiser(gateset).synthesise(U, eps);
        check_stage_0(a, U, eps, gateset);
        CHECK(a.cost.value == d["cost"].get<int64_t>());
        const json& recorded = run["result"];
        CHECK(a.cost.value <= recorded["cost"].get<int64_t>());
        if (gateset.step != 1 && recorded["certified"].get<bool>())
            CHECK(a.cost.value == recorded["cost"].get<int64_t>());
    }
}

TEST_CASE("stage 0 answers every ZZ rotation of the benchmark workloads" *
          doctest::test_suite("slow")) {
    // trotter, qaoa and adiabatic (golden_rotation.json): Clifford+T certified,
    // Clifford+sqrt(T) default and fast, each at its recorded cost
    const QuickSynth clifford_t = synthesiser(GateSet::clifford_t());
    const QuickSynth sqrt_t = synthesiser(GateSet::sqrt_t());
    const QuickSynth fast = synthesiser(GateSet::sqrt_t(), true);
    std::map<std::string, int> rows;
    for (const json& row : table_rows(golden("rotation")["zz"])) {
        const double eps = row["epsilon"];
        INFO(row["workload"].get<std::string>(), " ", row["label"].get<std::string>());
        const Unitary U = z_rotation(row["theta"]);
        const Synthesised t = clifford_t.synthesise(U, eps);
        const Synthesised a = sqrt_t.synthesise(U, eps);
        const Synthesised f = fast.synthesise(U, eps);
        check_stage_0(t, U, eps, GateSet::clifford_t());
        check_stage_0(a, U, eps, GateSet::sqrt_t());
        check_stage_0(f, U, eps, GateSet::sqrt_t());
        CHECK(t.cost.value == row["t_count"].get<int64_t>());
        CHECK(a.cost.value == row["cost"].get<int64_t>());
        CHECK(f.cost.value == row["fast_cost"].get<int64_t>());
        ++rows[row["workload"]];
    }
    CHECK(rows == std::map<std::string, int>{{"adiabatic", 80}, {"qaoa", 200}, {"trotter", 112}});
}

}  // namespace quicksynth::test
