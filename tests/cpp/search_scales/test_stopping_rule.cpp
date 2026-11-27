// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Model stopping rule against golden_window.json (tests/fixtures/make_window.py).

#include "fixtures.hpp"
#include "oracles.hpp"
#include "search_scales/replay.hpp"

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/search_scales/scale_bound.hpp"
#include "quicksynth/search_scales/scan.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace quicksynth::test {

namespace {

GateSet gateset_of(int step) { return step == 1 ? GateSet::sqrt_t() : GateSet::clifford_t(); }

// Shared model for a setting at the scan's first scale.
std::shared_ptr<const StoppingModel> model_of(int step, int64_t t, int64_t q, double eps) {
    const GateSet g = gateset_of(step);
    return stopping_model(g, {t, q}, eps, SearchRegion(eps, g).first_scale());
}

// Seconds of scale start + i; infinite if absent.
std::function<double(int)> timer(const StoppingModel& m,
                                 const std::vector<std::optional<double>>& s) {
    return [start = m.start(), s](int g) {
        const int i = g - start;
        if (i < 0 || i >= static_cast<int>(s.size()) || !s[static_cast<size_t>(i)])
            return std::numeric_limits<double>::infinity();
        return *s[static_cast<size_t>(i)];
    };
}

std::vector<std::optional<double>> seconds_of(const json& list) {
    std::vector<std::optional<double>> out;
    for (const json& x : list) out.push_back(x.is_null() ? std::nullopt : std::optional(number(x)));
    return out;
}

std::array<double, 3> phi_of(const json& phi) {
    return {number(phi[0]), number(phi[1]), number(phi[2])};
}

ScaleScan scan_of(const GateSet& g, const Price& price, StoppingRule stopping_rule,
                  std::optional<double> throughput = std::nullopt) {
    ScaleScan::Settings s;
    s.stopping_rule = stopping_rule;
    s.throughput = throughput;
    return ScaleScan(1e-3, g, price, 1, s);
}

}  // namespace

// --- census and bound -----------------------------------------------------------------

TEST_CASE("the census equals the reference on the overlap") {
    for (const json& ref : golden("window")["census"]) {
        INFO(ref["step"].get<int>(), " ", ref["t"].get<int>(), ",", ref["q"].get<int>());
        const json& N = ref["N"];
        const Census census(gateset_of(ref["step"]), {ref["t"], ref["q"]},
                            static_cast<int>(N.size()) - 1);
        const size_t width = std::min<size_t>(N[0].size(), census.cost_cap() + 1);
        for (size_t K = 0; K < N.size(); ++K)
            for (size_t c = 0; c < width; ++c) {
                const double exact = number(N[K][c]);
                const double ours = census.count(static_cast<int>(K), static_cast<int64_t>(c));
                if (exact < 0x1p53)
                    CHECK(ours == exact);
                else
                    CHECK(close(ours, exact, 1e-12));
            }
    }
}

TEST_CASE("the integer bound is the census minimum") {
    for (const json& ref : golden("window")["bounds"]) {
        const int step = ref["step"];
        const Price price{ref["t"], ref["q"]};
        const GateSet g = gateset_of(step);
        for (size_t i = 0; i < ref["K"].size(); ++i)
            CHECK(integer_scale_bound(g, price, ref["K"][i]) == ref["c_min"][i]);
        // core census: least cost at scales >= K
        const int kmax = step == 1 ? 36 : 48;
        const Census census(g, price, kmax);
        std::vector<int64_t> least;
        for (int K = 0; K <= kmax; ++K) {
            int64_t c = 0;
            while (census.count(K, c) == 0) ++c;
            least.push_back(c);
        }
        for (int K = 0; K < kmax - 6; ++K)
            CHECK(integer_scale_bound(g, price, K) ==
                  *std::min_element(least.begin() + K, least.end()));
    }
}

TEST_CASE("the integer bound is never weaker than the fixed one") {
    for (const Price& price : {Price{1, 1}, Price{1, 3}, Price{2, 3}, Price{2, 5}, Price{1, 2}}) {
        const ScaleScan fixed = scan_of(GateSet::sqrt_t(), price, StoppingRule::fixed);
        const ScaleScan model = scan_of(GateSet::sqrt_t(), price, StoppingRule::model);
        for (int K = 2; K < 30; ++K) {
            CHECK(model.scale_lower_bound(K) >= fixed.scale_lower_bound(K));
            CHECK(model.scale_lower_bound(K) ==
                  static_cast<double>(integer_scale_bound(GateSet::sqrt_t(), price, K)));
        }
    }
}

TEST_CASE("a T price above two q counts as two q") {
    // QQ = T, so (3, 1) acts as (2, 1)
    const Census a(GateSet::sqrt_t(), {3, 1}, 14), b(GateSet::sqrt_t(), {2, 1}, 14);
    REQUIRE(a.cost_cap() == b.cost_cap());
    for (int K = 0; K <= 14; ++K)
        for (int64_t c = 0; c <= a.cost_cap(); ++c) CHECK(a.count(K, c) == b.count(K, c));
    for (int K = 0; K < 20; ++K)
        CHECK(integer_scale_bound(GateSet::sqrt_t(), {3, 1}, K) ==
              integer_scale_bound(GateSet::sqrt_t(), {2, 1}, K));
}

// --- model tests ----------------------------------------------------------------------

TEST_CASE("the model equals the reference on a grid") {
    for (const json& s : golden("window")["settings"]) {
        INFO(s["step"].get<int>(), " ", s["t"].get<int>(), ",", s["q"].get<int>(), " ",
             number(s["epsilon"]));
        const auto m = model_of(s["step"], s["t"], s["q"], number(s["epsilon"]));
        CHECK(m->start() == s["start"]);
        CHECK(close(m->mu(), number(s["mu"]), 1e-14));
        const auto sec = timer(*m, seconds_of(s["seconds"]));
        const double rate = number(s["rate"]);
        for (const json& g : s["grid"]) {
            INFO(g.dump());
            const int K = g["K"];
            const int64_t b = g["b"];
            CHECK(close(m->saving(K, b), number(g["saving"])));
            // no per-scale overhead
            CHECK_FALSE(m->dp_needed(K, sec));
            const double one_step = number(g["one_step"]), want = number(g["dp"]);
            CHECK(m->walk(K, b, rate, sec) == (one_step > 0));
            const double dp = m->dp_value(K, b, rate, sec);
            if (std::isinf(want)) {
                CHECK(dp == want);
            } else {
                CHECK(close(dp, want, 1e-9, 1e-12));
                CHECK((dp > 0) == (one_step > 0));  // DP stops where the test does
            }
        }
        for (const json& sp : s["split"])
            CHECK(close(m->split_saving(sp["Kc"], phi_of(sp["phi"]), sp["b"], sp["L"]),
                        number(sp["saving"])));
    }
}

TEST_CASE("the dynamic programme replaces the one-step test under a large overhead") {
    const auto m = model_of(1, 1, 3, 1e-4);
    const Enumerator enumerator(1e-4, GateSet::sqrt_t(), {1, 3}, 1, {});
    std::vector<std::optional<double>> overhead;  // o = 0.2 s per scale
    for (int i = 0; i < 12; ++i) {
        const double base = 2 * enumerator.estimate_points(m->start() + i) / 7.9e6;
        overhead.push_back(base > 0.5 ? std::nullopt : std::optional(0.2 + base));
    }
    const auto sec = timer(*m, overhead);
    int checked = 0;
    for (int K = m->start(); K < m->start() + 3; ++K) {
        if (!overhead[static_cast<size_t>(K + 2 - m->start())]) continue;
        CHECK(m->dp_needed(K, sec));
        for (int64_t b = 15; b < 40; ++b) {
            CHECK(m->walk(K, b, 0.3, sec) == (m->dp_value(K, b, 0.3, sec) > 0));
            ++checked;
        }
    }
    CHECK(checked > 0);
}

TEST_CASE("a split of complete coverage is valued as complete enumeration") {
    const auto m = model_of(1, 1, 3, 1e-4);
    for (int L = m->start(); L < m->start() + 4; ++L) {
        const int Kc = L + 1;
        for (int64_t b = 0; b < 60; b += 5) {
            INFO(L, " ", b);
            // phi = (1, 0, 0): next group only
            CHECK(close(m->split_saving(Kc, {1.0, 0.0, 0.0}, b, L), m->saving(L, b)));
            CHECK(close(m->split_p_cheaper(Kc, {1.0, 0.0, 0.0}, b, L),
                        -std::expm1(-m->intensity_sums(Kc - m->start())[static_cast<size_t>(b)])));
            // phi = 1 two composites on: both groups in full
            const std::vector<double>& one = m->intensity_sums(Kc - m->start());
            const std::vector<double>& two = m->intensity_sums(Kc + 1 - m->start());
            long double want = 0;
            for (int64_t c = 1; c <= b; ++c)
                want += -std::expm1(-(static_cast<long double>(one[static_cast<size_t>(c)]) +
                                      two[static_cast<size_t>(c)]));
            CHECK(close(m->split_saving(Kc + 1, {1.0, 1.0, 1.0}, b, L), static_cast<double>(want)));
        }
    }
}

TEST_CASE("the split's coverage keeps the plus-two offset") {
    // K(P) + K(R) = K(PR) + 2, after a T syllable with l(P) = l(R) = 1 mod 4
    const std::string P = "SHQSHQSHQSHT", R = "SHQSHTSHTSHTH";
    CHECK(word_matrix(P).denominator_exponent == 4);
    CHECK(word_matrix(R).denominator_exponent == 4);
    CHECK(word_matrix(P + R).denominator_exponent == 6);
    for (const Price& price : {Price{1, 1}, Price{1, 3}, Price{2, 3}, Price{1, 2}, Price{3, 6}}) {
        const std::array<double, 3> phi = split_coverage(GateSet::sqrt_t(), price, 3, 600);
        CHECK(phi[2] > 0);
        CHECK(phi[0] + phi[1] + phi[2] <= 1);
    }
    CHECK(split_coverage(GateSet::sqrt_t(), {3, 6}, 3, 600) ==
          split_coverage(GateSet::sqrt_t(), {1, 2}, 3, 600));
    CHECK(split_coverage(GateSet::sqrt_t(), {2, 6}, 2, 128) ==
          split_coverage(GateSet::sqrt_t(), {1, 3}, 2, 128));
}

TEST_CASE("the Clifford+T coverage is offset zero for its sets") {
    // 24, 96, 384 prefixes of scales 2, 3, 4, and larger sets
    for (const auto& [k, size] :
         std::vector<std::pair<int, int64_t>>{{2, 24}, {3, 96}, {4, 384}, {3, 600}})
        CHECK(split_coverage(GateSet::clifford_t(), {1, 3}, k, size) ==
              std::array<double, 3>{1.0, 0.0, 0.0});
    CHECK(split_coverage(GateSet::clifford_t(), {1, 3}, 3, 95) ==
          split_coverage(GateSet::sqrt_t(), {5, 7}, 3, 1));
}

TEST_CASE("the model is built once per setting and outlives the cache") {
    const auto a = model_of(1, 1, 3, 2e-4), b = model_of(1, 1, 3, 2e-4);
    CHECK(a == b);
    CHECK(model_of(1, 1, 3, 3e-4)->mu() != a->mu());
    const auto held = model_of(1, 1, 3, 1.23e-4);
    std::vector<double> want;
    for (int64_t c = 0; c < 60; c += 7) want.push_back(held->saving(held->start(), c));
    for (int i = 0; i < 40; ++i) model_of(1, 1, 3, 1e-4 * (1 + i / 100.0));  // > 32 cached
    const auto again = model_of(1, 1, 3, 1.23e-4);
    for (size_t i = 0; i < want.size(); ++i) {
        const int64_t c = static_cast<int64_t>(7 * i);
        CHECK(held->saving(held->start(), c) == want[i]);
        CHECK(again->saving(again->start(), c) == want[i]);
    }
}

TEST_CASE("Clifford+T ignores the price's q") {
    const auto a = model_of(2, 1, 0, 1e-3), b = model_of(2, 1, 3, 1e-3);
    CHECK(a->start() == b->start());
    CHECK(a->groups() == b->groups());
    CHECK(a->cost_cap() == b->cost_cap());
    for (int K = a->start(); K < a->start() + 5; ++K)
        for (int64_t c = 0; c < 40; c += 3) CHECK(a->saving(K, c) == b->saving(K, c));
    const Unitary& V = calibration_targets()[57];
    std::set<std::optional<int64_t>> costs;
    for (int64_t q : {0, 3})
        for (StoppingRule w : {StoppingRule::fixed, StoppingRule::model}) {
            const ScaleScan::Result r = scan_of(GateSet::clifford_t(), {1, q}, w).run(V);
            costs.insert(best_cost(r.best));
        }
    CHECK(costs.size() == 1);
}

// --- fallback -------------------------------------------------------------------------

TEST_CASE("the model falls back on targets with a stabiliser") {
    for (double theta : {0.3, 1.1, 2.9}) {
        for (char axis : {'x', 'y', 'z'})
            CHECK(stopping_fallback(unitary_of(rotation(theta, axis))));
        // Clifford conjugate of a diagonal, with a phase
        const Matrix D = rotation(theta, 'z');
        CHECK(stopping_fallback(unitary_of(mul(mul(mul(mul(kH, kS), D), dagger_of(kS)), kH))));
        CHECK(stopping_fallback(unitary_of(scaled(D, std::polar(1.0L, 0.7L)))));
    }
    // symmetric up to phase: V^T = V
    const Matrix A = mul(mul(rotation(0.4, 'z'), rotation(0.9, 'x')), rotation(0.4, 'z'));
    for (int i = 0; i < 4; ++i) CHECK(std::abs(transpose_of(A)[i] - A[i]) < 1e-15L);
    CHECK(stopping_fallback(unitary_of(A)));
    for (size_t i = 58; i < 78; ++i) CHECK_FALSE(stopping_fallback(calibration_targets()[i]));
}

TEST_CASE("the model falls back on an axis on a mirror plane") {
    // U -> C U^T C^-1 preserves cost and distance: hits come in pairs
    struct Case {
        double nx, ny, nz;
        Matrix C;
    };
    for (const Case& k : {Case{0.3, 0.0, 0.8, kI},                 // n_y = 0: V^T = V
                          Case{0.0, 0.6, 0.5, kZ},                 // n_x = 0
                          Case{0.7, 0.4, 0.0, kX},                 // n_z = 0
                          Case{0.5, 0.3, 0.5, kH},                 // |n_x| = |n_z|
                          Case{0.6, 0.6, 0.2, kS},                 // n_x = n_y
                          Case{0.6, -0.6, 0.2, dagger_of(kS)}}) {  // n_x = -n_y
        INFO(k.nx, " ", k.ny, " ", k.nz);
        const Matrix V = rotation(1.3, k.nx, k.ny, k.nz);
        CHECK(stopping_fallback(unitary_of(V)));
        Rng rng(2);
        for (int i = 0; i < 20; ++i) {
            const Matrix U = haar_unitary(rng);
            const Matrix W = mul(mul(k.C, transpose_of(U)), dagger_of(k.C));
            CHECK(std::abs(std::abs(trace_of(mul(dagger_of(W), V))) -
                           std::abs(trace_of(mul(dagger_of(U), V)))) < 1e-12L);
        }
    }
    CHECK_FALSE(stopping_fallback(unitary_of(rotation(1.3, 0.2, 0.5, 0.8))));
}

// --- stopping rule in the scan --------------------------------------------------------

TEST_CASE("the scan's own default is the fixed stopping rule") {
    // QuickSynth's default is the model stopping rule (test_policy.cpp)
    CHECK(ScaleScan::Settings{}.stopping_rule == StoppingRule::fixed);
    CHECK(scan_of(GateSet::sqrt_t(), {1, 3}, StoppingRule::model).model() != nullptr);
    CHECK(scan_of(GateSet::sqrt_t(), {1, 3}, StoppingRule::fixed).model() == nullptr);
}

TEST_CASE("a price beyond the model's cost cap runs the fixed stopping rule") {
    // no model at sqrt(T) price 1e9 or T price 1e6; always one at the catalysis
    // price
    CHECK(model_of(1, 1, 1000000000, 1e-3) == nullptr);
    CHECK(model_of(2, 1000000, 3, 1e-3) == nullptr);
    CHECK_FALSE(has_stopping_model(GateSet::sqrt_t(), Price{1, 1000000000}, 1e-3));
    CHECK_FALSE(has_stopping_model(GateSet::clifford_t(), Price{1000000, 3}, 1e-3));
    for (double eps : {1e-2, 1e-8, 1e-13}) {
        INFO(eps);
        const auto m = model_of(1, 1, 3, eps);
        REQUIRE(m != nullptr);
        CHECK(has_stopping_model(GateSet::sqrt_t(), Price{1, 3}, eps));
        CHECK(m->cost_cap() <= kMaxModelCosts);
    }
    const Unitary& V = calibration_targets()[78];
    for (const GateSet& g : {GateSet::sqrt_t(), GateSet::clifford_t()}) {
        INFO(g.name());
        const Price price = g.step == 1 ? Price{1, 1000000000} : Price{1000000, 3};
        const ScaleScan::Result fixed = scan_of(g, price, StoppingRule::fixed).run(V);
        const ScaleScan::Result model = scan_of(g, price, StoppingRule::model).run(V);
        REQUIRE(fixed.best.has_value());
        REQUIRE(model.best.has_value());
        CHECK(model.best->cost.value == fixed.best->cost.value);
        REQUIRE(model.log.size() == fixed.log.size());
        for (size_t i = 0; i < fixed.log.size(); ++i)
            CHECK(model.log[i].scale == fixed.log[i].scale);
    }
}

TEST_CASE("the log records every walk of the scan") {
    const Unitary& V = calibration_targets()[79];
    for (const GateSet& g : {GateSet::sqrt_t(), GateSet::clifford_t()})
        for (StoppingRule w : {StoppingRule::fixed, StoppingRule::model}) {
            INFO(g.name(), " ", w == StoppingRule::model ? "model" : "fixed");
            const Price price = g.step == 1 ? Price{1, 3} : Price{1, 1};
            std::vector<Enumerator::Result> walks;
            const ScaleScan::Result r =
                scan_of(g, price, w).run(V, [&](const Enumerator::Result& walk, double) {
                    walks.push_back(walk);
                });
            REQUIRE(r.log.size() == walks.size());
            double estimated = 0;
            std::vector<int64_t> costs;
            for (size_t i = 0; i < r.log.size(); ++i) {
                const SearchLog& e = r.log[i];
                CHECK_FALSE(e.prefix_scale.has_value());
                CHECK(e.scale == walks[i].K);
                CHECK(e.det_class == walks[i].det_class);
                CHECK(e.points == walks[i].points);
                CHECK(e.seconds == walks[i].seconds);
                estimated += e.estimated;
                if (e.best_cost) costs.push_back(*e.best_cost);
            }
            CHECK(std::is_sorted(costs.rbegin(), costs.rend()));
            REQUIRE_FALSE(costs.empty());
            REQUIRE(r.best.has_value());
            CHECK(costs.back() == r.best->cost.value);
            CHECK(close(estimated, r.estimated_points, 1e-6));
        }
}

TEST_CASE("the model scan stops at the bound or the value test") {
    size_t i = 80;
    const std::vector<std::pair<GateSet, Price>> settings{
        {GateSet::sqrt_t(), {1, 1}}, {GateSet::sqrt_t(), {2, 3}}, {GateSet::clifford_t(), {1, 1}}};
    for (const auto& [g, price] : settings)
        for (int n = 0; n < 3; ++n, ++i) {
            const auto label = g.name();  // C++17: INFO cannot capture a binding
            INFO(label, " calibration target ", i);
            const Unitary& V = calibration_targets()[i];
            const ScaleScan::Result fixed = scan_of(g, price, StoppingRule::fixed, 7.9e6).run(V);
            const ScaleScan::Result model = scan_of(g, price, StoppingRule::model, 7.9e6).run(V);
            CHECK(model.first_hit_scale == fixed.first_hit_scale);
            REQUIRE(fixed.best.has_value());
            REQUIRE(model.best.has_value());
            if (model.certified) CHECK(model.best->cost.value <= fixed.best->cost.value);
            if (g.step == 2) {
                // stopped by the integer bound, as the fixed stopping rule
                CHECK(model.certified);
                CHECK(fixed.certified);
                CHECK(model.best->cost.value == fixed.best->cost.value);
            }
        }
}

// --- finer levels ---------------------------------------------------------------------

TEST_CASE("the census by level sums to the census by scale") {
    // no structured bindings: INFO's lambda cannot capture them in C++17 (Clang)
    for (const auto& entry :
         {std::pair{GateSet::sqrt_t(), Price{1, 3}}, std::pair{GateSet::sqrt_t(), Price{2, 3}},
          std::pair{GateSet::clifford_t(), Price{1, 1}}}) {
        const GateSet& g = entry.first;
        const Price& price = entry.second;
        INFO(g.name(), " ", price.t, ",", price.q);
        const int kmax = 6;
        const int L = levels_per_scale(g);
        const Census scales(g, price, kmax);
        const Census levels(g, price, kmax, true);
        REQUIRE(levels.rows() == L * kmax + 1);
        CHECK(levels.cost_cap() == scales.cost_cap());
        for (int64_t c = 0; c <= scales.cost_cap(); ++c) {
            CHECK(levels.count(0, c) == scales.count(0, c));
            for (int K = 1; K <= kmax; ++K) {
                double sum = 0.0;
                for (int l = L * (K - 1) + 1; l <= L * K; ++l) sum += levels.count(l, c);
                CHECK(close(sum, scales.count(K, c), 1e-12, 0.0));
            }
        }
    }
}

TEST_CASE("without a split the census gain's scan equals the replay's per-level census") {
    // census.py of the fine-level replay at the catalysis price: (epsilon, cov, a,
    // best) -> scan saving
    struct Case {
        double eps;
        int cov, a;
        int64_t best;
        double scan;
    };
    const std::vector<Case> cases{
        {1e-4, 36, 40, 27, 0.42958256312432297},
        {1e-8, 76, 79, 54, 0.03775769332367615},
    };
    const GateSet g = GateSet::sqrt_t();
    for (const Case& t : cases) {
        INFO(t.eps);
        const CensusGain gain(g, {1, 3}, t.eps, SearchRegion(t.eps, g).first_scale());
        CHECK(close(gain.scan(t.cov, t.a, {}, t.best), t.scan, 1e-9));
    }
}

TEST_CASE("the census gain prices each prefix by the gates it reaches first") {
    // A model built from the prefixes' words (cost ball grown in Python, children
    // found by word): at 1e-6, the scan at level 60 with best cost 43, each prefix
    // of the radius-6 scale-2 set walked from 60 to 66.
    const GateSet g = GateSet::sqrt_t();
    const double eps = 1e-6;
    const PrefixSet set = collect_prefixes(*recorded_prefix_ball("sqrtT"), {2}, std::nullopt)[0];
    REQUIRE(set.prefixes.size() == 128);
    const CensusGain gain(g, {1, 3}, eps, SearchRegion(eps, g).first_scale(), &set);
    const std::vector<std::pair<double, int>> want{
        {0.04957893939051, 8},   {0.023357889697294, 8},  {0.019180543106759, 16},
        {0.013814884546511, 16}, {0.008989246984042, 16}, {0.005645040275914, 16},
        {0.004825602036229, 32}, {0.001677496643311, 8},  {0.0, 8}};
    std::vector<PrefixWalk> every;
    std::vector<int> classes;
    std::vector<std::pair<double, int>> got;
    for (size_t p = 0; p < set.prefixes.size(); ++p) {
        const double v = gain.split(60, {{p, 60, 66}}, 43);
        const auto same = [&](const std::pair<double, int>& w) { return close(w.first, v, 1e-9); };
        const auto it = std::find_if(got.begin(), got.end(), same);
        if (it == got.end())
            got.emplace_back(v, 1);
        else
            ++it->second;
        every.push_back({p, 60, 66});
        classes.push_back(gain.prefix_class(p));
    }
    std::sort(got.rbegin(), got.rend());
    REQUIRE(got.size() == want.size());
    for (size_t i = 0; i < want.size(); ++i) {
        CHECK(close(got[i].first, want[i].first, 1e-9, 1e-15));
        CHECK(got[i].second == want[i].second);
    }
    std::sort(classes.begin(), classes.end());
    CHECK(std::unique(classes.begin(), classes.end()) - classes.begin() == 10);
    // together they save less than the sum of their savings
    CHECK(close(gain.split(60, every, 43), 1.1573630689145498, 1e-9));
    // the scan of level 61 saves what no prefix reached
    CHECK(close(gain.scan(60, 61, std::vector<int>(128, 60), 43), 0.33311420172232875, 1e-9));
    CHECK(close(gain.scan(60, 61, std::vector<int>(128, 66), 43), 0.0252741403236901, 1e-9));
    CHECK(gain.batch(128) == kPrefixBatch);
}

TEST_CASE("a level gain grows with the levels and shrinks with the split's coverage") {
    const GateSet g = GateSet::sqrt_t();
    const double eps = 1e-6;
    const PrefixSet set = collect_prefixes(*recorded_prefix_ball("sqrtT"), {2}, std::nullopt)[0];
    const size_t n = set.prefixes.size();
    const CensusGain census(g, {1, 3}, eps, SearchRegion(eps, g).first_scale(), &set);
    const DecayGain decay(56, n);
    // every prefix from composite level `from` to `to`
    const auto all = [&](int from, int to) {
        std::vector<PrefixWalk> walks;
        for (size_t p = 0; p < n; ++p) walks.push_back({p, from, to});
        return walks;
    };
    const std::vector<int> none(n, 56), deep(n, 64);
    for (const LevelGain* gain : std::initializer_list<const LevelGain*>{&census, &decay}) {
        const int64_t best = 41;
        CHECK(gain->scan(56, 56, none, best) == 0.0);
        CHECK(gain->split(56, all(58, 58), best) == 0.0);
        double before = 0.0;
        for (int a = 57; a < 70; ++a) {
            const double s = gain->scan(56, a, none, best);
            CHECK(s > before);
            CHECK(s < static_cast<double>(best));
            CHECK(gain->scan(56, a, deep, best) < s);
            CHECK(gain->split(56, all(56, a), best) < s);
            before = s;
        }
        // a cheaper answer leaves less to save
        CHECK(gain->scan(56, 62, none, best - 5) < gain->scan(56, 62, none, best));
    }
}

TEST_CASE("the decay gain saves a share of the best cost that decays per level") {
    const DecayGain gain(10, 4);
    const std::vector<int> none(4, 10), next(4, 11);
    CHECK(gain.scan(4, 10, none, 100) == 0.0);
    CHECK(close(gain.scan(10, 11, none, 100), 100 * kDecayFirstSaving));
    CHECK(close(gain.scan(11, 12, none, 100), 100 * kDecayFirstSaving * kDecayPerLevel));
    std::vector<PrefixWalk> all;
    for (size_t p = 0; p < 4; ++p) all.push_back({p, 10, 11});
    CHECK(close(gain.split(10, all, 100), 100 * kDecaySplitCoverage * kDecayFirstSaving));
    // one prefix saves its share of the level; the split walks them all together
    CHECK(close(gain.split(10, {{2, 10, 11}}, 100), 25 * kDecaySplitCoverage * kDecayFirstSaving));
    CHECK(gain.batch(4) == 4);
    // the split's coverage of the scan's next level is deducted
    CHECK(close(gain.scan(10, 11, next, 100), 100 * (1 - kDecaySplitCoverage) * kDecayFirstSaving));
    CHECK(gain.rates(6.0).scan == 3.0);
    CHECK(gain.rates(6.0).split == 1.0);
}

TEST_CASE("targets outside the census get the decay gain") {
    const GateSet g = GateSet::sqrt_t();
    const int start = SearchRegion(1e-4, g).first_scale();
    const auto is_census = [](const std::unique_ptr<LevelGain>& gain) {
        return dynamic_cast<const CensusGain*>(gain.get()) != nullptr;
    };
    CHECK(is_census(level_gain(g, {1, 3}, 1e-4, calibration_targets()[58], start, 30, nullptr)));
    CHECK_FALSE(
        is_census(level_gain(g, {1, 3}, 1e-4, unitary_of(rotation(0.3, 'z')), start, 30, nullptr)));
    CHECK_FALSE(is_census(
        level_gain(g, {1, 1000000000}, 1e-4, calibration_targets()[58], start, 30, nullptr)));
}

TEST_CASE("the walk-time curve passes through its knots and grows with the points") {
    CHECK(walk_seconds(0.0) == 21e-6);
    CHECK(walk_seconds(0.05) == 21e-6);
    CHECK(close(walk_seconds(1e4), 3.5e-3));
    CHECK(close(walk_seconds(1e6), 0.123));
    // values of the replay's curve between and past the knots
    CHECK(close(walk_seconds(3.0), 4.040666185385785e-05));
    CHECK(close(walk_seconds(2e4), 0.0057775530676363865));
    CHECK(close(walk_seconds(4e6), 0.3478965363437814));
    double before = 0.0;
    for (double est = 0.1; est < 1e9; est *= 1.7) {
        CHECK(walk_seconds(est) > before);
        before = walk_seconds(est);
    }
}

}  // namespace quicksynth::test
