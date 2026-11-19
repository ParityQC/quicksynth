// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The z-rotation search against the joint walk's hits at every level it
// walked, its Clifford+T level bound against priced completions, exact angles, the Clifford+T
// optimum of Ross-Selinger with the level order corrected (pygridsynth, tolerance 2 eps, as
// golden_rotation.json's clifford_t table), the log's one entry per (K, sub, class), cut walks
// and the stops that clear the certificate, the Clifford+sqrt(T) count rule and caps, and the
// lazy phase sweep against the full one.

#include "quicksynth/synthesise_rotation/rotation_search.hpp"

#include <doctest/doctest.h>

#include "quicksynth/synthesise_rotation/rotation_body.hpp"
#include "quicksynth/solve_norm_equation/relative_norm.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/number_field.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

namespace quicksynth::test {

namespace {

Unitary z_rotation(double theta) {
    return {std::cos(theta / 2), -std::sin(theta / 2), 0, 0, 0, 0,
            std::cos(theta / 2), std::sin(theta / 2)};
}

const Price kCliffordT = Price::clifford_t(1);
const Price kCatalysis{1, 3};

RotationSearch::Result search(const GateSet& g, const Price& price, double theta, double eps,
                              const RotationSearch::Settings& settings = {}) {
    return RotationSearch(g, price, settings).run(z_rotation(theta), eps);
}

// Clifford+sqrt(T) settings that stop by the bound alone.
RotationSearch::Settings without_count() {
    RotationSearch::Settings s;
    s.solvable = INT32_MAX;
    s.levels_past_hit = 1000;
    return s;
}

int level_of(const SearchLog& e, int L) { return L * e.scale - e.sub; }

// Level of the last log entry.
int last_level(const GateSet& g, const RotationSearch::Result& r) {
    const int L = g.degree() / 2;
    return L * r.log.back().scale - r.log.back().sub;
}

// The least cost of the joint walk's hits over every (K, sub, class) up to
// level `last`.
std::optional<int64_t> walk_minimum(const GateSet& g, const Price& price, double theta, double eps,
                                    int last) {
    Enumerator enumerator(eps, g, price, 1, {});
    const int L = g.degree() / 2;
    std::optional<int64_t> best;
    for (int l = 0; l <= last; ++l) {
        const int K = (l + L - 1) / L, sub = L * K - l;
        for (int j = 0; j < GateSet::kDeterminantClasses; ++j)
            for (const Enumerator::Hit& h : enumerator.enumerate(z_rotation(theta), K, j, sub).hits)
                if (!best || h.cost.value < *best) best = h.cost.value;
    }
    return best;
}

void check_against_walk(const GateSet& g, const Price& price, double eps, int angles,
                        uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> angle(0, 2 * M_PI);
    for (int trial = 0; trial < angles; ++trial) {
        const double theta = angle(rng);
        const RotationSearch::Result r = search(g, price, theta, eps, without_count());
        INFO(g.name(), " eps ", eps, " theta ", theta);
        REQUIRE(r.best);
        CHECK(r.stop == RotationStop::bound);
        CHECK(r.certified);
        CHECK(chord_distance(matrix_direction(r.best->element), su2_vector(z_rotation(theta))) <=
              eps);
        CHECK(walk_minimum(g, price, theta, eps, last_level(g, r)) == r.best->cost.value);
    }
}

}  // namespace

TEST_CASE("the rotation search finds the joint walk's least cost on Clifford+T") {
    check_against_walk(GateSet::clifford_t(), kCliffordT, 1e-2, 6, 11);
    check_against_walk(GateSet::clifford_t(), kCliffordT, 1e-3, 3, 12);
}

TEST_CASE("the rotation search finds the joint walk's least cost on Clifford+sqrt(T)") {
    check_against_walk(GateSet::sqrt_t(), kCatalysis, 3e-2, 3, 13);
}

TEST_CASE("Clifford+T completions at level l cost l - 2 or more with the class parity") {
    const GateSet g = GateSet::clifford_t();
    const Pricer pricer = Pricer::for_gateset(g, kCliffordT);
    PrimeIdeals<4> primes;
    std::mt19937_64 rng(21);
    std::uniform_real_distribution<double> angle(0, 2 * M_PI);
    int64_t priced = 0, at_bound = 0;
    for (const double eps : {1e-2, 1e-3}) {
        const int last = RotationSearch::first_useful_level(g, eps) + 3;
        for (int trial = 0; trial < 4; ++trial) {
            const Unitary D = z_rotation(angle(rng));
            const Vec4 target = su2_vector(D);
            Enumerator completer(eps, g, kCliffordT, 1, {});
            for (int l = 0; l <= last; ++l) {
                const int K = (l + 1) / 2, sub = 2 * K - l;
                for (int j = 0; j < GateSet::kDeterminantClasses; ++j) {
                    const RotationBody<4> body(search_direction(D, j, g), sub, eps);
                    for (const auto& p : body.points(K)) {
                        if (!p.new_at_level || !body.contains(p.u, K)) continue;
                        RealElement<2, __int128> power{};
                        power[0] = static_cast<__int128>(1) << K;
                        const NormSolutions<4> s = solve_relative_norm<4>(
                            real_sub(power, relative_norm(p.u)), primes, kChoiceCap, true);
                        if (s.outcome != NormOutcome::solved) continue;
                        std::vector<int64_t> rows;
                        for (const auto& t : s.choices)
                            for (int i = 0; i < 2; ++i) {
                                const auto phased = ring_times_zeta(t, i);
                                for (int k = 0; k < 4; ++k)
                                    rows.push_back(static_cast<int64_t>(p.u[k]));
                                for (int k = 0; k < 4; ++k)
                                    rows.push_back(static_cast<int64_t>(phased[k]));
                            }
                        for (const CyclotomicMatrix& e :
                             completer.matrices(rows.data(), rows.size() / 8, K, j)) {
                            if (chord_distance(matrix_direction(e), target) > eps) continue;
                            const int64_t cost = pricer.price(e).value;
                            INFO("l ", l, " class ", j, " cost ", cost);
                            CHECK(cost >= l - 2);
                            CHECK(cost % 2 == j);
                            ++priced;
                            at_bound += cost <= l - 1;
                        }
                    }
                }
            }
        }
    }
    CHECK(priced > 100);
    CHECK(at_bound > 0);
}

TEST_CASE("exact z-rotations cost their magic gates") {
    for (int k = -8; k <= 8; ++k) {
        INFO("k ", k);
        const RotationSearch::Result t =
            search(GateSet::clifford_t(), kCliffordT, k * M_PI / 4, 1e-6);
        REQUIRE(t.best);
        CHECK(t.best->cost.value == (k % 2 != 0));
        CHECK(t.certified);
        // pi/8 is one sqrt(T) at price 3, below any three T
        const RotationSearch::Result q = search(GateSet::sqrt_t(), kCatalysis, k * M_PI / 8, 1e-6);
        REQUIRE(q.best);
        CHECK(q.best->cost.value == (k % 2 != 0 ? 3 : k % 4 != 0 ? 1 : 0));
        CHECK(q.best->cost.n_q == (k % 2 != 0));
        CHECK(q.certified);
    }
}

TEST_CASE("the rotation search equals the corrected Ross-Selinger T-count") {
    struct Case {
        double theta, eps;
        int64_t t_count;
    };
    // golden_rotation.json's clifford_t angles 0, 1, 2 and 500, 501 and
    // near-exact angles, pygridsynth up to phase at tolerance 2 eps, levels in
    // T order.
    const Case cases[] = {
        {1.3217781747557364, 1e-3, 28},
        {1.3443949205279058, 1e-3, 29},
        {2.541722128811673, 1e-3, 28},
        {1.2e-5, 1.7e-6, 65},
        {0.01, 2.5e-8, 76},
        {M_PI / 4 + 3e-6, 1e-6, 72},
        {1e-3, 1e-4, 42},
    };
    for (const Case& c : cases) {
        INFO("theta ", c.theta, " eps ", c.eps);
        const RotationSearch::Result r = search(GateSet::clifford_t(), kCliffordT, c.theta, c.eps);
        REQUIRE(r.best);
        CHECK(r.best->cost.value == c.t_count);
        CHECK(r.best->cost.n_t == c.t_count);
        CHECK(r.certified);
    }
}

TEST_CASE("the rotation log holds one walk per level and class") {
    const GateSet g = GateSet::clifford_t();
    const RotationSearch::Result r = search(g, kCliffordT, 0.7, 1e-4);
    REQUIRE(r.best);
    RotationWalk total;
    for (size_t i = 0; i < r.log.size(); ++i) {
        const SearchLog& e = r.log[i];
        const int l = 2 * e.scale - e.sub;
        REQUIRE(e.det_class);
        // walk order; before the first hit no walk is skipped
        if (!r.log[i].best_cost) {
            CHECK(l == static_cast<int>(i / 2));
            CHECK(e.det_class == static_cast<int>(i % 2));
        }
        if (i > 0) {
            const SearchLog& p = r.log[i - 1];
            CHECK(2 * (2 * p.scale - p.sub) + *p.det_class < 2 * l + *e.det_class);
        }
        CHECK((e.sub == 0 || e.sub == 1));
        REQUIRE(e.rotation);
        CHECK(e.rotation->candidates <= e.points);
        CHECK(e.rotation->solvable <= e.rotation->candidates);
        total.candidates += e.rotation->candidates;
        total.solvable += e.rotation->solvable;
        total.pricings += e.rotation->pricings;
    }
    CHECK(total.candidates == r.counts.candidates);
    CHECK(total.solvable == r.counts.solvable);
    CHECK(total.pricings == r.counts.pricings);
    CHECK(r.counts.solvable > 0);
    CHECK(r.log.back().best_cost == r.best->cost.value);
    // stopped once the next level's bound, T-count >= l - 2, reaches the best
    // cost, not before
    const int last = last_level(g, r);
    CHECK(last + 1 - 2 >= r.best->cost.value);
    CHECK(last - 2 < r.best->cost.value);
    CHECK(r.summary.source == Source::rotation);
    CHECK(r.summary.scan_cost == r.best->cost.value);
    CHECK(r.summary.scan_certified);
    CHECK(r.summary.scales.back() == r.log.back().scale);
}

TEST_CASE("a cut walk or the level cap clears the certificate") {
    const GateSet g = GateSet::clifford_t();
    const Unitary D = z_rotation(0.7);
    const RotationSearch::Result whole = RotationSearch(g, kCliffordT, {}).run(D, 1e-4);
    REQUIRE(whole.best);
    // the search goes on past a cut walk, to the same least cost here
    RotationSearch::Settings few;
    few.point_limit = 2;
    const RotationSearch::Result cut = RotationSearch(g, kCliffordT, few).run(D, 1e-4);
    CHECK(cut.stop == RotationStop::bound);
    CHECK_FALSE(cut.certified);
    REQUIRE(cut.best);
    CHECK(cut.best->cost.value == whole.best->cost.value);
    CHECK(cut.summary.last_full_scale < whole.summary.last_full_scale);
    int64_t cuts = 0;
    for (const SearchLog& e : cut.log) {
        CHECK(e.points <= 3);
        cuts += e.points == 3;
    }
    CHECK(cuts > 1);
    RotationSearch::Settings shallow;
    shallow.levels_past_useful = -40;
    const RotationSearch::Result capped = RotationSearch(g, kCliffordT, shallow).run(D, 1e-4);
    CHECK(capped.stop == RotationStop::level_cap);
    CHECK_FALSE(capped.certified);
    CHECK_FALSE(capped.best);
}

TEST_CASE("the Clifford+sqrt(T) count rule finishes whole levels until enough u are solvable") {
    const GateSet g = GateSet::sqrt_t();
    for (const double theta : {0.7, 2.3, 4.1}) {
        std::optional<int64_t> previous;
        for (const int S : {4, 16, 64}) {
            INFO("theta ", theta, " S ", S);
            RotationSearch::Settings settings;
            settings.solvable = S;
            const RotationSearch::Result r = search(g, kCatalysis, theta, 1e-3, settings);
            REQUIRE(r.best);
            CHECK(r.stop == RotationStop::count);
            CHECK_FALSE(r.certified);
            const int last = level_of(r.log.back(), 4);
            int64_t before = 0, total = 0, walks_at_last = 0;
            for (const SearchLog& e : r.log) {
                total += e.rotation->solvable;
                if (level_of(e, 4) < last) before += e.rotation->solvable;
                walks_at_last += level_of(e, 4) == last;
            }
            CHECK(total == r.counts.solvable);
            CHECK(before < S);
            CHECK(total >= S);
            CHECK(walks_at_last == GateSet::kDeterminantClasses);
            // a larger budget walks the same levels further, so its best is no worse
            if (previous) CHECK(r.best->cost.value <= *previous);
            previous = r.best->cost.value;
            // deterministic
            const RotationSearch::Result again = search(g, kCatalysis, theta, 1e-3, settings);
            CHECK(again.best->element.matrix == r.best->element.matrix);
            CHECK(again.best->element.denominator_exponent == r.best->element.denominator_exponent);
            CHECK(again.log.size() == r.log.size());
        }
    }
}

// Within a few eps of a multiple of pi/4 the body of some levels is a slab of
// millions of points (rotation_body.hpp).
TEST_CASE("near-exact angles answer past cut walks at the reference T-count") {
    const GateSet g = GateSet::clifford_t();
    struct Case {
        double theta, eps;
        int64_t t_count;  // pygridsynth, as above
        bool certified;
    };
    const Case cases[] = {
        {3e-8, 1e-8, 99, false},
        {-1e-9, 1e-10, 121, true},
        {M_PI / 4 - 1e-9, 1e-10, 120, true},
    };
    for (const Case& c : cases) {
        INFO("theta ", c.theta, " eps ", c.eps);
        const RotationSearch::Result r = search(g, kCliffordT, c.theta, c.eps);
        REQUIRE(r.best);
        CHECK(r.stop == RotationStop::bound);
        CHECK(r.best->cost.value == c.t_count);
        CHECK(r.certified == c.certified);
        bool cut = false;
        for (const SearchLog& e : r.log) cut = cut || e.points > kRotationPointLimit;
        CHECK(cut == !c.certified);
    }
}

TEST_CASE("a Clifford+sqrt(T) walk is cut at a multiple of the count rule's u") {
    const GateSet g = GateSet::sqrt_t();
    RotationSearch::Settings settings;
    settings.solvable = 2;
    const RotationSearch::Result r = search(g, kCatalysis, 3e-6, 1e-6, settings);
    REQUIRE(r.best);
    CHECK(r.stop == RotationStop::count);
    int64_t cut = 0;
    for (const SearchLog& e : r.log) {
        CHECK(e.rotation->solvable <= kRotationWalkSolvable * settings.solvable);
        cut += e.rotation->solvable == kRotationWalkSolvable * settings.solvable;
    }
    CHECK(cut >= 1);
}

TEST_CASE("the Clifford+sqrt(T) caps leave no answer to the count rule") {
    const GateSet g = GateSet::sqrt_t();
    const double eps = 1e-3;
    RotationSearch::Settings settings;
    settings.solvable = INT32_MAX;
    settings.levels_past_hit = 0;
    const RotationSearch::Result r = search(g, kCatalysis, 0.7, eps, settings);
    CHECK(r.stop == RotationStop::hit_cap);
    CHECK_FALSE(r.certified);
    // the cap counts from the first solvable level, or the first useful one if later
    int first = -1;
    for (const SearchLog& e : r.log)
        if (first < 0 && e.rotation->solvable > 0) first = level_of(e, 4);
    REQUIRE(first >= 0);
    CHECK(level_of(r.log.back(), 4) == std::max(first, RotationSearch::first_useful_level(g, eps)));

    RotationSearch::Settings late;
    late.deadline = 0.0;
    const RotationSearch::Result timed = search(g, kCatalysis, 0.7, eps, late);
    CHECK(timed.stop == RotationStop::time_limit);
    CHECK(timed.log.empty());
    CHECK_FALSE(timed.best);
}

// Phase t -> t zeta^i is conjugation by diag(1, zeta)^i, whose least weight is
// q_cost for odd i and min(t_cost, 2 q_cost) for i = 2; the lazy sweep relies on
// the phases' costs differing by at most twice that.
TEST_CASE("Clifford+sqrt(T) phases of a completion differ by at most twice the phase's weight") {
    const GateSet g = GateSet::sqrt_t();
    PrimeIdeals<8> primes;
    std::mt19937_64 rng(41);
    std::uniform_real_distribution<double> angle(0, 2 * M_PI);
    int64_t completions = 0, spread = 0;
    for (const Price price : {Price{1, 3}, Price{2, 2}, Price{3, 2}}) {
        const Pricer pricer = Pricer::for_gateset(g, price);
        const int64_t weight[4] = {0, price.q, std::min(price.t, 2 * price.q), price.q};
        const double eps = 1e-2;
        const int last = RotationSearch::first_useful_level(g, eps) + 4;
        for (int trial = 0; trial < 2; ++trial) {
            const Unitary D = z_rotation(angle(rng));
            Enumerator completer(eps, g, price, 1, {});
            for (int l = 0; l <= last; ++l) {
                const int K = (l + 3) / 4, sub = 4 * K - l;
                for (int j = 0; j < GateSet::kDeterminantClasses; ++j) {
                    const RotationBody<8> body(search_direction(D, j, g), sub, eps);
                    for (const auto& p : body.points(K)) {
                        if (!p.new_at_level || !body.contains(p.u, K)) continue;
                        RealElement<4, __int128> power{};
                        power[0] = static_cast<__int128>(1) << K;
                        const NormSolutions<8> s = solve_relative_norm<8>(
                            real_sub(power, relative_norm(p.u)), primes, kChoiceCap, false);
                        if (s.outcome != NormOutcome::solved) continue;
                        for (const auto& t : s.choices) {
                            int64_t cost[4];
                            for (int i = 0; i < 4; ++i) {
                                const auto phased = ring_times_zeta(t, i);
                                std::vector<int64_t> row;
                                for (int k = 0; k < 8; ++k)
                                    row.push_back(static_cast<int64_t>(p.u[k]));
                                for (int k = 0; k < 8; ++k)
                                    row.push_back(static_cast<int64_t>(phased[k]));
                                const auto e = completer.matrices(row.data(), 1, K, j);
                                REQUIRE(e.size() == 1);
                                cost[i] = pricer.price(e[0]).value;
                            }
                            for (int i = 1; i < 4; ++i) {
                                INFO("price ", price.t, ",", price.q, " phase ", i);
                                CHECK(std::abs(cost[i] - cost[0]) <= 2 * weight[i]);
                                spread += cost[i] != cost[0];
                            }
                            ++completions;
                        }
                    }
                }
            }
        }
    }
    CHECK(completions > 100);
    CHECK(spread > 0);
}

TEST_CASE("the lazy phase sweep finds the full sweep's element with or without halving") {
    const GateSet g = GateSet::sqrt_t();
    std::mt19937_64 rng(43);
    std::uniform_real_distribution<double> angle(0, 2 * M_PI);
    int64_t lazy_pricings = 0, full_pricings = 0;
    for (const double eps : {1e-3, 1e-4}) {
        for (int trial = 0; trial < 3; ++trial) {
            const double theta = angle(rng);
            INFO("theta ", theta, " eps ", eps);
            RotationSearch::Settings full;
            full.lazy_phases = false;
            RotationSearch::Settings whole = full;
            whole.transpose_halving = false;
            RotationSearch::Settings lazy_whole;
            lazy_whole.transpose_halving = false;
            const RotationSearch::Result a = search(g, kCatalysis, theta, eps);
            const RotationSearch::Result b = search(g, kCatalysis, theta, eps, full);
            const RotationSearch::Result c = search(g, kCatalysis, theta, eps, whole);
            const RotationSearch::Result d = search(g, kCatalysis, theta, eps, lazy_whole);
            REQUIRE(a.best);
            REQUIRE(b.best);
            REQUIRE(c.best);
            REQUIRE(d.best);
            // the same walk order with fewer pricings: the same element
            CHECK(a.best->element.matrix == b.best->element.matrix);
            CHECK(a.best->element.denominator_exponent == b.best->element.denominator_exponent);
            CHECK(d.best->element.matrix == c.best->element.matrix);
            // halving drops transposes, which cost the same
            CHECK(a.best->cost.value == c.best->cost.value);
            CHECK(a.log.size() == c.log.size());
            for (size_t i = 0; i < a.log.size(); ++i)
                CHECK(a.log[i].best_cost == c.log[i].best_cost);
            lazy_pricings += a.counts.pricings;
            full_pricings += b.counts.pricings;
        }
    }
    CHECK(2 * lazy_pricings < full_pricings);
}

TEST_CASE("the rotation search checks its target") {
    const RotationSearch s(GateSet::clifford_t(), kCliffordT, {});
    CHECK_THROWS_AS(s.run(Unitary{0.6, 0, 0.8, 0, -0.8, 0, 0.6, 0}, 1e-3), std::invalid_argument);
    CHECK_THROWS_AS(s.run(z_rotation(0.3), 0.0), std::invalid_argument);
    CHECK_THROWS_AS(s.run(z_rotation(0.3), 1.0), std::invalid_argument);
}

}  // namespace quicksynth::test
