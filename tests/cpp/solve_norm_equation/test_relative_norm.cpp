// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The relative norm solver against three oracles: every t of small coefficient
// norm (t t* has constant coefficient sum c_i^2, so a box search is
// exhaustive), PARI/GP solution counts at N_{K+/Q}(xi) of 2^20 .. 2^45
// (tests/fixtures/golden_relative_norm.json), and round trips xi = t t* of
// random t at the same sizes. The transpose halving is checked by pricing every
// phase of both halves on elements of random words.

#include "quicksynth/solve_norm_equation/relative_norm.hpp"

#include <doctest/doctest.h>

#include "fixtures.hpp"
#include "oracles.hpp"
#include "quicksynth/data_structures/cost.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/price_candidates/price.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace quicksynth::test {

namespace {

template <size_t N>
using Ring = RingElement<N, __int128>;
template <size_t N>
using Real = RealElement<N / 2, __int128>;

constexpr int kNoCap = 1 << 20;

// The least of zeta^j t over j, as a key for t up to roots of unity.
template <size_t N>
Ring<N> canonical(const Ring<N>& t) {
    Ring<N> best = t;
    for (int j = 1; j < 2 * static_cast<int>(N); ++j) best = std::min(best, ring_times_zeta(t, j));
    return best;
}

template <size_t N>
std::set<Ring<N>> canonical_set(const std::vector<Ring<N>>& ts) {
    std::set<Ring<N>> out;
    for (const Ring<N>& t : ts) out.insert(canonical<N>(t));
    return out;
}

// Every solution is exact, the choices are distinct up to roots of unity and
// there is one per choice.
template <size_t N>
void check_solved(const Real<N>& xi, const NormSolutions<N>& s) {
    REQUIRE(s.outcome == NormOutcome::solved);
    CHECK_FALSE(s.capped);
    CHECK(static_cast<int64_t>(s.choices.size()) == s.count);
    for (const Ring<N>& t : s.choices) CHECK(relative_norm(t) == xi);
    CHECK(canonical_set<N>(s.choices).size() == s.choices.size());
}

// t* = zeta^m times the choice with every exponent mirrored.
template <size_t N>
void check_mirror(const NormSolutions<N>& s) {
    const size_t n = s.choices.size();
    for (size_t c = 0; c < n; ++c)
        CHECK(canonical<N>(conjugate(s.choices[c])) == canonical<N>(s.choices[n - 1 - c]));
}

template <size_t N>
using Oracle = std::map<std::array<__int128, N / 2>, std::set<Ring<N>>>;

// All t with sum c_i^2 <= bound, by the tower coordinates of t t*.
template <size_t N>
Oracle<N> small_solutions(int64_t bound) {
    Oracle<N> out;
    Ring<N> t{};
    std::function<void(size_t, int64_t)> fill = [&](size_t i, int64_t left) {
        if (i == N) {
            out[relative_norm(t).c].insert(canonical<N>(t));
            return;
        }
        const auto r = static_cast<int64_t>(std::sqrt(static_cast<double>(left)));
        for (int64_t c = -r; c <= r; ++c) {
            if (c * c > left) continue;
            t[i] = c;
            fill(i + 1, left - c * c);
        }
        t[i] = 0;
    };
    fill(0, bound);
    return out;
}

template <size_t N>
void check_exhaustive(int64_t bound, const std::vector<int64_t>& box) {
    const Oracle<N> oracle = small_solutions<N>(bound);
    PrimeIdeals<N> primes;
    int64_t solved = 0, none = 0;
    // Every totally positive xi with constant coefficient up to `bound`; the
    // other tower coordinates are bounded by box[i] times it.
    Real<N> xi;
    std::function<void(size_t)> fill = [&](size_t i) {
        if (i == N / 2) {
            if (!totally_positive(xi)) return;
            CAPTURE(static_cast<int64_t>(xi[0]));
            CAPTURE(static_cast<int64_t>(xi[1]));
            const NormSolutions<N> s = solve_relative_norm<N>(xi, primes, kNoCap);
            const auto it = oracle.find(xi.c);
            if (it == oracle.end()) {
                CHECK(s.outcome == NormOutcome::none);
                ++none;
            } else {
                check_solved<N>(xi, s);
                CHECK(canonical_set<N>(s.choices) == it->second);
                ++solved;
            }
            return;
        }
        const int64_t r = box[i] * static_cast<int64_t>(xi[0]);
        for (int64_t c = -r; c <= r; ++c) {
            xi[i] = c;
            fill(i + 1);
        }
    };
    for (int64_t a0 = 1; a0 <= bound; ++a0) {
        xi = Real<N>{};
        xi[0] = a0;
        fill(1);
    }
    MESSAGE("N = ", N, ": ", solved, " solvable and ", none, " unsolvable xi");
    CHECK(solved > 0);
    CHECK(none > 0);
}

int64_t uniform_in(Rng& rng, int64_t bound) { return rng.below(2 * bound + 1) - bound; }

template <size_t N>
double log2_norm(const Ring<N>& t) {
    const int_mp n = ring_norm(t);
    return n.zero() ? -1 : std::log2(n.to_double(0));
}

// A random t with log2 N(t t*) in [20, 45]: coefficient bound from the
// target, or a product of small elements to small powers (repeated primes).
template <size_t N>
Ring<N> random_t(Rng& rng, bool products) {
    for (;;) {
        Ring<N> t{};
        if (products) {
            t = ring_times_zeta(ring_zeta_power<N, __int128>(0),
                                static_cast<int>(rng.below(2 * N)));
            const int factors = 1 + static_cast<int>(rng.below(4));
            for (int f = 0; f < factors; ++f) {
                Ring<N> s{};
                for (__int128& c : s) c = uniform_in(rng, 2);
                if (log2_norm<N>(s) <= 1) continue;
                const int e = 1 + static_cast<int>(rng.below(3));
                for (int i = 0; i < e; ++i) t = ring_mul(t, s);
            }
        } else {
            const double bits = 20 + 25 * rng.uniform();
            // N(t t*) ~ (sum c_i^2)^(N/2) with sum c_i^2 ~ N B^2 / 3
            const double b = std::sqrt(3.0 / N * std::pow(2.0, bits / (N / 2)));
            for (__int128& c : t)
                c = uniform_in(rng, std::max<int64_t>(1, static_cast<int64_t>(b)));
        }
        const double bits = log2_norm<N>(t);
        if (bits >= 20 && bits <= 45) return t;
    }
}

template <size_t N>
void check_round_trips(uint64_t seed, int count, bool products) {
    Rng rng(seed);
    PrimeIdeals<N> primes;
    int64_t largest = 0;
    for (int i = 0; i < count; ++i) {
        const Ring<N> t = random_t<N>(rng, products);
        const Real<N> xi = relative_norm(t);
        const NormSolutions<N> s = solve_relative_norm<N>(xi, primes, kNoCap);
        check_solved<N>(xi, s);
        CHECK(canonical_set<N>(s.choices).count(canonical<N>(t)) == 1);
        check_mirror<N>(s);
        // the halved list is the first half of the full one
        const NormSolutions<N> half = solve_relative_norm<N>(xi, primes, kNoCap, true);
        REQUIRE(half.count == (s.count + 1) / 2);
        for (int64_t c = 0; c < half.count; ++c) CHECK(half.choices[c] == s.choices[c]);
        largest = std::max(largest, s.count);
    }
    MESSAGE("N = ", N, ": up to ", largest, " choices per xi");
}

template <size_t N>
Real<N> real_of(const json& coefficients) {
    Ring<N> x{};
    for (size_t i = 0; i < N; ++i) x[i] = coefficients[i].get<int64_t>();
    return to_real(x);
}

template <size_t N>
void check_pari(const json& c, PrimeIdeals<N>& primes) {
    const std::string family = c[1];
    const int64_t count = c[2];
    const Real<N> xi = real_of<N>(c[3]);
    CAPTURE(family);
    CAPTURE(c[3].dump());
    const NormSolutions<N> s = solve_relative_norm<N>(xi, primes, kNoCap);
    if (count == 0) {
        CHECK(s.outcome == NormOutcome::none);
    } else {
        check_solved<N>(xi, s);
        CHECK(s.count == count);
    }
}

// Z[zeta_2N] coefficients into Z[zeta16] and back.
template <size_t N>
Z16 to_z16(const Ring<N>& x) {
    Z16 out{};
    for (size_t i = 0; i < N; ++i) out[i * (8 / N)] = static_cast<int64_t>(x[i]);
    return out;
}

template <size_t N>
Ring<N> from_z16(const Z16& x) {
    Ring<N> out{};
    for (size_t i = 0; i < 8; ++i) {
        if (i % (8 / N) != 0) {
            REQUIRE(x[i] == 0);
            continue;
        }
        out[i / (8 / N)] = x[i];
    }
    return out;
}

bool same_cost(const Cost& a, const Cost& b) {
    return a.value == b.value && a.n_t == b.n_t && a.n_q == b.n_q;
}

// Prices [[u, -t'* w], [t', u* w]] / sqrt2^k for t' = zeta^i t_c, every
// phase i and every choice c of t t* = xi, with u, t and w = det / 2^k from
// the element of a random word; the multisets of costs of c and its mirror
// count - 1 - c agree.
template <size_t N>
void check_transpose_costs(const GateSet& gateset, const std::string& letters, uint64_t seed) {
    const Pricer pricer = Pricer::for_gateset(gateset, Price{1, 3});
    PrimeIdeals<N> primes;
    Rng rng(seed);
    // The premise: an element and its transpose cost the same.
    for (int i = 0; i < 200; ++i) {
        const std::string word = random_word(rng, letters, 48);
        INFO(word);
        const CyclotomicMatrix e = word_matrix(word);
        CyclotomicMatrix transposed = e;
        std::swap(transposed.matrix[0][1], transposed.matrix[1][0]);
        CHECK(same_cost(pricer.price(transposed), pricer.price(e)));
    }
    int tested = 0, pairs = 0;
    for (int attempt = 0; attempt < 4000 && tested < 40; ++attempt) {
        const std::string word = random_word(rng, letters, 48);
        const CyclotomicMatrix e = word_matrix(word);
        const Ring<N> u = from_z16<N>(e.matrix[0][0]), t = from_z16<N>(e.matrix[1][0]);
        const Ring<N> det = from_z16<N>(zeta_sub(zeta_mul(e.matrix[0][0], e.matrix[1][1]),
                                                 zeta_mul(e.matrix[0][1], e.matrix[1][0])));
        Ring<N> w = det;
        for (__int128& c : w) {
            REQUIRE(c % (__int128{1} << e.denominator_exponent) == 0);
            c >>= e.denominator_exponent;
        }
        const Real<N> xi = relative_norm(t);
        if (xi == Real<N>{} || real_norm(xi).bits() > 62) continue;
        const NormSolutions<N> s = solve_relative_norm<N>(xi, primes, kNoCap);
        REQUIRE(s.outcome == NormOutcome::solved);
        if (s.count < 2) continue;
        INFO(word);
        ++tested;
        const Cost own = pricer.price(e);
        bool found = false;
        std::vector<std::vector<int64_t>> costs(s.choices.size());
        for (size_t c = 0; c < s.choices.size(); ++c)
            for (int i = 0; i < 2 * static_cast<int>(N); ++i) {
                const Ring<N> ti = ring_times_zeta(s.choices[c], i);
                CyclotomicMatrix m;
                m.matrix[0][0] = to_z16<N>(u);
                m.matrix[1][0] = to_z16<N>(ti);
                m.matrix[0][1] = to_z16<N>(ring_neg(ring_mul(conjugate(ti), w)));
                m.matrix[1][1] = to_z16<N>(ring_mul(conjugate(u), w));
                m.denominator_exponent = e.denominator_exponent;
                reduce_sqrt2(m.matrix, m.denominator_exponent);
                const Cost cost = pricer.price(m);
                costs[c].push_back(cost.value * 1000000 + cost.n_t * 1000 + cost.n_q);
                found |= ti == t && same_cost(cost, own);
            }
        CHECK(found);
        for (size_t c = 0; c < costs.size(); ++c) {
            std::vector<int64_t> a = costs[c], b = costs[costs.size() - 1 - c];
            std::sort(a.begin(), a.end());
            std::sort(b.begin(), b.end());
            CHECK(a == b);
            ++pairs;
        }
    }
    MESSAGE(std::string(gateset.name()), ": ", tested, " elements, ", pairs,
            " choices priced at every phase");
    CHECK(tested >= 20);
}

}  // namespace

TEST_CASE("relative norm solutions equal the exhaustive search at small norms") {
    check_exhaustive<4>(60, {0, 1});
    check_exhaustive<8>(10, {0, 1, 2, 2});
}

TEST_CASE("relative norm solutions match PARI counts at norms 2^20 to 2^45") {
    PrimeIdeals<4> omega;
    PrimeIdeals<8> zeta16;
    int checked = 0;
    for (const json& c : golden("relative_norm")["cases"]) {
        if (c[0] == 4)
            check_pari<4>(c, omega);
        else
            check_pari<8>(c, zeta16);
        ++checked;
    }
    CHECK(checked == 2696);
}

TEST_CASE("relative norm round trips recover random t at norms 2^20 to 2^45") {
    check_round_trips<4>(41, 200, false);
    check_round_trips<4>(42, 120, true);
    check_round_trips<8>(81, 200, false);
    check_round_trips<8>(82, 120, true);
}

TEST_CASE("an inert prime to an odd power has no solution") {
    PrimeIdeals<4> omega;
    PrimeIdeals<8> zeta16;
    for (int64_t p : {7, 23, 31, 47, 71}) {
        Real<4> xi;
        xi[0] = p;
        CHECK(solve_relative_norm<4>(xi, omega).outcome == NormOutcome::none);
        xi[0] = p * p * p;
        CHECK(solve_relative_norm<4>(xi, omega).outcome == NormOutcome::none);
        xi[0] = p * p;
        const NormSolutions<4> s = solve_relative_norm<4>(xi, omega);
        check_solved<4>(xi, s);
        CHECK(s.count == 1);
    }
    for (int64_t p : {31, 47, 79, 127}) {
        Real<8> xi;
        xi[0] = p;
        CHECK(solve_relative_norm<8>(xi, zeta16).outcome == NormOutcome::none);
        xi[0] = p * p;
        check_solved<8>(xi, solve_relative_norm<8>(xi, zeta16));
    }
}

TEST_CASE("the relative norm solver handles zero and signs and the cap and large norms") {
    PrimeIdeals<8> primes;
    const NormSolutions<8> zero = solve_relative_norm<8>(Real<8>{}, primes);
    REQUIRE(zero.outcome == NormOutcome::solved);
    CHECK(zero.choices == std::vector<Ring<8>>{Ring<8>{}});
    Real<8> negative;
    negative[0] = -5;
    CHECK(solve_relative_norm<8>(negative, primes).outcome == NormOutcome::none);
    // 1 - sqrt2 is negative at the identity place, its conjugate positive
    Real<4> mixed;
    mixed[0] = 1;
    mixed[1] = -1;
    PrimeIdeals<4> omega;
    CHECK(solve_relative_norm<4>(mixed, omega).outcome == NormOutcome::none);
    // 17 splits into 4 pairs of primes: 17^3 has 4^4 choices
    Real<8> many;
    many[0] = 4913;
    const NormSolutions<8> capped = solve_relative_norm<8>(many, primes);
    REQUIRE(capped.outcome == NormOutcome::solved);
    CHECK(capped.count == 256);
    CHECK(capped.capped);
    CHECK(capped.choices.size() == static_cast<size_t>(kChoiceCap));
    for (const Ring<8>& t : capped.choices) CHECK(relative_norm(t) == many);
    // N(x) = x^2 past 2^63: undecided, never asserted
    Real<4> big;
    big[0] = (__int128{1} << 32) + 15;
    const NormSolutions<4> u = solve_relative_norm<4>(big, omega);
    CHECK(u.outcome == NormOutcome::undecided);
    CHECK(u.undecided == NormUndecided::norm_too_large);
}

TEST_CASE("transpose halving keeps every cost") {
    check_transpose_costs<4>(GateSet::clifford_t(), "HST", 4);
    check_transpose_costs<8>(GateSet::sqrt_t(), "HSQ", 8);
}

}  // namespace quicksynth::test
