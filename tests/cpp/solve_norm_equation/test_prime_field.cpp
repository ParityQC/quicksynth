// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// F_p polynomial arithmetic against its defining identities, the oracle's
// factorisation steps (fp_poly_oracle.hpp) on products of irreducibles found by
// root search, and the closed-form factors of x^4 + 1 and x^8 + 1 mod primes of
// every residue class against their product, the splitting law of cyclotomic
// fields and Cantor-Zassenhaus.

#include "quicksynth/solve_norm_equation/prime_field.hpp"

#include <doctest/doctest.h>

#include "fp_poly_oracle.hpp"
#include "oracles.hpp"
#include "quicksynth/solve_norm_equation/integer_factor.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace quicksynth::test {

namespace {

uint64_t random_below(Rng& rng, uint64_t p) { return rng.engine()() % p; }

// Uniform coefficients up to `degree`, leading one nonzero.
FpPoly random_poly(Rng& rng, int degree, uint64_t p) {
    FpPoly a;
    for (int i = 0; i < degree; ++i) a.c[i] = random_below(rng, p);
    a.c[degree] = 1 + random_below(rng, p - 1);
    a.degree = degree;
    return a;
}

FpPoly random_monic(Rng& rng, int degree, uint64_t p) {
    FpPoly a = random_poly(rng, degree, p);
    a.c[degree] = 1;
    return a;
}

bool has_root(const FpPoly& a, uint64_t p) {
    for (uint64_t x = 0; x < p; ++x) {
        uint64_t v = 0;
        for (int i = a.degree; i >= 0; --i) v = (mul_mod(v, x, p) + a.c[i]) % p;
        if (v == 0) return true;
    }
    return false;
}

// A monic irreducible of degree 2 or 3 (irreducible iff rootless) not in `avoid`.
FpPoly random_irreducible(Rng& rng, int degree, uint64_t p, const std::vector<FpPoly>& avoid) {
    for (;;) {
        const FpPoly a = random_monic(rng, degree, p);
        if (!has_root(a, p) && std::find(avoid.begin(), avoid.end(), a) == avoid.end()) return a;
    }
}

FpPoly product(const std::vector<FpPoly>& factors, uint64_t p) {
    FpPoly r = fp_poly({1}, p);
    for (const FpPoly& f : factors) r = poly_mul(r, f, p);
    return r;
}

bool divides(const FpPoly& d, const FpPoly& a, uint64_t p) {
    return poly_divmod(a, d, p).second.degree < 0;
}

FpPoly phi(int n) {
    FpPoly f;
    f.c[0] = 1;
    f.c[n] = 1;
    f.degree = n;
    return f;
}

const std::vector<uint64_t> kFieldPrimes = {
    3, 13, 101, 65537, (uint64_t{1} << 61) - 1, (uint64_t{1} << 63) - 25};

// Odd primes below 2000, for each odd class mod 16 the first primes in it
// from 2^40 and from 2^62, and random primes of 11 to 63 bits.
std::vector<uint64_t> odd_primes() {
    std::vector<uint64_t> primes;
    Rng rng(46);
    while (primes.size() < 100) {
        const int bits = 11 + static_cast<int>(rng.below(53));
        const uint64_t p = (rng.engine()() >> (64 - bits)) | (uint64_t{1} << (bits - 1)) | 1;
        if (is_prime_u64(p)) primes.push_back(p);
    }
    for (uint64_t p = 3; p < 2000; p += 2)
        if (is_prime_u64(p)) primes.push_back(p);
    for (uint64_t base : {uint64_t{1} << 40, uint64_t{1} << 62})
        for (uint64_t r = 1; r < 16; r += 2) {
            uint64_t p = base + r;
            for (int found = 0; found < 3; p += 16)
                if (is_prime_u64(p)) {
                    primes.push_back(p);
                    ++found;
                }
        }
    return primes;
}

void check_cyclotomic(int n) {
    for (uint64_t p : odd_primes()) {
        CAPTURE(p);
        const int f = cyclotomic_residue_degree(p, n);
        const std::vector<FpPoly> factors = cyclotomic_factors(n, p);
        REQUIRE(static_cast<int>(factors.size()) * f == n);
        for (size_t i = 0; i < factors.size(); ++i) {
            // Every irreducible factor of Phi_{2n} mod p has degree f, so a
            // factor of degree f is irreducible.
            CHECK(factors[i].degree == f);
            CHECK(factors[i].c[f] == 1);
            if (i > 0) CHECK(factors[i - 1] < factors[i]);
        }
        CHECK(product(factors, p) == phi(n));
        // The closed forms against the general random splitting.
        if (f < n) CHECK(equal_degree_factors(phi(n), f, p) == factors);
        // The distinct-degree step sees one degree.
        const auto groups = distinct_degree_factors(phi(n), p);
        REQUIRE(groups.size() == 1);
        CHECK(groups[0].first == f);
        CHECK(groups[0].second == phi(n));
    }
}

}  // namespace

TEST_CASE("F_p division leaves a remainder of lower degree") {
    Rng rng(41);
    for (uint64_t p : kFieldPrimes)
        for (int i = 0; i < 200; ++i) {
            const FpPoly a = random_poly(rng, static_cast<int>(rng.below(9)), p);
            const FpPoly b = random_poly(rng, static_cast<int>(rng.below(9)), p);
            const auto [q, r] = poly_divmod(a, b, p);
            CHECK(r.degree < b.degree);
            CHECK(poly_add(poly_mul(q, b, p), r, p) == a);
            CHECK(poly_sub(poly_add(a, b, p), b, p) == a);
        }
    CHECK_THROWS_AS(poly_divmod(fp_poly({1, 1}, 7), FpPoly{}, 7), std::domain_error);
    CHECK_THROWS_AS(poly_mul(phi(8), fp_poly({0, 1}, 7), 7), std::length_error);
    CHECK_THROWS_AS(fp_poly({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 11), std::length_error);
    CHECK(fp_poly({7, 14, 0}, 7).degree == -1);
}

TEST_CASE("F_p gcds are monic common divisors") {
    Rng rng(42);
    for (uint64_t p : kFieldPrimes)
        for (int i = 0; i < 200; ++i) {
            const FpPoly g = random_poly(rng, 1 + static_cast<int>(rng.below(3)), p);
            const FpPoly a = poly_mul(g, random_poly(rng, static_cast<int>(rng.below(6)), p), p);
            const FpPoly b = poly_mul(g, random_poly(rng, static_cast<int>(rng.below(6)), p), p);
            const FpPoly d = poly_gcd(a, b, p);
            CHECK(d.c[d.degree] == 1);
            CHECK(divides(d, a, p));
            CHECK(divides(d, b, p));
            CHECK(divides(g, d, p));
        }
    const FpPoly a = fp_poly({4, 2}, 7);
    CHECK(poly_gcd(a, FpPoly{}, 7) == fp_poly({2, 1}, 7));
    CHECK(poly_gcd(FpPoly{}, FpPoly{}, 7).degree == -1);
}

TEST_CASE("F_p powers agree with repeated multiplication") {
    Rng rng(43);
    for (uint64_t p : kFieldPrimes)
        for (int i = 0; i < 50; ++i) {
            const FpPoly m = random_poly(rng, 1 + static_cast<int>(rng.below(8)), p);
            const FpPoly a = random_poly(rng, static_cast<int>(rng.below(9)), p);
            FpPoly power = poly_divmod(fp_poly({1}, p), m, p).second;
            for (uint64_t e = 0; e < 40; ++e) {
                CHECK(poly_powmod(a, e, m, p) == power);
                power = poly_mulmod(power, a, m, p);
            }
        }
}

TEST_CASE("distinct-degree factorisation groups the factors by degree") {
    Rng rng(44);
    for (uint64_t p : {uint64_t{3}, uint64_t{13}, uint64_t{101}, uint64_t{1009}}) {
        CAPTURE(p);
        const uint64_t r1 = random_below(rng, p);
        const uint64_t r2 = (r1 + 1 + random_below(rng, p - 1)) % p;
        const FpPoly linear = poly_mul(fp_poly({p - r1, 1}, p), fp_poly({p - r2, 1}, p), p);
        const FpPoly quadratic = random_irreducible(rng, 2, p, {});
        const FpPoly cubic = random_irreducible(rng, 3, p, {});
        const auto groups = distinct_degree_factors(product({linear, quadratic, cubic}, p), p);
        REQUIRE(groups.size() == 3);
        CHECK(groups[0] == std::make_pair(1, linear));
        CHECK(groups[1] == std::make_pair(2, quadratic));
        CHECK(groups[2] == std::make_pair(3, cubic));
        const auto cubics = distinct_degree_factors(cubic, p);
        REQUIRE(cubics.size() == 1);
        CHECK(cubics[0] == std::make_pair(3, cubic));
    }
}

TEST_CASE("Cantor-Zassenhaus recovers distinct irreducibles of one degree") {
    Rng rng(45);
    for (uint64_t p : {uint64_t{3}, uint64_t{13}, uint64_t{101}, uint64_t{1009}}) {
        CAPTURE(p);
        for (int degree : {1, 2, 3}) {
            // As many as fit degree 8, at most the number of monic irreducibles.
            const uint64_t available = degree == 1   ? p
                                       : degree == 2 ? (p * p - p) / 2
                                                     : (p * p * p - p) / 3;
            const int count = static_cast<int>(std::min<uint64_t>(8 / degree, available));
            std::vector<FpPoly> factors;
            while (static_cast<int>(factors.size()) < count) {
                const FpPoly h = degree == 1 ? fp_poly({random_below(rng, p), 1}, p)
                                             : random_irreducible(rng, degree, p, factors);
                if (std::find(factors.begin(), factors.end(), h) == factors.end())
                    factors.push_back(h);
            }
            std::sort(factors.begin(), factors.end());
            CHECK(equal_degree_factors(product(factors, p), degree, p) == factors);
        }
    }
    CHECK_THROWS_AS(equal_degree_factors(phi(8), 3, 17), std::invalid_argument);
    CHECK_THROWS_AS(equal_degree_factors(phi(4), 2, 2), std::invalid_argument);
}

TEST_CASE("the residue degree of p follows its class mod 8 and mod 16") {
    // Z[omega]: 1 splits completely, 3 5 7 into two primes of degree 2.
    // Z[zeta16]: 1 completely, 7 9 15 into four of degree 2, 3 5 11 13 into two of degree 4.
    const int deg4[8] = {0, 1, 0, 2, 0, 2, 0, 2};
    const int deg8[16] = {0, 1, 0, 4, 0, 4, 0, 2, 0, 2, 0, 4, 0, 4, 0, 2};
    for (uint64_t p : odd_primes()) {
        CHECK(cyclotomic_residue_degree(p, 4) == deg4[p % 8]);
        CHECK(cyclotomic_residue_degree(p, 8) == deg8[p % 16]);
    }
    CHECK_THROWS_AS(cyclotomic_residue_degree(2, 4), std::invalid_argument);
    CHECK_THROWS_AS(cyclotomic_residue_degree(3, 6), std::invalid_argument);
}

TEST_CASE("x^4 + 1 factors mod p into factors of the residue degree") { check_cyclotomic(4); }

TEST_CASE("x^8 + 1 factors mod p into factors of the residue degree") { check_cyclotomic(8); }

TEST_CASE("cyclotomic_factors rejects other n and even p") {
    CHECK_THROWS_AS(cyclotomic_factors(2, 17), std::invalid_argument);
    CHECK_THROWS_AS(cyclotomic_factors(8, 2), std::invalid_argument);
}

}  // namespace quicksynth::test
