// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Prime generators for every odd p < 2000 and 40-bit p of every residue class:
// exact norm p^f, membership in (p, h(zeta)) by reduction mod (p, h), and the
// decomposition of p as the product of the Galois images, split or inert by the
// class of p, each K+ prime once.

#include "quicksynth/solve_norm_equation/prime_ideals.hpp"

#include <doctest/doctest.h>

#include "quicksynth/solve_norm_equation/integer_factor.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace quicksynth::test {

namespace {

template <size_t N>
using Ring = RingElement<N, __int128>;
template <size_t N>
using RingMp = RingElement<N, int_mp>;

int_mp power(uint64_t p, int f) {
    int_mp out(1);
    for (int i = 0; i < f; ++i) out = out * int_mp(static_cast<__int128>(p));
    return out;
}

// x lies in (p, h(zeta)): x mod p is a multiple of h.
template <size_t N>
bool in_prime(const Ring<N>& x, uint64_t p, const FpPoly& h) {
    FpPoly a;
    for (size_t i = 0; i < N; ++i) {
        __int128 c = x[i] % static_cast<__int128>(p);
        if (c < 0) c += p;
        a.c[i] = static_cast<uint64_t>(c);
    }
    a.degree = static_cast<int>(N) - 1;
    while (a.degree >= 0 && a.c[a.degree] == 0) --a.degree;
    return poly_divmod(a, h, p).second.degree < 0;
}

// a = unit * b.
template <size_t N>
bool associate(const RingMp<N>& a, const RingMp<N>& b) {
    const std::optional<RingMp<N>> q = ring_divide(a, b);
    return q && ring_norm(*q) == int_mp(1);
}

template <size_t N>
RingMp<N> rational(uint64_t p) {
    RingMp<N> out{};
    out[0] = int_mp(static_cast<__int128>(p));
    return out;
}

template <size_t N>
void check_generators(uint64_t p) {
    CAPTURE(p);
    CAPTURE(N);
    const int f = cyclotomic_residue_degree(p, static_cast<int>(N));
    for (const FpPoly& h : cyclotomic_factors(static_cast<int>(N), p)) {
        const std::optional<Ring<N>> s = short_generator<N>(p, h);
        const std::optional<Ring<N>> e = euclid_generator<N>(p, h);
        REQUIRE(s);
        CHECK(ring_norm(*s) == power(p, f));
        CHECK(in_prime<N>(*s, p, h));
        // Euclid's gcd is the same ideal; it may carry a different unit.
        if (e) {
            CHECK(ring_norm(*e) == power(p, f));
            CHECK(in_prime<N>(*e, p, h));
            CHECK(associate<N>(to_int_mp(*e), to_int_mp(*s)));
        }
        CHECK(prime_generator<N>(p, h) == s);
    }
}

template <size_t N>
void check_decomposition(uint64_t p, PrimeIdeals<N>& ideals) {
    CAPTURE(p);
    CAPTURE(N);
    const std::optional<PrimeDecomposition<N>> d = ideals.above(p);
    REQUIRE(d);
    const int f = cyclotomic_residue_degree(p, static_cast<int>(N));
    const bool inert = p % (2 * N) == 2 * N - 1;
    CHECK(d->f == f);
    CHECK(d->primes.size() == (inert ? N / f : N / (2 * f)));
    // The primes of K above p, each once: P, and P* for a split pair.
    std::vector<RingMp<N>> all;
    for (const PrimeAbove<N>& q : d->primes) {
        CHECK(q.inert == inert);
        CHECK(ring_norm(q.pi) == power(p, f));
        CHECK(q.norm == relative_norm(q.pi));
        const RingMp<N> pi = to_int_mp(q.pi);
        all.push_back(pi);
        if (inert)
            CHECK(associate<N>(conjugate(pi), pi));
        else
            all.push_back(conjugate(pi));
        // pi pi* generates the K+ prime below (split) or its square (inert).
        const RealElement<N / 2, int_mp> below = relative_norm(pi);
        CHECK(real_norm(below) == power(p, f));
        RealElement<N / 2, int_mp> rational_p;
        rational_p[0] = int_mp(static_cast<__int128>(p));
        CHECK(real_divide(inert ? real_mul(rational_p, rational_p) : rational_p, below));
        if (inert) CHECK_FALSE(real_divide(rational_p, below));
    }
    for (size_t i = 0; i < all.size(); ++i)
        for (size_t j = 0; j < i; ++j) CHECK_FALSE(associate<N>(all[i], all[j]));
    // p is unramified: the product of the primes above it.
    RingMp<N> product = ring_zeta_power<N, int_mp>(0);
    for (const RingMp<N>& pi : all) product = ring_mul(product, pi);
    CHECK(associate<N>(rational<N>(p), product));
}

std::vector<uint64_t> odd_primes_below(uint64_t bound) {
    std::vector<uint64_t> out;
    for (uint64_t p = 3; p < bound; p += 2)
        if (is_prime_u64(p)) out.push_back(p);
    return out;
}

// The first `count` primes above 2^40 in each class mod 16.
std::vector<uint64_t> primes_40_bits(int count) {
    std::vector<uint64_t> out;
    for (uint64_t r = 1; r < 16; r += 2) {
        int found = 0;
        for (uint64_t p = (uint64_t{1} << 40) + r; found < count; p += 16)
            if (is_prime_u64(p)) {
                out.push_back(p);
                ++found;
            }
    }
    return out;
}

}  // namespace

TEST_CASE("prime generators have norm p^f and lie in their prime for p below 2000") {
    for (uint64_t p : odd_primes_below(2000)) {
        check_generators<4>(p);
        check_generators<8>(p);
    }
}

TEST_CASE("prime generators of 40-bit primes in every class mod 16") {
    for (uint64_t p : primes_40_bits(3)) {
        check_generators<4>(p);
        check_generators<8>(p);
    }
}

TEST_CASE("the Galois images of one generator decompose p") {
    PrimeIdeals<4> omega;
    PrimeIdeals<8> zeta16;
    std::vector<uint64_t> primes = odd_primes_below(2000);
    for (uint64_t p : primes_40_bits(2)) primes.push_back(p);
    for (uint64_t p : primes) {
        check_decomposition<4>(p, omega);
        check_decomposition<8>(p, zeta16);
    }
}

TEST_CASE("the prime cache returns what it computed and survives eviction") {
    PrimeIdeals<8> cached(2);
    const std::optional<PrimeDecomposition<8>> first = cached.above(17);
    REQUIRE(first);
    for (uint64_t p : {3, 5, 7, 11, 13}) REQUIRE(cached.above(p));
    for (uint64_t p : {17, 17}) {
        const std::optional<PrimeDecomposition<8>> again = cached.above(p);
        REQUIRE(again);
        REQUIRE(again->primes.size() == first->primes.size());
        for (size_t i = 0; i < first->primes.size(); ++i)
            CHECK(again->primes[i].pi == first->primes[i].pi);
    }
    CHECK_THROWS_AS(cached.above(2), std::invalid_argument);
}

TEST_CASE("prime generators reject inputs outside their contract") {
    const FpPoly h = cyclotomic_factors(4, 17).front();
    CHECK_THROWS_AS(short_generator<4>(16, h), std::invalid_argument);
    CHECK_THROWS_AS(euclid_generator<4>(uint64_t{1} << 63 | 1, h), std::invalid_argument);
    FpPoly whole;
    whole.c = {1, 0, 0, 0, 1};
    whole.degree = 4;
    CHECK_THROWS_AS(prime_generator<4>(17, whole), std::invalid_argument);
}

}  // namespace quicksynth::test
