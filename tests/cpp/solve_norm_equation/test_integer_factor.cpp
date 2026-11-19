// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Primality, factorisation and square roots mod p against a sieve, trial
// division and known strong pseudoprimes; every factorisation is multiplied
// back.

#include "quicksynth/solve_norm_equation/integer_factor.hpp"

#include <doctest/doctest.h>

#include "oracles.hpp"

#include <algorithm>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace quicksynth::test {

namespace {

// Smallest prime factor of every n below `limit` (0 and 1 map to 0).
std::vector<uint32_t> smallest_factors(uint32_t limit) {
    std::vector<uint32_t> spf(limit, 0);
    for (uint32_t i = 2; i < limit; ++i) {
        if (spf[i] != 0) continue;
        for (uint64_t j = i; j < limit; j += i)
            if (spf[j] == 0) spf[j] = i;
    }
    return spf;
}

// Primality by trial division; for the oracle's primes below 2^32.
bool prime_by_trial(uint64_t n) {
    if (n < 2) return false;
    for (uint64_t d = 2; d * d <= n; ++d)
        if (n % d == 0) return false;
    return true;
}

// A random prime of `bits` bits, below 2^32.
uint64_t random_prime(Rng& rng, int bits) {
    for (;;) {
        const uint64_t n = (uint64_t{1} << (bits - 1)) |
                           static_cast<uint64_t>(rng.below(int64_t{1} << (bits - 1))) | 1;
        if (prime_by_trial(n)) return n;
    }
}

// Checks that `f` is a valid factorisation of n: increasing primes with
// positive exponents whose product is n.
void check_factorisation(uint64_t n, const std::optional<std::vector<PrimePower>>& f) {
    REQUIRE(f.has_value());
    unsigned __int128 product = 1;
    uint64_t last = 0;
    for (const PrimePower& pe : *f) {
        CHECK(pe.p > last);
        CHECK(pe.e >= 1);
        CHECK(is_prime_u64(pe.p));
        last = pe.p;
        for (int i = 0; i < pe.e; ++i) product *= pe.p;
    }
    CHECK(product == n);
}

}  // namespace

TEST_CASE("is_prime_u64 agrees with a sieve below 10^7") {
    constexpr uint32_t kLimit = 10'000'000;
    const std::vector<uint32_t> spf = smallest_factors(kLimit);
    uint32_t wrong = 0, primes = 0;
    for (uint32_t n = 0; n < kLimit; ++n) {
        const bool prime = n >= 2 && spf[n] == n;
        primes += prime;
        wrong += prime != is_prime_u64(n);
    }
    CHECK(primes == 664579);
    CHECK(wrong == 0);
}

TEST_CASE("is_prime_u64 rejects the strong pseudoprimes to the first prime bases") {
    // The least strong pseudoprimes to the bases 2; 2, 3; ...; 2 .. 23, and
    // Carmichael numbers.
    for (uint64_t n : {uint64_t{2047}, uint64_t{1373653}, uint64_t{25326001}, uint64_t{3215031751},
                       uint64_t{2152302898747}, uint64_t{3474749660383}, uint64_t{341550071728321},
                       uint64_t{3825123056546413051}, uint64_t{561}, uint64_t{41041},
                       uint64_t{825265}, uint64_t{321197185}})
        CHECK_FALSE(is_prime_u64(n));
    CHECK(uint64_t{3825123056546413051} == uint64_t{149491} * 747451 * 34233211);
}

TEST_CASE("is_prime_u64 accepts the largest primes below powers of two") {
    for (uint64_t p : {uint64_t{65537}, uint64_t{2147483647}, uint64_t{998244353},
                       uint64_t{1000000007}, (uint64_t{1} << 61) - 1, (uint64_t{1} << 62) - 57,
                       (uint64_t{1} << 63) - 25, uint64_t{18446744073709551557u}}) {
        CHECK(is_prime_u64(p));
    }
    const uint64_t m31 = 2147483647, m61 = (uint64_t{1} << 61) - 1;
    CHECK_FALSE(is_prime_u64(m31 * m31));
    CHECK_FALSE(is_prime_u64(3 * m61));
    CHECK_FALSE(is_prime_u64(UINT64_MAX));
}

TEST_CASE("factor_u64 of 1 is empty and of 0 throws") {
    const auto one = factor_u64(1);
    REQUIRE(one.has_value());
    CHECK(one->empty());
    CHECK_THROWS_AS(factor_u64(0), std::invalid_argument);
}

TEST_CASE("factor_u64 agrees with a sieve below 2^18") {
    constexpr uint32_t kLimit = 1 << 18;
    const std::vector<uint32_t> spf = smallest_factors(kLimit);
    for (uint32_t n = 2; n < kLimit; ++n) {
        std::vector<PrimePower> expected;
        for (uint32_t m = n; m > 1; m /= spf[m]) {
            if (!expected.empty() && expected.back().p == spf[m])
                ++expected.back().e;
            else
                expected.push_back({spf[m], 1});
        }
        const auto f = factor_u64(n);
        REQUIRE(f.has_value());
        if (*f != expected) FAIL_CHECK("n = " << n);
    }
}

TEST_CASE("factor_u64 splits 50-bit and 63-bit semiprimes") {
    Rng rng(31);
    for (int i = 0; i < 200; ++i) {
        const uint64_t p = random_prime(rng, 25), q = random_prime(rng, 25);
        const auto f = factor_u64(p * q);
        check_factorisation(p * q, f);
        CHECK(f->size() == (p == q ? 1u : 2u));
    }
    // Two factors near 2^31.5: rho's longest case below 2^63.
    for (int i = 0; i < 20; ++i) {
        const uint64_t p = random_prime(rng, 32), q = random_prime(rng, 31);
        const auto f = factor_u64(p * q);
        check_factorisation(p * q, f);
        CHECK(f->size() == (p == q ? 1u : 2u));
    }
}

TEST_CASE("factor_u64 splits prime powers and products of many primes") {
    Rng rng(32);
    const uint64_t m31 = 2147483647;
    check_factorisation(m31 * m31, factor_u64(m31 * m31));
    CHECK(factor_u64(m31 * m31)->front() == PrimePower{m31, 2});
    const uint64_t p = 1031;  // the first prime above the trial bound
    uint64_t power = 1;
    for (int e = 1; e <= 6; ++e) {
        power *= p;
        const auto f = factor_u64(power);
        REQUIRE(f.has_value());
        CHECK(*f == std::vector<PrimePower>{{p, e}});
    }
    for (int i = 0; i < 300; ++i) {
        // Up to five primes above 2^10 with repeats, below 2^63.
        uint64_t n = 1;
        const int count = 1 + static_cast<int>(rng.below(5));
        for (int k = 0; k < count; ++k) {
            const int bits = 11 + static_cast<int>(rng.below(std::min(32, 63 / count) - 10));
            const uint64_t q = random_prime(rng, bits);
            n *= q;
            if (rng.below(3) == 0 && k + 1 < count && n < (uint64_t{1} << 62) / q) n *= q;
        }
        // Small factors too.
        n *= uint64_t{1} << rng.below(2);
        check_factorisation(n, factor_u64(n));
    }
}

TEST_CASE("factor_u64 gives nullopt when rho has no attempts left") {
    const uint64_t p = 1000003, q = 1000033;
    CHECK_FALSE(factor_u64(p * q, 0).has_value());
    // Trial division and primality need no attempt.
    const auto prime = factor_u64(uint64_t{2} * 3 * 1000003, 0);
    REQUIRE(prime.has_value());
    CHECK(*prime == std::vector<PrimePower>{{2, 1}, {3, 1}, {1000003, 1}});
}

TEST_CASE("sqrt_mod finds every square mod small primes") {
    const std::vector<uint32_t> spf = smallest_factors(2000);
    for (uint64_t p = 2; p < 2000; ++p) {
        if (spf[p] != p) continue;
        std::vector<bool> square(p, false);
        for (uint64_t x = 0; x < p; ++x) square[x * x % p] = true;
        for (uint64_t a = 0; a < p; ++a) {
            const auto r = sqrt_mod(a, p);
            REQUIRE(r.has_value() == square[a]);
            if (r) {
                CHECK(*r * *r % p == a);
                CHECK(*r <= p - *r);
            }
        }
    }
}

TEST_CASE("sqrt_mod is exact mod large primes of every residue class mod 16") {
    Rng rng(33);
    // 998244353 = 119 2^23 + 1 and 65537 = 2^16 + 1 exercise Tonelli-Shanks'
    // longest chains; the others cover each odd class mod 16 near 2^40 and 2^62.
    std::vector<uint64_t> primes = {65537, 998244353, (uint64_t{1} << 61) - 1,
                                    (uint64_t{1} << 63) - 25, uint64_t{18446744073709551557u}};
    for (uint64_t base : {uint64_t{1} << 40, uint64_t{1} << 62})
        for (uint64_t r = 1; r < 16; r += 2) {
            uint64_t p = base + r;
            while (!is_prime_u64(p)) p += 16;
            primes.push_back(p);
        }
    for (uint64_t p : primes) {
        int squares = 0;
        for (int i = 0; i < 200; ++i) {
            const uint64_t x = rng.engine()() % p;
            const uint64_t a = rng.below(2) == 0 ? mul_mod(x, x, p) : rng.engine()() % p;
            const auto r = sqrt_mod(a, p);
            const bool residue = a == 0 || pow_mod(a, (p - 1) / 2, p) == 1;
            REQUIRE(r.has_value() == residue);
            if (r) {
                ++squares;
                CHECK(mul_mod(*r, *r, p) == a);
                CHECK(*r <= p - *r);
            }
        }
        CHECK(squares > 0);
    }
}

TEST_CASE("sqrt_mod rejects a composite modulus") {
    CHECK_THROWS_AS(sqrt_mod(1, 65), std::invalid_argument);
}

}  // namespace quicksynth::test
