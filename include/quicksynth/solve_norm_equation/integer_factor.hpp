// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Arithmetic on unsigned 64-bit integers for the norm equation: primality,
// factorisation of the norm N_{K+/Q}(xi) and square roots mod p. Products go
// through unsigned __int128, so every modulus below 2^64 is exact.
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

namespace quicksynth {

#pragma GCC visibility push(hidden)

// a + b and a - b mod m, for a, b < m.
inline uint64_t add_mod(uint64_t a, uint64_t b, uint64_t m) {
    return a >= m - b ? a - (m - b) : a + b;
}
inline uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t m) { return a >= b ? a - b : a + (m - b); }

// a b mod m and a^e mod m, for m >= 1.
uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t m);
uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t m);

// Deterministic Miller-Rabin: the twelve prime bases 2 .. 37 have no common
// strong pseudoprime below 3.1e23 [Sorenson, Webster, Math. Comp. 86 (2017)],
// so the answer is exact for every n < 2^64.
bool is_prime_u64(uint64_t n);

struct PrimePower {
    uint64_t p;
    int e;
    bool operator==(const PrimePower& o) const { return p == o.p && e == o.e; }
};

// Pollard-Brent attempts before factor_u64 gives up on a cofactor.
inline constexpr int kRhoAttempts = 32;

// The prime factorisation of n >= 1, by increasing p (empty for n = 1); throws
// std::invalid_argument for n = 0. Trial division to 2^10, then Pollard's rho
// with Brent's cycle detection and batched gcds [Brent, BIT 20 (1980)]: attempt
// i iterates x -> x^2 + (i + 1) from x = 2, so the result does not depend on
// any random state. A cofactor that `attempts` attempts leave composite gives
// nullopt.
std::optional<std::vector<PrimePower>> factor_u64(uint64_t n, int attempts = kRhoAttempts);

// x with x^2 = a mod p and x <= p - x, or nullopt if a is not a square mod p;
// p prime. Tonelli-Shanks [Cohen, GTM 138, Algorithm 1.5.1].
std::optional<uint64_t> sqrt_mod(uint64_t a, uint64_t p);

#pragma GCC visibility pop

}  // namespace quicksynth
