// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/solve_norm_equation/integer_factor.hpp"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

namespace quicksynth {

namespace {

// Trial division bound; every cofactor left has only prime factors above it.
constexpr uint64_t kTrialBound = 1024;
// Brent's cycle length cap per attempt. A cofactor below 2^64 with no factor
// below 2^10 has one below 2^32, which rho finds after about 2^16 steps.
constexpr uint64_t kRhoCycleLimit = uint64_t{1} << 22;
// Products of |x - y| between gcds [Brent 1980, section 7].
constexpr uint64_t kRhoBatch = 128;

uint64_t distance(uint64_t a, uint64_t b) { return a > b ? a - b : b - a; }

// A nontrivial factor of the odd composite n > 4, or 0 if this attempt fails.
uint64_t brent_rho(uint64_t n, uint64_t c) {
    const auto f = [n, c](uint64_t y) { return add_mod(mul_mod(y, y, n), c, n); };
    uint64_t y = 2, x = y, ys = y, q = 1, g = 1;
    for (uint64_t r = 1; g == 1; r *= 2) {
        if (r > kRhoCycleLimit) return 0;
        x = y;
        for (uint64_t i = 0; i < r; ++i) y = f(y);
        for (uint64_t k = 0; k < r && g == 1; k += kRhoBatch) {
            ys = y;
            for (uint64_t i = 0; i < std::min(kRhoBatch, r - k); ++i) {
                y = f(y);
                q = mul_mod(q, distance(x, y), n);
            }
            g = std::gcd(q, n);
        }
    }
    // The batch overshot to a multiple of n: step from its start one at a time.
    if (g == n) {
        do {
            ys = f(ys);
            g = std::gcd(distance(x, ys), n);
        } while (g == 1);
    }
    return g == n ? 0 : g;
}

}  // namespace

uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t m) {
    return static_cast<uint64_t>(static_cast<unsigned __int128>(a) * b % m);
}

uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t m) {
    uint64_t r = 1 % m;
    a %= m;
    for (; e != 0; e >>= 1) {
        if (e & 1) r = mul_mod(r, a, m);
        a = mul_mod(a, a, m);
    }
    return r;
}

bool is_prime_u64(uint64_t n) {
    static constexpr uint64_t kBases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return false;
    for (uint64_t p : kBases)
        if (n % p == 0) return n == p;
    uint64_t d = n - 1;
    int s = 0;
    for (; (d & 1) == 0; d >>= 1) ++s;
    for (uint64_t a : kBases) {
        uint64_t x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (int i = 1; i < s && witness; ++i) {
            x = mul_mod(x, x, n);
            if (x == n - 1) witness = false;
        }
        if (witness) return false;
    }
    return true;
}

std::optional<std::vector<PrimePower>> factor_u64(uint64_t n, int attempts) {
    if (n == 0) throw std::invalid_argument("factor_u64: n is 0");
    std::vector<uint64_t> primes;
    for (uint64_t p = 2; p < kTrialBound && p * p <= n; p += (p == 2 ? 1 : 2))
        for (; n % p == 0; n /= p) primes.push_back(p);
    std::vector<uint64_t> pending;
    if (n > 1) pending.push_back(n);
    while (!pending.empty()) {
        const uint64_t m = pending.back();
        pending.pop_back();
        if (is_prime_u64(m)) {
            primes.push_back(m);
            continue;
        }
        uint64_t d = 0;
        for (int i = 0; i < attempts && d == 0; ++i) d = brent_rho(m, static_cast<uint64_t>(i) + 1);
        if (d == 0) return std::nullopt;
        pending.push_back(d);
        pending.push_back(m / d);
    }
    std::sort(primes.begin(), primes.end());
    std::vector<PrimePower> out;
    for (uint64_t p : primes) {
        if (!out.empty() && out.back().p == p)
            ++out.back().e;
        else
            out.push_back({p, 1});
    }
    return out;
}

std::optional<uint64_t> sqrt_mod(uint64_t a, uint64_t p) {
    a %= p;
    if (p == 2 || a == 0) return a;
    if (pow_mod(a, (p - 1) / 2, p) != 1) return std::nullopt;
    uint64_t r;
    if (p % 4 == 3) {
        r = pow_mod(a, (p + 1) / 4, p);
    } else {
        // p - 1 = q 2^s, q odd; z a non-square, so c = z^q has order 2^s.
        uint64_t q = p - 1;
        int s = 0;
        for (; (q & 1) == 0; q >>= 1) ++s;
        uint64_t z = 2;
        while (pow_mod(z, (p - 1) / 2, p) != p - 1)
            if (++z == p) throw std::invalid_argument("sqrt_mod: p is not prime");
        uint64_t c = pow_mod(z, q, p);
        uint64_t t = pow_mod(a, q, p);
        r = pow_mod(a, (q + 1) / 2, p);
        // Invariant r^2 = a t, with the order of t a power of two below 2^m.
        for (int m = s; t != 1;) {
            int i = 0;
            for (uint64_t u = t; u != 1; u = mul_mod(u, u, p))
                if (++i >= m) throw std::invalid_argument("sqrt_mod: p is not prime");
            uint64_t b = c;
            for (int j = 0; j < m - i - 1; ++j) b = mul_mod(b, b, p);
            m = i;
            c = mul_mod(b, b, p);
            t = mul_mod(t, c, p);
            r = mul_mod(r, b, p);
        }
    }
    return std::min(r, p - r);
}

}  // namespace quicksynth
