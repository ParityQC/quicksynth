// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/solve_norm_equation/prime_field.hpp"

#include "quicksynth/solve_norm_equation/integer_factor.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace quicksynth {

namespace {

uint64_t inv_mod(uint64_t a, uint64_t p) { return pow_mod(a, p - 2, p); }

void trim(FpPoly& a) {
    while (a.degree >= 0 && a.c[a.degree] == 0) --a.degree;
}

// a mod m, quotient into q.
FpPoly reduce(FpPoly a, const FpPoly& m, uint64_t p, FpPoly& q) {
    if (m.degree < 0) throw std::domain_error("prime_field: division by the zero polynomial");
    const uint64_t lead_inv = inv_mod(m.c[m.degree], p);
    q = FpPoly{};
    for (int i = a.degree; i >= m.degree; --i) {
        const uint64_t f = mul_mod(a.c[i], lead_inv, p);
        if (f == 0) continue;
        q.c[i - m.degree] = f;
        q.degree = std::max(q.degree, i - m.degree);
        for (int j = 0; j <= m.degree; ++j)
            a.c[i - m.degree + j] = sub_mod(a.c[i - m.degree + j], mul_mod(f, m.c[j], p), p);
    }
    FpPoly r;
    for (int i = 0; i < std::min(m.degree, a.degree + 1); ++i) r.c[i] = a.c[i];
    r.degree = std::min(m.degree, a.degree + 1) - 1;
    trim(r);
    return r;
}

// A primitive m-th root of unity mod p, m a power of two dividing p - 1: z^((p - 1) / m)
// for a non-square z has order m.
uint64_t root_of_unity(uint64_t m, uint64_t p) {
    uint64_t z = 2;
    while (pow_mod(z, (p - 1) / 2, p) != p - 1)
        if (++z == p) throw std::invalid_argument("prime_field: p is not prime");
    return pow_mod(z, (p - 1) / m, p);
}

// A square root of a that the residue class of p guarantees.
uint64_t square_root(uint64_t a, uint64_t p) {
    const std::optional<uint64_t> r = sqrt_mod(a, p);
    if (!r) throw std::invalid_argument("prime_field: p is not prime");
    return *r;
}

}  // namespace

bool FpPoly::operator<(const FpPoly& o) const {
    if (degree != o.degree) return degree < o.degree;
    for (int i = degree; i >= 0; --i)
        if (c[i] != o.c[i]) return c[i] < o.c[i];
    return false;
}

std::pair<FpPoly, FpPoly> poly_divmod(const FpPoly& a, const FpPoly& b, uint64_t p) {
    FpPoly q;
    FpPoly r = reduce(a, b, p, q);
    return {q, r};
}

int cyclotomic_residue_degree(uint64_t p, int n) {
    if ((p & 1) == 0 || n < 1 || (n & (n - 1)) != 0)
        throw std::invalid_argument("cyclotomic_residue_degree: needs p odd and n a power of two");
    const uint64_t m = 2 * static_cast<uint64_t>(n);
    const uint64_t q = p % m;
    int f = 1;
    for (uint64_t t = q; t != 1; t = t * q % m) ++f;
    return f;
}

std::vector<FpPoly> cyclotomic_factors(int n, uint64_t p) {
    if (n != 4 && n != 8) throw std::invalid_argument("cyclotomic_factors: n is 4 or 8");
    if ((p & 1) == 0) throw std::invalid_argument("cyclotomic_factors: p is odd");
    std::vector<FpPoly> out;
    // Monic x^k + a x^j + b, k > j.
    const auto trinomial = [p](int k, int j, uint64_t a, uint64_t b) {
        FpPoly h;
        h.c[k] = 1;
        h.c[j] = a % p;
        h.c[0] = add_mod(h.c[0], b % p, p);
        h.degree = k;
        return h;
    };
    const uint64_t minus_one = p - 1;
    if (p % (2 * n) == 1) {
        // Split completely: x - r^j, r of order 2n, j odd.
        const uint64_t r = root_of_unity(2 * n, p), r2 = mul_mod(r, r, p);
        for (uint64_t x = r, j = 0; j < static_cast<uint64_t>(n); ++j, x = mul_mod(x, r2, p))
            out.push_back(trinomial(1, 0, 0, p - x));
    } else if (n == 8 && p % 16 == 9) {
        // x^2 - w, w^4 = -1.
        const uint64_t w = root_of_unity(8, p), w2 = mul_mod(w, w, p);
        for (uint64_t x = w, j = 0; j < 4; ++j, x = mul_mod(x, w2, p))
            out.push_back(trinomial(2, 0, 0, p - x));
    } else if (p % 8 == 5) {
        // x^(n/2) -+ i, i^2 = -1.
        const uint64_t i = root_of_unity(4, p);
        out = {trinomial(n / 2, 0, 0, i), trinomial(n / 2, 0, 0, p - i)};
    } else if (p % 8 == 3) {
        // y^2 -+ a y - 1, a^2 = -2, y = x^(n/4).
        const uint64_t a = square_root(p - 2, p);
        out = {trinomial(n / 2, n / 4, a, minus_one), trinomial(n / 2, n / 4, p - a, minus_one)};
    } else if (n == 4) {
        // p = 7 mod 8: x^2 -+ s x + 1, s^2 = 2.
        const uint64_t s = square_root(2, p);
        out = {trinomial(2, 1, s, 1), trinomial(2, 1, p - s, 1)};
    } else {
        // p = 7 or 15 mod 16: x^8 + 1 = (x^4 - s x^2 + 1)(x^4 + s x^2 + 1), s^2 = 2,
        // and x^4 -+ s x^2 + 1 = (x^2 + b x + c)(x^2 - b x + c) with c = 1 for
        // p = 15 (b^2 = 2 +- s) and c = -1 for p = 7 (b^2 = -2 +- s) mod 16.
        const uint64_t s = square_root(2, p);
        const uint64_t c = p % 16 == 15 ? 1 : minus_one;
        const uint64_t two_c = add_mod(c, c, p);
        for (uint64_t t : {s, p - s}) {
            const uint64_t b = square_root(add_mod(two_c, t, p), p);
            out.push_back(trinomial(2, 1, b, c));
            out.push_back(trinomial(2, 1, p - b, c));
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace quicksynth
