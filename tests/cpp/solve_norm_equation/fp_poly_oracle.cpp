// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "fp_poly_oracle.hpp"

#include "quicksynth/solve_norm_equation/integer_factor.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <vector>

namespace quicksynth::test {

namespace {

constexpr int kMax = FpPoly::kMaxDegree;
// Random trials per split. One succeeds with probability at least 1/2.
constexpr int kSplitAttempts = 64;

uint64_t inv_mod(uint64_t a, uint64_t p) { return pow_mod(a, p - 2, p); }

// A product before reduction, up to degree 2 kMax.
struct Wide {
    std::array<uint64_t, 2 * kMax + 1> c{};
    int degree = -1;
};

template <class Poly>
void trim(Poly& a) {
    while (a.degree >= 0 && a.c[a.degree] == 0) --a.degree;
}

Wide multiply(const FpPoly& a, const FpPoly& b, uint64_t p) {
    Wide w;
    if (a.degree < 0 || b.degree < 0) return w;
    w.degree = a.degree + b.degree;
    for (int i = 0; i <= a.degree; ++i)
        for (int j = 0; j <= b.degree; ++j)
            w.c[i + j] = add_mod(w.c[i + j], mul_mod(a.c[i], b.c[j], p), p);
    trim(w);
    return w;
}

// a mod m, quotient into q when given.
template <class Poly>
FpPoly reduce(Poly a, const FpPoly& m, uint64_t p, FpPoly* q = nullptr) {
    if (m.degree < 0) throw std::domain_error("prime_field: division by the zero polynomial");
    const uint64_t lead_inv = inv_mod(m.c[m.degree], p);
    if (q != nullptr) *q = FpPoly{};
    for (int i = a.degree; i >= m.degree; --i) {
        const uint64_t f = mul_mod(a.c[i], lead_inv, p);
        if (f == 0) continue;
        if (q != nullptr) {
            if (i - m.degree > kMax) throw std::length_error("prime_field: quotient past degree 8");
            q->c[i - m.degree] = f;
            q->degree = std::max(q->degree, i - m.degree);
        }
        for (int j = 0; j <= m.degree; ++j)
            a.c[i - m.degree + j] = sub_mod(a.c[i - m.degree + j], mul_mod(f, m.c[j], p), p);
    }
    FpPoly r;
    for (int i = 0; i < std::min(m.degree, a.degree + 1); ++i) r.c[i] = a.c[i];
    r.degree = std::min(m.degree, a.degree + 1) - 1;
    trim(r);
    return r;
}

FpPoly x_poly() {
    FpPoly x;
    x.c[1] = 1;
    x.degree = 1;
    return x;
}

FpPoly one_poly() {
    FpPoly one;
    one.c[0] = 1;
    one.degree = 0;
    return one;
}

// SplitMix64 [Steele, Lea, Flood, OOPSLA 2014]: a fixed stream for the splits.
struct SplitMix {
    uint64_t state;
    uint64_t next() {
        uint64_t z = (state += 0x9e3779b97f4a7c15u);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9u;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebu;
        return z ^ (z >> 31);
    }
};

// a^((p^d - 1) / 2) mod f, as (a a^p ... a^(p^(d-1)))^((p - 1) / 2) so that
// every exponent fits 64 bits.
FpPoly half_power(const FpPoly& a, int d, const FpPoly& f, uint64_t p) {
    FpPoly conjugate = a, product = a;
    for (int i = 1; i < d; ++i) {
        conjugate = poly_powmod(conjugate, p, f, p);
        product = poly_mulmod(product, conjugate, f, p);
    }
    return poly_powmod(product, (p - 1) / 2, f, p);
}

// A proper monic factor of f, whose irreducible factors all have degree d < deg f.
FpPoly split(const FpPoly& f, int d, uint64_t p, SplitMix& rng) {
    for (int attempt = 0; attempt < kSplitAttempts; ++attempt) {
        FpPoly a;
        for (int i = 0; i < f.degree; ++i) a.c[i] = rng.next() % p;
        a.degree = f.degree - 1;
        trim(a);
        if (a.degree < 1) continue;
        FpPoly g = poly_gcd(a, f, p);
        if (g.degree > 0 && g.degree < f.degree) return g;
        g = poly_gcd(poly_sub(half_power(a, d, f, p), one_poly(), p), f, p);
        if (g.degree > 0 && g.degree < f.degree) return g;
    }
    throw std::runtime_error("prime_field: Cantor-Zassenhaus found no split");
}

void split_into(const FpPoly& f, int d, uint64_t p, SplitMix& rng, std::vector<FpPoly>& out) {
    if (f.degree == d) {
        out.push_back(f);
        return;
    }
    const FpPoly g = split(f, d, p, rng);
    split_into(g, d, p, rng, out);
    split_into(poly_divmod(f, g, p).first, d, p, rng, out);
}

}  // namespace

FpPoly fp_poly(std::initializer_list<uint64_t> coefficients, uint64_t p) {
    if (coefficients.size() > kMax + 1) throw std::length_error("prime_field: degree past 8");
    FpPoly a;
    int i = 0;
    for (uint64_t x : coefficients) a.c[i++] = x % p;
    a.degree = i - 1;
    trim(a);
    return a;
}

FpPoly poly_add(const FpPoly& a, const FpPoly& b, uint64_t p) {
    FpPoly r;
    r.degree = std::max(a.degree, b.degree);
    for (int i = 0; i <= r.degree; ++i) r.c[i] = add_mod(a.c[i], b.c[i], p);
    trim(r);
    return r;
}

FpPoly poly_sub(const FpPoly& a, const FpPoly& b, uint64_t p) {
    FpPoly r;
    r.degree = std::max(a.degree, b.degree);
    for (int i = 0; i <= r.degree; ++i) r.c[i] = sub_mod(a.c[i], b.c[i], p);
    trim(r);
    return r;
}

FpPoly poly_mul(const FpPoly& a, const FpPoly& b, uint64_t p) {
    const Wide w = multiply(a, b, p);
    if (w.degree > kMax) throw std::length_error("prime_field: product past degree 8");
    FpPoly r;
    std::copy_n(w.c.begin(), kMax + 1, r.c.begin());
    r.degree = w.degree;
    return r;
}

FpPoly poly_mulmod(const FpPoly& a, const FpPoly& b, const FpPoly& m, uint64_t p) {
    return reduce(multiply(a, b, p), m, p);
}

FpPoly poly_powmod(const FpPoly& a, uint64_t e, const FpPoly& m, uint64_t p) {
    FpPoly base = reduce(a, m, p), r = reduce(one_poly(), m, p);
    for (; e != 0; e >>= 1) {
        if (e & 1) r = poly_mulmod(r, base, m, p);
        base = poly_mulmod(base, base, m, p);
    }
    return r;
}

FpPoly poly_monic(const FpPoly& a, uint64_t p) {
    if (a.degree < 0) return a;
    const uint64_t inv = inv_mod(a.c[a.degree], p);
    FpPoly r = a;
    for (int i = 0; i <= r.degree; ++i) r.c[i] = mul_mod(r.c[i], inv, p);
    return r;
}

FpPoly poly_gcd(const FpPoly& a, const FpPoly& b, uint64_t p) {
    FpPoly x = a, y = b;
    while (y.degree >= 0) {
        FpPoly r = reduce(x, y, p);
        x = y;
        y = r;
    }
    return poly_monic(x, p);
}

std::vector<std::pair<int, FpPoly>> distinct_degree_factors(const FpPoly& f, uint64_t p) {
    std::vector<std::pair<int, FpPoly>> out;
    const FpPoly x = x_poly();
    FpPoly rest = f;
    // w = x^(p^d) mod rest.
    FpPoly w = reduce(x, rest, p);
    for (int d = 1; 2 * d <= rest.degree; ++d) {
        w = poly_powmod(w, p, rest, p);
        const FpPoly g = poly_gcd(poly_sub(w, x, p), rest, p);
        if (g.degree > 0) {
            out.emplace_back(d, g);
            rest = poly_divmod(rest, g, p).first;
            w = reduce(w, rest, p);
        }
    }
    if (rest.degree > 0) out.emplace_back(rest.degree, rest);
    return out;
}

std::vector<FpPoly> equal_degree_factors(const FpPoly& f, int d, uint64_t p) {
    if (p < 3 || d < 1 || f.degree < d || f.degree % d != 0)
        throw std::invalid_argument("equal_degree_factors: needs p odd and d dividing deg f");
    SplitMix rng{p};
    std::vector<FpPoly> out;
    split_into(f, d, p, rng, out);
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace quicksynth::test
