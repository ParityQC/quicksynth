// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Z[omega] and Z[zeta16] and their real subrings against brute force: products
// of all Galois conjugates for norms, fixed-point radicals of 600 bits for
// signs, and the cyclotomic units for units and square roots.

#include "quicksynth/data_structures/number_field.hpp"

#include <doctest/doctest.h>

#include "oracles.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace quicksynth::test {

namespace {

template <size_t N>
using Ring = RingElement<N, int_mp>;
template <size_t P>
using Real = RealElement<P, int_mp>;

int_mp abs_of(const int_mp& x) { return x.sign() < 0 ? -x : x; }

int_mp isqrt(const int_mp& x) {
    if (x.zero()) return x;
    int_mp y = shift_left(int_mp(1), (x.bits() + 1) / 2 + 1);
    for (;;) {
        const int_mp next = half_of(y + x / y);
        if (compare(next, y) >= 0) return y;
        y = next;
    }
}

int64_t uniform_in(Rng& rng, int64_t bound) { return rng.below(2 * bound + 1) - bound; }

template <size_t N>
Ring<N> random_ring(Rng& rng, int64_t bound) {
    Ring<N> x;
    for (int_mp& c : x) c = int_mp(uniform_in(rng, bound));
    return x;
}

template <size_t P>
Real<P> random_real(Rng& rng, int64_t bound) {
    Real<P> x;
    for (size_t i = 0; i < P; ++i) x[i] = int_mp(uniform_in(rng, bound));
    return x;
}

template <size_t N>
Ring<N> one() {
    return ring_zeta_power<N, int_mp>(0);
}

template <size_t P>
Real<P> real_one() {
    Real<P> r;
    r[0] = int_mp(1);
    return r;
}

template <size_t P>
Real<P> real_neg(const Real<P>& r) {
    return real_sub(Real<P>{}, r);
}

// The product of x over every automorphism, which is N_{K/Q}(x).
template <size_t N>
Ring<N> product_of_conjugates(const Ring<N>& x) {
    Ring<N> p = one<N>();
    for (int j = 1; j < 2 * static_cast<int>(N); j += 2) p = ring_mul(p, galois(x, j));
    return p;
}

// Fixed-point values 2^kS sqrt2 and 2^kS gamma_+-, rounded down.
constexpr int kS = 600;
struct Radicals {
    int_mp unit, sqrt2, gamma_plus, gamma_minus;
    Radicals() {
        unit = shift_left(int_mp(1), kS);
        sqrt2 = isqrt(shift_left(int_mp(2), 2 * kS));
        gamma_plus = isqrt(shift_left(shift_left(int_mp(2), kS) + sqrt2, kS));
        gamma_minus = isqrt(shift_left(shift_left(int_mp(2), kS) - sqrt2, kS));
    }
};
const Radicals& radicals() {
    static const Radicals r;
    return r;
}

// The sign of r at place p from the fixed-point radicals; nullopt if their
// rounding could change it.
template <size_t P>
std::optional<int> oracle_sign(const Real<P>& r, int place) {
    const Radicals& k = radicals();
    const int s = P == 2 ? (place == 0 ? 1 : -1) : std::array<int, 4>{1, -1, -1, 1}[place];
    const int t = std::array<int, 4>{1, 1, -1, -1}[place];
    // value 2^(2 kS), and a bound on its error
    int_mp value = shift_left(r[0] * k.unit + r[1] * k.sqrt2 * s, kS);
    int_mp error = shift_left(abs_of(r[1]) + int_mp(1), kS);
    if constexpr (P == 4) {
        const int_mp& gamma = s > 0 ? k.gamma_plus : k.gamma_minus;
        const int_mp b = r[2] * k.unit + r[3] * k.sqrt2 * s;
        value = value + b * gamma * t;
        error = error + shift_left(abs_of(r[2]) + abs_of(r[3]) * 3 + int_mp(1), kS + 2);
    }
    if (compare(abs_of(value), shift_left(error, 2)) <= 0) {
        if (value.zero() && r == Real<P>{}) return 0;
        return std::nullopt;
    }
    return value.sign();
}

template <size_t P>
void check_signs(const Real<P>& r) {
    bool positive = true;
    for (int p = 0; p < static_cast<int>(P); ++p) {
        const std::optional<int> want = oracle_sign(r, p);
        REQUIRE(want.has_value());
        CHECK(real_sign(r, p) == *want);
        positive = positive && *want > 0;
        // __int128 digits agree, through the int_mp fallback when squares overflow
        bool fits = true;
        for (size_t i = 0; i < P; ++i) fits = fits && r[i].bits() <= 126;
        if (fits) CHECK(real_sign(RealElement<P, __int128>{to_int128(r.c)}, p) == *want);
    }
    CHECK(totally_positive(r) == positive);
}

// The cyclotomic units zeta^((1-a)/2) (1 - zeta^a) / (1 - zeta) of Z[zeta16],
// a = 3, 5, 7, real and with -1 generating the units of Z[g]; 1 + sqrt2 for Z[sqrt2].
std::vector<Real<4>> units_of_z_g() {
    std::vector<Real<4>> out;
    const Ring<8> denominator = ring_sub(one<8>(), ring_zeta_power<8, int_mp>(1));
    for (int a : {3, 5, 7}) {
        const Ring<8> numerator = ring_sub(one<8>(), ring_zeta_power<8, int_mp>(a));
        const std::optional<Ring<8>> q = ring_divide(numerator, denominator);
        REQUIRE(q.has_value());
        const Ring<8> e = ring_times_zeta(*q, (1 - a) / 2);
        REQUIRE(is_real(e));
        out.push_back(to_real(e));
    }
    return out;
}

template <size_t P>
Real<P> real_pow(const Real<P>& x, int n) {
    Real<P> base = x;
    if (n < 0) {
        const std::optional<Real<P>> inverse = real_divide(real_one<P>(), x);
        REQUIRE(inverse.has_value());
        base = *inverse;
        n = -n;
    }
    Real<P> out = real_one<P>();
    for (int i = 0; i < n; ++i) out = real_mul(out, base);
    return out;
}

// A random unit of Z[g] (P = 4) or Z[sqrt2] (P = 2) with exponents in [-e, e].
template <size_t P>
Real<P> random_unit(Rng& rng, int e) {
    Real<P> u = real_one<P>();
    if constexpr (P == 2) {
        u = real_pow(Real<2>{{int_mp(1), int_mp(1)}}, static_cast<int>(uniform_in(rng, e)));
    } else {
        static const std::vector<Real<4>> units = units_of_z_g();
        for (const Real<4>& g : units)
            u = real_mul(u, real_pow(g, static_cast<int>(uniform_in(rng, e))));
    }
    return rng.below(2) != 0 ? real_neg(u) : u;
}

template <size_t N>
void check_ring_axioms(uint64_t seed) {
    Rng rng(seed);
    for (int i = 0; i < 300; ++i) {
        const Ring<N> a = random_ring<N>(rng, 1000), b = random_ring<N>(rng, 1000),
                      c = random_ring<N>(rng, 1000);
        CHECK(ring_mul(ring_mul(a, b), c) == ring_mul(a, ring_mul(b, c)));
        CHECK(ring_mul(a, b) == ring_mul(b, a));
        CHECK(ring_mul(a, ring_add(b, c)) == ring_add(ring_mul(a, b), ring_mul(a, c)));
        CHECK(ring_mul(a, one<N>()) == a);
        CHECK(ring_add(a, ring_neg(a)) == Ring<N>{});
        CHECK(ring_sub(a, b) == ring_add(a, ring_neg(b)));
        const int j = static_cast<int>(uniform_in(rng, 40));
        CHECK(ring_times_zeta(a, j) == ring_mul(a, ring_zeta_power<N, int_mp>(j)));
        // __int128 digits give the same values
        const RingElement<N, __int128> a128 = to_int128(a), b128 = to_int128(b);
        CHECK(to_int_mp(ring_mul(a128, b128)) == ring_mul(a, b));
        CHECK(to_int_mp(ring_sub(a128, b128)) == ring_sub(a, b));

        constexpr size_t P = N / 2;
        const Real<P> r = random_real<P>(rng, 1000), s = random_real<P>(rng, 1000);
        CHECK(to_real(from_real(r)) == r);
        CHECK(from_real(real_mul(r, s)) == ring_mul(from_real(r), from_real(s)));
        CHECK(from_real(real_add(r, s)) == ring_add(from_real(r), from_real(s)));
        CHECK(real_mul(r, s) == real_mul(s, r));
        CHECK(to_int_mp(real_mul(RealElement<P, __int128>{to_int128(r.c)},
                                 RealElement<P, __int128>{to_int128(s.c)})
                            .c) == real_mul(r, s).c);
    }
    CHECK(ring_zeta_power<N, int_mp>(static_cast<int>(N)) == ring_neg(one<N>()));
    CHECK(ring_zeta_power<N, int_mp>(2 * static_cast<int>(N)) == one<N>());
    CHECK(ring_zeta_power<N, int_mp>(-1) == conjugate(ring_zeta_power<N, int_mp>(1)));
    CHECK_FALSE(is_real(ring_zeta_power<N, int_mp>(1)));
    CHECK_THROWS_AS(to_real(ring_zeta_power<N, int_mp>(1)), std::invalid_argument);
}

template <size_t N>
void check_galois(uint64_t seed) {
    Rng rng(seed);
    const int period = 2 * static_cast<int>(N);
    for (int i = 0; i < 100; ++i) {
        const Ring<N> a = random_ring<N>(rng, 1000), b = random_ring<N>(rng, 1000);
        for (int j = 1; j < period; j += 2) {
            CHECK(galois(ring_mul(a, b), j) == ring_mul(galois(a, j), galois(b, j)));
            CHECK(galois(ring_add(a, b), j) == ring_add(galois(a, j), galois(b, j)));
            for (int k = 1; k < period; k += 2)
                CHECK(galois(galois(a, j), k) == galois(a, (j * k) % period));
            // embedding after sigma_j is the embedding at j
            const std::complex<double> x = embed(galois(a, j), 1), y = embed(a, j);
            CHECK(std::abs(x - y) <= 1e-9 * (1 + std::abs(y)));
        }
        CHECK(galois(a, 1) == a);
        CHECK(galois(a, period - 1) == conjugate(a));
        CHECK(conjugate(conjugate(a)) == a);
    }
    CHECK_THROWS_AS(galois(one<N>(), 2), std::invalid_argument);
}

template <size_t N>
void check_norms(uint64_t seed) {
    Rng rng(seed);
    constexpr size_t P = N / 2;
    for (int i = 0; i < 100; ++i) {
        // coefficients to 2^40: norms of 330 bits for N = 8
        const int64_t bound = i % 2 == 0 ? 1000 : int64_t{1} << 40;
        const Ring<N> a = random_ring<N>(rng, bound), b = random_ring<N>(rng, bound);
        CHECK(relative_norm(ring_mul(a, b)) == real_mul(relative_norm(a), relative_norm(b)));
        CHECK(ring_norm(ring_mul(a, b)) == ring_norm(a) * ring_norm(b));
        CHECK(ring_norm(a).sign() >= 0);

        // the tower norm is the brute-force product of all conjugates
        const Ring<N> brute = product_of_conjugates(a);
        Ring<N> constant{};
        constant[0] = ring_norm(a);
        CHECK(brute == constant);

        const Real<P> r = random_real<P>(rng, bound), s = random_real<P>(rng, bound);
        CHECK(real_norm(real_mul(r, s)) == real_norm(r) * real_norm(s));
        // over K+: the conjugates zeta -> zeta^(2p+1), p < P, one per place
        Ring<N> real_brute = one<N>();
        for (int p = 0; p < static_cast<int>(P); ++p)
            real_brute = ring_mul(real_brute, galois(from_real(r), 2 * p + 1));
        Ring<N> real_constant{};
        real_constant[0] = real_norm(r);
        CHECK(real_brute == real_constant);
        // __int128 digits widen to the same norm
        CHECK(ring_norm(RingElement<N, __int128>(to_int128(a))) == ring_norm(a));
        CHECK(real_norm(RealElement<P, __int128>{to_int128(r.c)}) == real_norm(r));
    }
}

template <size_t N>
void check_division(uint64_t seed) {
    Rng rng(seed);
    constexpr size_t P = N / 2;
    for (int i = 0; i < 200; ++i) {
        const Ring<N> a = random_ring<N>(rng, 1000), b = random_ring<N>(rng, 50);
        if (b == Ring<N>{}) continue;
        const Ring<N> ab = ring_mul(a, b);
        const std::optional<Ring<N>> q = ring_divide(ab, b);
        REQUIRE(q.has_value());
        CHECK(*q == a);
        // __int128 digits
        const auto q128 = ring_divide(RingElement<N, __int128>(to_int128(ab)),
                                      RingElement<N, __int128>(to_int128(b)));
        REQUIRE(q128.has_value());
        CHECK(to_int_mp(*q128) == a);
        // off by one is not divisible unless b is a unit
        if (!(ring_norm(b) == 1)) CHECK_FALSE(ring_divide(ring_add(ab, one<N>()), b).has_value());

        const Real<P> r = random_real<P>(rng, 1000), s = random_real<P>(rng, 50);
        if (s == Real<P>{}) continue;
        const std::optional<Real<P>> rq = real_divide(real_mul(r, s), s);
        REQUIRE(rq.has_value());
        CHECK(*rq == r);
        if (!(abs_of(real_norm(s)) == 1))
            CHECK_FALSE(real_divide(real_add(real_mul(r, s), real_one<P>()), s).has_value());
    }
    CHECK_THROWS_AS(ring_divide(one<N>(), Ring<N>{}), std::domain_error);
    CHECK_THROWS_AS(real_divide(real_one<P>(), Real<P>{}), std::domain_error);
    // __int128 arithmetic refuses to wrap
    RingElement<N, __int128> big{};
    big[0] = static_cast<__int128>(1) << 100;
    CHECK_THROWS_AS(ring_mul(big, big), std::overflow_error);
}

}  // namespace

TEST_CASE("Z[omega] and Z[sqrt2] satisfy the ring axioms") { check_ring_axioms<4>(21); }
TEST_CASE("Z[zeta16] and Z[g] satisfy the ring axioms") { check_ring_axioms<8>(22); }

TEST_CASE("the Galois automorphisms of Z[omega] are homomorphisms") { check_galois<4>(23); }
TEST_CASE("the Galois automorphisms of Z[zeta16] are homomorphisms") { check_galois<8>(24); }

TEST_CASE("the tower norms of Z[omega] are multiplicative and the conjugate product") {
    check_norms<4>(25);
}
TEST_CASE("the tower norms of Z[zeta16] are multiplicative and the conjugate product") {
    check_norms<8>(26);
}

TEST_CASE("exact division in Z[omega] and Z[sqrt2] checks the remainder") { check_division<4>(27); }
TEST_CASE("exact division in Z[zeta16] and Z[g] checks the remainder") { check_division<8>(28); }

TEST_CASE("real signs are exact on random values") {
    Rng rng(29);
    for (int i = 0; i < 500; ++i) {
        check_signs(random_real<2>(rng, int64_t{1} << 30));
        check_signs(random_real<4>(rng, int64_t{1} << 30));
    }
    check_signs(Real<2>{});
    check_signs(Real<4>{});
}

// Units with a conjugate near zero: (sqrt2 - 1)^n, Pell's p - q sqrt2, and
// products of the cyclotomic units of Z[g], whose place values reach 2^-100
// while the coefficients reach 2^100 and their squares pass 128 bits.
TEST_CASE("real signs are exact on values near zero at a place") {
    for (int n = 1; n <= 80; ++n) {
        const Real<2> small = real_pow(Real<2>{{int_mp(-1), int_mp(1)}}, n);
        check_signs(small);
        check_signs(real_neg(small));
        check_signs(real_add(small, real_one<2>()));
    }
    Rng rng(30);
    const std::vector<Real<4>> units = units_of_z_g();
    for (const Real<4>& e : units) CHECK(abs_of(real_norm(e)) == int_mp(1));
    for (int i = 0; i < 300; ++i) {
        const Real<4> u = random_unit<4>(rng, 12);
        check_signs(u);
        // a unit plus an integer: the small conjugate moves off zero
        check_signs(real_add(u, real_one<4>()));
    }
    for (int n = 1; n <= 30; ++n)
        for (const Real<4>& e : units) check_signs(real_pow(e, n));
}

TEST_CASE("embeddings agree with each other and with the exact signs") {
    Rng rng(31);
    for (int i = 0; i < 300; ++i) {
        const Ring<8> x = random_ring<8>(rng, 1 << 20);
        const Ring<4> w = random_ring<4>(rng, 1 << 20);
        for (int k = 1; k < 16; k += 2) {
            const ComplexDD v = embed_dd(x, k);
            const std::complex<double> d = embed(x, k);
            CHECK(std::abs(std::complex<double>(to_double(v.re), to_double(v.im)) - d) <= 1e-7);
        }
        // |x_k|^2 = (x x*)_k to double-double precision
        for (int p = 0; p < 4; ++p) {
            const ComplexDD v = embed_dd(x, 2 * p + 1);
            const DD square = v.re * v.re + v.im * v.im;
            const DD norm = real_embed_dd(relative_norm(x), p);
            CHECK(std::abs(to_double(square - norm)) <= 1e-26 * (to_double(norm) + 1e12));
        }
        for (int p = 0; p < 2; ++p) {
            const ComplexDD v = embed_dd(w, 2 * p + 1);
            const DD norm = real_embed_dd(relative_norm(w), p);
            CHECK(std::abs(to_double(v.re * v.re + v.im * v.im - norm)) <=
                  1e-26 * (to_double(norm) + 1e12));
        }
        const Real<4> r = random_real<4>(rng, 1 << 20);
        for (int p = 0; p < 4; ++p) {
            const double z = embed(from_real(r), 2 * p + 1).real();
            CHECK(std::abs(to_double(real_embed_dd(r, p)) - z) <= 1e-7);
        }
    }
    // a unit whose place values span about 2^-25 .. 2^25: the double-double
    // values carry the exact signs
    const Real<4> u = real_pow(units_of_z_g()[0], 12);
    for (int p = 0; p < 4; ++p) CHECK(dd_sign(real_embed_dd(u, p)) == real_sign(u, p));
}

TEST_CASE("unit_sqrt recovers square roots of products of the cyclotomic units") {
    Rng rng(32);
    for (int i = 0; i < 300; ++i) {
        const Real<4> eta = random_unit<4>(rng, 6);
        const Real<4> u = real_mul(eta, eta);
        CHECK(totally_positive(u));
        const std::optional<Real<4>> root = unit_sqrt(u);
        REQUIRE(root.has_value());
        CHECK((*root == eta || *root == real_neg(eta)));
        CHECK(real_sign(*root, 0) == 1);
        // __int128 digits
        const auto root128 = unit_sqrt(RealElement<4, __int128>{to_int128(u.c)});
        REQUIRE(root128.has_value());
        CHECK(to_int_mp(root128->c) == root->c);

        const Real<2> lambda_power = random_unit<2>(rng, 30);
        const std::optional<Real<2>> root2 = unit_sqrt(real_mul(lambda_power, lambda_power));
        REQUIRE(root2.has_value());
        CHECK((*root2 == lambda_power || *root2 == real_neg(lambda_power)));
    }
    // Weber: a unit is a square exactly when it is totally positive
    for (int i = 0; i < 300; ++i) {
        const Real<4> u = random_unit<4>(rng, 6);
        CHECK(unit_sqrt(u).has_value() == totally_positive(u));
        const Real<2> v = random_unit<2>(rng, 30);
        CHECK(unit_sqrt(v).has_value() == totally_positive(v));
    }
    // squares that are not units, and totally positive non-squares
    CHECK(unit_sqrt(Real<4>{{int_mp(9), int_mp(), int_mp(), int_mp()}}) ==
          Real<4>{{int_mp(3), int_mp(), int_mp(), int_mp()}});
    CHECK(unit_sqrt(Real<4>{{int_mp(2), int_mp(1), int_mp(), int_mp()}}) ==
          Real<4>{{int_mp(), int_mp(), int_mp(1), int_mp()}});  // sqrt(2 + sqrt2) = g
    CHECK_FALSE(unit_sqrt(Real<2>{{int_mp(2), int_mp(1)}}).has_value());
    CHECK_FALSE(unit_sqrt(Real<4>{{int_mp(3), int_mp(), int_mp(), int_mp()}}).has_value());
}

}  // namespace quicksynth::test
