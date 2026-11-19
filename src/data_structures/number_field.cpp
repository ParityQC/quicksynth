// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/number_field.hpp"

#include <array>
#include <complex>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace quicksynth {

namespace {

// Digit operations: checked in __int128, exact in int_mp.
__int128 dadd(__int128 a, __int128 b) { return checked_add(a, b); }
__int128 dsub(__int128 a, __int128 b) { return checked_sub(a, b); }
__int128 dmul(__int128 a, __int128 b) { return checked_mul(a, b); }
int dsign(__int128 a) { return (a > 0) - (a < 0); }
double ddouble(__int128 a) { return static_cast<double>(a); }

int_mp dadd(const int_mp& a, const int_mp& b) { return a + b; }
int_mp dsub(const int_mp& a, const int_mp& b) { return a - b; }
int_mp dmul(const int_mp& a, const int_mp& b) { return a * b; }
int dsign(const int_mp& a) { return a.sign(); }
double ddouble(const int_mp& a) { return a.to_double(0); }

template <class D>
D dneg(const D& a) {
    return dsub(D(0), a);
}

template <class D>
D dtimes(const D& a, int s) {
    return dmul(a, D(s));
}

template <class D, size_t M>
std::array<int_mp, M> widened(const std::array<D, M>& x) {
    if constexpr (std::is_same_v<D, int_mp>)
        return x;
    else
        return to_int_mp(x);
}

template <class D, size_t M>
std::array<D, M> narrowed(const std::array<int_mp, M>& x) {
    if constexpr (std::is_same_v<D, int_mp>)
        return x;
    else
        return to_int128(x);
}

template <size_t P, class D>
RealElement<P, int_mp> widened(const RealElement<P, D>& r) {
    return {widened(r.c)};
}

template <class D, size_t P>
RealElement<P, D> narrowed(const RealElement<P, int_mp>& r) {
    return {narrowed<D>(r.c)};
}

// Z[sqrt2] as (a, b) = a + b sqrt2.
template <class D>
struct Z2 {
    D a, b;
};

template <class D>
Z2<D> z2_add(const Z2<D>& x, const Z2<D>& y) {
    return {dadd(x.a, y.a), dadd(x.b, y.b)};
}

template <class D>
Z2<D> z2_sub(const Z2<D>& x, const Z2<D>& y) {
    return {dsub(x.a, y.a), dsub(x.b, y.b)};
}

template <class D>
Z2<D> z2_mul(const Z2<D>& x, const Z2<D>& y) {
    return {dadd(dmul(x.a, y.a), dtimes(dmul(x.b, y.b), 2)), dadd(dmul(x.a, y.b), dmul(x.b, y.a))};
}

// x (2 + sqrt2) = g^2 x.
template <class D>
Z2<D> z2_times_g2(const Z2<D>& x) {
    return {dtimes(dadd(x.a, x.b), 2), dadd(x.a, dtimes(x.b, 2))};
}

// Sign of a + s b sqrt2, s = +-1.
template <class D>
int z2_sign(const Z2<D>& x, int s) {
    const int sa = dsign(x.a), sb = s * dsign(x.b);
    if (sa == 0) return sb;
    if (sb == 0 || sa == sb) return sa;
    const int c = dsign(dsub(dmul(x.a, x.a), dtimes(dmul(x.b, x.b), 2)));
    return c > 0 ? sa : sb;  // a^2 = 2 b^2 only at 0
}

// The tower parts A, B of r = A + B g (P = 4).
template <class D>
Z2<D> part_a(const RealElement<4, D>& r) {
    return {r[0], r[1]};
}
template <class D>
Z2<D> part_b(const RealElement<4, D>& r) {
    return {r[2], r[3]};
}

// N_{K+/Q(sqrt2)}(A + B g) = (A + B g)(A - B g) = A^2 - B^2 (2 + sqrt2).
template <class D>
Z2<D> relative_tower_norm(const RealElement<4, D>& r) {
    const Z2<D> a = part_a(r), b = part_b(r);
    return z2_sub(z2_mul(a, a), z2_times_g2(z2_mul(b, b)));
}

// Signs of sqrt2 and g at place p (zeta -> zeta^(2p+1)): 2 cos((2p+1) pi/4)
// and 2 cos((2p+1) pi/8).
constexpr int kSqrt2Sign[4] = {1, -1, -1, 1};
constexpr int kGSign[4] = {1, 1, -1, -1};

template <size_t P>
void check_place(int place) {
    if (place < 0 || place >= static_cast<int>(P))
        throw std::invalid_argument("number_field: no such place");
}

template <size_t P, class W>
int sign_in(const RealElement<P, W>& r, int place) {
    const int s = P == 2 ? (place == 0 ? 1 : -1) : kSqrt2Sign[place];
    if constexpr (P == 2) {
        return z2_sign(Z2<W>{r[0], r[1]}, s);
    } else {
        const int sa = z2_sign(part_a(r), s);
        const int sb = kGSign[place] * z2_sign(part_b(r), s);
        if (sa == 0) return sb;
        if (sb == 0 || sa == sb) return sa;
        // |A| against |B gamma|: A^2 - B^2 gamma^2 is the relative norm at s
        return z2_sign(relative_tower_norm(r), s) > 0 ? sa : sb;
    }
}

// c with r c = real_norm(r): the product of r's other conjugates.
template <size_t P>
RealElement<P, int_mp> cofactor(const RealElement<P, int_mp>& r) {
    if constexpr (P == 2) {
        return {{r[0], -r[1]}};
    } else {
        const Z2<int_mp> n = relative_tower_norm(r);
        const RealElement<4, int_mp> conj_g{{r[0], r[1], -r[2], -r[3]}};
        const RealElement<4, int_mp> conj_sqrt2{{n.a, -n.b, int_mp(), int_mp()}};
        return real_mul(conj_g, conj_sqrt2);
    }
}

// q = x / n coefficientwise if n divides every coefficient.
template <size_t M>
bool divide_exactly(std::array<int_mp, M>& x, const int_mp& n) {
    for (int_mp& c : x) {
        int_mp q = c / n;
        if (q * n != c) return false;
        c = std::move(q);
    }
    return true;
}

template <class D>
DD to_dd(const D& x) {
    return dd(x);
}

// sqrt2 and gamma = 2 cos(pi/8), 2 cos(3 pi/8) in double-double.
DD sqrt2_value() { return sqrt2_dd(); }
DD gamma_value(int s) { return dd_ldexp(cos_pi16(s > 0 ? 2 : 6), 1); }

}  // namespace

void throw_int128_overflow() {
    throw std::overflow_error("number_field: an __int128 value passes 127 bits");
}

template <size_t N, class D>
RingElement<N, D> ring_zeta_power(int j) {
    RingElement<N, D> out{};
    const int period = 2 * static_cast<int>(N);
    j = ((j % period) + period) % period;
    if (j < static_cast<int>(N))
        out[j] = D(1);
    else
        out[j - N] = D(-1);
    return out;
}

template <size_t N, class D>
RingElement<N, D> ring_add(const RingElement<N, D>& a, const RingElement<N, D>& b) {
    RingElement<N, D> out;
    for (size_t i = 0; i < N; ++i) out[i] = dadd(a[i], b[i]);
    return out;
}

template <size_t N, class D>
RingElement<N, D> ring_sub(const RingElement<N, D>& a, const RingElement<N, D>& b) {
    RingElement<N, D> out;
    for (size_t i = 0; i < N; ++i) out[i] = dsub(a[i], b[i]);
    return out;
}

template <size_t N, class D>
RingElement<N, D> ring_neg(const RingElement<N, D>& a) {
    RingElement<N, D> out;
    for (size_t i = 0; i < N; ++i) out[i] = dneg(a[i]);
    return out;
}

template <size_t N, class D>
RingElement<N, D> ring_mul(const RingElement<N, D>& a, const RingElement<N, D>& b) {
    RingElement<N, D> out{};
    for (size_t i = 0; i < N; ++i) {
        if (dsign(a[i]) == 0) continue;
        for (size_t j = 0; j < N; ++j) {
            if (dsign(b[j]) == 0) continue;
            const D p = dmul(a[i], b[j]);
            if (i + j < N)
                out[i + j] = dadd(out[i + j], p);
            else
                out[i + j - N] = dsub(out[i + j - N], p);  // zeta^N = -1
        }
    }
    return out;
}

template <size_t N, class D>
RingElement<N, D> ring_times_zeta(const RingElement<N, D>& x, int j) {
    const int period = 2 * static_cast<int>(N);
    j = ((j % period) + period) % period;
    RingElement<N, D> out;
    for (size_t i = 0; i < N; ++i) {
        const size_t k = (i + static_cast<size_t>(j)) % (2 * N);
        if (k < N)
            out[k] = x[i];
        else
            out[k - N] = dneg(x[i]);
    }
    return out;
}

template <size_t N, class D>
RingElement<N, D> galois(const RingElement<N, D>& x, int j) {
    const int period = 2 * static_cast<int>(N);
    j = ((j % period) + period) % period;
    if (j % 2 == 0) throw std::invalid_argument("number_field: galois needs an odd exponent");
    RingElement<N, D> out{};
    for (size_t i = 0; i < N; ++i) {
        const size_t k = (i * static_cast<size_t>(j)) % (2 * N);
        if (k < N)
            out[k] = dadd(out[k], x[i]);
        else
            out[k - N] = dsub(out[k - N], x[i]);
    }
    return out;
}

template <size_t N, class D>
RingElement<N, D> conjugate(const RingElement<N, D>& x) {
    // zeta^-i = -zeta^(N-i)
    RingElement<N, D> out{};
    out[0] = x[0];
    for (size_t i = 1; i < N; ++i) out[N - i] = dneg(x[i]);
    return out;
}

template <size_t N, class D>
bool is_real(const RingElement<N, D>& x) {
    return conjugate(x) == x;
}

template <size_t N, class D>
RealElement<N / 2, D> to_real(const RingElement<N, D>& x) {
    if (!is_real(x)) throw std::invalid_argument("number_field: the element is not real");
    if constexpr (N == 4) {
        // c0 + c1 (omega - omega^3) = c0 + c1 sqrt2
        return {{x[0], x[1]}};
    } else {
        // zeta^2 - zeta^6 = sqrt2, zeta - zeta^7 = g, zeta^3 - zeta^5 = g (sqrt2 - 1)
        return {{x[0], x[2], dsub(x[1], x[3]), x[3]}};
    }
}

template <size_t P, class D>
RingElement<2 * P, D> from_real(const RealElement<P, D>& r) {
    RingElement<2 * P, D> out{};
    if constexpr (P == 2) {
        out[0] = r[0];
        out[1] = r[1];
        out[3] = dneg(r[1]);
    } else {
        // sqrt2 g = zeta + zeta^3 - zeta^5 - zeta^7
        out[0] = r[0];
        out[2] = r[1];
        out[6] = dneg(r[1]);
        out[1] = dadd(r[2], r[3]);
        out[7] = dneg(out[1]);
        out[3] = r[3];
        out[5] = dneg(r[3]);
    }
    return out;
}

template <size_t N, class D>
RealElement<N / 2, D> relative_norm(const RingElement<N, D>& x) {
    return to_real(ring_mul(x, conjugate(x)));
}

template <size_t N, class D>
int_mp ring_norm(const RingElement<N, D>& x) {
    return real_norm(relative_norm(widened(x)));
}

template <size_t P, class D>
RealElement<P, D> real_add(const RealElement<P, D>& a, const RealElement<P, D>& b) {
    RealElement<P, D> out;
    for (size_t i = 0; i < P; ++i) out[i] = dadd(a[i], b[i]);
    return out;
}

template <size_t P, class D>
RealElement<P, D> real_sub(const RealElement<P, D>& a, const RealElement<P, D>& b) {
    RealElement<P, D> out;
    for (size_t i = 0; i < P; ++i) out[i] = dsub(a[i], b[i]);
    return out;
}

template <size_t P, class D>
RealElement<P, D> real_mul(const RealElement<P, D>& a, const RealElement<P, D>& b) {
    if constexpr (P == 2) {
        const Z2<D> p = z2_mul(Z2<D>{a[0], a[1]}, Z2<D>{b[0], b[1]});
        return {{p.a, p.b}};
    } else {
        // (A + B g)(C + E g) = (AC + BE g^2) + (AE + BC) g
        const Z2<D> A = part_a(a), B = part_b(a), C = part_a(b), E = part_b(b);
        const Z2<D> first = z2_add(z2_mul(A, C), z2_times_g2(z2_mul(B, E)));
        const Z2<D> second = z2_add(z2_mul(A, E), z2_mul(B, C));
        return {{first.a, first.b, second.a, second.b}};
    }
}

template <size_t P, class D>
int_mp real_norm(const RealElement<P, D>& r) {
    const RealElement<P, int_mp> w = widened(r);
    Z2<int_mp> z;
    if constexpr (P == 2)
        z = {w[0], w[1]};
    else
        z = relative_tower_norm(w);
    return z.a * z.a - z.b * z.b * 2;
}

template <size_t P, class D>
int real_sign(const RealElement<P, D>& r, int place) {
    check_place<P>(place);
    if constexpr (std::is_same_v<D, int_mp>) {
        return sign_in(r, place);
    } else {
        // Squares pass 127 bits only for digits past about 2^62.
        try {
            return sign_in(r, place);
        } catch (const std::overflow_error&) {
            return sign_in(widened(r), place);
        }
    }
}

template <size_t P, class D>
bool totally_positive(const RealElement<P, D>& r) {
    for (size_t p = 0; p < P; ++p)
        if (real_sign(r, static_cast<int>(p)) <= 0) return false;
    return true;
}

template <size_t N, class D>
std::optional<RingElement<N, D>> ring_divide(const RingElement<N, D>& a,
                                             const RingElement<N, D>& b) {
    const RingElement<N, int_mp> wb = widened(b);
    const RealElement<N / 2, int_mp> r = relative_norm(wb);
    const int_mp n = real_norm(r);
    if (n.zero()) throw std::domain_error("number_field: division by zero");
    RingElement<N, int_mp> q =
        ring_mul(ring_mul(widened(a), conjugate(wb)), from_real(cofactor(r)));
    if (!divide_exactly(q, n)) return std::nullopt;
    return narrowed<D>(q);
}

template <size_t P, class D>
std::optional<RealElement<P, D>> real_divide(const RealElement<P, D>& a,
                                             const RealElement<P, D>& b) {
    const RealElement<P, int_mp> wb = widened(b);
    const int_mp n = real_norm(wb);
    if (n.zero()) throw std::domain_error("number_field: division by zero");
    RealElement<P, int_mp> q = real_mul(widened(a), cofactor(wb));
    if (!divide_exactly(q.c, n)) return std::nullopt;
    return narrowed<D>(q);
}

template <size_t N, class D>
ComplexDD embed_dd(const RingElement<N, D>& x, int k) {
    // zeta^(i k) = exp(i pi (16 i k / N) / 16)
    ComplexDD out{};
    for (size_t i = 0; i < N; ++i) {
        if (dsign(x[i]) == 0) continue;
        const int m = static_cast<int>(i) * k * (16 / static_cast<int>(N));
        const DD c = to_dd(x[i]);
        out.re = out.re + c * cos_pi16(m);
        out.im = out.im + c * sin_pi16(m);
    }
    return out;
}

template <size_t N, class D>
std::complex<double> embed(const RingElement<N, D>& x, int k) {
    std::complex<double> out = 0;
    for (size_t i = 0; i < N; ++i) {
        if (dsign(x[i]) == 0) continue;
        const int m = static_cast<int>(i) * k * (16 / static_cast<int>(N));
        out += ddouble(x[i]) * std::complex<double>(cos_pi16(m).hi, sin_pi16(m).hi);
    }
    return out;
}

template <size_t P, class D>
DD real_embed_dd(const RealElement<P, D>& r, int place) {
    check_place<P>(place);
    const int s = P == 2 ? (place == 0 ? 1 : -1) : kSqrt2Sign[place];
    const DD root2 = sqrt2_value() * static_cast<double>(s);
    const DD a = to_dd(r[0]) + to_dd(r[1]) * root2;
    if constexpr (P == 2) {
        return a;
    } else {
        const DD b = to_dd(r[2]) + to_dd(r[3]) * root2;
        return a + b * gamma_value(s) * static_cast<double>(kGSign[place]);
    }
}

template <size_t P, class D>
std::optional<RealElement<P, D>> unit_sqrt(const RealElement<P, D>& u) {
    if (!totally_positive(u)) return std::nullopt;
    std::array<DD, P> root;
    for (size_t p = 0; p < P; ++p) root[p] = dd_sqrt(real_embed_dd(u, static_cast<int>(p)));
    const RealElement<P, int_mp> want = widened(u);
    const DD two_root2 = dd_ldexp(sqrt2_value(), 1);
    // One sign per place but place 0's: eta and -eta give the same square.
    for (unsigned pattern = 0; pattern < (1u << (P - 1)); ++pattern) {
        std::array<DD, P> v = root;
        for (size_t p = 1; p < P; ++p)
            if ((pattern >> (p - 1)) & 1u) v[p] = -v[p];
        // Invert the place map: places with sqrt2 -> +sqrt2 and -sqrt2 give the
        // Q(sqrt2) parts at each sign, whose sum and difference give a and b.
        RealElement<P, int_mp> eta;
        if constexpr (P == 2) {
            eta[0] = nearest_integer(dd_ldexp(v[0] + v[1], -1));
            eta[1] = nearest_integer((v[0] - v[1]) / two_root2);
        } else {
            // places 0, 3 have sqrt2 -> +sqrt2, g -> +-gamma+; places 1, 2 the minus sign
            const DD a_plus = dd_ldexp(v[0] + v[3], -1);
            const DD b_plus = (v[0] - v[3]) / dd_ldexp(gamma_value(1), 1);
            const DD a_minus = dd_ldexp(v[1] + v[2], -1);
            const DD b_minus = (v[1] - v[2]) / dd_ldexp(gamma_value(-1), 1);
            eta[0] = nearest_integer(dd_ldexp(a_plus + a_minus, -1));
            eta[1] = nearest_integer((a_plus - a_minus) / two_root2);
            eta[2] = nearest_integer(dd_ldexp(b_plus + b_minus, -1));
            eta[3] = nearest_integer((b_plus - b_minus) / two_root2);
        }
        if (real_mul(eta, eta) != want) continue;
        // A place value near zero does not move the rounding, so the pattern
        // that matched may have either sign there; fix the sign exactly.
        if (real_sign(eta, 0) < 0)
            for (int_mp& c : eta.c) c = -c;
        return narrowed<D>(eta);
    }
    return std::nullopt;
}

#define QUICKSYNTH_RING(N, D)                                                                \
    template RingElement<N, D> ring_zeta_power<N, D>(int);                                   \
    template RingElement<N, D> ring_add(const RingElement<N, D>&, const RingElement<N, D>&); \
    template RingElement<N, D> ring_sub(const RingElement<N, D>&, const RingElement<N, D>&); \
    template RingElement<N, D> ring_neg(const RingElement<N, D>&);                           \
    template RingElement<N, D> ring_mul(const RingElement<N, D>&, const RingElement<N, D>&); \
    template RingElement<N, D> ring_times_zeta(const RingElement<N, D>&, int);               \
    template RingElement<N, D> galois(const RingElement<N, D>&, int);                        \
    template RingElement<N, D> conjugate(const RingElement<N, D>&);                          \
    template bool is_real(const RingElement<N, D>&);                                         \
    template RealElement<N / 2, D> to_real(const RingElement<N, D>&);                        \
    template RealElement<N / 2, D> relative_norm(const RingElement<N, D>&);                  \
    template int_mp ring_norm(const RingElement<N, D>&);                                     \
    template std::optional<RingElement<N, D>> ring_divide(const RingElement<N, D>&,          \
                                                          const RingElement<N, D>&);         \
    template std::complex<double> embed(const RingElement<N, D>&, int);                      \
    template ComplexDD embed_dd(const RingElement<N, D>&, int);

#define QUICKSYNTH_REAL(P, D)                                                                \
    template RingElement<2 * P, D> from_real(const RealElement<P, D>&);                      \
    template RealElement<P, D> real_add(const RealElement<P, D>&, const RealElement<P, D>&); \
    template RealElement<P, D> real_sub(const RealElement<P, D>&, const RealElement<P, D>&); \
    template RealElement<P, D> real_mul(const RealElement<P, D>&, const RealElement<P, D>&); \
    template int_mp real_norm(const RealElement<P, D>&);                                     \
    template int real_sign(const RealElement<P, D>&, int);                                   \
    template bool totally_positive(const RealElement<P, D>&);                                \
    template std::optional<RealElement<P, D>> real_divide(const RealElement<P, D>&,          \
                                                          const RealElement<P, D>&);         \
    template DD real_embed_dd(const RealElement<P, D>&, int);                                \
    template std::optional<RealElement<P, D>> unit_sqrt(const RealElement<P, D>&);

QUICKSYNTH_RING(4, __int128)
QUICKSYNTH_RING(8, __int128)
QUICKSYNTH_RING(4, int_mp)
QUICKSYNTH_RING(8, int_mp)
QUICKSYNTH_REAL(2, __int128)
QUICKSYNTH_REAL(4, __int128)
QUICKSYNTH_REAL(2, int_mp)
QUICKSYNTH_REAL(4, int_mp)

#undef QUICKSYNTH_RING
#undef QUICKSYNTH_REAL

}  // namespace quicksynth
