// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The rings of integers O_K = Z[zeta], zeta = exp(i pi / N), of the two gate
// sets: N = 4 (Z[omega], Clifford+T) and N = 8 (Z[zeta16], Clifford+sqrt(T)),
// and their real subrings O_K+ = Z[sqrt2] and Z[g], g = zeta + zeta^-1 =
// 2 cos(pi/8), of degree P = N / 2. The z-rotation path solves relative norm
// equations t t* = xi in them [Kliuchnikov, Bocharov, Roetteler, Yard,
// arXiv:1510.03888; Washington, GTM 83].
//
// A ring element is its coefficients over zeta^0 .. zeta^(N-1). A real
// element is held in tower coordinates, which make norms and signs exact:
//   P = 2: (a, b) = a + b sqrt2;
//   P = 4: (a0, a1, b0, b1) = (a0 + a1 sqrt2) + (b0 + b1 sqrt2) g, g^2 = 2 + sqrt2.
// The places of K+ are numbered p = 0 .. P-1, place p being zeta -> zeta^(2p+1);
// place 0 is the identity embedding.
//
// Digits are __int128 or int_mp. __int128 arithmetic is overflow-checked and
// throws std::overflow_error past 127 bits; norms, signs and divisions widen
// to int_mp internally where their intermediates need it.
#pragma once

#include "quicksynth/data_structures/double_double.hpp"
#include "quicksynth/data_structures/int_mp.hpp"

#include <array>
#include <complex>
#include <cstddef>
#include <optional>
#include <stdexcept>

namespace quicksynth {

#pragma GCC visibility push(hidden)

template <size_t N, class D>
using RingElement = std::array<D, N>;

template <size_t P, class D>
struct RealElement {
    std::array<D, P> c{};  // tower coordinates
    D& operator[](size_t i) { return c[i]; }
    const D& operator[](size_t i) const { return c[i]; }
    bool operator==(const RealElement& o) const { return c == o.c; }
    bool operator!=(const RealElement& o) const { return c != o.c; }
};

struct ComplexDD {
    DD re, im;
};

// zeta^j for any integer j.
template <size_t N, class D>
RingElement<N, D> ring_zeta_power(int j);

template <size_t N, class D>
RingElement<N, D> ring_add(const RingElement<N, D>& a, const RingElement<N, D>& b);
template <size_t N, class D>
RingElement<N, D> ring_sub(const RingElement<N, D>& a, const RingElement<N, D>& b);
template <size_t N, class D>
RingElement<N, D> ring_neg(const RingElement<N, D>& a);
template <size_t N, class D>
RingElement<N, D> ring_mul(const RingElement<N, D>& a, const RingElement<N, D>& b);
// zeta^j x.
template <size_t N, class D>
RingElement<N, D> ring_times_zeta(const RingElement<N, D>& x, int j);

// The automorphism zeta -> zeta^j, j odd; std::invalid_argument for even j.
template <size_t N, class D>
RingElement<N, D> galois(const RingElement<N, D>& x, int j);
// Complex conjugation, zeta -> zeta^-1.
template <size_t N, class D>
RingElement<N, D> conjugate(const RingElement<N, D>& x);

// x = x*.
template <size_t N, class D>
bool is_real(const RingElement<N, D>& x);
// Tower coordinates of a real x; std::invalid_argument if x is not real.
template <size_t N, class D>
RealElement<N / 2, D> to_real(const RingElement<N, D>& x);
template <size_t P, class D>
RingElement<2 * P, D> from_real(const RealElement<P, D>& r);

// x x*, in K+.
template <size_t N, class D>
RealElement<N / 2, D> relative_norm(const RingElement<N, D>& x);
// N_{K/Q}(x) = N_{K+/Q}(x x*) >= 0.
template <size_t N, class D>
int_mp ring_norm(const RingElement<N, D>& x);

template <size_t P, class D>
RealElement<P, D> real_add(const RealElement<P, D>& a, const RealElement<P, D>& b);
template <size_t P, class D>
RealElement<P, D> real_sub(const RealElement<P, D>& a, const RealElement<P, D>& b);
template <size_t P, class D>
RealElement<P, D> real_mul(const RealElement<P, D>& a, const RealElement<P, D>& b);

// N_{K+/Q}(r) through the tower: P = 2: a^2 - 2 b^2; P = 4: N_{Q(sqrt2)/Q} of
// N_{K+/Q(sqrt2)}(A + B g) = A^2 - B^2 (2 + sqrt2).
template <size_t P, class D>
int_mp real_norm(const RealElement<P, D>& r);

// Exact sign (-1, 0, 1) of r at place p. a + b sqrt2 with a, b of opposite
// signs has the sign of the larger of a^2 and 2 b^2; A + B gamma, gamma =
// +-sqrt(2 +- sqrt2), likewise from A^2 - B^2 gamma^2 in Z[sqrt2].
template <size_t P, class D>
int real_sign(const RealElement<P, D>& r, int place);
template <size_t P, class D>
bool totally_positive(const RealElement<P, D>& r);

// a / b if b divides a, else nullopt; b nonzero (std::domain_error). The
// quotient is a b* c / n with b b* c = n = N(b) an integer, c the tower
// cofactor; it exists iff n divides every coefficient. For __int128 digits a
// quotient past 127 bits throws std::overflow_error.
template <size_t N, class D>
std::optional<RingElement<N, D>> ring_divide(const RingElement<N, D>& a,
                                             const RingElement<N, D>& b);
template <size_t P, class D>
std::optional<RealElement<P, D>> real_divide(const RealElement<P, D>& a,
                                             const RealElement<P, D>& b);

// x at zeta -> exp(i pi k / N).
template <size_t N, class D>
std::complex<double> embed(const RingElement<N, D>& x, int k);
template <size_t N, class D>
ComplexDD embed_dd(const RingElement<N, D>& x, int k);
// r at place p.
template <size_t P, class D>
DD real_embed_dd(const RealElement<P, D>& r, int place);

// eta with eta^2 = u, positive at place 0, if u is a square; u is meant to be
// a totally positive unit, which always is one in these rings (Weber; see
// arXiv:1510.03888 and Washington, GTM 83). The roots of the place
// values are taken in double-double, mapped back to coordinates over each
// sign pattern by the tower's explicit inverse, rounded and checked exactly.
// The rounding is exact while the place values stay below about 2^90; past
// that a square may be missed (nullopt).
template <size_t P, class D>
std::optional<RealElement<P, D>> unit_sqrt(const RealElement<P, D>& u);

// a + b, a - b and a b, std::overflow_error past 127 bits.
[[noreturn]] void throw_int128_overflow();
inline __int128 checked_add(__int128 a, __int128 b) {
    __int128 r;
    if (__builtin_add_overflow(a, b, &r)) throw_int128_overflow();
    return r;
}
inline __int128 checked_sub(__int128 a, __int128 b) {
    __int128 r;
    if (__builtin_sub_overflow(a, b, &r)) throw_int128_overflow();
    return r;
}
inline __int128 checked_mul(__int128 a, __int128 b) {
    __int128 r;
    if (__builtin_mul_overflow(a, b, &r)) throw_int128_overflow();
    return r;
}

// Digit conversions; to_int128 throws std::overflow_error past 127 bits.
template <size_t M>
std::array<int_mp, M> to_int_mp(const std::array<__int128, M>& x) {
    std::array<int_mp, M> out;
    for (size_t i = 0; i < M; ++i) out[i] = int_mp(x[i]);
    return out;
}

template <size_t M>
std::array<__int128, M> to_int128(const std::array<int_mp, M>& x) {
    std::array<__int128, M> out;
    for (size_t i = 0; i < M; ++i) {
        if (x[i].bits() > 127) throw std::overflow_error("number_field: a digit passes 127 bits");
        out[i] = x[i].to_int128();
    }
    return out;
}

#pragma GCC visibility pop

}  // namespace quicksynth
