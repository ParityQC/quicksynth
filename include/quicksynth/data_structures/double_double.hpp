// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Double-double numbers: an unevaluated sum hi + lo of two doubles, about 106
// significant bits [Dekker, Numer. Math. 18 (1971); Shewchuk, DCG 18 (1997)].
// The z-rotation path needs them where a thin circular segment of depth
// R eps^2 / 2 sits at distance R from the origin.
//
// The exact products use Veltkamp's split, not a fused multiply-add: the build
// compiles with -ffp-contract=off, and the error-free transformations need
// every operation rounded on its own. Arguments stay below about 2^996 in
// magnitude, where the split cannot overflow.
#pragma once

#include <cmath>
#include <cstdint>

namespace quicksynth {

class int_mp;  // data_structures/int_mp.hpp

#pragma GCC visibility push(hidden)

// hi + lo with hi = fl(hi + lo).
struct DD {
    double hi = 0;
    double lo = 0;
};

// a + b exactly.
inline DD two_sum(double a, double b) {
    const double s = a + b;
    const double bb = s - a;
    return {s, (a - (s - bb)) + (b - bb)};
}

// a + b exactly, for |a| >= |b| or a = 0.
inline DD quick_two_sum(double a, double b) {
    const double s = a + b;
    return {s, b - (s - a)};
}

// a = hi + lo with hi and lo of at most 26 significant bits each [Veltkamp].
inline void split(double a, double& hi, double& lo) {
    const double c = 134217729.0 * a;  // 2^27 + 1
    hi = c - (c - a);
    lo = a - hi;
}

// a b exactly [Dekker 1971].
inline DD two_prod(double a, double b) {
    double ah, al, bh, bl;
    split(a, ah, al);
    split(b, bh, bl);
    const double p = a * b;
    return {p, ((ah * bh - p) + ah * bl + al * bh) + al * bl};
}

inline DD operator-(const DD& a) { return {-a.hi, -a.lo}; }

inline DD operator+(const DD& a, const DD& b) {
    DD s = two_sum(a.hi, b.hi);
    const DD t = two_sum(a.lo, b.lo);
    s = quick_two_sum(s.hi, s.lo + t.hi);
    return quick_two_sum(s.hi, s.lo + t.lo);
}

inline DD operator+(const DD& a, double b) {
    const DD s = two_sum(a.hi, b);
    return quick_two_sum(s.hi, s.lo + a.lo);
}

inline DD operator-(const DD& a, const DD& b) { return a + (-b); }
inline DD operator-(const DD& a, double b) { return a + (-b); }

inline DD operator*(const DD& a, const DD& b) {
    const DD p = two_prod(a.hi, b.hi);
    return quick_two_sum(p.hi, p.lo + (a.hi * b.lo + a.lo * b.hi));
}

inline DD operator*(const DD& a, double b) {
    const DD p = two_prod(a.hi, b);
    return quick_two_sum(p.hi, p.lo + a.lo * b);
}

// a / b, b nonzero: three quotient digits, each from the remainder so far.
DD operator/(const DD& a, const DD& b);

// sqrt(a), a >= 0: one Newton step from the double root. The dd_ prefix keeps
// std's double functions visible to unqualified calls in this namespace.
DD dd_sqrt(const DD& a);

// a 2^e, exact unless it leaves the normal range.
inline DD dd_ldexp(const DD& a, int e) { return {std::ldexp(a.hi, e), std::ldexp(a.lo, e)}; }

inline double to_double(const DD& a) { return a.hi + a.lo; }

// -1, 0 or 1 as a < b, a == b or a > b; exact for normalised operands.
inline int compare(const DD& a, const DD& b) {
    if (a.hi != b.hi) return a.hi < b.hi ? -1 : 1;
    return (a.lo > b.lo) - (a.lo < b.lo);
}
inline bool operator<(const DD& a, const DD& b) { return compare(a, b) < 0; }
inline bool operator<=(const DD& a, const DD& b) { return compare(a, b) <= 0; }
inline bool operator>(const DD& a, const DD& b) { return compare(a, b) > 0; }
inline bool operator>=(const DD& a, const DD& b) { return compare(a, b) >= 0; }
inline bool operator==(const DD& a, const DD& b) { return compare(a, b) == 0; }
inline bool operator!=(const DD& a, const DD& b) { return compare(a, b) != 0; }

inline int dd_sign(const DD& a) { return (a.hi > 0) - (a.hi < 0); }
inline DD dd_abs(const DD& a) { return a.hi < 0 ? -a : a; }

// Integers to double-double: hi is x rounded to a double, lo the rest rounded;
// exact up to 106 bits.
DD dd(int64_t x);
DD dd(__int128 x);
// x 2^exponent.
DD dd(const int_mp& x, int exponent = 0);

// The integer nearest a; a halfway case rounds either way.
int_mp nearest_integer(const DD& a);

// cos(k pi / 16) and sin(k pi / 16), correctly rounded, for any integer k.
DD cos_pi16(int k);
DD sin_pi16(int k);
// sqrt(2), correctly rounded.
DD sqrt2_dd();

#pragma GCC visibility pop

}  // namespace quicksynth
