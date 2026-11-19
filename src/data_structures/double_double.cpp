// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/double_double.hpp"
#include "quicksynth/data_structures/int_mp.hpp"

#include <cmath>
#include <cstdint>

namespace quicksynth {

DD operator/(const DD& a, const DD& b) {
    const double q1 = a.hi / b.hi;
    DD r = a - b * q1;
    const double q2 = r.hi / b.hi;
    r = r - b * q2;
    const double q3 = r.hi / b.hi;
    return quick_two_sum(q1, q2) + q3;
}

DD dd_sqrt(const DD& a) {
    if (a.hi <= 0) return {};
    const double y = std::sqrt(a.hi);
    // a - y^2 is exact to double-double precision since y^2 = two_prod(y, y).
    const DD r = a - two_prod(y, y);
    return quick_two_sum(y, r.hi / (2 * y));
}

DD dd(int64_t x) {
    const double hi = static_cast<double>(x);
    // In 128 bits: hi may be 2^63, outside int64.
    const __int128 rest = static_cast<__int128>(x) - static_cast<__int128>(hi);
    return quick_two_sum(hi, static_cast<double>(rest));
}

DD dd(__int128 x) {
    const double hi = static_cast<double>(x);
    __int128 rest;
    if (std::fabs(hi) < 0x1p126) {
        rest = x - static_cast<__int128>(hi);
    } else {
        // hi may be +-2^127, outside __int128; its half is not.
        const __int128 half = static_cast<__int128>(hi * 0.5);
        rest = (x - half) - half;
    }
    return quick_two_sum(hi, static_cast<double>(rest));
}

DD dd(const int_mp& x, int exponent) {
    const double hi = x.to_double(exponent);
    if (hi == 0) return {};
    // hi 2^-exponent is x rounded to 53 bits, an integer.
    const int_mp rest = x - int_mp::from_double(std::ldexp(hi, -exponent));
    return quick_two_sum(hi, rest.to_double(exponent));
}

int_mp nearest_integer(const DD& a) {
    const double whole = std::nearbyint(a.hi);
    // a.hi - whole is exact: both are multiples of ulp(a.hi) within 1/2.
    const double fraction = std::nearbyint((a.hi - whole) + a.lo);
    return int_mp::from_double(whole) + int_mp::from_double(fraction);
}

namespace {

// cos(k pi / 16), k = 0..8, rounded to double-double from exact nested radicals
// (cos(pi/8) = sqrt(2 + sqrt2) / 2, cos(pi/16) = sqrt(2 + sqrt(2 + sqrt2)) / 2).
constexpr DD kCosPi16[9] = {
    {0x1p+0, 0},
    {0x1.f6297cff75cb0p-1, 0x1.562172a361fd3p-56},
    {0x1.d906bcf328d46p-1, 0x1.457e610231ac2p-56},
    {0x1.a9b66290ea1a3p-1, 0x1.9f630e8b6dac8p-60},
    {0x1.6a09e667f3bcdp-1, -0x1.bdd3413b26456p-55},
    {0x1.1c73b39ae68c8p-1, 0x1.b25dd267f6600p-55},
    {0x1.87de2a6aea963p-2, -0x1.72cedd3d5a610p-57},
    {0x1.8f8b83c69a60bp-3, -0x1.26d19b9ff8d82p-57},
    {0, 0},
};

}  // namespace

DD cos_pi16(int k) {
    k = ((k % 32) + 32) % 32;
    if (k > 16) k = 32 - k;  // cos is even
    return k <= 8 ? kCosPi16[k] : -kCosPi16[16 - k];
}

DD sin_pi16(int k) { return cos_pi16(8 - k); }

DD sqrt2_dd() { return dd_ldexp(kCosPi16[4], 1); }

}  // namespace quicksynth
