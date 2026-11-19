// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Double-double arithmetic against exact values: every double is a dyadic
// rational, so values scaled by a power of two compare exactly in int_mp.

#include "quicksynth/data_structures/double_double.hpp"
#include "quicksynth/data_structures/int_mp.hpp"

#include <doctest/doctest.h>

#include "oracles.hpp"

#include <cmath>
#include <cstdint>
#include <limits>

namespace quicksynth::test {

namespace {

constexpr int kScale = 400;  // fixed-point bits of the exact references

int_mp abs_of(const int_mp& x) { return x.sign() < 0 ? -x : x; }

// x 2^s exactly; s large enough that it is an integer.
int_mp exact(double x, int s = kScale) {
    const double scaled = std::ldexp(x, s);
    REQUIRE(std::trunc(scaled) == scaled);
    return int_mp::from_double(scaled);
}

int_mp exact(const DD& x, int s = kScale) { return exact(x.hi, s) + exact(x.lo, s); }

// |got - want| <= 2^-bits |want|.
bool within(const int_mp& got, const int_mp& want, int bits) {
    return compare(shift_left(abs_of(got - want), bits), abs_of(want)) <= 0;
}

// floor(sqrt(x)), x >= 0, by Newton's iteration from above.
int_mp isqrt(const int_mp& x) {
    if (x.zero()) return x;
    int_mp y = shift_left(int_mp(1), (x.bits() + 1) / 2 + 1);
    for (;;) {
        const int_mp next = half_of(y + x / y);
        if (compare(next, y) >= 0) return y;
        y = next;
    }
}

// A normalised double-double near 2^exponent with a random sign and lo.
DD random_dd(Rng& rng, int exponent) {
    const double hi = std::ldexp(1 + rng.uniform(), exponent) * (rng.below(2) != 0 ? -1 : 1);
    const double lo = std::ldexp(rng.uniform() - 0.5, exponent - 53);
    return quick_two_sum(hi, lo);
}

bool normalised(const DD& x) { return x.hi + x.lo == x.hi; }

}  // namespace

TEST_CASE("two_sum and two_prod are exact") {
    Rng rng(11);
    for (int i = 0; i < 20000; ++i) {
        const double a = std::ldexp(rng.uniform() - 0.5, static_cast<int>(rng.below(120)) - 60);
        const double b = std::ldexp(rng.uniform() - 0.5, static_cast<int>(rng.below(120)) - 60);
        const DD s = two_sum(a, b);
        CHECK(exact(s) == exact(a) + exact(b));
        CHECK(normalised(s));
        const DD p = two_prod(a, b);
        CHECK(exact(p, 2 * kScale - 100) == exact(a, kScale - 50) * exact(b, kScale - 50));
        CHECK(normalised(p));
    }
    // Operands whose 26-bit halves carry: all 53 bits set.
    const double full = std::ldexp(static_cast<double>((int64_t{1} << 53) - 1), -53);
    const DD p = two_prod(full, -full);
    CHECK(exact(p, 2 * kScale - 100) == exact(full, kScale - 50) * exact(-full, kScale - 50));
}

TEST_CASE("double-double operations are accurate to about 2^-104") {
    Rng rng(12);
    for (int i = 0; i < 20000; ++i) {
        const DD a = random_dd(rng, static_cast<int>(rng.below(60)) - 30);
        const DD b = random_dd(rng, static_cast<int>(rng.below(60)) - 30);
        const double c = std::ldexp(rng.uniform() + 0.5, static_cast<int>(rng.below(40)) - 20);
        const int_mp ea = exact(a), eb = exact(b);
        // scale 2 kScale for the products
        const int_mp ea2 = shift_left(ea, kScale);

        const DD sum = a + b;
        CHECK(normalised(sum));
        CHECK(within(exact(sum), ea + eb, 104));
        CHECK(within(exact(a - b), ea - eb, 104));
        CHECK(within(exact(a + c), ea + exact(c), 104));

        const DD product = a * b;
        CHECK(normalised(product));
        CHECK(within(shift_left(exact(product), kScale), ea * eb, 102));
        CHECK(within(shift_left(exact(a * c), kScale), ea * exact(c), 103));

        // q b against a, and r^2 against |a|
        const DD q = a / b;
        CHECK(normalised(q));
        CHECK(within(exact(q) * eb, ea2, 102));
        const DD r = dd_sqrt(dd_abs(a));
        CHECK(normalised(r));
        CHECK(within(exact(r) * exact(r), abs_of(ea2), 102));
    }
    CHECK(dd_sqrt(DD{}) == DD{});
}

TEST_CASE("double-double comparison is exact") {
    Rng rng(13);
    for (int i = 0; i < 5000; ++i) {
        const DD a = random_dd(rng, 0);
        DD b = random_dd(rng, 0);
        if (rng.below(2) != 0) b = quick_two_sum(a.hi, std::ldexp(rng.uniform() - 0.5, -53));
        const int want = compare(exact(a), exact(b));
        CHECK(compare(a, b) == want);
        CHECK((a < b) == (want < 0));
        CHECK((a >= b) == (want >= 0));
    }
    CHECK(dd_sign(DD{-1, 0}) == -1);
    CHECK(dd_sign(DD{}) == 0);
}

TEST_CASE("integers convert to double-double exactly up to 106 bits") {
    const int64_t small[] = {0,
                             1,
                             -1,
                             (int64_t{1} << 53) + 1,
                             std::numeric_limits<int64_t>::max(),
                             std::numeric_limits<int64_t>::min(),
                             -987654321987654321};
    for (int64_t x : small) {
        CHECK(exact(dd(x)) == shift_left(int_mp(x), kScale));
        CHECK(normalised(dd(x)));
    }
    Rng rng(14);
    for (int i = 0; i < 2000; ++i) {
        const int bits = 1 + static_cast<int>(rng.below(127));
        __int128 x = 0;
        for (int b = 0; b < bits; ++b) x = (x << 1) | static_cast<__int128>(rng.below(2));
        if (rng.below(2) != 0) x = -x;
        const int_mp want = shift_left(int_mp(x), kScale);
        const DD got = dd(x);
        CHECK(normalised(got));
        if (bits <= 106)
            CHECK(exact(got) == want);
        else
            CHECK(within(exact(got), want, 105));
        const DD via_mp = dd(int_mp(x));
        CHECK(via_mp == got);
    }
    const __int128 max = ~(static_cast<__int128>(1) << 127);
    for (__int128 x : {max, -max, -max - 1}) {
        CHECK(within(exact(dd(x)), shift_left(int_mp(x), kScale), 105));
    }
    // int_mp past 128 bits, with an exponent
    const int_mp big = shift_left(int_mp(int64_t{0x123456789abcdef}), 150) + int_mp(int64_t{77});
    const DD scaled = dd(big, -200);
    CHECK(within(exact(scaled, kScale + 200), shift_left(big, kScale), 105));
    CHECK(dd(int_mp(), 30) == DD{});
}

TEST_CASE("nearest_integer is within one half") {
    Rng rng(15);
    for (int i = 0; i < 5000; ++i) {
        const DD a = random_dd(rng, static_cast<int>(rng.below(140)) - 10);
        const int_mp n = nearest_integer(a);
        const int_mp twice_gap = shift_left(abs_of(exact(a) - shift_left(n, kScale)), 1);
        CHECK(compare(twice_gap, shift_left(int_mp(1), kScale)) <= 0);
    }
    // a fraction in lo below a hi past 2^53
    CHECK(nearest_integer(quick_two_sum(0x1p60, 0.75)) == shift_left(int_mp(1), 60) + int_mp(1));
    CHECK(nearest_integer(quick_two_sum(-0x1p60, -0.75)) ==
          -(shift_left(int_mp(1), 60) + int_mp(1)));
    CHECK(nearest_integer(DD{2.5, -0x1p-60}) == int_mp(2));
}

TEST_CASE("cos and sin of multiples of pi/16 are correctly rounded") {
    // Exact references 2^S cos(k pi/16) from nested radicals, rounded down.
    constexpr int kS = 300;
    const int_mp one = shift_left(int_mp(1), kS);
    const auto root = [&](const int_mp& x) { return isqrt(shift_left(x, kS)); };  // sqrt(x/2^S) 2^S
    const int_mp two = shift_left(one, 1);
    const int_mp s2 = root(two);
    const int_mp a = root(two + s2), b = root(two - s2);  // 2 cos(pi/8), 2 cos(3 pi/8)
    const int_mp want[9] = {one,         half_of(root(two + a)), half_of(a), half_of(root(two + b)),
                            half_of(s2), half_of(root(two - b)), half_of(b), half_of(root(two - a)),
                            int_mp()};
    for (int k = 0; k <= 8; ++k) {
        const DD c = cos_pi16(k);
        CHECK(normalised(c));
        // within one ulp of lo, plus the reference's rounding
        const int_mp gap = abs_of(exact(c, kS) - want[k]);
        CHECK(compare(gap, shift_left(int_mp(1), kS - 106)) <= 0);
    }
    CHECK(within(exact(sqrt2_dd(), kS), s2, 106));
    for (int k = -40; k <= 40; ++k) {
        const double angle = k * M_PI / 16;
        CHECK(std::abs(to_double(cos_pi16(k)) - std::cos(angle)) < 1e-15);
        CHECK(std::abs(to_double(sin_pi16(k)) - std::sin(angle)) < 1e-15);
        const DD c = cos_pi16(k), s = sin_pi16(k);
        CHECK(within(exact(c * c + s * s), shift_left(int_mp(1), kScale), 103));
    }
}

// The rotation body is a circular segment |u| <= R, Re u >= R sqrt(1 - eps^2)
// of depth R eps^2 / 2; at eps = 1e-13 and R = 2^33 that is 4e-17, below a
// double's resolution at R. Points at eighths of the depth either side of the
// chord must be told apart, and the depth must be free of cancellation.
TEST_CASE("the segment's chord and depth are exact enough at 1e-10 and 1e-13") {
    for (double eps : {1e-10, 1e-13}) {
        for (int K : {40, 41, 66, 67}) {
            CAPTURE(eps);
            CAPTURE(K);
            const DD e2 = two_prod(eps, eps);
            const DD c = dd_sqrt(DD{1, 0} - e2);
            const DD R =
                K % 2 == 0 ? DD{std::ldexp(1.0, K / 2), 0} : dd_ldexp(sqrt2_dd(), (K - 1) / 2);
            const DD chord = R * c;
            const DD depth_ratio = e2 / (c + 1.0);  // 1 - sqrt(1 - eps^2)
            const DD depth = R * depth_ratio;

            // (1 - h')^2 = 1 - eps^2: 2 h' - h'^2 = eps^2 to 2^-100 of eps^2
            const int_mp h = exact(depth_ratio);
            const int_mp lhs = shift_left(h, kScale + 1) - h * h;
            CHECK(within(lhs, exact(eps) * exact(eps), 100));

            // x^2 against R^2 (1 - eps^2) = 2^K (1 - eps^2), exactly
            const int_mp rhs =
                shift_left(shift_left(int_mp(1), 2 * kScale) - exact(eps) * exact(eps), K);
            int double_wrong = 0;
            for (int j = -12; j <= 12; ++j) {
                if (j == 0) continue;
                const DD x = chord + depth * (j / 8.0);
                const int_mp ex = exact(x);
                const int want = compare(ex * ex, rhs);
                CHECK(want == (j > 0 ? 1 : -1));
                CHECK(compare(x, chord) == want);
                const double chord_double = to_double(R) * std::sqrt(1 - eps * eps);
                if ((to_double(x) >= chord_double) != (want > 0)) ++double_wrong;
            }
            // doubles alone cannot place these points
            CHECK(double_wrong > 0);
        }
    }
}

}  // namespace quicksynth::test
