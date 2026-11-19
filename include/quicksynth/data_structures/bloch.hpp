// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Bloch image of a Clifford+sqrt(T) unitary over Z[lambda], lambda = zeta + zeta^-1.
//
// For U = M / sqrt2^k the image is N u^(k+1) / lambda^(4(k+1)), with real
// N[i][j] = tr(sigma_i M sigma_j M^dagger) and the unit u = 2 lambda^2 - 1
// (1/2 = u / lambda^4). Common powers of lambda are stripped from N before
// multiplying by u^(k+1), which keeps digits small.
//
// Templated over the digit type: __int128 for the scan (try_*_128). The type int_mp
// is used for elements of any size.
#pragma once

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"

#include <array>
#include <cstdint>
#include <optional>
#include <string>

namespace quicksynth {

class int_mp;               // data_structures/int_mp.hpp
struct CyclotomicMatrixMp;  // data_structures/int_mp.hpp

// Element of Z[lambda] as coefficients of 1, lambda, lambda^2, lambda^3, and a
// row-major 3x3 matrix of them.
template <class D>
using Lambda = std::array<D, 4>;
template <class D>
using LambdaMatrix = std::array<Lambda<D>, 9>;

// Bloch image m / lambda^denominator_exponent with the least exponent.
template <class D>
struct BlochImage {
    LambdaMatrix<D> m{};
    int64_t denominator_exponent = 0;
};
using BlochImage128 = BlochImage<__int128>;

// Digit operations for __int128; data_structures/int_mp.hpp defines int_mp's.
inline bool is_even(__int128 x) { return (x & 1) == 0; }
inline __int128 half_of(__int128 x) { return x / 2; }

// Hidden, as data_structures/cyclotomic_matrix.hpp's templates are.
#pragma GCC visibility push(hidden)

// Entry (row, col) of a row-major matrix.
template <class D>
const Lambda<D>& at(const LambdaMatrix<D>& m, int row, int col) {
    return m[3 * row + col];
}

template <class D>
Lambda<D>& at(LambdaMatrix<D>& m, int row, int col) {
    return m[3 * row + col];
}

template <class D>
bool is_zero(const Lambda<D>& x) {
    return x[0] == 0 && x[1] == 0 && x[2] == 0 && x[3] == 0;
}

template <class D>
Lambda<D> add(const Lambda<D>& a, const Lambda<D>& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]};
}

template <class D>
Lambda<D> sub(const Lambda<D>& a, const Lambda<D>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]};
}

template <class D>
Lambda<D> scale(const Lambda<D>& a, int s) {
    return {a[0] * s, a[1] * s, a[2] * s, a[3] * s};
}

// s * x for a sign s = +-1.
template <class D>
Lambda<D> times_sign(const Lambda<D>& x, int s) {
    return s == 1 ? x : scale(x, -1);
}

// x * lambda, reduced by lambda^4 = 4 lambda^2 - 2.
template <class D>
Lambda<D> times_lambda(const Lambda<D>& x) {
    return {x[3] * -2, x[0], x[1] + x[3] * 4, x[2]};
}

template <class D>
Lambda<D> times_lambda_pow(Lambda<D> x, int n) {
    for (int i = 0; i < n; ++i) x = times_lambda(x);
    return x;
}

// x * u, u = 2 lambda^2 - 1.
template <class D>
Lambda<D> times_u(const Lambda<D>& x) {
    return {-x[0] - x[2] * 4, -x[1] - x[3] * 4, x[0] * 2 + x[2] * 7, x[1] * 2 + x[3] * 7};
}

// x * theta, theta = lambda^2 - 3.
template <class D>
Lambda<D> times_theta(const Lambda<D>& x) {
    return sub(times_lambda_pow(x, 2), scale(x, 3));
}

// True if lambda divides x, i.e. its constant digit is even
// (1/lambda = (4 lambda - lambda^3) / 2). Zero is divisible by every power.
template <class D>
bool divisible(const Lambda<D>& x) {
    return is_even(x[0]);
}

// x / lambda; requires divisible(x).
template <class D>
Lambda<D> over_lambda(const Lambda<D>& x) {
    const D half = half_of(x[0]);
    return {x[1], x[2] + half * 4, x[3], -half};
}

// Divides out the largest common power of lambda in place and returns its
// exponent; a zero matrix is unchanged.
template <class D>
int64_t strip_lambda(LambdaMatrix<D>& n) {
    int64_t stripped = 0;
    for (;;) {
        bool all = true, any = false;
        for (const Lambda<D>& cell : n) {
            if (!divisible(cell)) all = false;
            if (!is_zero(cell)) any = true;
        }
        if (!all || !any) return stripped;
        for (Lambda<D>& cell : n)
            if (!is_zero(cell)) cell = over_lambda(cell);
        ++stripped;
    }
}

// m X in place; X's image is diag(1, -1, -1), so columns 1 and 2 are negated.
template <class D>
void apply_x(LambdaMatrix<D>& m) {
    for (int r = 0; r < 3; ++r) {
        at(m, r, 1) = scale(at(m, r, 1), -1);
        at(m, r, 2) = scale(at(m, r, 2), -1);
    }
}

#pragma GCC visibility pop

// Scan pricing in 128-bit digits. None if a coefficient of u reaches 2^59 or an
// image digit would pass 2^100; the caller then uses the int_mp versions below.

// Least exponent only.
std::optional<int64_t> try_bloch_exponent_128(const M2& u, int k);

// Image (matrix, least exponent).
std::optional<BlochImage128> try_bloch_image_128(const M2& u, int k);

// Image of u / sqrt2^k in int_mp digits: 36 decimal digits (row, column, power of
// lambda) and the least exponent. Defined for any 2x2 u; the descent checks
// unitarity.
struct BlochDigits {
    std::array<std::string, 36> m;
    int64_t exponent = 0;
};

// Digit of lambda^power in entry (row, col).
inline const std::string& at(const BlochDigits& b, int row, int col, int power) {
    return b.m[static_cast<size_t>(12 * row + 4 * col + power)];
}
inline std::string& at(BlochDigits& b, int row, int col, int power) {
    return b.m[static_cast<size_t>(12 * row + 4 * col + power)];
}
BlochDigits bloch_digits(const M2& u, int k);
BlochDigits bloch_digits(const CyclotomicMatrixDigits& u);

// Least exponent in int_mp digits; half of it is a Clifford+T element's T-count.
int64_t bloch_exponent(const M2& u, int k);
int64_t bloch_exponent(const CyclotomicMatrixDigits& u);

// Image of u in int_mp digits.
BlochImage<int_mp> bloch_image(const CyclotomicMatrixMp& u);

}  // namespace quicksynth
