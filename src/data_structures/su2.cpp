// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/su2.hpp"

#include <algorithm>
#include <cmath>

namespace quicksynth {

namespace {

void complex_product(double xr, double xi, double yr, double yi, double& zr, double& zi) {
    zr = xr * yr - xi * yi;
    zi = xr * yi + xi * yr;
}

// The first column (a, b) of V = (a c; b d) as (Re a, Im a, Re b, Im b), with V
// first multiplied by the global phase that makes det V = 1 / phase:
// (a, b) / sqrt(phase det V), principal root.
// With its determinant fixed, a unitary is determined by its first column up to
// sign, so this vector stands for V in the lattice search. phase = 1 gives V's
// SU(2) vector; search_direction() passes zeta_n^-j for determinant class j.
Vec4 su2_direction(const Unitary& V, double phase_re, double phase_im) {
    const double ar = V[0], ai = V[1], cr = V[2], ci = V[3];
    const double br = V[4], bi = V[5], dr = V[6], di = V[7];
    // x = phase det V, det V = a d - c b
    double pr, pi, qr, qi, xr, xi;
    complex_product(ar, ai, dr, di, pr, pi);
    complex_product(cr, ci, br, bi, qr, qi);
    complex_product(pr - qr, pi - qi, phase_re, phase_im, xr, xi);
    // r = sqrt(x), principal root, each branch free of cancellation
    const double modulus = std::sqrt(xr * xr + xi * xi);
    double rr, ri;
    if (xr >= 0.0) {
        rr = std::sqrt((modulus + xr) / 2);
        ri = xi / (2 * rr);
    } else {
        ri = std::copysign(std::sqrt((modulus - xr) / 2), xi);
        rr = xi / (2 * ri);
    }
    // (a, b) / r = (a, b) conj(r) / |r|^2
    const double norm = rr * rr + ri * ri;
    return {(ar * rr + ai * ri) / norm, (ai * rr - ar * ri) / norm, (br * rr + bi * ri) / norm,
            (bi * rr - br * ri) / norm};
}

}  // namespace

std::array<double, 2> zeta_power(int j) {
    const double sqrt2 = std::sqrt(2.0);
    const double c1 = std::sqrt(2.0 + sqrt2) / 2;  // cos(pi/8)
    const double s1 = std::sqrt(2.0 - sqrt2) / 2;  // sin(pi/8)
    const double h = sqrt2 / 2;                    // cos(pi/4)
    const std::array<double, 2> half_turn[8] = {{1.0, 0.0}, {c1, s1},  {h, h},  {s1, c1},
                                                {0.0, 1.0}, {-s1, c1}, {-h, h}, {-c1, s1}};
    j = ((j % 16) + 16) % 16;
    if (j < 8) return half_turn[j];
    return {-half_turn[j - 8][0], -half_turn[j - 8][1]};
}

double sqrt2_power(int n) {
    if (n < 0) {
        const int K = -n;
        return std::ldexp(K % 2 == 0 ? 1.0 : std::sqrt(0.5), -(K / 2));
    }
    return std::ldexp(n % 2 == 0 ? 1.0 : std::sqrt(2.0), n / 2);
}

Vec4 search_direction(const Unitary& V, int det_class, const GateSet& gateset) {
    const auto phase = zeta_power(-gateset.step * det_class);  // det Q^-j
    return su2_direction(V, phase[0], phase[1]);
}

Vec4 su2_vector(const Unitary& V) { return su2_direction(V, 1.0, 0.0); }

double chord_distance(const Vec4& a, const Vec4& b) {
    double minus = 0.0, plus = 0.0;
    for (int k = 0; k < 4; ++k) {
        minus += (a[k] - b[k]) * (a[k] - b[k]);
        plus += (a[k] + b[k]) * (a[k] + b[k]);
    }
    const double s = std::min(minus, plus);
    return std::sqrt(std::max(s - s * s / 4, 0.0));
}

Unitary matrix_unitary(int k, const std::array<double, 32>& coefficients) {
    const double sc = sqrt2_power(-k);
    Unitary V;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            double re = 0.0, im = 0.0;
            for (int t = 0; t < kZetaDegree; ++t) {
                const auto z = zeta_power(t);
                re += coefficients[16 * i + 8 * j + t] * z[0];
                im += coefficients[16 * i + 8 * j + t] * z[1];
            }
            V[4 * i + 2 * j] = re * sc;
            V[4 * i + 2 * j + 1] = im * sc;
        }
    return V;
}

Unitary matrix_unitary(const CyclotomicMatrix& e) {
    std::array<double, 32> coefficients;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int t = 0; t < kZetaDegree; ++t)
                coefficients[16 * i + 8 * j + t] = static_cast<double>(e.matrix[i][j][t]);
    return matrix_unitary(e.denominator_exponent, coefficients);
}

Vec4 matrix_direction(const CyclotomicMatrix& e) { return su2_vector(matrix_unitary(e)); }

Unitary dagger(const Unitary& M) { return {M[0], -M[1], M[4], -M[5], M[2], -M[3], M[6], -M[7]}; }

Unitary matmul(const Unitary& A, const Unitary& B) {
    Unitary C;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            const double ar0 = A[4 * i], ai0 = A[4 * i + 1], ar1 = A[4 * i + 2], ai1 = A[4 * i + 3];
            const double br0 = B[2 * j], bi0 = B[2 * j + 1], br1 = B[4 + 2 * j],
                         bi1 = B[4 + 2 * j + 1];
            C[4 * i + 2 * j] = (ar0 * br0 - ai0 * bi0) + (ar1 * br1 - ai1 * bi1);
            C[4 * i + 2 * j + 1] = (ar0 * bi0 + ai0 * br0) + (ar1 * bi1 + ai1 * br1);
        }
    return C;
}

double ball_measure(double epsilon) {
    const double a = std::asin(epsilon);
    if (a < 1e-3) {  // series, to avoid cancellation
        const double x = 2 * a;
        return (std::pow(x, 3) / 6 - std::pow(x, 5) / 120 + std::pow(x, 7) / 5040) / M_PI;
    }
    return (2 * a - std::sin(2 * a)) / M_PI;
}

}  // namespace quicksynth
