// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// LLL reduction of the walk's basis, with the unimodular transform and its
// exact inverse. A Schnorr-Euchner floating-point LLL with L2's tolerances.
// Cohen's integral LLL on int_mp is its fallback.
// Also: exact_product, the reduced basis A T (A the ellipsoid's basis, T the
// integer LLL transform) with each entry rounded once; ordered_matvec, which maps
// the centre c to the reduced coordinates T^-1 c; product, the plain matrix product;
// and qr_positive, the QR factors the walk's forms are built from.
#pragma once

#include "quicksynth/data_structures/int_mp.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace quicksynth {

// A unimodular n x n matrix T and its inverse, row-major.
struct LllTransform {
    std::vector<int64_t> T, T_inv;
};

// The integer basis rows of A (n x n, row-major). Row j is column j of A times
// 2^precision_bits, rounded half to even.
std::vector<double> lll_rows(int n, const double* A, int precision_bits);

// T with B T LLL-reduced, where row j of `rows` is column j of B, and T^-1.
// False if an entry of T is not below `limit` in absolute value or one of T^-1
// is not below 2^63 - 1. Throws std::invalid_argument on linearly dependent rows.
bool lll_transform(int n, const std::vector<double>& rows, int64_t limit, LllTransform& out);
bool lll_transform(int n, const std::vector<int_mp>& rows, int64_t limit, LllTransform& out);

// The same with T^-1 of any size, for bases skewed enough that it passes 64
// bits (the z-rotation body below eps = 1e-12).
#pragma GCC visibility push(hidden)
struct LllTransformMp {
    std::vector<int64_t> T;
    std::vector<int_mp> T_inv;
};
#pragma GCC visibility pop
bool lll_transform(int n, const std::vector<int_mp>& rows, int64_t limit, LllTransformMp& out);

// First `rows` rows of A T, each entry the exact product rounded once. A spans 1/eps^2 while
// A T is O(1), so doubles would cancel. Uses 128-bit integers, or int_mp for a row that
// spans too many binades or whose sums could overflow 128 bits.
void exact_product(int n, int rows, const double* A, const int64_t* T, double* out);

// M x for an n x n matrix M of integers or doubles, each row summed left to right.
template <typename T>
std::vector<double> ordered_matvec(int n, const std::vector<T>& M, const std::vector<double>& x) {
    std::vector<double> out(n);
    for (int i = 0; i < n; ++i) {
        const T* row = M.data() + static_cast<size_t>(i) * n;
        double s = static_cast<double>(row[0]) * x[0];
        for (int j = 1; j < n; ++j) s += static_cast<double>(row[j]) * x[j];
        out[i] = s;
    }
    return out;
}

// The product A B of two n x n matrices, each entry summed in the order of j.
template <typename T>
std::vector<T> product(int n, const std::vector<T>& A, const std::vector<T>& B) {
    std::vector<T> AB(static_cast<size_t>(n) * n, T{0});
    for (int r = 0; r < n; ++r)
        for (int j = 0; j < n; ++j)
            for (int c = 0; c < n; ++c)
                AB[static_cast<size_t>(r) * n + c] +=
                    A[static_cast<size_t>(r) * n + j] * B[static_cast<size_t>(j) * n + c];
    return AB;
}

// The factors of A = Q R, row-major.
struct QrFactors {
    std::vector<double> Q;  // m x m, orthogonal
    std::vector<double> R;  // n x n, upper triangular, its rows m to n-1 zero
};

// Householder QR of the m x n row-major A, m <= n (a block of rows of a basis),
// every diagonal entry of R made non-negative by flipping the sign of its row of
// R and its column of Q. R is padded with zero rows to n x n, the shape of a
// walk's form: a block of m rows constrains only the last m coordinates the walk
// fixes.
QrFactors qr_positive(int m, int n, const double* A);

}  // namespace quicksynth
