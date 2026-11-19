// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Exact Clifford+sqrt(T) elements M / sqrt2^k, k least [paper: def:exact]. M is
// 2x2 over Z[zeta16], zeta = exp(i pi/8); arithmetic is templated on the digit
// type (int64, __int128, int_mp).
#pragma once

#include "quicksynth/data_structures/gate_set.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace quicksynth {

constexpr int kZetaDegree = 8;  // Z[zeta16] as Z^8, zeta^8 = -1

// Z[zeta16] element and 2x2 matrix in digits D.
template <class D>
using Z16Of = std::array<D, kZetaDegree>;
template <class D>
using M2Of = std::array<std::array<Z16Of<D>, 2>, 2>;

using Z16 = Z16Of<int64_t>;
using M2 = M2Of<int64_t>;

// Hidden: exported calls go through the symbol table, too slow here.
#pragma GCC visibility push(hidden)

template <class D>
Z16Of<D> zeta_mul(const Z16Of<D>& a, const Z16Of<D>& b) {
    Z16Of<D> out{};
    for (int i = 0; i < kZetaDegree; ++i) {
        if (a[i] == 0) continue;
        for (int j = 0; j < kZetaDegree; ++j) {
            const int s = i + j;
            if (s < kZetaDegree)
                out[s] += a[i] * b[j];
            else
                out[s - kZetaDegree] -= a[i] * b[j];
        }
    }
    return out;
}

template <class D>
Z16Of<D> zeta_add(const Z16Of<D>& a, const Z16Of<D>& b) {
    Z16Of<D> out;
    for (int t = 0; t < kZetaDegree; ++t) out[t] = a[t] + b[t];
    return out;
}

template <class D>
Z16Of<D> zeta_sub(const Z16Of<D>& a, const Z16Of<D>& b) {
    Z16Of<D> out;
    for (int t = 0; t < kZetaDegree; ++t) out[t] = a[t] - b[t];
    return out;
}

// -a
template <class D>
Z16Of<D> zeta_neg(const Z16Of<D>& a) {
    Z16Of<D> out;
    for (int t = 0; t < kZetaDegree; ++t) out[t] = -a[t];
    return out;
}

// zeta^i -> zeta^-i = -zeta^(8-i)
template <class D>
Z16Of<D> zeta_conj(const Z16Of<D>& a) {
    Z16Of<D> out{};
    out[0] = a[0];
    for (int t = 1; t < kZetaDegree; ++t) out[kZetaDegree - t] = -a[t];
    return out;
}

// zeta^j * x for any integer j (zeta^16 = 1): a shift with a sign at each wrap,
// zeta^8 = -1.
template <class D>
Z16Of<D> times_zeta(const Z16Of<D>& x, int j) {
    j = ((j % 16) + 16) % 16;
    Z16Of<D> out;
    for (int t = 0; t < kZetaDegree; ++t) {
        const int s = t + j;
        out[s % kZetaDegree] = (s / kZetaDegree) % 2 != 0 ? -x[t] : x[t];
    }
    return out;
}

template <class D>
M2Of<D> matrix_mul(const M2Of<D>& a, const M2Of<D>& b) {
    M2Of<D> out;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            out[i][j] = zeta_add(zeta_mul(a[i][0], b[0][j]), zeta_mul(a[i][1], b[1][j]));
    return out;
}

// Conjugate transpose.
template <class D>
M2Of<D> matrix_dagger(const M2Of<D>& a) {
    M2Of<D> out;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) out[i][j] = zeta_conj(a[j][i]);
    return out;
}

// An int64 matrix in digits D.
template <class D>
M2Of<D> widened(const M2& m) {
    M2Of<D> out;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int t = 0; t < kZetaDegree; ++t) out[i][j][t] = D(m[i][j][t]);
    return out;
}
#pragma GCC visibility pop

// An exact element M / sqrt2^k with k least.
struct CyclotomicMatrix {
    M2 matrix{};                   // M
    int denominator_exponent = 0;  // k
};

// k and M's 32 row-major coefficients of any size, as hex text.
struct CyclotomicMatrixDigits {
    int denominator_exponent = 0;  // k
    std::array<std::string, 32> hex;
};

// Divides out common sqrt2 factors (k least).
void reduce_sqrt2(M2& m, int& k);

// The product a b, reduced.
CyclotomicMatrix compose(const CyclotomicMatrix& a, const CyclotomicMatrix& b);

// Whether u1 u1* + u2 u2* is exactly 2^K: (u1, u2) / sqrt2^K is a unit vector, the
// first column of an exact unitary. Every coefficient below 2^60 in absolute value,
// so that the norm is exact in 128 bits.
bool unit_norm(const Z16& u1, const Z16& u2, int K);

// S = [[u1, -u2*], [u2, u1*]] / sqrt2^K, K = denominator_exponent, reduced: the
// exact element with first column (u1, u2) / sqrt2^K, in SU(2) exactly when
// unit_norm(u1, u2, K). The S of U = zeta_n^p S Q^j (su2.hpp).
CyclotomicMatrix su2_element(const Z16& u1, const Z16& u2, int denominator_exponent);

// S Q^j for the gate set's Q = diag(1, zeta_n): S's second column times
// zeta_n^j. A unit factor, so S stays reduced.
CyclotomicMatrix times_q_power(const CyclotomicMatrix& S, int j, const GateSet& gateset);

// Reduced element of a word. Throws std::invalid_argument on an unknown letter,
// std::overflow_error past 64 bits.
CyclotomicMatrix word_matrix(const std::string& word);

// Without the 64-bit limit.
CyclotomicMatrixDigits word_matrix_digits(const std::string& word);

// std::invalid_argument naming the first unknown letter of `word`.
std::invalid_argument unknown_letter(const std::string& word);

}  // namespace quicksynth
