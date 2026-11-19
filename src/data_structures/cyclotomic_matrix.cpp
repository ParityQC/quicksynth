// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"

#include "quicksynth/data_structures/int_mp.hpp"

#include <optional>
#include <utility>

namespace quicksynth {

namespace {

// --- dividing by sqrt2, in int64 and in digits of any size -----------------------

// The digit operations the templates below need; int_mp.hpp has int_mp's.
bool is_even(int64_t x) { return (x & 1) == 0; }
int64_t half_of(int64_t x) { return x / 2; }

// x / sqrt2 is in Z[zeta] iff x_i and x_(i+4) agree in parity
// (sqrt2 = zeta^2 - zeta^6).
template <class D>
bool divisible_by_sqrt2(const Z16Of<D>& x) {
    for (int i = 0; i < 4; ++i)
        if (is_even(x[i]) != is_even(x[i + 4])) return false;
    return true;
}

// x / sqrt2, for an x that divisible_by_sqrt2() passed.
template <class D>
Z16Of<D> divide_by_sqrt2(const Z16Of<D>& x) {
    Z16Of<D> y = zeta_sub(times_zeta(x, 2), times_zeta(x, 6));
    for (D& v : y) v = half_of(v);
    return y;
}

// M / sqrt2^k with every common factor sqrt2 stripped, so that k is least.
template <class D>
void strip_sqrt2(M2Of<D>& m, int& k) {
    while (k > 0) {
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                if (!divisible_by_sqrt2(m[i][j])) return;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) m[i][j] = divide_by_sqrt2(m[i][j]);
        --k;
    }
}

// u <- u g for the gate g a letter names: a column operation, as each gate is
// monomial or H.
bool apply_gate(CyclotomicMatrixMp& u, char letter) {
    for (int r = 0; r < 2; ++r) {
        Z16Mp& a = u.matrix[r][0];
        Z16Mp& b = u.matrix[r][1];
        switch (letter) {
            case 'H': {
                Z16Mp s = zeta_add(a, b), d = zeta_sub(a, b);
                a = std::move(s);
                b = std::move(d);
                break;
            }
            case 'S': b = times_zeta(b, 4); break;
            case 's': b = times_zeta(b, 12); break;
            case 'T': b = times_zeta(b, 2); break;
            case 't': b = times_zeta(b, 14); break;
            case 'Q': b = times_zeta(b, 1); break;
            case 'q': b = times_zeta(b, 15); break;
            case 'X': std::swap(a, b); break;
            case 'Y': {  // [[0, -i], [i, 0]]
                Z16Mp c0 = times_zeta(b, 4), c1 = times_zeta(a, 12);
                a = std::move(c0);
                b = std::move(c1);
                break;
            }
            case 'Z': b = times_zeta(b, 8); break;
            default: return false;
        }
    }
    if (letter == 'H') ++u.denominator_exponent;
    return true;
}

}  // namespace

void reduce_sqrt2(M2& m, int& k) { strip_sqrt2(m, k); }

CyclotomicMatrix compose(const CyclotomicMatrix& a, const CyclotomicMatrix& b) {
    CyclotomicMatrix out;
    out.matrix = matrix_mul(a.matrix, b.matrix);
    out.denominator_exponent = a.denominator_exponent + b.denominator_exponent;
    reduce_sqrt2(out.matrix, out.denominator_exponent);
    return out;
}

bool unit_norm(const Z16& u1, const Z16& u2, int K) {
    // Each product of two int64 coefficients is taken in 128 bits, below 2^120,
    // and a coefficient of the norm, a sum of 16 of them, stays below 2^124.
    // (zeta_mul over __int128 digits would multiply 128 x 128 bits, twice as slow.)
    __int128 norm[kZetaDegree] = {};
    for (const Z16* u : {&u1, &u2}) {
        const Z16 conjugate = zeta_conj(*u);
        for (int i = 0; i < kZetaDegree; ++i) {
            if ((*u)[i] == 0) continue;
            for (int j = 0; j < kZetaDegree; ++j) {
                const __int128 product = static_cast<__int128>((*u)[i]) * conjugate[j];
                if (i + j < kZetaDegree)
                    norm[i + j] += product;
                else
                    norm[i + j - kZetaDegree] -= product;  // zeta^8 = -1
            }
        }
    }
    // 2^K is 2^K times 1 and nothing at zeta, ..., zeta^7
    if (norm[0] != (static_cast<__int128>(1) << K)) return false;
    for (int t = 1; t < kZetaDegree; ++t)
        if (norm[t] != 0) return false;
    return true;
}

CyclotomicMatrix su2_element(const Z16& u1, const Z16& u2, int denominator_exponent) {
    CyclotomicMatrix out;
    out.matrix[0][0] = u1;
    out.matrix[0][1] = zeta_neg(zeta_conj(u2));
    out.matrix[1][0] = u2;
    out.matrix[1][1] = zeta_conj(u1);
    out.denominator_exponent = denominator_exponent;
    reduce_sqrt2(out.matrix, out.denominator_exponent);
    return out;
}

CyclotomicMatrix times_q_power(const CyclotomicMatrix& S, int j, const GateSet& gateset) {
    CyclotomicMatrix out = S;
    for (int r = 0; r < 2; ++r) out.matrix[r][1] = times_zeta(S.matrix[r][1], gateset.step * j);
    return out;
}

// --- words -----------------------------------------------------------------------

CyclotomicMatrix word_matrix(const std::string& word) {
    const CyclotomicMatrixMp u = word_matrix_mp(word);
    CyclotomicMatrix e;
    e.denominator_exponent = u.denominator_exponent;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int t = 0; t < kZetaDegree; ++t) {
                if (!u.matrix[i][j][t].fits_int64())
                    throw std::overflow_error("word_matrix: a coefficient passes 64 bits");
                e.matrix[i][j][t] = u.matrix[i][j][t].to_int64();
            }
    return e;
}

CyclotomicMatrixDigits word_matrix_digits(const std::string& word) {
    const CyclotomicMatrixMp u = word_matrix_mp(word);
    CyclotomicMatrixDigits out;
    out.denominator_exponent = u.denominator_exponent;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int t = 0; t < kZetaDegree; ++t)
                out.hex[16 * i + 8 * j + t] = u.matrix[i][j][t].str(16);
    return out;
}

std::invalid_argument unknown_letter(const std::string& word) {
    const std::string alphabet = "HSXYZTQstq";
    for (char letter : word)
        if (alphabet.find(letter) == std::string::npos)
            return std::invalid_argument("gate word: unknown letter '" + std::string(1, letter) +
                                         "' (the alphabet is " + alphabet + ")");
    return std::invalid_argument("gate word: unknown letter");
}

// --- CyclotomicMatrixMp ----------------------------------------------------------

bool word_matrix_mp(const std::string& word, CyclotomicMatrixMp& out) {
    CyclotomicMatrixMp u;
    u.matrix[0][0][0] = int_mp(int64_t{1});
    u.matrix[1][1][0] = int_mp(int64_t{1});
    for (char letter : word) {
        if (!apply_gate(u, letter)) return false;
        strip_sqrt2(u.matrix, u.denominator_exponent);
    }
    out = std::move(u);
    return true;
}

CyclotomicMatrixMp word_matrix_mp(const std::string& word) {
    CyclotomicMatrixMp u;
    if (!word_matrix_mp(word, u)) throw unknown_letter(word);
    return u;
}

CyclotomicMatrixMp widen(const M2& u, int k) {
    CyclotomicMatrixMp out;
    out.matrix = widened<int_mp>(u);
    out.denominator_exponent = k;
    strip_sqrt2(out.matrix, out.denominator_exponent);
    return out;
}

CyclotomicMatrixMp from_digits(const CyclotomicMatrixDigits& u) {
    CyclotomicMatrixMp out;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int t = 0; t < kZetaDegree; ++t)
                if (!int_mp::parse(u.hex[16 * i + 8 * j + t].c_str(), 16, out.matrix[i][j][t]))
                    throw std::invalid_argument("element: a coefficient that is not hexadecimal");
    out.denominator_exponent = u.denominator_exponent;
    strip_sqrt2(out.matrix, out.denominator_exponent);
    return out;
}

std::optional<int> phase_between(const CyclotomicMatrixMp& a, const CyclotomicMatrixMp& b) {
    if (a.denominator_exponent != b.denominator_exponent) return std::nullopt;
    for (int j = 0; j < 16; ++j) {
        bool same = true;
        for (int r = 0; r < 2 && same; ++r)
            for (int c = 0; c < 2 && same; ++c)
                same = a.matrix[r][c] == times_zeta(b.matrix[r][c], j);
        if (same) return j;
    }
    return std::nullopt;
}

bool equal_up_to_phase(const CyclotomicMatrixMp& a, const CyclotomicMatrixMp& b) {
    return phase_between(a, b).has_value();
}

}  // namespace quicksynth
