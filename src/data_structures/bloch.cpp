// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/bloch.hpp"

#include "quicksynth/data_structures/int_mp.hpp"

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>

namespace quicksynth {

namespace {

M2 pauli(int which) {
    M2 m{};
    if (which == 0) {  // X
        m[0][1][0] = 1;
        m[1][0][0] = 1;
    } else if (which == 1) {  // Y = [[0, -i], [i, 0]], i = zeta^4
        m[0][1][4] = -1;
        m[1][0][4] = 1;
    } else {  // Z
        m[0][0][0] = 1;
        m[1][1][0] = -1;
    }
    return m;
}

const M2 kPauli[3] = {pauli(0), pauli(1), pauli(2)};

// Numerator N of the SO(3) image N / 2^(k+1) of U = M / sqrt2^k, row-major
// over Z[lambda]. Entries are real for any M; std::logic_error otherwise.
template <class D>
LambdaMatrix<D> so3_numerator(const M2Of<D>& M) {
    const M2Of<D> Md = matrix_dagger(M);
    LambdaMatrix<D> n;
    for (int j = 0; j < 3; ++j) {
        const M2Of<D> xj = matrix_mul(matrix_mul(M, widened<D>(kPauli[j])), Md);
        for (int i = 0; i < 3; ++i) {
            const M2Of<D> sx = matrix_mul(widened<D>(kPauli[i]), xj);
            const Z16Of<D> tr = zeta_add(sx[0][0], sx[1][1]);
            if (zeta_conj(tr) != tr)
                throw std::logic_error("Bloch image: a trace that is not real");
            // to basis 1, lambda, lambda^2, lambda^3, peeling from the top power:
            // lambda = (0,1,0,0,0,0,0,-1), lambda^2 = (2,0,1,0,0,0,-1,0),
            // lambda^3 = (0,3,0,1,0,-1,0,-3) in powers of zeta
            at(n, i, j) = {tr[0] - tr[2] * 2, tr[1] - tr[3] * 3, tr[2], tr[3]};
        }
    }
    return n;
}

// Least exponent of the image N u^(k+1) / lambda^(4(k+1) - stripped).
int64_t image_exponent(int k, int64_t stripped) {
    return 4 * (static_cast<int64_t>(k) + 1) - stripped;
}

// N in __int128 digits; none if a coefficient of M reaches 2^59, which keeps N
// below 2^127. int64 would overflow near Bloch exponent 48 (Clifford+T at 1e-5),
// as each multiplication by u grows a digit up to ninefold.
std::optional<LambdaMatrix<__int128>> so3_numerator_128(const M2& M) {
    constexpr int64_t kLimit = int64_t{1} << 59;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int64_t c : M[i][j])
                if (c >= kLimit || c <= -kLimit) return std::nullopt;
    return so3_numerator(widened<__int128>(M));
}

constexpr __int128 kWideGuard = __int128{1} << 100;

bool wide_too_big(const Lambda<__int128>& d) {
    for (int t = 0; t < 4; ++t)
        if (d[t] > kWideGuard || d[t] < -kWideGuard) return true;
    return false;
}

int64_t least_exponent(const CyclotomicMatrixMp& u) {
    LambdaMatrix<int_mp> n = so3_numerator(u.matrix);
    return image_exponent(u.denominator_exponent, strip_lambda(n));
}

BlochDigits digits_of(const CyclotomicMatrixMp& u) {
    const BlochImage<int_mp> image = bloch_image(u);
    BlochDigits out;
    out.exponent = image.denominator_exponent;
    for (int r = 0; r < 3; ++r)
        for (int col = 0; col < 3; ++col)
            for (int p = 0; p < 4; ++p) at(out, r, col, p) = at(image.m, r, col)[p].str(10);
    return out;
}

}  // namespace

std::optional<int64_t> try_bloch_exponent_128(const M2& u, int k) {
    std::optional<LambdaMatrix<__int128>> n = so3_numerator_128(u);
    if (!n) return std::nullopt;
    return image_exponent(k, strip_lambda(*n));
}

std::optional<BlochImage128> try_bloch_image_128(const M2& u, int k) {
    std::optional<LambdaMatrix<__int128>> n = so3_numerator_128(u);
    if (!n) return std::nullopt;
    BlochImage128 out;
    out.denominator_exponent = image_exponent(k, strip_lambda(*n));
    out.m = std::move(*n);
    for (int step = 0; step <= k; ++step)
        for (auto& cell : out.m) {
            cell = times_u(cell);
            if (wide_too_big(cell)) return std::nullopt;
        }
    return out;
}

BlochDigits bloch_digits(const M2& u, int k) { return digits_of(widen(u, k)); }

BlochDigits bloch_digits(const CyclotomicMatrixDigits& u) { return digits_of(from_digits(u)); }

int64_t bloch_exponent(const M2& u, int k) { return least_exponent(widen(u, k)); }

int64_t bloch_exponent(const CyclotomicMatrixDigits& u) { return least_exponent(from_digits(u)); }

BlochImage<int_mp> bloch_image(const CyclotomicMatrixMp& u) {
    BlochImage<int_mp> out;
    out.m = so3_numerator(u.matrix);
    out.denominator_exponent = image_exponent(u.denominator_exponent, strip_lambda(out.m));
    // times u^(k+1); see image_exponent
    for (int step = 0; step <= u.denominator_exponent; ++step)
        for (Lambda<int_mp>& d : out.m) d = times_u(d);
    return out;
}

}  // namespace quicksynth
