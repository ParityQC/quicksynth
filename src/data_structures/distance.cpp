// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/distance.hpp"

#include "quicksynth/data_structures/int_mp.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace quicksynth {

namespace {

// The widest coefficient summed unscaled: eight terms below 2^1000 stay far
// from a double's overflow at 2^1024.
constexpr int kWidestUnscaled = 1000;

// The matrix of an element of any size: its coefficients as doubles, summed as
// matrix_unitary sums those of an int64 element (the same bits where they fit).
// A coefficient past 64 bits is truncated to a double (int_mp::truncated_frexp).
// Past kWidestUnscaled bits every coefficient is scaled by 2^-shift and the
// scale 2^(-k/2) by 2^shift, by taking k - 2 shift. Both are exact, so the
// bits are those of the unscaled sum wherever that is finite, and a word of
// any length gives a finite matrix rather than inf and then NaN.
Unitary matrix_of(const CyclotomicMatrixMp& u) {
    int widest = 0;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int t = 0; t < kZetaDegree; ++t)
                widest = std::max(widest, u.matrix[i][j][t].bits());
    const int shift = widest > kWidestUnscaled ? widest - kWidestUnscaled : 0;
    std::array<double, 32> coefficients;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int t = 0; t < kZetaDegree; ++t) {
                const int_mp& c = u.matrix[i][j][t];
                double x;
                if (c.fits_int64()) {
                    x = std::ldexp(static_cast<double>(c.to_int64()), -shift);
                } else {
                    int exponent;
                    const double mantissa = c.truncated_frexp(exponent);
                    x = std::ldexp(mantissa, exponent - shift);
                }
                coefficients[16 * i + 8 * j + t] = x;
            }
    return matrix_unitary(u.denominator_exponent - 2 * shift, coefficients);
}

double distance(const CyclotomicMatrixMp& w, const Unitary& target) {
    for (double x : target)
        if (!std::isfinite(x))
            throw std::invalid_argument("precise_distance: a target entry is not finite");
    // ||U^dagger U - I||, Frobenius
    const Unitary gram = matmul(dagger(target), target);
    double defect = 0.0;
    for (int i = 0; i < 8; ++i) {
        const double x = gram[i] - (i == 0 || i == 6 ? 1.0 : 0.0);
        defect += x * x;
    }
    defect = std::sqrt(defect);
    if (!(defect <= kTargetUnitarity)) {
        char text[160];
        std::snprintf(text, sizeof text,
                      "precise_distance: the target is not unitary: ||U^dagger U - I|| = %.1e "
                      "exceeds %g; take it to its nearest unitary first",
                      defect, kTargetUnitarity);
        throw std::invalid_argument(text);
    }
    return chord_distance(su2_vector(matrix_of(w)), su2_vector(target));
}

}  // namespace

double precise_distance(const std::string& word, const Unitary& target) {
    return distance(word_matrix_mp(word), target);
}

}  // namespace quicksynth
