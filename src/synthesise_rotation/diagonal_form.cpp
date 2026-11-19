// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/synthesise_rotation/diagonal_form.hpp"

#include "quicksynth/data_structures/clifford.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace quicksynth {

namespace {

// X, Y, Z as complex matrices (su2.hpp layout)
const std::array<Unitary, 3> kPaulis{{
    {0, 0, 1, 0, 1, 0, 0, 0},
    {0, 0, 0, -1, 0, 1, 0, 0},
    {1, 0, 0, 0, 0, 0, -1, 0},
}};

// The 24 Cliffords of spelling_tables().words as unitaries, in that order.
const std::vector<Unitary>& clifford_unitaries() {
    static const std::vector<Unitary> out = [] {
        std::vector<Unitary> u;
        for (const auto& entry : spelling_tables().words)
            u.push_back(matrix_unitary(word_matrix(entry.second)));
        return u;
    }();
    return out;
}

// First Clifford in word order with Bloch-image entry (row, column) equal to
// `value`.
size_t first_clifford(int row, int column, int value) {
    const auto& words = spelling_tables().words;
    for (size_t c = 0; c < words.size(); ++c)
        if (words[c].first[row][column] == value) return c;
    return words.size();  // unreachable: the Cliffords act transitively on signed axes
}

}  // namespace

BlochRotation bloch_rotation(const Unitary& U) {
    const Unitary ud = dagger(U);
    BlochRotation r{};
    for (int k = 0; k < 3; ++k) {
        const Unitary image = matmul(matmul(U, kPaulis[k]), ud);
        for (int j = 0; j < 3; ++j) {
            const Unitary m = matmul(kPaulis[j], image);
            r[j][k] = 0.5 * (m[0] + m[6]);
        }
    }
    return r;
}

std::optional<DiagonalForm> diagonal_form(const Unitary& U, double tolerance) {
    const BlochRotation r = bloch_rotation(U);
    const Vec4 target = su2_vector(U);
    const auto& words = spelling_tables().words;
    const auto& unitaries = clifford_unitaries();
    std::optional<DiagonalForm> best;
    for (int i : {2, 0, 1})
        for (int j : {2, 0, 1}) {
            if (!(std::fabs(r[j][i]) >= 0.5)) continue;
            const int s = r[j][i] > 0 ? 1 : -1;
            const size_t c2 = first_clifford(2, i, 1);  // R(C2) e_i = e_z
            const size_t c1 = first_clifford(j, 2, s);  // R(C1) e_z = s e_j
            const Unitary inner = matmul(matmul(dagger(unitaries[c1]), U), dagger(unitaries[c2]));
            const Vec4 v = su2_vector(inner);
            const double norm = std::hypot(v[0], v[1]);
            if (!(norm > 0)) continue;
            const double zr = v[0] / norm, zi = v[1] / norm;
            const Unitary d{zr, zi, 0, 0, 0, 0, zr, -zi};
            const Unitary w = matmul(matmul(unitaries[c1], d), unitaries[c2]);
            const double residual = chord_distance(target, su2_vector(w));
            if (best && !(residual < best->residual)) continue;
            best = DiagonalForm{words[c1].second, words[c2].second, d, residual};
        }
    if (!best || !(best->residual <= tolerance)) return std::nullopt;
    return best;
}

}  // namespace quicksynth
