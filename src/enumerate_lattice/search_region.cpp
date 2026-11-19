// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/enumerate_lattice/search_region.hpp"

#include "quicksynth/enumerate_lattice/lll.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

namespace quicksynth {

namespace {

// Minimum cap-box half-thickness (~11 ulps of 1). Below eps = 1e-7 the cap is
// too thin for a safe double-rounding margin.
constexpr double kHalfAlongFloor = 2.4e-15;

// How far |y_p|^2 of a leaf may be from 1 at any place.
constexpr double kUnitTolerance = 1e-6;

// pi^(dim/2) / (dim/2)!, the volume of the unit ball in dim dimensions.
double unit_ball(const GateSet& gateset) {
    double power = 1.0, factorial = 1.0;
    for (int i = 1; i <= gateset.dim() / 2; ++i) {
        power *= M_PI;
        factorial *= i;
    }
    return power / factorial;
}

// Volume ratio of the Loewner ellipsoid to the product region.
double ellipsoid_over_region(const GateSet& gateset) {
    // P^(2P) 2 (2P-2)! / ((2P)! (P-1)^(2P-2)), every factor exact as a double
    const int P = gateset.num_places();
    double p2p = 1.0, q = 1.0, f_small = 1.0, f_big = 1.0;
    for (int i = 0; i < 2 * P; ++i) p2p *= P;
    for (int i = 0; i < 2 * P - 2; ++i) q *= P - 1;
    for (int i = 2; i <= 2 * P - 2; ++i) f_small *= i;
    for (int i = 2; i <= 2 * P; ++i) f_big *= i;
    return p2p * 2.0 * f_small / f_big / q;
}

// Estimated points at which a scale is worth walking.
double first_scale_points(const GateSet& gateset) { return 1e-2 / ellipsoid_over_region(gateset); }

// Embedding at K = 0, from the powers of zeta.
std::vector<double> unscaled_embedding(const GateSet& gateset) {
    const int d = gateset.degree(), n = gateset.dim();
    std::vector<double> E(static_cast<size_t>(n) * n, 0.0);
    for (int p = 0; p < gateset.num_places(); ++p) {
        const int a = gateset.places[p];
        for (int i = 0; i < d; ++i) {
            const auto z = zeta_power(a * gateset.step * i);
            E[static_cast<size_t>(4 * p) * n + i] = z[0];
            E[static_cast<size_t>(4 * p + 1) * n + i] = z[1];
            E[static_cast<size_t>(4 * p + 2) * n + d + i] = z[0];
            E[static_cast<size_t>(4 * p + 3) * n + d + i] = z[1];
        }
    }
    return E;
}

// Multiplication by delta = 1 - zeta on the power-basis coefficients (u1, u2):
// (delta x)_i = x_i - x_{i-1}, and (delta x)_0 = x_0 + x_{d-1} since zeta^d = -1.
std::vector<int64_t> delta_matrix(const GateSet& gateset) {
    const int d = gateset.degree(), n = gateset.dim();
    std::vector<int64_t> M(static_cast<size_t>(n) * n, 0);
    for (int r = 0; r < n; ++r) {
        const int i = r % d, prev = r - i + (i + d - 1) % d;
        M[static_cast<size_t>(r) * n + r] = 1;
        M[static_cast<size_t>(r) * n + prev] = i == 0 ? 1 : -1;
    }
    return M;
}

// Multiplication by 1/delta on the power-basis coefficients: 1/delta = (1 + zeta
// + ... + zeta^(d-1)) / 2, since (1 - zeta)(1 + ... + zeta^(d-1)) = 1 - zeta^d = 2.
// Per block 1/2 on and below the diagonal, -1/2 above (zeta^d = -1).
std::vector<double> inverse_delta_matrix(const GateSet& gateset) {
    const int d = gateset.degree(), n = gateset.dim();
    std::vector<double> M(static_cast<size_t>(n) * n, 0.0);
    for (int r = 0; r < n; ++r) {
        const int i = r % d, block = r - i;
        for (int j = 0; j < d; ++j) M[static_cast<size_t>(r) * n + block + j] = j <= i ? 0.5 : -0.5;
    }
    return M;
}

// M^power, the identity at power 0 (sub = 0, the whole lattice). Exact for
// delta_matrix and inverse_delta_matrix, whose powers have entries k / 2^power.
template <typename T>
std::vector<T> matrix_power(int n, const std::vector<T>& M, int power) {
    std::vector<T> out(static_cast<size_t>(n) * n, T{0});
    for (int i = 0; i < n; ++i) out[static_cast<size_t>(i) * n + i] = T{1};
    for (int p = 0; p < power; ++p) out = product(n, M, out);
    return out;
}

}  // namespace

SearchRegion::SearchRegion(double epsilon, GateSet gateset)
    : epsilon_(epsilon), gateset_(std::move(gateset)) {
    s_ = std::sqrt(std::max(0.0, 1.0 - epsilon * epsilon));
    half_along_ = std::max((1.0 - s_) / 2, kHalfAlongFloor);
    // the volume of the 4-dimensional box that encloses the cap
    // 2 * half_along_ along the target, 2 * epsilon across it,
    // and 2 * epsilon in the other two directions
    cap_volume_ = (2 * half_along_) * (2 * epsilon) * (2 * epsilon) * (2 * epsilon);
    // 5% of the cap thickness, at least 1e-13: the leaf dot product rounds by
    // <= ~1.6e-14, so no element within epsilon is dropped
    s_min_ = s_ - std::max(0.05 * (1.0 - s_), 1e-13);
}

double SearchRegion::estimate_points(int K) const {
    // Ellipsoid volume / covolume = unit_ball prod(2h) 2^(K dim/2), since
    // M^t M = P I (unscaled embedding) and det S_cap = 1 / prod(2h).
    return std::ldexp(unit_ball(gateset_) * cap_volume_, K * gateset_.dim() / 2) /
           ellipsoid_over_region(gateset_);
}

int SearchRegion::first_scale() const {
    int K = 0;
    while (estimate_points(K) < first_scale_points(gateset_)) ++K;
    return K;
}

int SearchRegion::first_level(int L) const {
    int l = 0;
    for (;; ++l) {
        const LatticeLevel at = level_at(l, L);
        if (estimate_points(at.K, at.sub) >= first_scale_points(gateset_)) return l;
    }
}

int levels_per_scale(const GateSet& gateset) { return gateset.degree() / 2; }

LatticeLevel level_at(int l, int L) {
    const int K = (l + L - 1) / L;
    return {K, L * K - l};
}

LoewnerEllipsoid SearchRegion::ellipsoid(const Vec4& target, int denominator_exponent,
                                         const std::vector<double>& embedding) const {
    const int n = gateset_.dim();
    const std::vector<double> E_inv = embedding_inverse(gateset_, denominator_exponent);
    // S_cap is chosen such that the enclosed area is the smallest possible elipsoid that encloses
    // the cap box: || S_cap (y- y0) ||^2 <= 1. The right side is 1 so S_cap has to contain the
    // dimensions of the box: 2eps perpendicular to v (the target vector) and 2h along v, where
    // h = (1-s)/2 is the half-thickness of the cap box along v :
    // S_cap = I/(2 eps) + (1/(2 h0) - 1/(2 eps)) v v^t of the region
    // and its centre y0 = (1 + s)/2 v with the target vector v.
    const double across = 1.0 / (2.0 * epsilon_);
    const double along = 1.0 / (2.0 * half_along_) - across;
    // set the elements of S_cap
    double S[4][4];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            S[i][j] = (along * target[i]) * target[j] + (i == j ? across : 0.0);
    const double c = (1.0 + s_) / 2;
    const double y0[4] = {c * target[0], c * target[1], c * target[2], c * target[3]};
    // A = S E and c0 = E^-1 y0, the sums in a fixed order: S is the identity for all places but the
    // the first place (first four rows)
    std::vector<double> A = embedding;
    // Just update the first four rows of A
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < n; ++j) {
            double a = S[i][0] * embedding[j];
            a += S[i][1] * embedding[static_cast<size_t>(n) + j];
            a += S[i][2] * embedding[static_cast<size_t>(2 * n) + j];
            a += S[i][3] * embedding[static_cast<size_t>(3 * n) + j];
            A[static_cast<size_t>(i) * n + j] = a;
        }
    // The centers are 0 everywhere but for the first four rows (first place).
    std::vector<double> c0(n, 0.0);
    for (int j = 0; j < n; ++j) {
        const double* row = E_inv.data() + static_cast<size_t>(j) * n;
        double x = row[0] * y0[0] + row[1] * y0[1];
        x += row[2] * y0[2];
        x += row[3] * y0[3];
        c0[j] = x;
    }
    return {std::move(A), std::move(c0)};
}

WalkLattice SearchRegion::walk_lattice(const std::vector<double>& E) const {
    return {gateset_.dim(), E.data(), s_min_, kUnitTolerance};
}

std::vector<double> embedding(const GateSet& gateset, int K) {
    std::vector<double> E = unscaled_embedding(gateset);
    const double s = sqrt2_power(-K);
    for (double& x : E) x *= s;
    return E;
}

std::vector<double> embedding_inverse(const GateSet& gateset, int K) {
    // uses MᵗM = P·I (the columns are orthogonal; this is the orthogonality of the characters).
    // So E⁻¹ = Mᵗ·√2^K / P, with no matrix inversion.
    const int n = gateset.dim();
    const std::vector<double> M = unscaled_embedding(gateset);
    const double s = sqrt2_power(K);
    const double P = gateset.num_places();
    std::vector<double> out(static_cast<size_t>(n) * n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            out[static_cast<size_t>(j) * n + i] = (M[static_cast<size_t>(i) * n + j] * s) / P;
    return out;
}

std::vector<int64_t> delta_power(const GateSet& gateset, int sub) {
    return matrix_power(gateset.dim(), delta_matrix(gateset), sub);
}

std::vector<double> inverse_delta_power(const GateSet& gateset, int sub) {
    return matrix_power(gateset.dim(), inverse_delta_matrix(gateset), sub);
}

}  // namespace quicksynth
