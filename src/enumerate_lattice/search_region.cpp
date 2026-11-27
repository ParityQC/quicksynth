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
#include "quicksynth/data_structures/double_double.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace quicksynth {

namespace {

// Minimum cap-box half-thickness (~11 ulps of 1). Below eps = 1e-7 the cap is
// too thin for a safe double-rounding margin.
constexpr double kHalfAlongFloor = 2.4e-15;

// How far |y_p|^2 of a leaf may be from 1 at any place.
constexpr double kUnitTolerance = 1e-6;

// The ellipsoids below enclose the cap at the identity place: the points y of the unit 4-sphere
// with t = y . v >= s = sqrt(1 - eps^2), where t is y's coordinate along the target direction v
// and r = |y - t v| its distance across. Each is given by its semi-axes along v and across and
// its centre's coordinate t (CapEllipsoid).

// Allowance for the rounding of the walk's coordinates at the identity place, along v: 64 ulps
// of 1. The cylinder and the dome touch the cap, so a solution can lie on their boundary; without
// this pad the walk dropped planted solutions within 4 to 16 ulps of it.
constexpr double kCapPad = 64 * std::numeric_limits<double>::epsilon();

// The rounding allowance across v: kCapPad, and at least 1e-9 eps.
double across_pad(double epsilon) { return std::max(1e-9 * epsilon, kCapPad); }

// The box's ellipsoid: the Loewner ellipsoid of the box [s, 1] x [-eps, eps]^3 around the cap,
// centred mid-box, with semi-axes 2 h along v and 2 eps across, where h = half_along is the box's
// half-width along v (floored at kHalfAlongFloor). A solution is at most 1/2 inside it, so it
// needs no pad.
CapEllipsoid box_ellipsoid(double epsilon, double s, double half_along) {
    return {2.0 * half_along, 2.0 * epsilon, (1.0 + s) / 2};
}

// The cylinder's ellipsoid: the Loewner ellipsoid of the cylinder {s <= t <= 1, r <= eps} around
// the cap, centred mid-cylinder, with semi-axes 2 h along v and 2 eps / sqrt(3) across;
// 3^(3/2) = 5.2 times less volume than the box's. A solution can lie on its boundary, so it is
// padded: h by kCapPad, eps by across_pad().
CapEllipsoid cylinder_ellipsoid(double epsilon, double s, double half_along) {
    return {2.0 * (half_along + kCapPad), 2.0 / std::sqrt(3.0) * (epsilon + across_pad(epsilon)),
            (1.0 + s) / 2};
}

// The dome's ellipsoid. At depth u = 1 - t below the tip v, the cap lies in the paraboloid
// {r^2 / 2 <= u <= U}, U = 1 - s, since 1 - sqrt(1 - r^2) >= r^2 / 2. The paraboloid's Loewner
// ellipsoid is centred at depth 0.8 U, with semi-axes 0.8 U along v and (4 / sqrt(15)) sqrt(2 U)
// across, 9.1 times less volume than the box's. It touches the cap at the tip and along the rim
// (t = s, r = eps), so it is padded: the paraboloid starts kCapPad above the tip and is 2 kCapPad
// thicker, and the semi-axis across grows by the factor 1 + across_pad() / eps.
CapEllipsoid dome_ellipsoid(double epsilon, double s) {
    // U = 1 - s without cancellation, padded at both ends
    const double U = epsilon * epsilon / (1.0 + s) + 2.0 * kCapPad;
    return {0.8 * U,
            4.0 / std::sqrt(15.0) * std::sqrt(2.0 * U) * (1.0 + across_pad(epsilon) / epsilon),
            1.0 + kCapPad - 0.8 * U};
}

// Volume of a cap ellipsoid, in units of the unit 4-ball's: one semi-axis along v, three across.
double volume(const CapEllipsoid& cap) { return cap.along * cap.across * cap.across * cap.across; }

// The cap's ellipsoid the walk uses: the smallest of the box's, the cylinder's and the dome's,
// the box's on ties. That is the dome's down to eps ~ 1e-6, the cylinder's near 3e-7 and the box's
// below, where the pads outgrow the cap's thickness.
CapEllipsoid smallest_cap_ellipsoid(double epsilon, double s, double half_along) {
    CapEllipsoid best = box_ellipsoid(epsilon, s, half_along);
    for (const CapEllipsoid& cap :
         {cylinder_ellipsoid(epsilon, s, half_along), dome_ellipsoid(epsilon, s)})
        if (volume(cap) < volume(best)) best = cap;
    return best;
}

// The form of a cap ellipsoid around the target direction v, |S (y - y0)|^2 <= 1:
// S = I / across + (1 / along - 1 / across) v v^t and y0 = centre v.
struct CapForm {
    double S[4][4];
    double y0[4];
};

CapForm cap_form(const CapEllipsoid& cap, const Vec4& v) {
    const double across = 1.0 / cap.across;
    const double along = 1.0 / cap.along - across;
    CapForm out{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) out.S[i][j] = (along * v[i]) * v[j] + (i == j ? across : 0.0);
    for (int i = 0; i < 4; ++i) out.y0[i] = cap.centre * v[i];
    return out;
}

// The double-double set-up. Below eps ~ 3.5e-7 the pads above, ulps of 1 along v, are no longer
// small against the cap's thickness, and the dome stops being the smallest. They cover the
// set-up's rounding: S's entries reach 1 / along and the walk's centre B^-1 E^-1 y0 is about as
// large, both formed in doubles. Formed in double-double (precise_ellipsoid(), finish_setup()),
// the walk's coordinates are within about 1e-11 of exact and the walk encloses the dome of the
// cap's true thickness. The walk's own form has no slack (kFormSlack pads only the cap and
// conjugate forms), so the kPreciseCapPad pads below, about 1.5e-6 of its squared radius at the
// tip and 3.4e-6 at the rim, are the margin for that rounding.

// The precise dome's pads. eps grows by 1e-14, more than the double acceptance test
// (chord_distance) misjudges a distance by. Along v the paraboloid starts 1e-6 U above the tip
// and is 2e-6 U thicker, across the semi-axis grows by 1e-6: a solution exactly at the tip (an
// exact target) lies on the unpadded boundary, and the walk dropped some there.
constexpr double kAcceptanceFuzz = 1e-14;
constexpr double kPreciseCapPad = 1e-6;

// The least semi-axis along v of the precise dome, 2^-60; it binds below eps ~ 1.5e-9. to_ball's
// rows at the identity place and the walk's centre reach about 1 / along: at 2^-60 the LLL rows
// (48 fractional bits) stay below the 2^120 of FloatLll, past which the reduction falls back to
// exact integers, and the centre's double-double error, about 1e-32 / along ~ 1e-14, far below
// the kPreciseCapPad margin.
constexpr double kMinAlong = 0x1p-60;

// The dome's ellipsoid, as dome_ellipsoid(), for the double-double set-up: of the cap at
// eps + kAcceptanceFuzz, its thickness padded by kPreciseCapPad at both ends, and its semi-axis
// along v at least kMinAlong; a longer semi-axis about the same centre still holds the cap.
PreciseCap precise_dome(double epsilon) {
    const double e = epsilon + kAcceptanceFuzz;
    // U = 1 - s without cancellation
    const double U = e * e / (1.0 + std::sqrt(1.0 - e * e));
    const double pad = kPreciseCapPad * U, thickness = U + 2.0 * pad;
    return {std::max(0.8 * thickness, kMinAlong),
            4.0 / std::sqrt(15.0) * std::sqrt(2.0 * thickness) * (1.0 + kPreciseCapPad),
            0.8 * thickness - pad};
}

// sqrt2^K in double-double.
DD sqrt2_power_dd(int K) {
    const int whole = K >= 0 ? K / 2 : -((1 - K) / 2);  // floor(K / 2)
    return dd_ldexp(K - 2 * whole == 1 ? sqrt2_dd() : DD{1.0, 0.0}, whole);
}

// The Loewner ellipsoid of the double-double set-up, around `cap`: the identity place's rows of
// to_ball, S E with S = I / across + (1 / along - 1 / across) v v^t, and the centre E^-1 y0 with
// y0 = (1 - depth) v, formed in double-double for v normalised. Only E's identity rows
// M / sqrt2^K enter, M's entries correctly rounded: the other places' rows of to_ball are
// `embedding`'s (their forms are the unit balls), and E^-1 = M^t sqrt2^K / P as M^t M = P I.
LoewnerEllipsoid precise_ellipsoid(const GateSet& gateset, const PreciseCap& cap,
                                   const Vec4& target, int denominator_exponent,
                                   const std::vector<double>& embedding) {
    const int n = gateset.dim(), d = gateset.degree();
    // search_direction() leaves |v|^2 - 1 at a few ulps
    DD norm_sq{};
    for (const double x : target) norm_sq = norm_sq + two_prod(x, x);
    const DD norm = dd_sqrt(norm_sq);
    DD v[4];
    for (int i = 0; i < 4; ++i) v[i] = DD{target[i], 0.0} / norm;
    const DD across = DD{1.0, 0.0} / DD{cap.across, 0.0};
    const DD along = DD{1.0, 0.0} / DD{cap.along, 0.0} - across;
    DD S[4][4];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) S[i][j] = along * v[i] * v[j] + (i == j ? across : DD{});
    const DD height = DD{1.0, 0.0} - cap.depth;
    // M's rows at the identity place: Re and Im of zeta^(a step i) = exp(i pi 2 a step i / 16) in
    // rows 0, 1 at column i and in rows 2, 3 at column d + i
    std::vector<DD> M(static_cast<size_t>(kRowsPerPlace) * n);
    for (int i = 0; i < d; ++i) {
        const int k = 2 * gateset.places[0] * gateset.step * i;
        M[i] = M[static_cast<size_t>(2) * n + d + i] = cos_pi16(k);
        M[static_cast<size_t>(n) + i] = M[static_cast<size_t>(3) * n + d + i] = sin_pi16(k);
    }
    const DD down = sqrt2_power_dd(-denominator_exponent);
    const DD up =
        sqrt2_power_dd(denominator_exponent) / DD{static_cast<double>(gateset.num_places()), 0.0};
    LoewnerEllipsoid out{embedding, std::vector<double>(n),
                         std::vector<double>(static_cast<size_t>(kRowsPerPlace) * n),
                         std::vector<double>(n)};
    for (int j = 0; j < n; ++j) {
        DD c{};
        for (int r = 0; r < kRowsPerPlace; ++r)
            c = c + M[static_cast<size_t>(r) * n + j] * (height * v[r]);
        c = c * up;
        out.centre[j] = c.hi;
        out.centre_low[j] = c.lo;
        for (int i = 0; i < kRowsPerPlace; ++i) {
            DD x{};
            for (int r = 0; r < kRowsPerPlace; ++r)
                x = x + S[i][r] * M[static_cast<size_t>(r) * n + j];
            x = x * down;
            out.to_ball[static_cast<size_t>(i) * n + j] = x.hi;
            out.to_ball_low[static_cast<size_t>(i) * n + j] = x.lo;
        }
    }
    return out;
}

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
    // the ellipsoid the walk encloses the cap in; the estimates count the box's
    cap_ = smallest_cap_ellipsoid(epsilon, s_, half_along_);
    // where the pads make the dome not the smallest, the set-up is in double-double instead
    if (volume(cap_) < volume(dome_ellipsoid(epsilon, s_))) precise_cap_ = precise_dome(epsilon);
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
    if (precise_cap_)
        return precise_ellipsoid(gateset_, *precise_cap_, target, denominator_exponent, embedding);
    const int n = gateset_.dim();
    const std::vector<double> E_inv = embedding_inverse(gateset_, denominator_exponent);
    // S_cap is chosen such that the enclosed area is the cap's ellipsoid
    // (smallest_cap_ellipsoid()):
    // || S_cap (y- y0) ||^2 <= 1, with its centre y0 on the target vector v (cap_form()).
    const CapForm form = cap_form(cap_, target);
    const auto& S = form.S;
    const auto& y0 = form.y0;
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
    return {std::move(A), std::move(c0), {}, {}};
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
