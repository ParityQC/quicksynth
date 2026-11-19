// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The diagonal entries u of the candidates near a z-rotation, one level at a
// time [Ross, Selinger, arXiv:1403.2975, sections 5 and 7].
//
// A candidate [[u, -t* w], [t, u* w]] / R, R = sqrt2^K, u, t in O_K = Z[zeta]
// (N = 4: Clifford+T, N = 8: Clifford+sqrt(T)), has SU(2) vector (u, t) / R
// up to the determinant phase w that the class j absorbs. Its distance to
// diag(a, a*) depends on u alone: d^2 R^2 = R^2 - Re(u a*)^2 with
// t t* = R^2 - u u*. So d <= eps with a solvable t needs u in the body
//   Re(u a*) >= R sqrt(1 - eps^2) at the identity place, a circular segment,
//   xi = 2^K - u u* totally >= 0: |sigma(u)| <= R at every place.
// The level (K, sub) holds the u of delta^sub O_K, delta = 1 - zeta.
//
// The segment, of depth h = R (1 - sqrt(1 - eps^2)) and half-chord R eps, lies
// in the ellipse with centre kSegmentCentre h above the chord and semi-axes
// kSegmentRadial h and kSegmentTangential R eps: its form is at most 0.99974
// on the segment for every eps in (0, 1). With the conjugate disks the body
// lies in an ellipsoid of form at most P = N / 2, enumerated at kBodyInflation
// P. The whitened basis is K-independent up to the factor 1 / R, so one exact
// LLL serves every K. Points come out of a Fincke-Pohst enumeration
// [Fincke, Pohst, Math. Comp. 44 (1985)] on the QR factor of the reduced
// form, with an inclusive margin, so rounding can only add points.
#pragma once

#include "quicksynth/data_structures/double_double.hpp"
#include "quicksynth/data_structures/number_field.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace quicksynth {

#pragma GCC visibility push(hidden)

inline constexpr double kSegmentCentre = 0.336;
inline constexpr double kSegmentRadial = 0.665;
inline constexpr double kSegmentTangential = 1.159;
inline constexpr double kBodyInflation = 1.01;
inline constexpr int64_t kNoPointLimit = INT64_MAX - 1;

template <size_t N>
class RotationBody {
  public:
    using Ring = RingElement<N, __int128>;
    struct Point {
        Ring u;
        // delta does not divide u / delta^sub: u is not in a coarser level.
        bool new_at_level = false;
    };

    // The body around a = (direction[0], direction[1]), normalised here; for
    // class j of a diagonal target D this is search_direction(D, j). `sub` in
    // [0, N / 2), eps in (0, 1); std::invalid_argument otherwise.
    // `reduction_scale` is the scale whose form is LLL-reduced; the points do
    // not depend on it.
    RotationBody(const Vec4& direction, int sub, double epsilon, int reduction_scale = 0);

    int sub() const { return sub_; }
    double epsilon() const { return epsilon_; }

    // Every u of delta^sub O_K in the enclosing ellipsoid at scale K, a superset
    // of the body, each once; the walk stops after limit + 1 points, so more
    // than `limit` means the list is cut. std::overflow_error if a coefficient
    // passes 127 bits.
    //
    // A direction within about eps of a field direction, such as a root of
    // unity, makes the radial form vanish on a sublattice of rank N / 2. The
    // ellipsoid then holds a slab of its cosets, about 2^K points per scale on
    // Clifford+T and 4^K on Clifford+sqrt(T) past a few scales below the first
    // useful level, while the body holds few of them.
    std::vector<Point> points(int K, int64_t limit = kNoPointLimit) const;

    // The points of points(K) in the same order, each new one (all with
    // `all`) passed to visit until it returns false or more than `limit`
    // points are walked. Returns the number walked, new or not.
    int64_t walk(int K, int64_t limit, bool all,
                 const std::function<bool(const Point&)>& visit) const;

    // Whether u is in the body at scale K: xi = 2^K - u u* is 0 or totally
    // positive (exact), and Re(u a*)^2 >= 2^K (1 - eps^2) with Re(u a*) > 0
    // (double-double).
    bool contains(const Ring& u, int K) const;

  private:
    int sub_;
    double epsilon_;
    ComplexDD a_;     // the direction, |a| = 1
    DD cos_squared_;  // 1 - eps^2
    int reduction_scale_;
    // Coefficients of delta^sub zeta^i, column i.
    std::array<std::array<int64_t, N>, N> basis_{};
    // x* / R: the ellipsoid's centre in the basis delta^sub zeta^i.
    std::array<DD, N> centre_{};
    // T with W T LLL-reduced, W the whitened basis, and T^-1, which passes 64
    // bits below eps = 1e-12; row-major.
    std::vector<int64_t> transform_;
    std::vector<int_mp> inverse_;
    // Upper-triangular factor of W T at reduction_scale, row-major.
    std::array<double, N * N> triangle_{};
};

#pragma GCC visibility pop

}  // namespace quicksynth
