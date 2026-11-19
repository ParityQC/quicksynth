// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The search region for one precision.
//
// A unitary (u1, u2) / sqrt2^K is a lattice point in R^dim, in the cap u . v >= sqrt(1 - eps^2)
// at the identity place and on the unit 4-sphere at the others. The cap is boxed (eps^2/2
// along v, eps across), the product region enclosed in its Loewner ellipsoid. LLL inputs use a
// fixed rounding order, so they are reproducible.
#pragma once

#include "quicksynth/enumerate_lattice/walk.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

namespace quicksynth {

// Loewner ellipsoid of the region around a direction v at a scale K:
// |to_ball (c - centre)|^2 <= walk_radius_sq() for the power-basis coefficients c.
struct LoewnerEllipsoid {
    // n x n, S E (E: embedding at K): power-basis coefficients to the coordinates in which
    // the ellipsoid is the ball of radius sqrt(walk_radius_sq()); the walk and LLL measure
    // lengths there
    std::vector<double> to_ball;
    std::vector<double> centre;  // n, E^-1 y0
};

class SearchRegion {
  public:
    SearchRegion(double epsilon, GateSet gateset);

    double epsilon() const { return epsilon_; }
    const GateSet& gateset() const { return gateset_; }

    // Expected lattice points at scale K (direction-independent), and the
    // first scale worth walking.
    double estimate_points(int K) const;
    int first_scale() const;
    // Finer levels: level l of L per scale is the sublattice delta^sub Z[zeta]^2 at
    // scale K = ceil(l / L), sub = L K - l, delta = 1 - zeta (index 4 per level).
    double estimate_points(int K, int sub) const {
        return std::ldexp(estimate_points(K), -2 * sub);
    }
    int first_level(int L) const;

    // Minimum leaf dot product with the target; a margin below the cap keeps
    // every element within epsilon.
    double s_min() const { return s_min_; }
    // Squared radii of the walk's ellipsoid and of the conjugate block.
    double walk_radius_sq() const { return gateset_.num_places(); }
    double conj_radius_sq() const { return gateset_.num_places() - 1; }
    // Per-place conjugate forms; 0 with a single conjugate place.
    int conjugate_forms() const {
        return gateset_.num_places() > 2 ? gateset_.num_places() - 1 : 0;
    }

    // Loewner ellipsoid around `target`, the target's direction v, at the scale
    // K = `denominator_exponent`, whose embedding is `embedding`.
    LoewnerEllipsoid ellipsoid(const Vec4& target, int denominator_exponent,
                               const std::vector<double>& embedding) const;

    // The lattice at the scale whose embedding is E, with this region's leaf
    // tests: the cap y_0 . v >= s_min() and |y_p|^2 = 1 at every place, up to a
    // tolerance. A view of E, which must outlive it.
    WalkLattice walk_lattice(const std::vector<double>& E) const;

  private:
    double epsilon_;
    GateSet gateset_;
    double s_;           // sqrt(1 - eps^2)
    double half_along_;  // box half-width along the target
    double cap_volume_;  // box volume prod(2 h)
    double s_min_;
};

// The embedding's rows per place: Re and Im of u1, then of u2.
constexpr int kRowsPerPlace = 4;

// Embedding at scale K (dim x dim, row-major, scaled by 2^(-K/2)) and its inverse. Rows:
// Re/Im of u1, u2 per place; columns: ring-basis coefficients of u1, u2.
std::vector<double> embedding(const GateSet& gateset, int K);
std::vector<double> embedding_inverse(const GateSet& gateset, int K);

// Levels per sqrt2-scale: degree / 2 (sqrt2 = delta^(degree/2) x unit).
int levels_per_scale(const GateSet& gateset);
struct LatticeLevel {
    int K = 0;    // sqrt2-scale
    int sub = 0;  // delta exponent of the sublattice, 0 .. L - 1
};
LatticeLevel level_at(int l, int L);

// The sublattice delta^sub Z[zeta]^2 of a level: M_sub, multiplication by delta^sub on
// the power-basis coefficients (u1, u2), and its inverse, whose entries k / 2^sub
// doubles hold exactly.
std::vector<int64_t> delta_power(const GateSet& gateset, int sub);
std::vector<double> inverse_delta_power(const GateSet& gateset, int sub);

}  // namespace quicksynth
