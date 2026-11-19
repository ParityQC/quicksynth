// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Streaming Fincke-Pohst enumeration with an in-loop unitarity/cap filter.
//
// n below is the lattice dimension (`lattice_dimension`): 16 on
// Clifford+sqrt(T), 8 on Clifford+T, the integer coefficients of a candidate's
// first column.
//
// Enumerates x in Z^n with |R (x - c)|^2 <= radius_sq, R upper triangular, from the last
// coordinate down (`lvl`). A leaf is kept if y = E (T x) has unit norm within tol per
// 4-block and y1 . v >= s_min. T x is exact in int64; forming E T in doubles would lose
// ~1e-7 at eps = 1e-5 (T entries ~1e10).
#pragma once

#include <cstdint>
#include <vector>

namespace quicksynth {

// The region |upper_triangular (x - center) + offset|^2 <= radius_sq of R^n,
// the Ellipsoid of the walk. `upper_triangular` is n x n, so row i of
// upper_triangular (x - center) involves only x_i, ..., x_{n-1}: fixing the
// last coordinates leaves an interval for the next one. A row that is
// identically zero constrains nothing at its coordinate, so a form of rank r
// only bites at the last r coordinates. `center` may be null, meaning the
// origin; `offset` may be null, meaning none.
//
// A walk takes several of them, its forms, and enumerates their intersection.
// The first is the walk's own ellipsoid: full rank, centred, without an offset,
// and its rows bound the coordinates the walk fixes. On the scan it is the
// Loewner relaxation of the search region -- a product of one cap and three unit
// balls -- and much larger than it; passing the cap block (rank 4, radius 1), the
// conjugate block (rank 12, radius 3) and one form per conjugate place (rank 4,
// radius 1) as further forms prunes the walk towards the search region.
//
// The value of row `i` is
// `sum_{j>=i} upper_triangular[i][j] (x_j - center_j) + offset[i]`. The
// offset `R c`, computed exactly, evaluates an origin-centred form about the centre c with
// O(1) operands: `|R (x - c) + R c|^2 = |R x|^2`. Otherwise coordinates reach 1e14 at 1e-8.
struct Ellipsoid {
    const double* upper_triangular = nullptr;
    const double* center = nullptr;
    const double* offset = nullptr;
    double radius_sq = 0.0;
};

// The lattice at one scale and the tests a leaf must pass, shared by every walk
// of a batch: the embedding y = E c of the power-basis coefficients c (n x n,
// row-major, 4 rows per place), the cap y1 . v >= s_min at the identity place
// and the unitarity test | |y_p|^2 - 1 | < tol at every place.
struct WalkLattice {
    int lattice_dimension = 0;
    const double* embedding = nullptr;
    double s_min = 0.0;
    double tol = 0.0;
};

// One walk on that lattice: its forms, the walk's ellipsoid first (at most
// eight), the transform from the reduced coordinates to the power basis (n x n
// int64, row-major), the cap's direction v (4) and `base` (n int64), T times the
// integer origin the walk is centred on: the coefficients reported are T x +
// base, so the caller can keep the enumerated coordinates of order one. A null
// base means the origin.
struct WalkInputs {
    std::vector<Ellipsoid> forms;
    const int64_t* transform = nullptr;
    const double* v = nullptr;
    const int64_t* base = nullptr;
};

struct WalkResult {
    std::vector<int64_t> coefficients;  // flat, n per survivor, power basis
    int64_t visited = 0;                // leaves of the enumeration tree
};

// Every point of the lattice in all the forms of `inputs`, filtered at the
// leaves as above. `max_points` > 0 caps the leaves.
//
// With threads > 1, the top one to six coordinates are split into work items; output is in
// item order. A positive max_points is shared, so truncation is not deterministic.
WalkResult walk(const WalkLattice& lattice, const WalkInputs& inputs, int64_t max_points,
                int threads);

// Independent walks on one lattice run over `threads` workers, one whole walk
// at a time: the tree of a small walk holds less work per thread than handing
// it out costs, so a batch of small walks is split by walk instead.
// `max_points` caps each walk's leaves (0 = no cap). It bounds the time; the volume
// estimate can be off by 10^4 on a rotated region.
std::vector<WalkResult> walk_batch(const WalkLattice& lattice, const std::vector<WalkInputs>& walks,
                                   int64_t max_points, int threads);

}  // namespace quicksynth
