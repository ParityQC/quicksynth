// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/enumerate_lattice/walk.hpp"

#include "quicksynth/support/parallel.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

namespace quicksynth {

namespace {

constexpr double kSlack = 1e-9;      // a coordinate's term may exceed its budget by this
constexpr double kRound = 1e-12;     // widen an interval before rounding to integers
constexpr int64_t kLeafBlock = 256;  // most leaves a walker claims from a shared budget at once

// Index range of the walk. Bounds are clamped: casting a larger double to int64 is undefined
// (INT64_MIN on x86 empties the interval). Round-off QR diagonals produce such bounds; NaN
// means no bound.
constexpr int64_t kIndexMax = std::numeric_limits<int64_t>::max() / 4;

int64_t lower_index(double v) {
    if (!(v > -static_cast<double>(kIndexMax))) return -kIndexMax;
    return v < static_cast<double>(kIndexMax) ? static_cast<int64_t>(v) : kIndexMax;
}

int64_t upper_index(double v) {
    if (!(v < static_cast<double>(kIndexMax))) return kIndexMax;
    return v > -static_cast<double>(kIndexMax) ? static_cast<int64_t>(v) : -kIndexMax;
}

// a + b * c, wrapping at 2^64. set_coordinate multiplies x by a column of T, and
// x lies near kIndexMax only at a level the forms leave unbounded (the clamp in
// enter_level()): the product then overflows, which in int64 is undefined.
// Unsigned arithmetic wraps by definition, and the conversion back is two's
// complement.
inline int64_t wrapping_mul_add(int64_t a, int64_t b, int64_t c) {
    return static_cast<int64_t>(static_cast<uint64_t>(a) +
                                static_cast<uint64_t>(b) * static_cast<uint64_t>(c));
}

// T by column, out[lvl * n + r] = T[r][lvl], so a column update is contiguous.
std::vector<int64_t> transpose(int lattice_dimension, const int64_t* T) {
    std::vector<int64_t> out(static_cast<size_t>(lattice_dimension) * lattice_dimension);
    for (int r = 0; r < lattice_dimension; ++r)
        for (int l = 0; l < lattice_dimension; ++l)
            out[static_cast<size_t>(l) * lattice_dimension + r] =
                T[static_cast<size_t>(r) * lattice_dimension + l];
    return out;
}

// The walk's quadratic forms |R (x - c) + g|^2 <= radius_sq (the Ellipsoids
// of walk.hpp, R upper triangular) side by side in `NL` lanes: form f in lane f,
// the unused lanes zero, and so is a null centre or offset. A coordinate's
// partial sums are then one loop over j with the forms in parallel lanes, each
// lane adding its products in the order a single form would, so the sums are
// the same numbers.
//
// The layout is this class's alone: the walker asks for the NL lanes of one
// entry of R, or of one element of a vector, and for single values by form.
// Indices are those of the matrices: R by row and column, the offset and the
// diagonal by row, the centre by column (c_j is subtracted from x_j, which
// column j multiplies).
//
// What the forms are, as the enumerator builds them (walk_of in
// enumerator.cpp). With P places and y_p the candidate's four reals (Re u,
// Im u, Re t, Im t) at place p, the search region is the cap at the identity place
// times the unit ball at every other place. It is no ellipsoid, so the walk
// carries ellipsoids around it:
//
//   lane 0     the walk's ellipsoid, |S (y_0 - y_cap)|^2 + sum_{p >= 1} |y_p|^2
//              <= P, S taking the cap to the unit ball (SearchRegion::ellipsoid). Full rank,
//              so it gives every coordinate an interval: what makes the walk
//              possible, and much larger than the search region.
//   lane 1     the cap: its first term alone, <= 1. Rank 4, highest row 3, the
//              same centre.
//   lane 2     the conjugate block: the sum alone, <= P - 1. Rank n - 4, highest
//              row n - 5, centred at the origin and so with an offset.
//   lanes 3..  one per conjugate place when P > 2: |y_p|^2 <= 1. Rank 4,
//              highest row 3, with an offset.
//
// The walk's ellipsoid gives the cap and the conjugate places one budget,
// which a point can shift from one to the other; the others give each its
// own. Every element of the scale in the cap satisfies them all (each place
// has norm exactly 1), so they cut only points the leaf test would reject;
// each radius has 1e-9 of slack for rounding. Clifford+T (P = 2, n = 8) has
// lanes 0-2 on NL = 4, a per-place form there would repeat the conjugate block;
// Clifford+sqrt(T) (P = 4, n = 16) has lanes 0-5 on NL = 8.
template <int NL>
class Forms {
  public:
    Forms(int lattice_dimension, const std::vector<Ellipsoid>& ellipsoids)
        : lattice_dimension_(lattice_dimension),
          count_(static_cast<int>(ellipsoids.size())),
          upper_triangular_(static_cast<size_t>(lattice_dimension) * lattice_dimension * NL, 0.0),
          center_(static_cast<size_t>(lattice_dimension) * NL, 0.0),
          offset_(static_cast<size_t>(lattice_dimension) * NL, 0.0),
          inverse_diagonal_(static_cast<size_t>(lattice_dimension) * NL, 0.0) {
        if (count_ < 1 || ellipsoids[0].offset != nullptr)
            throw std::invalid_argument("walk: the walk's ellipsoid first, no offset");
        for (int f = 0; f < count_; ++f) set_lane(f, ellipsoids[f]);
    }

    int count() const { return count_; }
    double radius_sq(int f) const { return radius_sq_[f]; }  // 0 in an unused lane

    // The NL lanes of R_f[row][column], of c_f[column], of g_f[row] and of
    // 1 / R_f[row][row] -- 0 where the diagonal is not positive: the form then
    // gives no interval at that coordinate, only the test of inside_every_form().
    // Elsewhere a coordinate costs multiplications and no divisions.
    const double* upper_triangular(int row, int column) const {
        return upper_triangular_.data() + entry_index(row, column);
    }
    const double* centers(int column) const { return center_.data() + vector_index(column); }
    const double* offsets(int row) const { return offset_.data() + vector_index(row); }
    const double* inverse_diagonals(int row) const {
        return inverse_diagonal_.data() + vector_index(row);
    }

    // The NL lanes of R_f[row][row].
    const double* diagonals(int row) const { return upper_triangular(row, row); }

  private:
    // Where the lanes of entry (row, column) of R begin, and those of element
    // `index` of a vector.
    size_t entry_index(int row, int column) const {
        // Each matrix position takes up NL slots, so the group for a position starts at position *
        // NL
        return (static_cast<size_t>(row) * lattice_dimension_ + column) * NL;
    }
    // Each vector position takes up NL slots, so the group for a position starts at position * NL
    static size_t vector_index(int index) { return static_cast<size_t>(index) * NL; }

    // Form f into lane f.
    void set_lane(int f, const Ellipsoid& ellipsoid) {
        radius_sq_[f] = ellipsoid.radius_sq;
        for (int row = 0; row < lattice_dimension_; ++row) {
            const double* row_entries =
                ellipsoid.upper_triangular + static_cast<size_t>(row) * lattice_dimension_;
            for (int column = 0; column < lattice_dimension_; ++column)
                upper_triangular_[entry_index(row, column) + f] = row_entries[column];
            const size_t slot = vector_index(row) + f;
            if (ellipsoid.offset != nullptr) offset_[slot] = ellipsoid.offset[row];
            if (row_entries[row] > 0.0) inverse_diagonal_[slot] = 1.0 / row_entries[row];
        }
        if (ellipsoid.center != nullptr)
            for (int column = 0; column < lattice_dimension_; ++column)
                center_[vector_index(column) + f] = ellipsoid.center[column];
    }

    int lattice_dimension_;
    int count_;
    std::vector<double> upper_triangular_, center_, offset_, inverse_diagonal_;
    std::array<double, NL> radius_sq_{};
};

// The candidate at one place: (Re u, Im u, Re t, Im t).
using PlaceVector = std::array<double, 4>;

// The test a complete point must pass, on its power-basis coefficients c.
//
// The embedding E is lattice_dimension x lattice_dimension, row-major. Its
// columns are the coefficients of the candidate's first column (u, t): the
// first half those of u, the second half those of t. Its rows come in blocks
// of four, one block per place p (a Galois conjugate zeta -> zeta^a; p = 0 is
// the identity):
//
//     row 4p + 0: Re u at p      row 4p + 2: Re t at p
//     row 4p + 1: Im u at p      row 4p + 3: Im t at p
//
// so y = E c lists, place by place, the candidate's four real components
// there. A point passes when they form a unit vector at every place, within
// tol, and the identity place's lies in the cap y . v >= s_min.
struct LeafTest {
    int lattice_dimension;
    const double* embedding;
    const double* v;
    double s_min;
    double tol;

    // The conjugate places first: the walk's forms already hold its leaves
    // near the cap, so the identity place passes almost all of them and the
    // rejections come from the conjugates. The test is a conjunction of
    // independent checks, so the order changes only the time.
    bool passes(const double* coefficients) const {
        for (int place = 1; place < lattice_dimension / 4; ++place)
            if (!unit_vector_at(place, coefficients)) return false;
        const std::optional<PlaceVector> y = unit_vector_at(0, coefficients);
        return y && in_cap(*y);
    }

    // Rows 4 place, ..., 4 place + 3 of y, if they form a unit vector within
    // tol. The squared norm only rises as the components come in, so the test
    // stops at the first one that takes it past 1 + tol.
    std::optional<PlaceVector> unit_vector_at(int place, const double* coefficients) const {
        PlaceVector y;
        double norm_sq = 0.0;
        for (int component = 0; component < 4; ++component) {
            y[component] = embedded(4 * place + component, coefficients);
            norm_sq += y[component] * y[component];
            if (norm_sq > 1.0 + tol) return std::nullopt;
        }
        if (std::fabs(norm_sq - 1.0) >= tol) return std::nullopt;
        return y;
    }

    // Row `row` of y = embedding c.
    double embedded(int row, const double* coefficients) const {
        const double* entries = &embedding[static_cast<size_t>(row) * lattice_dimension];
        double value = 0.0;
        for (int column = 0; column < lattice_dimension; ++column)
            value += entries[column] * coefficients[column];
        return value;
    }

    // The cap at the identity place: y . v >= s_min.
    bool in_cap(const PlaceVector& y) const {
        double dot = 0.0;
        for (int component = 0; component < 4; ++component) dot += y[component] * v[component];
        return dot >= s_min;
    }
};

// The walk of the forms below coordinate `top`, with x[top..n) fixed to a
// prefix (prefix[k] is x[top + k]; `top == n` fixes nothing). A coordinate's
// interval is the intersection of the
// intervals the forms allow.
//
// At coordinate `lvl` the walk works on row `lvl` of every form: its entries
// at the columns j > lvl, the coordinates already fixed, make the partial sum,
// and its diagonal turns what is left of the radius into an interval for
// x[lvl].
template <int NL>
class Walker {
  public:
    Walker(int lattice_dimension, const Forms<NL>& forms, const int64_t* transform_columns,
           const LeafTest& test, const int64_t* base)
        : lattice_dimension_(lattice_dimension),
          form_count_(forms.count()),
          forms_(forms),
          transform_columns_(transform_columns),
          test_(test),
          x_(lattice_dimension, 0),
          upper_bound_(lattice_dimension, 0),
          remaining_(static_cast<size_t>(lattice_dimension + 1) * NL, 0.0),
          displacement_(static_cast<size_t>(lattice_dimension) * NL, 0.0),
          cached_sums_stride_(static_cast<size_t>(lattice_dimension + 1) * NL),
          cached_sums_(static_cast<size_t>(lattice_dimension) * cached_sums_stride_, 0.0),
          stale_from_(lattice_dimension + 1, lattice_dimension - 1),
          image_(static_cast<size_t>(lattice_dimension + 1) * lattice_dimension, 0),
          coefficients_as_double_(lattice_dimension, 0.0) {
        // S(n) = g in every row of the cache: the sum before any column.
        for (int row = 0; row < lattice_dimension; ++row) {
            double* sum = cached_sum_from(cached_sums_of_row(row), lattice_dimension);
            for (int f = 0; f < NL; ++f) sum[f] = forms.offsets(row)[f];
        }
        if (base != nullptr)
            for (int coefficient = 0; coefficient < lattice_dimension; ++coefficient)
                image_[static_cast<size_t>(lattice_dimension) * lattice_dimension + coefficient] =
                    base[coefficient];
    }

    // A depth-first walk of the tree below `top`: each level lvl fixes one variable x[lvl],
    // and level 0 completes a point. Every step does one of four things:
    //
    //   x[lvl] past the upper end of its interval      -> up a level, next value there
    //   x[lvl] outside a form (inside_every_form() fails) -> next value at this level
    //   lvl == 0, a complete point                     -> count it, test it, next value
    //   otherwise                                      -> down a level, enter_level()
    //
    // A step costs little because the walk only ever changes the lowest
    // coordinate it is working on, x[lvl], and keeps two running sums per
    // level, so that what was computed for the coordinates above stays valid:
    //
    //   image[lvl] = T[:, lvl:] x[lvl:] + base, with T the transformation
    //   matrix from the walk's coordinates to the power basis: the exact
    //   coefficients of the coordinates fixed so far, and image[0] those of the
    //   complete point that the leaf test needs. Going down a level sets
    //   image[lvl] from image[lvl + 1] (set_coordinate()), and the next value
    //   of x[lvl] adds one column of T to it (next_value()): n additions instead
    //   of the n x n product T x + base per leaf, and in integers, so exactly
    //   the same numbers.
    //
    //   The partial sums g_l + sum_{j > l} R[l][j] (x_j - c_j) of every form's
    //   row l (cached_sums_). A change of x[lvl] leaves the rows below stale only
    //   from column lvl down, and update_partial_sums() adds just those terms -- a
    //   single one for row lvl - 1 after the next value of x[lvl] -- in the
    //   order a full sum from the top would, so the result is the same double.
    //
    // `budget` is a shared count of the leaves the `walkers` walkers have
    // claimed; with max_points > 0 every walker stops once max_points are
    // claimed (take_leaf()).
    WalkResult run(int top, const int64_t* prefix, int64_t max_points, int walkers,
                   std::atomic<int64_t>& budget) {
        WalkResult result;
        if (!build_prefix(top, prefix)) return result;
        if (top == 0) {  // a complete point
            ++result.visited;
            if (is_valid_leaf()) keep(result);
            return result;
        }
        int level = top - 1;
        enter_level(level);
        // A leaf budget cannot bound the time on its own: a tree can have any
        // number of interior nodes that reach no leaf at all, as a thin rotated
        // region's can. The node cap, 64 nodes a leaf of the budget, is far above
        // what a tree with leaves needs, so it fires only on such trees. Without
        // a budget it is out of reach, so that a step costs one comparison.
        const int64_t node_cap =
            max_points > 0 ? max_points * 64 : std::numeric_limits<int64_t>::max();
        int64_t nodes = 0;
        while (true) {
            if (++nodes > node_cap) break;
            if (x_[level] > upper_bound_[level]) {
                // level used up. Go one level back (i.e. higher) if possible
                if (level == top - 1) break;  // we cannot go higher in level (top -1 was the start)
                ++level;
                next_value(level);
                continue;
            }
            if (!inside_every_form(level)) {
                next_value(level);
                continue;
            }
            if (level == 0) {
                ++result.visited;
                if (max_points > 0 && !take_leaf(max_points, walkers, budget)) break;
                if (is_valid_leaf()) keep(result);
                next_value(level);
                continue;
            }
            --level;
            enter_level(level);
        }
        // leaves claimed and not used go back to the other walkers
        if (claimed_leaves_ > 0) budget.fetch_sub(claimed_leaves_, std::memory_order_relaxed);
        return result;
    }

  private:
    // What each form has left of its radius when the walk comes to coordinate
    // `level`, in NL lanes. Level -1 is a spare row that level 0 writes into, so
    // that inside_every_form() needs no branch for it.
    double* remaining_at(int level) {
        return remaining_.data() + static_cast<size_t>(level + 1) * NL;
    }
    // x[coordinate] - c_f[coordinate] for every form f, in NL lanes.
    double* displacements_from_centers(int coordinate) {
        return displacement_.data() + static_cast<size_t>(coordinate) * NL;
    }
    // Row `row` of the cache: S(j) for j = row + 1, ..., n, each in NL lanes.
    double* cached_sums_of_row(int row) {
        return cached_sums_.data() + static_cast<size_t>(row) * cached_sums_stride_;
    }
    // The partial sum of row `row` for every form, in NL lanes: S(row + 1) in the
    // row's cache, once update_partial_sums() has brought it up to date. A form
    // that does not constrain the row has zero entries and offset there, so its
    // lane is exactly 0.
    const double* partial_sums_of_row(int row) {
        return cached_sum_from(cached_sums_of_row(row), row + 1);
    }
    // S(column) in a row of the cache, `row_sums`: the row's sum from `column`
    // up to the top coordinate, in NL lanes.
    static double* cached_sum_from(double* row_sums, int column) {
        return row_sums + static_cast<size_t>(column) * NL;
    }

    // x[top..n) = prefix, consumed from the top down by the same steps as the
    // walk's way down; false when it lies outside a form. What each form has
    // left then belongs to coordinate top - 1. Every form starts with its whole
    // radius (an unused lane with 0).
    bool build_prefix(int top, const int64_t* prefix) {
        double* budget = remaining_at(lattice_dimension_ - 1);
        for (int lane = 0; lane < NL; ++lane) budget[lane] = forms_.radius_sq(lane);
        for (int level = lattice_dimension_ - 1; level >= top; --level) {
            x_[level] = prefix[level - top];
            set_coordinate(level);
            update_partial_sums(level);
            if (!inside_every_form(level)) return false;
        }
        return true;
    }

    // x[level] - c_f[level] for every form.
    void set_displacement(int level) {
        const double value = static_cast<double>(x_[level]);
        double* displacement = displacements_from_centers(level);
        const double* center = forms_.centers(level);
        for (int f = 0; f < NL; ++f) displacement[f] = value - center[f];
    }

    // x[level] newly set: image[level] = image[level + 1] + T[:, level] x[level].
    void set_coordinate(int level) {
        const int64_t* above = image_.data() + static_cast<size_t>(level + 1) * lattice_dimension_;
        int64_t* here = image_.data() + static_cast<size_t>(level) * lattice_dimension_;
        const int64_t* column =
            transform_columns_ + static_cast<size_t>(level) * lattice_dimension_;
        const int64_t value = x_[level];
        for (int coefficient = 0; coefficient < lattice_dimension_; ++coefficient)
            here[coefficient] = wrapping_mul_add(above[coefficient], column[coefficient], value);
        set_displacement(level);
    }

    // ++x[level]: image[level] += T[:, level]. A step past the upper end of the
    // interval (once per interval, the most common step on Clifford+sqrt(T))
    // stops there: run() goes back up at once, and nothing reads this level's
    // image or displacements before the next enter_level(level).
    void next_value(int level) {
        if (++x_[level] > upper_bound_[level]) return;
        int64_t* here = image_.data() + static_cast<size_t>(level) * lattice_dimension_;
        const int64_t* column =
            transform_columns_ + static_cast<size_t>(level) * lattice_dimension_;
        for (int coefficient = 0; coefficient < lattice_dimension_; ++coefficient)
            here[coefficient] += column[coefficient];
        set_displacement(level);
    }

    // Brings the partial sum of row `level` up to date for every form:
    // g[level] + sum_{j > level} R[level][j] (x_j - c_j), the part of the row
    // already fixed, read through partial_sums_of_row(). Called on entering
    // coordinate level, with x[level + 1..n) fixed. The top row has no columns
    // above it: its partial sum is g.
    //
    // Row `level` of the cache, `sums`, holds S(j) for j = level + 1, ..., n,
    // each in NL lanes (cached_sum_from()).
    void update_partial_sums(int level) {
        double* sums = cached_sums_of_row(level);
        if (level + 1 < lattice_dimension_) {
            // 1. Bring row `level` of the cache up to date: S(j) = S(j + 1) +
            //    R[level][j] (x_j - c_j) from its first stale column down to
            //    level + 1. The sums above that column are still valid.
            for (int j = stale_from_[level + 1]; j > level; --j) {
                const double* upper_triangular_entry = forms_.upper_triangular(level, j);
                const double* displacement = displacements_from_centers(j);
                const double* sum_from_next = cached_sum_from(sums, j + 1);  // S(j + 1)
                double* sum_from_j = cached_sum_from(sums, j);               // S(j)
                for (int f = 0; f < NL; ++f)
                    sum_from_j[f] = sum_from_next[f] + upper_triangular_entry[f] * displacement[f];
            }
            // 2. Pass the staleness on (stale_from_[l + 1] is row l's first
            //    stale column). Whatever changed above row `level` has also
            //    changed for row level - 1, which the walk sums next; and the
            //    walk comes back to row `level` only after x[level + 1] has
            //    moved, so from column level + 1 at least.
            if (stale_from_[level + 1] > stale_from_[level])
                stale_from_[level] = stale_from_[level + 1];
            stale_from_[level + 1] = level + 1;
        }
    }

    // Going down to coordinate `level`, with x[level + 1..n) fixed: finds the
    // integers x[level] can take in every form and starts at the smallest.
    //
    // Row `level` of a form reads (R_ll (x - c) + p)^2 <= left: p is the row's
    // partial sum over the fixed coordinates (update_partial_sums()), `left` what is
    // left of the form's radius after the rows above. Solved for x, that is the
    // interval
    //
    //     middle -+ half_width,  middle = c - p / R_ll,  half_width = sqrt(left) / R_ll,
    //
    // and x[level] must lie in the intersection over the forms. A form whose
    // diagonal is not positive gives no interval (inverse 0): one that does not
    // constrain this row, whose row is zero, and one with a zero diagonal under
    // nonzero entries, whose budget inside_every_form() still checks.
    //
    // Each interval is widened by kRound, so that rounding cannot drop an
    // integer at its edge; inside_every_form() then checks each value exactly. The
    // intersection is rounded once -- ceil of the largest lower bound is the
    // largest ceil, so the integers are the same as rounding every form -- and
    // converted to int64 with a clamp at kIndexMax, past which a double has no
    // int64 value. A NaN bound compares false and drops out, the permissive
    // side.
    //
    // Ends with x[level] at the lower end and upper_bound_[level] at the upper
    // end. An empty interval (on Clifford+sqrt(T) over a third of the calls)
    // leaves x[level] > upper_bound_[level] and stops there: run() goes back up
    // at once. Otherwise the image and the displacements are brought up to date
    // for the new x[level] (set_coordinate()).
    void enter_level(int level) {
        constexpr double kEdge = static_cast<double>(kIndexMax);
        double lower = -kEdge;
        double upper = kEdge;
        update_partial_sums(level);
        const double* center = forms_.centers(level);
        const double* inverse_diagonal = forms_.inverse_diagonals(level);
        const double* partial = partial_sums_of_row(level);
        const double* remaining = remaining_at(level);
        for (int form = 0; form < form_count_; ++form) {
            const double inverse = inverse_diagonal[form];
            if (inverse <= 0.0) continue;
            const double left = remaining[form];
            const double half_width = std::sqrt(left < 0.0 ? 0.0 : left) * inverse;
            const double middle = center[form] - partial[form] * inverse;
            const double form_lower = middle - half_width - kRound;
            const double form_upper = middle + half_width + kRound;
            if (form_lower > lower) lower = form_lower;
            if (form_upper < upper) upper = form_upper;
        }
        // lower >= -kEdge and upper <= kEdge already; clamp the other side
        lower = std::ceil(lower);
        upper = std::floor(upper);
        const int64_t first = lower < kEdge ? static_cast<int64_t>(lower) : kIndexMax;
        const int64_t last = upper > -kEdge ? static_cast<int64_t>(upper) : -kIndexMax;
        upper_bound_[level] = last;
        x_[level] = first;
        // An empty interval: run() goes back up at once, and nothing reads this
        // level's image or displacements before the next enter_level(level).
        if (first > last) return;
        set_coordinate(level);
    }

    // Whether x[level..n) is inside every form; what each has left then belongs to
    // coordinate level - 1.
    //
    // One pass over all NL lanes, without a branch per form, so that it
    // vectorises. A form that does not constrain this row, and an unused lane,
    // has a diagonal and a partial sum of exactly 0 there: its term is 0 and its
    // budget passes on unchanged. The pass writes every lane even when a form
    // fails; enter_level(level - 1) reads that row only after a pass that
    // succeeds.
    bool inside_every_form(int level) {
        const double* diagonal = forms_.diagonals(level);
        const double* displacement = displacements_from_centers(level);
        const double* partial = partial_sums_of_row(level);
        const double* remaining = remaining_at(level);
        double* remaining_below = remaining_at(level - 1);
        bool outside = false;
        for (int lane = 0; lane < NL; ++lane) {
            const double term = diagonal[lane] * displacement[lane] + partial[lane];
            const double left = remaining[lane] - term * term;
            remaining_below[lane] = left;
            outside |= left < -kSlack;
        }
        return !outside;
    }

    // One leaf of the shared budget, false once max_points are used up. A walker
    // claims leaves from the shared count in blocks: an add on it moves its
    // cache line between the cores, and with one add a leaf a capped walk on
    // Clifford+T ran 3 times slower on 4 threads and 6 times on 8. A block is at
    // most kLeafBlock and at most 1 / (2 walkers) of what is left, so near the
    // end no walker holds many unused leaves while another runs out, and the
    // walkers together test close to max_points leaves. A claim is cut to what
    // is left, so they never test more, and a walk on one thread stops at the
    // same leaf as with one add a leaf.
    bool take_leaf(int64_t max_points, int walkers, std::atomic<int64_t>& budget) {
        if (claimed_leaves_ == 0) {
            const int64_t left = max_points - budget.load(std::memory_order_relaxed);
            const int64_t block = std::clamp<int64_t>(left / (2 * walkers), 1, kLeafBlock);
            const int64_t claimed_before = budget.fetch_add(block, std::memory_order_relaxed);
            claimed_leaves_ = std::clamp<int64_t>(max_points - claimed_before, 0, block);
            if (claimed_leaves_ == 0) return false;
        }
        --claimed_leaves_;
        return true;
    }

    // The leaf test on the exact coefficients image[0].
    bool is_valid_leaf() {
        for (int coefficient = 0; coefficient < lattice_dimension_; ++coefficient)
            coefficients_as_double_[coefficient] = static_cast<double>(image_[coefficient]);
        return test_.passes(coefficients_as_double_.data());
    }

    // Survivors live outside the walker; a member would make the compiler assume aliasing
    // and stop vectorising the lane loops.
    void keep(WalkResult& result) const {
        result.coefficients.insert(result.coefficients.end(), image_.data(),
                                   image_.data() + lattice_dimension_);
    }

    const int lattice_dimension_;
    const int form_count_;
    const Forms<NL>& forms_;
    const int64_t* transform_columns_;  // T by column, see transpose()
    const LeafTest test_;
    std::vector<int64_t> x_;
    std::vector<int64_t> upper_bound_;  // the upper end of x[level]'s interval, inclusive
    std::vector<double> remaining_;     // what is left of each radius, see remaining_at()
    // By coordinate, forms in lanes: x[j] - c_f[j].
    std::vector<double> displacement_;
    // Partial sums kept between visits: S(l, j) = g_l + sum_{k >= j} R[l][k]
    // (x_k - c_k) for l < j <= n, summed from the top coordinate down, row l
    // at cached_sums_of_row(l). Row l is then brought up to date from its
    // first stale column, the highest coordinate that changed since it was last summed, down to l +
    // 1 -- one term when only x[l + 1] moved, as when the walk steps to a sibling -- and S(l, l +
    // 1) is the same number a full sum in that order gives. stale_from_[l + 1] holds row l's first
    // stale column: the entry for row `level` is stale_from_[level + 1].
    const size_t cached_sums_stride_;
    std::vector<double> cached_sums_;
    std::vector<int> stale_from_;
    // image[lvl] = T[:, lvl:] x[lvl:] + base, the power-basis coefficients of the
    // coordinates fixed so far; image[n] = base. Updated by one column per step
    // instead of a full n x n product per leaf. Exact in int64. Starting at
    // base = T base_x, the integer origin the walk is centred on, keeps every
    // coordinate of order one where the centre would otherwise run to 1e14,
    // which is what the tight forms need to stay meaningful.
    std::vector<int64_t> image_;
    // The coefficients as doubles, converted once per leaf: int64 -> double has
    // no vector instruction before AVX-512, so converting inside the four dot
    // products would serialise them.
    std::vector<double> coefficients_as_double_;
    int64_t claimed_leaves_ = 0;  // taken from the shared budget, not used yet (take_leaf())
};

// The work items of a split walk: prefixes that fix x[top..n), each a block of
// `width` = n - top values (x[top], ..., x[n-1]), in the order the walk would
// visit them, and what the walk's ellipsoid has left after each.
struct Prefixes {
    int top = 0;
    int width = 0;
    std::vector<int64_t> values;
    std::vector<double> remaining;

    size_t count() const { return width > 0 ? values.size() / static_cast<size_t>(width) : 0; }
    // The k-th prefix: x[top], ..., x[n-1].
    const int64_t* prefix(size_t k) const { return values.data() + k * static_cast<size_t>(width); }

    // The items in the order to hand them out: most left of the ellipsoid
    // first. The subtree below a prefix grows with what is left, so this starts
    // the largest items first and keeps one of them from finishing the walk
    // alone while the other threads wait (longest processing time first).
    std::vector<size_t> largest_first() const {
        std::vector<size_t> order(count());
        std::iota(order.begin(), order.end(), size_t{0});
        std::stable_sort(order.begin(), order.end(),
                         [&](size_t a, size_t b) { return remaining[a] > remaining[b]; });
        return order;
    }
};

// Every prefix x[top..n) inside `ellipsoid`, the walk's own. This is the walk
// of Walker::run() on that one form, stopped at level `top`: the same four
// cases, without the lanes and the cached sums, since it runs only a few times
// per walk. Only the ellipsoid is used: a prefix the tighter forms reject is
// dropped when the worker consumes it (build_prefix()), so the split stays a
// split and the result does not depend on it.
Prefixes prefixes_inside(int lattice_dimension, const Ellipsoid& ellipsoid, int top) {
    const double* center = ellipsoid.center;
    auto upper_triangular = [&](int row, int column) {
        return ellipsoid.upper_triangular[static_cast<size_t>(row) * lattice_dimension + column];
    };
    std::vector<int64_t> x(lattice_dimension, 0), upper_bound(lattice_dimension, 0);
    std::vector<double> remaining(lattice_dimension, 0.0), partial(lattice_dimension, 0.0);
    // As Walker::enter_level(): the interval of x[level], x[level] at its lower end.
    auto enter_level = [&](int level) {
        double sum = 0.0;
        for (int j = level + 1; j < lattice_dimension; ++j)
            sum += upper_triangular(level, j) * (static_cast<double>(x[j]) - center[j]);
        partial[level] = sum;
        const double left = remaining[level] < 0 ? 0.0 : remaining[level];
        const double diagonal = upper_triangular(level, level);
        const double half_width = std::sqrt(left) / diagonal;
        const double middle = center[level] - sum / diagonal;
        x[level] = lower_index(std::ceil(middle - half_width - kRound));
        upper_bound[level] = upper_index(std::floor(middle + half_width + kRound));
    };
    Prefixes out{top, lattice_dimension - top, {}, {}};
    int level = lattice_dimension - 1;
    remaining[level] = ellipsoid.radius_sq;
    enter_level(level);
    while (true) {
        if (x[level] > upper_bound[level]) {  // level used up: up a level, next value there
            if (level == lattice_dimension - 1) break;
            ++level;
            ++x[level];
            continue;
        }
        const double term =
            upper_triangular(level, level) * (static_cast<double>(x[level]) - center[level]) +
            partial[level];
        const double left = remaining[level] - term * term;
        if (left < -kSlack) {  // outside the ellipsoid: next value
            ++x[level];
            continue;
        }
        if (level == top) {  // a whole prefix: keep it, next value
            out.values.insert(out.values.end(), x.begin() + top, x.end());
            out.remaining.push_back(left);
            ++x[level];
            continue;
        }
        remaining[level - 1] = left;  // otherwise: down a level
        --level;
        enter_level(level);
    }
    return out;
}

// The work items for a walk on `threads` threads: the prefixes of the fewest top
// coordinates (at most six) that give at least eight items per thread. Many
// more items than threads let the shared queue of parallel_for() even out their
// very different sizes. Each depth is enumerated from scratch, which costs
// microseconds next to the walk.
Prefixes split_for_threads(int lattice_dimension, const Ellipsoid& ellipsoid, int threads) {
    constexpr int kItemsPerThread = 8;
    constexpr int kMaxDepth = 6;
    Prefixes items;
    for (int depth = 1; depth <= std::min(kMaxDepth, lattice_dimension - 1); ++depth) {
        items = prefixes_inside(lattice_dimension, ellipsoid, lattice_dimension - depth);
        if (items.count() >= static_cast<size_t>(kItemsPerThread * threads)) break;
    }
    return items;
}

// One walk on NL lanes: whole on one thread, or split over `threads` by the
// prefixes of its top coordinates. The forms are laid out once and shared.
template <int NL>
WalkResult walk_on_lanes(const WalkLattice& lattice, const WalkInputs& inputs, int64_t max_points,
                         int threads) {
    const int lattice_dimension = lattice.lattice_dimension;
    const std::vector<int64_t> transform_columns = transpose(lattice_dimension, inputs.transform);
    const Forms<NL> forms(lattice_dimension, inputs.forms);
    const LeafTest test{lattice_dimension, lattice.embedding, inputs.v, lattice.s_min, lattice.tol};
    if (threads <= 1) {
        std::atomic<int64_t> budget{0};
        return Walker<NL>(lattice_dimension, forms, transform_columns.data(), test, inputs.base)
            .run(lattice_dimension, nullptr, max_points, 1, budget);
    }
    const Prefixes items = split_for_threads(lattice_dimension, inputs.forms[0], threads);
    // Handed out largest first, but each result goes to its item's place, so
    // that the concatenation below is in the walk's order.
    const std::vector<size_t> order = items.largest_first();
    std::vector<WalkResult> parts(items.count());
    std::atomic<int64_t> budget{0};
    parallel_for(parts.size(), threads, [&](size_t next) {
        if (max_points > 0 && budget.load(std::memory_order_relaxed) > max_points) return;
        const size_t k = order[next];
        parts[k] = Walker<NL>(lattice_dimension, forms, transform_columns.data(), test, inputs.base)
                       .run(items.top, items.prefix(k), max_points, threads, budget);
    });
    WalkResult result;  // concatenated in item order: deterministic
    for (auto& p : parts) {
        result.visited += p.visited;
        result.coefficients.insert(result.coefficients.end(), p.coefficients.begin(),
                                   p.coefficients.end());
    }
    return result;
}

}  // namespace

// On as few lanes as hold the forms.
WalkResult walk(const WalkLattice& lattice, const WalkInputs& inputs, int64_t max_points,
                int threads) {
    if (inputs.forms.size() <= 4) return walk_on_lanes<4>(lattice, inputs, max_points, threads);
    if (inputs.forms.size() <= 8) return walk_on_lanes<8>(lattice, inputs, max_points, threads);
    throw std::invalid_argument("walk: at most eight forms");
}

std::vector<WalkResult> walk_batch(const WalkLattice& lattice, const std::vector<WalkInputs>& walks,
                                   int64_t max_points, int threads) {
    std::vector<WalkResult> out(walks.size());
    parallel_for(out.size(), threads,
                 [&](size_t k) { out[k] = walk(lattice, walks[k], max_points, 1); });
    return out;
}

}  // namespace quicksynth
