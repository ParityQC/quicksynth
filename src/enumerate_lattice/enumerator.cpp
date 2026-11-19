// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/enumerate_lattice/enumerator.hpp"

#include "quicksynth/enumerate_lattice/lll.hpp"
#include "quicksynth/enumerate_lattice/setup.hpp"
#include "quicksynth/enumerate_lattice/walk.hpp"
#include "quicksynth/support/clock.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace quicksynth {

namespace {

constexpr double kFormSlack = 1e-9;  // rounding room of the cap and conjugate forms
constexpr int64_t kCoefficientLimit = int64_t{1} << 60;  // 128-bit norm check

// The forms of a walk, numbered as the walk's lanes (walk.cpp): its ellipsoid,
// then the forms that prune it.
constexpr int kEllipsoidForm = 0;
constexpr int kCapForm = 1;
constexpr int kConjugateForm = 2;   // the conjugate block
constexpr int kFirstPlaceForm = 3;  // + a: conjugate place a, where there are several

// Inputs of one walk at a scale, in the kernel's layout. Each form is
// |R (x - center) + offset|^2 <= radius_sq in the walk's coordinates x: the
// reduced ones less the integer origin the walk is centred on. Form 0, the
// walk's ellipsoid, has fields of its own; the pruning forms 1, 2, ... are
// stored back to back in the pruning_ arrays, form i in block i - 1.
struct Walk {
    int lattice_dimension = 0;  // n
    // form 0, the walk's ellipsoid (no offset)
    std::vector<double> ellipsoid_upper;   // n x n, its R at the walk's scale
    std::vector<double> ellipsoid_center;  // n, the reduced centre less the origin
    double ellipsoid_radius_sq = 0.0;      // P, the number of places
    // forms 1, 2, ...
    std::vector<double> pruning_upper;      // n x n each, R at the walk's scale, zero past its rank
    std::vector<double> pruning_center;     // n each: the cap's is ellipsoid_center, the others 0
    std::vector<double> pruning_offset;     // n each: Q^t E base for a conjugate form, else 0
    std::vector<double> pruning_radius_sq;  // each, its squared radius plus kFormSlack
    // T times the origin: the walk adds it back to every coefficient vector it
    // reports, so those are the power-basis coefficients
    std::vector<int64_t> base;
    const int64_t* transform = nullptr;  // n x n, the set-up's LLL transform T
    Vec4 direction{};                    // the target direction, for the leaves' cap test

    // n coordinates and `form_count` forms, every entry zero; the ellipsoid's
    // centre is left empty, for walk_of() to move in.
    Walk(int lattice_dimension, int form_count)
        : lattice_dimension(lattice_dimension),
          ellipsoid_upper(static_cast<size_t>(lattice_dimension) * lattice_dimension),
          pruning_upper(matrix_block(form_count)),
          pruning_center(vector_block(form_count)),
          pruning_offset(vector_block(form_count)),
          pruning_radius_sq(static_cast<size_t>(form_count - 1)) {}

    int form_count() const { return 1 + static_cast<int>(pruning_radius_sq.size()); }

    // Form i as the walk takes it, a view into this walk.
    Ellipsoid form(int i) const {
        if (i == kEllipsoidForm)
            return {ellipsoid_upper.data(), ellipsoid_center.data(), nullptr, ellipsoid_radius_sq};
        return {pruning_upper.data() + matrix_block(i), pruning_center.data() + vector_block(i),
                pruning_offset.data() + vector_block(i), pruning_radius_sq[i - 1]};
    }

    // The storage of form i, for walk_of() to fill; the ellipsoid has no offset.
    double* upper_of(int i) {
        return i == kEllipsoidForm ? ellipsoid_upper.data()
                                   : pruning_upper.data() + matrix_block(i);
    }
    double* center_of(int i) {
        return i == kEllipsoidForm ? ellipsoid_center.data()
                                   : pruning_center.data() + vector_block(i);
    }
    double* offset_of(int i) { return pruning_offset.data() + vector_block(i); }
    double& radius_sq_of(int i) {
        return i == kEllipsoidForm ? ellipsoid_radius_sq : pruning_radius_sq[i - 1];
    }

    // Form i's R: an R of the set-up, at its reference scale, times `power`.
    void set_upper(int i, const std::vector<double>& R, double power) {
        double* upper = upper_of(i);
        for (size_t k = 0; k < R.size(); ++k) upper[k] = R[k] * power;
    }

    // Pruning form i's offset Q^t y, which moves the form, centred at 0, to the
    // walk's origin: R x' + Q^t y = R x for x = x' + origin. `y` are the reals of
    // the origin's lattice point at the form's `rows` rows of the basis (see
    // embedded_rows()); computed this way every operand stays of order one,
    // where R origin would cancel huge products.
    void set_offset(int i, const QrFactors& form, int rows, const double* y) {
        double* offset = offset_of(i);
        for (int k = 0; k < rows; ++k) {
            double s = 0.0;
            for (int r = 0; r < rows; ++r) s += form.Q[static_cast<size_t>(r) * rows + k] * y[r];
            offset[k] = s;
        }
    }

    // Every form, the transform, the target direction and the base.
    WalkInputs inputs() const {
        WalkInputs out;
        out.forms.reserve(form_count());
        for (int i = 0; i < form_count(); ++i) out.forms.push_back(form(i));
        out.transform = transform;
        out.v = direction.data();
        out.base = base.data();
        return out;
    }

  private:
    // Where pruning form i >= 1 starts in the pruning_ arrays; for i = form_count
    // their size.
    size_t matrix_block(int i) const {
        return static_cast<size_t>(i - 1) * lattice_dimension * lattice_dimension;
    }
    size_t vector_block(int i) const { return static_cast<size_t>(i - 1) * lattice_dimension; }
};

// The factor 2^((reference_scale - scale) / 2) taking an R of a set-up to
// `scale`; exact, the two scales having the same parity.
double rescale_factor(int reference_scale, int scale) {
    return std::ldexp(1.0, (reference_scale - scale) / 2);
}

// The set-up's centre at the walk's scale, split into the integer point nearest
// it and the rest. The reduced centre is huge (about 1e14 to 1e15 at eps =
// 1e-8, where doubles are 0.01 to 0.1 apart), so the walk runs about that
// integer point and computes with the rest only: its operands stay of order
// one. Lattice points stay lattice points, the point being an integer one, and
// the subtraction is exact.
struct SplitCentre {
    std::vector<int64_t> origin;  // the nearest integer point
    std::vector<double> rest;     // the centre less the origin, each in [-1/2, 1/2]
};

SplitCentre split_centre(const std::vector<double>& reduced_center, double power) {
    const size_t n = reduced_center.size();
    SplitCentre out{std::vector<int64_t>(n), std::vector<double>(n)};
    for (size_t i = 0; i < n; ++i) {
        const double centre = reduced_center[i] / power;  // at the walk's scale
        out.origin[i] = static_cast<int64_t>(std::nearbyint(centre));
        out.rest[i] = centre - static_cast<double>(out.origin[i]);
    }
    return out;
}

// T x for the n x n transform T and an integer n-vector x, in wrapping
// arithmetic: the products may pass 64 bits, but the sum, a lattice point near
// the centre, does not, and wrapping keeps it exact.
std::vector<int64_t> apply_transform(const std::vector<int64_t>& T, const std::vector<int64_t>& x) {
    const size_t n = x.size();
    std::vector<int64_t> out(n);
    for (size_t r = 0; r < n; ++r) {
        uint64_t s = 0;
        for (size_t j = 0; j < n; ++j)
            s += static_cast<uint64_t>(T[r * n + j]) * static_cast<uint64_t>(x[j]);
        out[r] = static_cast<int64_t>(s);
    }
    return out;
}

// Rows `first_row` to n - 1 of the embedding times the power-basis coefficients
// x: the reals of that lattice point at the places those rows belong to.
std::vector<double> embedded_rows(const std::vector<double>& embedding,
                                  const std::vector<int64_t>& x, int first_row) {
    const int n = static_cast<int>(x.size());
    std::vector<double> out(static_cast<size_t>(n - first_row));
    for (int r = first_row; r < n; ++r) {
        double s = 0.0;
        for (int j = 0; j < n; ++j)
            s += embedding[static_cast<size_t>(r) * n + j] * static_cast<double>(x[j]);
        out[static_cast<size_t>(r - first_row)] = s;
    }
    return out;
}

// The walk of a set-up at `scale`, whose lattice is `embedding`: the set-up's
// forms rescaled from its reference scale, about the integer point nearest
// their centre (split_centre()).
Walk walk_of(const Setup& setup, int scale, const std::vector<double>& embedding) {
    const int n = setup.lattice_dimension;
    const int place_forms = static_cast<int>(setup.conjugate_places.size());
    const double power = rescale_factor(setup.reference_scale, scale);
    SplitCentre centre = split_centre(setup.reduced_center, power);
    Walk out(n, kFirstPlaceForm + place_forms);
    out.direction = setup.direction;
    out.transform = setup.basis->to_power_basis.data();
    // the reported coefficients are B x' + base for the walk's x' = x - origin
    out.base = apply_transform(setup.basis->to_power_basis, centre.origin);
    // the walk's ellipsoid and the cap, about the same centre
    out.ellipsoid_center = std::move(centre.rest);
    out.set_upper(kEllipsoidForm, setup.ellipsoid_upper, power);
    out.radius_sq_of(kEllipsoidForm) = setup.ellipsoid_radius_sq;
    out.set_upper(kCapForm, setup.cap_upper, power);
    std::copy(out.ellipsoid_center.begin(), out.ellipsoid_center.end(), out.center_of(kCapForm));
    out.radius_sq_of(kCapForm) = 1.0 + kFormSlack;
    // the conjugate forms, centred at 0; their offsets move them to the origin
    const std::vector<double> conjugate_reals = embedded_rows(embedding, out.base, kRowsPerPlace);
    out.set_upper(kConjugateForm, setup.conjugate_block.R, power);
    out.set_offset(kConjugateForm, setup.conjugate_block, n - kRowsPerPlace,
                   conjugate_reals.data());
    out.radius_sq_of(kConjugateForm) = setup.conjugate_radius_sq + kFormSlack;
    for (int a = 0; a < place_forms; ++a) {
        const int form = kFirstPlaceForm + a;
        out.set_upper(form, setup.conjugate_places[a].R, power);
        out.set_offset(form, setup.conjugate_places[a], kRowsPerPlace,
                       conjugate_reals.data() + kRowsPerPlace * a);
        out.radius_sq_of(form) = 1.0 + kFormSlack;
    }
    return out;
}

// The first column (u1, u2) of a walk survivor, its n power-basis coefficients
// (u1's, then u2's) placed at the powers zeta16^(step i) of Z[zeta16]. Throws
// std::overflow_error at 2^60, beyond which unit_norm()'s 128-bit sums overflow.
struct Column {
    Z16 u1, u2;
};

Column unpack_into_Z16(const int64_t* survivor, const GateSet& gateset) {
    const int d = gateset.degree();
    Column out{};
    for (int i = 0; i < d; ++i) {
        if (survivor[i] >= kCoefficientLimit || survivor[i] <= -kCoefficientLimit ||
            survivor[d + i] >= kCoefficientLimit || survivor[d + i] <= -kCoefficientLimit)
            throw std::overflow_error("scan: a survivor's coefficients exceed 60 bits");
        // Note: for Clifford + T step is 2 so every even place gets populated since we work in
        // Z[zeta16].
        out.u1[gateset.step * i] = survivor[i];
        out.u2[gateset.step * i] = survivor[d + i];
    }
    return out;
}

}  // namespace

Enumerator::Enumerator(double epsilon, GateSet gateset, Price price, int threads, Settings settings)
    : region_(epsilon, std::move(gateset)),
      pricer_(Pricer::for_gateset(region_.gateset(), price)),
      threads_(std::max(1, threads)),
      setups_(std::make_unique<SetupCache>(std::max<int64_t>(1, settings.transform_cache))) {}

Enumerator::~Enumerator() = default;
Enumerator::Enumerator(Enumerator&&) noexcept = default;

std::vector<CyclotomicMatrix> Enumerator::matrices(const int64_t* walk_survivors,
                                                   size_t num_survivors, int denominator_exponent,
                                                   int det_class) const {
    const GateSet& gateset = region_.gateset();
    const size_t n = static_cast<size_t>(gateset.dim());
    std::vector<CyclotomicMatrix> out;
    for (size_t s = 0; s < num_survivors; ++s) {
        const auto [u1, u2] = unpack_into_Z16(walk_survivors + s * n, gateset);
        if (!unit_norm(u1, u2, denominator_exponent)) continue;
        out.push_back(times_q_power(su2_element(u1, u2, denominator_exponent), det_class, gateset));
    }
    return out;
}

Enumerator::Result Enumerator::enumerate(const Unitary& V, int K, int det_class, int sub) {
    const double start = steady_seconds();
    const GateSet& gateset = region_.gateset();
    const auto setup = setups_->for_problems(region_, threads_, {{V, det_class}}, K, sub)[0];
    const std::vector<double>& E = setups_->embedding(gateset, K);
    const Walk w = walk_of(*setup, K, E);
    const int n = gateset.dim();
    const WalkResult result = walk(region_.walk_lattice(E), w.inputs(), 0, threads_);
    Result out;
    out.K = K;
    out.sub = sub;
    out.det_class = det_class;
    out.points = result.visited;
    const Vec4 target = su2_vector(V);
    for (const CyclotomicMatrix& e :
         matrices(result.coefficients.data(), result.coefficients.size() / n, K, det_class)) {
        const double distance = chord_distance(matrix_direction(e), target);
        if (distance > epsilon()) continue;  // admitted by the slack of the cap pre-filter only
        out.hits.push_back({Cost{}, e});
    }
    for (Hit& h : out.hits) h.cost = pricer_.price(h.element);
    out.seconds = steady_seconds() - start;
    return out;
}

std::vector<WalkResult> Enumerator::survivors_batch(const std::vector<Problem>& problems, int K,
                                                    int64_t walk_budget, int sub) {
    if (problems.empty()) return {};
    const GateSet& gateset = region_.gateset();
    const auto setups = setups_->for_problems(region_, threads_, problems, K, sub);
    const std::vector<double>& E = setups_->embedding(gateset, K);
    std::vector<Walk> walks;
    walks.reserve(setups.size());
    for (const auto& setup : setups) walks.push_back(walk_of(*setup, K, E));
    std::vector<WalkInputs> batch;
    batch.reserve(walks.size());
    for (const Walk& w : walks) batch.push_back(w.inputs());
    return walk_batch(region_.walk_lattice(E), batch, walk_budget, threads_);
}

}  // namespace quicksynth
