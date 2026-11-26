// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/enumerate_lattice/setup.hpp"

#include "quicksynth/enumerate_lattice/lll.hpp"
#include "quicksynth/data_structures/double_double.hpp"
#include "quicksynth/support/parallel.hpp"

#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

namespace quicksynth {

namespace {

constexpr int kLllBits = 48;                           // fractional bits of the LLL input
constexpr int64_t kTransformLimit = int64_t{1} << 40;  // the walk's int64 accumulator
// below 2^53: B^-1's entries k / 2^sub are then exact doubles
constexpr int64_t kInverseNumeratorLimit = int64_t{1} << 53;

DirectionBits direction_bits(const Vec4& v) {
    DirectionBits out;
    std::memcpy(out.data(), v.data(), sizeof(out));
    return out;
}

SubBits sub_bits(const Vec4& v, int sub) {
    const DirectionBits b = direction_bits(v);
    return {b[0], b[1], b[2], b[3], static_cast<uint64_t>(sub)};
}

// The sublattice delta^sub Z[zeta]^2: its delta matrix M_sub, whose columns delta^sub zeta^j
// span it, and the numerators 2^sub M_sub^-1 of its inverse.
struct Sublattice {
    int sub = 0;
    std::vector<int64_t> delta_matrix;
    std::vector<__int128> inverse_delta_matrix_numerators;
};

Sublattice sublattice(const GateSet& gateset, int sub) {
    Sublattice out{sub, delta_power(gateset, sub), {}};
    for (double x : inverse_delta_power(gateset, sub))
        out.inverse_delta_matrix_numerators.push_back(static_cast<__int128>(std::ldexp(x, sub)));
    return out;
}

std::vector<__int128> widened(const std::vector<int64_t>& v) { return {v.begin(), v.end()}; }

// k / 2^sub as a double, k in int64 range by a hardware conversion
double over_power_of_two(__int128 k, int sub) {
    const bool narrow =
        k >= std::numeric_limits<int64_t>::min() && k <= std::numeric_limits<int64_t>::max();
    const double x = narrow ? static_cast<double>(static_cast<int64_t>(k)) : static_cast<double>(k);
    return sub == 0 ? x : x * std::ldexp(1.0, -sub);
}

// The basis of `level` from a fresh LLL reduction of to_ball M_sub: B = M_sub T and
// B^-1 = T^-1 M_sub^-1.
std::shared_ptr<const SublatticeBasis> reduce(const GateSet& gateset,
                                              const std::vector<double>& to_ball,
                                              const Sublattice& level) {
    const int n = gateset.dim(), sub = level.sub;
    LllTransform transform;
    // to_ball M_sub; at sub 1 each entry is one rounded sum of two terms, as exact_product's
    const std::vector<double> rows = lll_rows(
        n,
        sub == 0
            ? to_ball.data()
            : product(n, to_ball,
                      std::vector<double>(level.delta_matrix.begin(), level.delta_matrix.end()))
                  .data(),
        kLllBits);
    if (!lll_transform(n, rows, kTransformLimit, transform))
        throw std::overflow_error("scan: an LLL transform or its inverse too wide for the walk");
    auto basis = std::make_shared<SublatticeBasis>();
    basis->sub = sub;
    if (sub == 0) {  // M_0 = 1: B = T
        basis->to_power_basis = std::move(transform.T);
        basis->from_power_basis.assign(transform.T_inv.begin(), transform.T_inv.end());
        return basis;
    }
    // below 2^44: M_sub's rows sum to at most 2^sub in absolute value, T's entries are below 2^40
    basis->to_power_basis = product(n, level.delta_matrix, transform.T);
    for (int64_t b : basis->to_power_basis)
        if (!(b < kTransformLimit && b > -kTransformLimit))
            throw std::overflow_error("scan: a lattice basis too wide for the walk");
    const std::vector<__int128> numerators =
        product(n, widened(transform.T_inv), level.inverse_delta_matrix_numerators);
    basis->from_power_basis.resize(numerators.size());
    for (size_t i = 0; i < numerators.size(); ++i)
        basis->from_power_basis[i] = over_power_of_two(numerators[i], sub);
    return basis;
}

// y / delta for y in delta Z[zeta]^2, per block by 2 / delta = 1 + zeta + ... + zeta^(d-1):
// coefficient 0 is (y_0 - y_1 - ... - y_(d-1)) / 2, each next one adds y_i. Empty when y is
// not divisible.
std::optional<std::vector<int64_t>> divided_by_delta(int d, const std::vector<int64_t>& y) {
    std::vector<int64_t> out(y.size());
    for (size_t o = 0; o < y.size(); o += static_cast<size_t>(d)) {
        int64_t twice = y[o];
        for (int t = 1; t < d; ++t) twice -= y[o + t];
        if (twice % 2 != 0) return std::nullopt;
        out[o] = twice / 2;
        for (int i = 1; i < d; ++i) out[o + i] = out[o + i - 1] + y[o + i];
    }
    return out;
}

// The label of each column of `basis`: bit 0 is the parity of the coefficient sum of the
// column's u1 / delta^sub, bit 1 that of its u2 / delta^sub. Zero exactly for the columns in
// delta^(sub+1) Z[zeta]^2, as zeta = 1 mod delta. Empty if a column is not divisible by
// delta^sub.
std::optional<std::vector<int>> column_labels(const GateSet& gateset,
                                              const SublatticeBasis& basis) {
    const int n = gateset.dim(), d = gateset.degree();
    std::vector<int> labels(n, 0);
    for (int c = 0; c < n; ++c) {
        std::vector<int64_t> y(n);
        for (int r = 0; r < n; ++r) y[r] = basis.to_power_basis[static_cast<size_t>(r) * n + c];
        for (int s = 0; s < basis.sub; ++s) {
            std::optional<std::vector<int64_t>> divided = divided_by_delta(d, y);
            if (!divided) return std::nullopt;
            y = std::move(*divided);
        }
        for (int r = 0; r < n; ++r) labels[c] ^= static_cast<int>(y[r] & 1) << (r < d ? 0 : 1);
    }
    return labels;
}

// Kernel, the n x n matrix of the combinations of columns whose labels cancel. Its pivots p
// and q are the first two columns with different nonzero labels; their labels and their sum
// are the three nonzero labels, so adding column p, column q or both cancels any label.
// Column j of Kernel is e_j plus e_p and e_q as column j's label needs, and columns p and q
// are 2 e_p and 2 e_q. Empty if no two columns have different nonzero labels.
std::optional<std::vector<int64_t>> parity_kernel(const std::vector<int>& labels) {
    const int n = static_cast<int>(labels.size());
    const auto at = [n](int r, int c) { return static_cast<size_t>(r) * n + c; };
    for (int p = 0; p < n; ++p)
        for (int q = p + 1; q < n; ++q) {
            if (labels[p] == 0 || labels[q] == 0 || labels[p] == labels[q]) continue;
            const int both = labels[p] ^ labels[q];
            std::vector<int64_t> kernel(static_cast<size_t>(n) * n, 0);
            for (int j = 0; j < n; ++j) {
                kernel[at(j, j)] = j == p || j == q ? 2 : 1;
                if (j == p || j == q) continue;
                kernel[at(p, j)] = labels[j] == labels[p] || labels[j] == both;
                kernel[at(q, j)] = labels[j] == labels[q] || labels[j] == both;
            }
            return kernel;
        }
    return std::nullopt;
}

// 2 Kernel^-1 = 3 - Kernel: Kernel = 1 + N, where N's only nonzero rows, p and q, have 1 on
// the diagonal and 0 in each other's column, so N^2 = N and
// Kernel (3 - Kernel) = 3 + 3N - (1 + 2N + N^2) = 2.
std::vector<int64_t> twice_kernel_inverse(int n, const std::vector<int64_t>& kernel) {
    std::vector<int64_t> out(kernel.size());
    for (size_t i = 0; i < out.size(); ++i)
        out[i] = (i % static_cast<size_t>(n + 1) == 0 ? 3 : 0) - kernel[i];
    return out;
}

// The numerators 2^sub B^-1 of `basis`'s inverse, integers. Empty if one is 2^53 or more,
// where its double may have been rounded.
std::optional<std::vector<int64_t>> inverse_numerators(const SublatticeBasis& basis) {
    const double two_to_sub = std::ldexp(1.0, basis.sub);
    std::vector<int64_t> out;
    out.reserve(basis.from_power_basis.size());
    for (double x : basis.from_power_basis) {
        const double k = x * two_to_sub;
        if (!(std::fabs(k) < static_cast<double>(kInverseNumeratorLimit))) return std::nullopt;
        out.push_back(static_cast<int64_t>(k));
    }
    return out;
}

// Computes the preconditioned inverse of the new basis:
// B_{s+1} = B_s·Kernel·U so B_{s+1}^-1 = U^{-1}·Kernel^{-1}·B_s^-1
// here we compute (2*Kernel^{-1}) · (2^s·B_s^-1). Both are integer valued since Kernel contains
// halves and (2^s·B_s^-1) has entires k/2^s. Eventually we have to divide by 2^{s+1} to get the
// actual inverse, but we can keep the numerator as an integer
std::optional<std::vector<int64_t>> preconditioned_inverse(int n, const SublatticeBasis& basis,
                                                           const std::vector<int64_t>& kernel) {
    const std::optional<std::vector<int64_t>> numerators = inverse_numerators(basis);
    if (!numerators) return std::nullopt;
    return product(n, twice_kernel_inverse(n, kernel), *numerators);
}

// A basis B of the sublattice `sub` before polishing, the parent's basis times Kernel, and the
// numerators 2^sub B^-1 of its inverse.
struct Preconditioned {
    int sub = 0;
    std::vector<int64_t> preconditioned_basis;
    std::vector<int64_t> preconditioned_inverse;
};

// The parent's columns combined by Kernel into ones divisible by delta once more. Empty if
// the parent's columns have no two different nonzero labels or its inverse is not exact.
std::optional<Preconditioned> precondition(const GateSet& gateset, const SublatticeBasis& parent) {
    const int n = gateset.dim();
    const std::optional<std::vector<int>> labels = column_labels(gateset, parent);
    if (!labels) return std::nullopt;
    const std::optional<std::vector<int64_t>> kernel = parity_kernel(*labels);
    if (!kernel) return std::nullopt;
    // below 2^57: a row of 2 Kernel^-1 has at most n entries of size 1, or a single 2
    const std::optional<std::vector<int64_t>> inverse = preconditioned_inverse(n, parent, *kernel);
    if (!inverse) return std::nullopt;
    // B Kernel below 3 * 2^40: a column of Kernel adds at most three of B's, whose entries are
    // below 2^40
    return Preconditioned{parent.sub + 1, product(n, parent.to_power_basis, *kernel), *inverse};
}

// Computes the next sublattice basis from the preconditioned input. The basis is tailored to the
// input `ball` which is why we have tp provide it here.
std::optional<SublatticeBasis> next_sublattice_from_preconditioned(
    int n, const std::vector<double>& to_ball, const Preconditioned& start) {
    std::vector<double> start_basis(static_cast<size_t>(n) * n);
    exact_product(n, n, to_ball.data(), start.preconditioned_basis.data(), start_basis.data());
    LllTransform U;
    if (!lll_transform(n, lll_rows(n, start_basis.data(), kLllBits), kTransformLimit, U))
        return std::nullopt;
    const std::vector<__int128> basis =
        product(n, widened(start.preconditioned_basis), widened(U.T));
    const std::vector<__int128> numerators =
        product(n, widened(U.T_inv), widened(start.preconditioned_inverse));

    SublatticeBasis out{start.sub, {}, {}};
    out.to_power_basis.reserve(basis.size());
    out.from_power_basis.reserve(numerators.size());
    for (__int128 b : basis) {
        if (!(b < kTransformLimit && b > -kTransformLimit)) return std::nullopt;
        out.to_power_basis.push_back(static_cast<int64_t>(b));
    }
    // final step to get the real inverse (this was left over from the preconditioned inverse)
    for (__int128 k : numerators) out.from_power_basis.push_back(over_power_of_two(k, out.sub));
    return out;
}

// The double-double set-up's part of finish_setup(): to_ball_low B added to the identity
// place's rows of the reduced basis to_ball B, each entry of order one then and rounded twice;
// and the centre B^-1 (centre + centre_low) in double-double, B^-1's entries k / 2^sub being
// exact doubles (SublatticeBasis).
void add_low_parts(int n, const LoewnerEllipsoid& ellipsoid, const SublatticeBasis& basis,
                   std::vector<double>& reduced_basis, Setup& setup) {
    std::vector<double> low(static_cast<size_t>(kRowsPerPlace) * n);
    exact_product(n, kRowsPerPlace, ellipsoid.to_ball_low.data(), basis.to_power_basis.data(),
                  low.data());
    for (size_t k = 0; k < low.size(); ++k) reduced_basis[k] += low[k];
    setup.reduced_center.resize(n);
    setup.reduced_center_low.resize(n);
    for (int i = 0; i < n; ++i) {
        DD c{};
        for (int j = 0; j < n; ++j)
            c = c + DD{ellipsoid.centre[j], ellipsoid.centre_low[j]} *
                        basis.from_power_basis[static_cast<size_t>(i) * n + j];
        setup.reduced_center[i] = c.hi;
        setup.reduced_center_low[i] = c.lo;
    }
}

// The set-up of the search direction `target` at the scale K =
// `denominator_exponent`, from the region's Loewner ellipsoid around it and the
// level's basis B: the reduced basis to_ball B, formed exactly; the centre in reduced
// coordinates, B^-1 centre; and the R factors of the walk's forms (QR of to_ball B and of
// its row blocks: the whole ellipsoid, the cap block, the conjugate block and
// one block per conjugate place where there are several), with the Q of the
// conjugate forms, which walk_of() needs for their offsets. The region's
// dimension and radii go along, so that a walk needs nothing else. Valid at
// every scale of the same parity, where walk_of() rescales it.
std::shared_ptr<const Setup> finish_setup(const SearchRegion& region, const Vec4& target,
                                          int denominator_exponent,
                                          std::shared_ptr<const SublatticeBasis> basis,
                                          const LoewnerEllipsoid& ellipsoid) {
    const GateSet& gateset = region.gateset();
    const int n = gateset.dim(), conjugate_rows = gateset.conj_dim();
    const size_t nn = static_cast<size_t>(n) * n;
    auto setup = std::make_shared<Setup>();
    setup->reference_scale = denominator_exponent;
    setup->direction = target;
    setup->lattice_dimension = n;
    setup->ellipsoid_radius_sq = region.walk_radius_sq();
    setup->conjugate_radius_sq = region.conj_radius_sq();
    setup->basis = std::move(basis);
    // B's entries exceed to_ball's precision, so to_ball B is formed exactly.
    std::vector<double> reduced_basis(nn);
    exact_product(n, n, ellipsoid.to_ball.data(), setup->basis->to_power_basis.data(),
                  reduced_basis.data());
    if (ellipsoid.centre_low.empty())
        setup->reduced_center = ordered_matvec(n, setup->basis->from_power_basis, ellipsoid.centre);
    else
        add_low_parts(n, ellipsoid, *setup->basis, reduced_basis, *setup);
    setup->ellipsoid_upper = qr_positive(n, n, reduced_basis.data()).R;
    // Pruning forms (squared radii): cap block (1), conjugate block (P - 1,
    // origin-centred), one per conjugate place if several (1). They cut only
    // points the leaf test rejects.
    // rows of the reduced basis by place: the identity place (the cap), then
    // the conjugate places
    const auto place_rows = [&](int place) {
        return reduced_basis.data() + static_cast<size_t>(kRowsPerPlace * place) * n;
    };
    setup->cap_upper = qr_positive(kRowsPerPlace, n, place_rows(0)).R;
    setup->conjugate_block = qr_positive(conjugate_rows, n, place_rows(1));
    for (int a = 0; a < region.conjugate_forms(); ++a)
        setup->conjugate_places.push_back(qr_positive(kRowsPerPlace, n, place_rows(1 + a)));
    return setup;
}

}  // namespace

// Compute the next sublattice using a preconditioned basis so that the LLL reduction is fast.
// B_{s+1} = B_s·Kernel·U where the preconditioning computes B_s·Kernel and U is the LLL transform.
std::optional<SublatticeBasis> next_sublattice_basis(const GateSet& gateset,
                                                     const SublatticeBasis& parent,
                                                     const std::vector<double>& to_ball) {
    const std::optional<Preconditioned> start = precondition(gateset, parent);
    if (!start) return std::nullopt;
    return next_sublattice_from_preconditioned(gateset.dim(), to_ball, *start);
}

const std::vector<double>& SetupCache::embedding(const GateSet& gateset, int denominator_exponent) {
    auto it = embeddings_.find(denominator_exponent);
    if (it == embeddings_.end())
        it =
            embeddings_
                .emplace(denominator_exponent, quicksynth::embedding(gateset, denominator_exponent))
                .first;
    return it->second;
}

std::vector<std::shared_ptr<const Setup>> SetupCache::for_problems(
    const SearchRegion& region, int threads, const std::vector<Enumerator::Problem>& problems,
    int denominator_exponent, int sub) {
    // Plan: each problem's set-up is stored, built for an earlier problem of
    // this batch, or built for this one.
    const size_t problem_count = problems.size();
    std::vector<Key> keys(problem_count);
    std::vector<std::shared_ptr<const Setup>> out(problem_count);
    std::vector<Vec4> to_build;                         // their directions, in order
    std::unordered_map<Key, size_t, KeyHash> built_at;  // a key's index in to_build
    auto batch = setups_.batch();
    for (size_t i = 0; i < problem_count; ++i) {
        const Enumerator::Problem& p = problems[i];
        const Vec4 direction = search_direction(p.target, p.det_class, region.gateset());
        keys[i] = Key{sub_bits(direction, sub), p.det_class, denominator_exponent & 1};
        if (const auto* stored = batch.stored(keys[i])) {
            out[i] = *stored;
            continue;
        }
        if (built_at.emplace(keys[i], to_build.size()).second) to_build.push_back(direction);
        if (!batch.inserted(keys[i])) batch.insert(keys[i]);
    }
    // Compute, then commit.
    const std::vector<std::shared_ptr<const Setup>> built =
        this->build(region, threads, to_build, denominator_exponent, sub);
    batch.commit([&](const Key& key) { return built[built_at.at(key)]; });
    for (size_t i = 0; i < problem_count; ++i)
        if (!out[i]) out[i] = built[built_at.at(keys[i])];
    return out;
}

// Build the Setup for each direction in `directions`, at the scale K = `denominator_exponent` in
// parallel
std::vector<std::shared_ptr<const Setup>> SetupCache::build(const SearchRegion& region, int threads,
                                                            const std::vector<Vec4>& directions,
                                                            int denominator_exponent, int sub) {
    const size_t count = directions.size();
    const GateSet& gateset = region.gateset();
    const std::vector<double>& E = this->embedding(gateset, denominator_exponent);
    // The ellipsoid is the same at every level; a level narrows the lattice to its basis.
    const auto ellipsoid_of = [&](const Vec4& direction) {
        return region.ellipsoid(direction, denominator_exponent, E);
    };
    // Plan: each direction's basis is stored, reduced for an earlier direction
    // of this batch, or reduced for this one.
    std::vector<SubBits> keys(count);
    std::vector<std::shared_ptr<const SublatticeBasis>> bases(count);
    std::vector<size_t> to_reduce;  // in order
    // the direction reducing it
    std::unordered_map<SubBits, size_t, DirectionBitsHash> reduced_at;
    // Finer levels on Clifford+sqrt(T) derive a sub's basis from the next
    // coarser one (next_sublattice_basis): one fresh LLL per direction, at sub 0,
    // and every sub on the way cached. Clifford+T reduces each sub afresh (its
    // pinned replays depend on it).
    const bool derive = sub > 0 && gateset.step == 1;
    // per derived direction: the finest stored coarser basis (null: reduce sub 0
    // afresh), then the bases of the subs between
    std::vector<std::shared_ptr<const SublatticeBasis>> chain_from(count);
    std::vector<std::vector<std::shared_ptr<const SublatticeBasis>>> chain(count);
    auto batch = bases_.batch();
    for (size_t b = 0; b < count; ++b) {
        keys[b] = sub_bits(directions[b], sub);
        if (const auto* stored = batch.stored(keys[b])) {
            bases[b] = *stored;
            continue;
        }
        if (reduced_at.emplace(keys[b], b).second) {
            to_reduce.push_back(b);
            for (int s = sub - 1; derive && s >= 0 && !chain_from[b]; --s)
                if (const auto* coarser = batch.stored(sub_bits(directions[b], s)))
                    chain_from[b] = *coarser;
            const int first = chain_from[b] ? chain_from[b]->sub + 1 : 0;
            for (int s = first; derive && s < sub; ++s) batch.insert(sub_bits(directions[b], s));
        }
        if (!batch.inserted(keys[b])) batch.insert(keys[b]);
    }
    // a fresh reduction's sublattice, built once: sub 0 where derived, else this sub
    const Sublattice level = sublattice(gateset, derive ? 0 : sub);
    const auto derive_chain = [&](size_t b, const std::vector<double>& to_ball) {
        std::shared_ptr<const SublatticeBasis> basis = chain_from[b];
        if (!basis) {
            basis = reduce(gateset, to_ball, level);
            chain[b].push_back(basis);
        }
        for (int s = basis->sub + 1; s <= sub; ++s) {
            std::optional<SublatticeBasis> next = next_sublattice_basis(gateset, *basis, to_ball);
            basis = next ? std::make_shared<const SublatticeBasis>(std::move(*next))
                         : reduce(gateset, to_ball, sublattice(gateset, s));
            if (s < sub) chain[b].push_back(basis);
        }
        return basis;
    };
    // Compute the reductions in parallel, and commit them.
    std::vector<LoewnerEllipsoid> ellipsoids(count);
    parallel_for(to_reduce.size(), threads, [&](size_t r) {
        const size_t b = to_reduce[r];
        ellipsoids[b] = ellipsoid_of(directions[b]);
        bases[b] = derive ? derive_chain(b, ellipsoids[b].to_ball)
                          : reduce(gateset, ellipsoids[b].to_ball, level);
    });
    // fill in the bases of the directions reduced by an earlier one
    for (size_t b = 0; b < count; ++b)
        if (!bases[b]) bases[b] = bases[reduced_at.at(keys[b])];
    // the bases of the subs on the way, by key
    std::unordered_map<SubBits, std::shared_ptr<const SublatticeBasis>, DirectionBitsHash> chained;
    for (size_t b : to_reduce)
        for (const auto& basis : chain[b])
            chained.emplace(sub_bits(directions[b], basis->sub), basis);
    batch.commit([&](const SubBits& key) {
        const auto c = chained.find(key);
        return c != chained.end() ? c->second : bases[reduced_at.at(key)];
    });
    // Then every set-up, in parallel.
    std::vector<std::shared_ptr<const Setup>> out(count);
    parallel_for(count, threads, [&](size_t b) {
        if (ellipsoids[b].to_ball.empty()) ellipsoids[b] = ellipsoid_of(directions[b]);
        out[b] = finish_setup(region, directions[b], denominator_exponent, bases[b], ellipsoids[b]);
    });
    return out;
}

}  // namespace quicksynth
