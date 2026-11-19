// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Single-scale enumerator: embedding, known elements, exact filter, threads.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/enumerate_lattice/lll.hpp"
#include "quicksynth/enumerate_lattice/setup.hpp"
#include "quicksynth/enumerate_lattice/walk.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace quicksynth::test {

namespace {

int parity_of(const std::string& word, char magic) {
    return static_cast<int>(std::count(word.begin(), word.end(), magic) % 2);
}

// distance compared against epsilon
double distance_of(const Enumerator::Hit& h, const Unitary& V) {
    return chord_distance(matrix_direction(h.element), su2_vector(V));
}

// The basis of sub 0 from a fresh reduction of A: T, with B^-1 = T^-1.
SublatticeBasis fresh_basis(int n, const std::vector<double>& A) {
    LllTransform t;
    REQUIRE(lll_transform(n, lll_rows(n, A.data(), 48), int64_t{1} << 40, t));
    return {0, t.T, std::vector<double>(t.T_inv.begin(), t.T_inv.end())};
}

// X Y for n x n integer matrices, exact.
std::vector<int64_t> integer_product(int n, const std::vector<int64_t>& X,
                                     const std::vector<int64_t>& Y) {
    std::vector<int64_t> out(static_cast<size_t>(n) * n, 0);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            for (int j = 0; j < n; ++j) out[r * n + c] += X[r * n + j] * Y[j * n + c];
    return out;
}

// Sum of log column norms of A T; with the volume fixed, larger is less reduced.
double log_norms(int n, const std::vector<double>& A, const std::vector<int64_t>& T) {
    std::vector<double> B(static_cast<size_t>(n) * n);
    exact_product(n, n, A.data(), T.data(), B.data());
    double sum = 0.0;
    for (int c = 0; c < n; ++c) {
        double sq = 0.0;
        for (int r = 0; r < n; ++r) sq += B[r * n + c] * B[r * n + c];
        sum += 0.5 * std::log(sq);
    }
    return sum;
}

}  // namespace

TEST_CASE("the embedding lists the conjugates") {
    // column i is zeta^(step i) / sqrt2^K at every place
    for (const GateSet& g : {GateSet::sqrt_t(), GateSet::clifford_t()})
        for (int K : {0, 1, 4}) {
            INFO(g.name(), " K ", K);
            const std::vector<double> E = embedding(g, K);
            const int n = g.dim(), d = g.degree();
            REQUIRE(E.size() == static_cast<size_t>(n * n));
            for (int p = 0; p < g.num_places(); ++p)
                for (int i = 0; i < d; ++i) {
                    const int a = g.places[static_cast<size_t>(p)];
                    const long double angle = std::acos(-1.0L) * ((a * g.step * i) % 16) / 8;
                    const long double scale = std::pow(2.0L, -K / 2.0L);
                    const double re = static_cast<double>(std::cos(angle) * scale);
                    const double im = static_cast<double>(std::sin(angle) * scale);
                    CHECK(std::abs(E[(4 * p) * n + i] - re) <= 1e-15);
                    CHECK(std::abs(E[(4 * p + 1) * n + i] - im) <= 1e-15);
                    CHECK(std::abs(E[(4 * p + 2) * n + d + i] - re) <= 1e-15);
                    CHECK(std::abs(E[(4 * p + 3) * n + d + i] - im) <= 1e-15);
                }
        }
}

TEST_CASE("the embedding of an element is its matrix at the identity place") {
    const CyclotomicMatrix U = word_matrix("HQSHTQHSQ");
    const std::vector<double> E = embedding(GateSet::sqrt_t(), U.denominator_exponent);
    std::vector<long double> coeff;
    for (int64_t c : U.matrix[0][0]) coeff.push_back(c);
    for (int64_t c : U.matrix[1][0]) coeff.push_back(c);
    std::vector<long double> y(16, 0);
    for (int r = 0; r < 16; ++r)
        for (int c = 0; c < 16; ++c) y[r] += E[r * 16 + c] * coeff[c];
    const Matrix M = matrix_value(U);
    // (Re u1, Im u1, Re u2, Im u2) of M at the identity place
    const std::vector<long double> want{M[0].real(), M[0].imag(), M[2].real(), M[2].imag()};
    for (int i = 0; i < 4; ++i) CHECK(std::abs(y[i] - want[i]) < 1e-12L);
    for (int p = 0; p < 4; ++p) {
        long double norm = 0;
        for (int i = 0; i < 4; ++i) norm += y[4 * p + i] * y[4 * p + i];
        CHECK(std::abs(norm - 1) < 1e-9L);  // U is unitary at every place
    }
}

TEST_CASE("a next sublattice basis spans the sublattice and is as reduced as a fresh one") {
    // log column-norm sums against a fresh reduction: equal on average, a single
    // basis up to 0.95 longer in 600 sampled directions
    double excess = 0.0;
    int count = 0;
    for (const GateSet& g : {GateSet::sqrt_t(), GateSet::clifford_t()})
        for (double eps : {1e-3, 1e-6, 1e-9}) {
            const SearchRegion region(eps, g);
            const int n = g.dim(), K = region.first_scale() + 1;
            for (size_t t = 0; t < 4; ++t)
                for (int det_class : {0, 1}) {
                    const Vec4 v = search_direction(calibration_targets()[t], det_class, g);
                    const std::vector<double> A = region.ellipsoid(v, K, embedding(g, K)).to_ball;
                    SublatticeBasis parent = fresh_basis(n, A);
                    for (int sub = 1; sub < levels_per_scale(g); ++sub) {
                        INFO(std::string(g.name()), " eps ", eps, " target ", t, " class ",
                             det_class, " sub ", sub);
                        const std::optional<SublatticeBasis> next =
                            next_sublattice_basis(g, parent, A);
                        REQUIRE(next);
                        REQUIRE(next->sub == sub);
                        // B = M_sub T with T = M_sub^-1 B and T^-1 = B^-1 M_sub integers: B
                        // spans what M_sub spans, delta^sub Z[zeta]^2; and B B^-1 = 1
                        const std::vector<int64_t> M = delta_power(g, sub);
                        const std::vector<double> M_inv = inverse_delta_power(g, sub);
                        for (int r = 0; r < n; ++r)
                            for (int c = 0; c < n; ++c) {
                                // with the numerators 2^sub M_sub^-1 and 2^sub B^-1, exact
                                __int128 t_rc = 0, t_inv_rc = 0, one = 0;
                                for (int j = 0; j < n; ++j) {
                                    const auto m_inv =
                                        static_cast<__int128>(std::ldexp(M_inv[r * n + j], sub));
                                    const auto b_inv = static_cast<__int128>(
                                        std::ldexp(next->from_power_basis[r * n + j], sub));
                                    t_rc += m_inv * next->to_power_basis[j * n + c];
                                    t_inv_rc += b_inv * M[j * n + c];
                                    one += next->to_power_basis[r * n + j] *
                                           static_cast<__int128>(
                                               std::ldexp(next->from_power_basis[j * n + c], sub));
                                }
                                CHECK(t_rc % (__int128{1} << sub) == 0);
                                CHECK(t_inv_rc % (__int128{1} << sub) == 0);
                                CHECK(one == (r == c ? __int128{1} << sub : 0));
                            }
                        std::vector<double> AM(static_cast<size_t>(n) * n);
                        exact_product(n, n, A.data(), M.data(), AM.data());
                        LllTransform fresh;
                        REQUIRE(
                            lll_transform(n, lll_rows(n, AM.data(), 48), int64_t{1} << 40, fresh));
                        const double more = log_norms(n, A, next->to_power_basis) -
                                            log_norms(n, A, integer_product(n, M, fresh.T));
                        CHECK(more < 1.4);
                        excess += more;
                        ++count;
                        parent = *next;
                    }
                }
        }
    CHECK(excess / count < 0.1);
}

TEST_CASE("a next sublattice basis rejects parents it cannot divide exactly") {
    const GateSet g = GateSet::sqrt_t();
    const SearchRegion region(1e-6, g);
    const int n = g.dim(), K = region.first_scale() + 1;
    const Vec4 v = search_direction(calibration_targets()[0], 0, g);
    const std::vector<double> A = region.ellipsoid(v, K, embedding(g, K)).to_ball;
    const SublatticeBasis parent = fresh_basis(n, A);
    REQUIRE(next_sublattice_basis(g, parent, A));
    // 2 B: every column has even parity, so no pivot pair
    SublatticeBasis doubled = parent;
    for (int64_t& b : doubled.to_power_basis) b *= 2;
    CHECK_FALSE(next_sublattice_basis(g, doubled, A));
    // B^-1 scaled to 2^62: the derived inverse's numerators pass 2^53
    SublatticeBasis huge = parent;
    for (double& b : huge.from_power_basis) b = b == 0 ? 0 : std::copysign(std::ldexp(1.0, 62), b);
    CHECK_FALSE(next_sublattice_basis(g, huge, A));
}

TEST_CASE("a transform cache cleared within a batch of finer levels keeps the survivors") {
    // capacity 1 clears the cache between the subs of one derived chain
    std::vector<Enumerator::Problem> problems;
    for (size_t t = 0; t < 6; ++t)
        problems.push_back({calibration_targets()[t], static_cast<int>(t % 2)});
    Enumerator large(1e-3, GateSet::sqrt_t(), {1, 3}, 1, {});
    Enumerator tiny(1e-3, GateSet::sqrt_t(), {1, 3}, 1, {1});
    const int K = large.first_scale() + 2;
    size_t survivors = 0;
    for (int sub : {3, 1, 2}) {
        INFO("sub ", sub);
        const std::vector<WalkResult> a = large.survivors_batch(problems, K, 0, sub);
        const std::vector<WalkResult> b = tiny.survivors_batch(problems, K, 0, sub);
        REQUIRE(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            const size_t n = static_cast<size_t>(GateSet::sqrt_t().dim());
            std::set<std::vector<int64_t>> rows_a, rows_b;
            for (size_t at = 0; at < a[i].coefficients.size(); at += n)
                rows_a.insert({a[i].coefficients.begin() + at, a[i].coefficients.begin() + at + n});
            for (size_t at = 0; at < b[i].coefficients.size(); at += n)
                rows_b.insert({b[i].coefficients.begin() + at, b[i].coefficients.begin() + at + n});
            CHECK(rows_a == rows_b);
            survivors += rows_a.size();
        }
    }
    CHECK(survivors > 0);
}

TEST_CASE("the walk finds a known element at its scale and determinant class") {
    const std::string word = "HQSHTQHSQ";
    const CyclotomicMatrix U = word_matrix(word);
    const Unitary V = matrix_unitary(U);
    const int det_class = parity_of(word, 'Q');
    Enumerator enumerator(1e-3, GateSet::sqrt_t(), {1, 3}, 1, {});
    const Enumerator::Result r = enumerator.enumerate(V, U.denominator_exponent, det_class);
    CHECK(r.points > 0);
    bool found = false;
    for (const Enumerator::Hit& h : r.hits) {
        found = found || distance_of(h, V) < 1e-12;
        CHECK(distance_of(h, V) < 1e-3);
        CHECK(h.cost.value >= 0);
    }
    CHECK(found);
    for (const Enumerator::Hit& h :
         enumerator.enumerate(V, U.denominator_exponent, 1 - det_class).hits)
        CHECK(distance_of(h, V) > 1e-9);
}

TEST_CASE("the batched walk finds the single walk's elements") {
    const std::string word = "HQSHTQHSQTHQ";
    const CyclotomicMatrix U = word_matrix(word);
    const Unitary V = matrix_unitary(U);
    const int det_class = parity_of(word, 'Q'), K = U.denominator_exponent;
    Enumerator single(1e-3, GateSet::sqrt_t(), {}, 1, {});
    const Enumerator::Result walked = single.enumerate(V, K, det_class);
    Enumerator batched(1e-3, GateSet::sqrt_t(), {}, 1, {});
    const std::vector<WalkResult> b = batched.survivors_batch({{V, det_class}}, K, 0);
    REQUIRE(b.size() == 1);
    CHECK(b[0].visited == walked.points);
    const size_t n = static_cast<size_t>(GateSet::sqrt_t().dim());
    const std::vector<CyclotomicMatrix> elements =
        batched.matrices(b[0].coefficients.data(), b[0].coefficients.size() / n, K, det_class);
    std::set<json> rows;
    for (const CyclotomicMatrix& e : elements) rows.insert(row_of(e));
    bool exact = false;
    for (const Enumerator::Hit& h : walked.hits) {
        CHECK(rows.count(row_of(h.element)) == 1);
        exact = exact || distance_of(h, V) < 1e-12;
    }
    CHECK(exact);
}

TEST_CASE("every hit is within epsilon") {
    // target 1 has a cost-27 element at 1.07e-4 that the slack pre-filter passes
    Enumerator enumerator(1e-4, GateSet::sqrt_t(), {}, 1, {});
    std::vector<Enumerator::Hit> hits;
    for (int K : {10, 11})
        for (int det_class : {0, 1}) {
            const Enumerator::Result r =
                enumerator.enumerate(calibration_targets()[1], K, det_class);
            hits.insert(hits.end(), r.hits.begin(), r.hits.end());
        }
    REQUIRE_FALSE(hits.empty());  // first-hit scale
    for (const Enumerator::Hit& h : hits) {
        CHECK(distance_of(h, calibration_targets()[1]) <= 1e-4);
        CHECK(h.cost.value >= 28);  // certified optimum at this epsilon
    }
}

TEST_CASE("the threaded walk returns the same survivors in the same order") {
    // scale 9: a few million points
    const Unitary& V = calibration_targets()[0];
    Enumerator serial(1e-3, GateSet::sqrt_t(), {}, 1, {});
    Enumerator parallel(1e-3, GateSet::sqrt_t(), {}, 5, {});
    bool any = false;
    for (int det_class : {0, 1}) {
        const Enumerator::Result a = serial.enumerate(V, 9, det_class);
        const Enumerator::Result b = parallel.enumerate(V, 9, det_class);
        CHECK(a.points == b.points);
        REQUIRE(a.hits.size() == b.hits.size());
        for (size_t i = 0; i < a.hits.size(); ++i) {
            CHECK(a.hits[i].cost.value == b.hits[i].cost.value);
            CHECK(row_of(a.hits[i].element) == row_of(b.hits[i].element));
        }
        any = any || !a.hits.empty();
    }
    CHECK(any);  // a scale with hits
}

TEST_CASE("on Clifford+T the walk finds a known element in its determinant class") {
    // T shifts the determinant exponent by 1, H and S by an even amount
    const std::string word = "HTSHTHTSHTHT";
    const CyclotomicMatrix U = word_matrix(word);
    const Unitary V = matrix_unitary(U);
    Enumerator enumerator(1e-2, GateSet::clifford_t(), {1, 3}, 1, {});
    std::set<int> exact;
    for (int det_class = 0; det_class < GateSet::kDeterminantClasses; ++det_class)
        for (const Enumerator::Hit& h :
             enumerator.enumerate(V, U.denominator_exponent, det_class).hits)
            if (distance_of(h, V) < 1e-12) exact.insert(det_class % 2);
    CHECK(exact == std::set<int>{parity_of(word, 'T')});
}

TEST_CASE("the walk honours its leaf budget") {
    Enumerator enumerator(0.05, GateSet::sqrt_t(), {1, 3}, 1, {});
    const std::vector<WalkResult> b =
        enumerator.survivors_batch({{calibration_targets()[37], 0}}, 6, 500);
    REQUIRE(b.size() == 1);
    CHECK(b[0].visited <= 501);
}

TEST_CASE("a capped walk on several threads uses its whole budget") {
    // a ball with every leaf kept: 938 leaves split into many small work items,
    // so a walker that holds a large block of the budget at the end would cut
    // the others short
    const int n = 12;
    std::vector<double> R(n * n, 0.0), center(n, 0.25), embedding(n * n, 0.0);
    std::vector<int64_t> T(n * n, 0);
    for (int i = 0; i < n; ++i) {
        R[i * n + i] = 1.0;
        T[i * n + i] = 1;
    }
    const double v[4] = {1.0, 0.0, 0.0, 0.0};
    const WalkLattice lattice{n, embedding.data(), -1e300, 1e300};
    WalkInputs inputs;
    inputs.forms.push_back({R.data(), center.data(), nullptr, 3.0});
    inputs.transform = T.data();
    inputs.v = v;
    const int64_t leaves = walk(lattice, inputs, 0, 1).visited;
    REQUIRE(leaves == 938);
    for (int64_t budget : {300, 1000})
        for (int threads : {1, 4, 8}) {
            const int64_t visited = walk(lattice, inputs, budget, threads).visited;
            // never more than the budget, plus the leaf each walker stops at; on
            // threads short of it only by blocks a preempted walker still holds
            CHECK(visited <= budget + threads);
            if (threads == 1)
                CHECK(visited == std::min(budget + 1, leaves));
            else
                CHECK(visited >= std::min(budget * 9 / 10, leaves));
        }
}

}  // namespace quicksynth::test
