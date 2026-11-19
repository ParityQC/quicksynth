// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// LLL reduction: T unimodular with its exact inverse, and the reduced basis
// checked in exact integers, for the float path and the exact fallback.
// The exact product A T of the walk basis and a transform.

#include "oracles.hpp"

#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/enumerate_lattice/lll.hpp"
#include "quicksynth/data_structures/int_mp.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <doctest/doctest.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

namespace quicksynth::test {

namespace {

constexpr int64_t kLimit = int64_t{1} << 40;
constexpr int kScanBits = 48;  // the scan's LLL input precision (enumerator.cpp)
// Input precision past 2^120, where the exact fallback reduces.
constexpr int kWideBits = kScanBits + 80;

// T X == I in 128-bit integers; n x n row-major.
bool is_identity_product(int n, const std::vector<int64_t>& T, const std::vector<int64_t>& X) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            __int128 s = 0;
            for (int l = 0; l < n; ++l) s += static_cast<__int128>(T[i * n + l]) * X[l * n + j];
            if (s != (i == j ? 1 : 0)) return false;
        }
    return true;
}

std::vector<int_mp> to_int_mp(const std::vector<double>& rows) {
    std::vector<int_mp> out;
    for (double x : rows) out.push_back(int_mp::from_double(x));
    return out;
}

// The rows of B T, where row j of `rows` is column j of B. Row j of the result
// is sum_i T[i][j] row i.
std::vector<int_mp> transformed_rows(int n, const std::vector<int_mp>& rows,
                                     const std::vector<int64_t>& T) {
    std::vector<int_mp> out(rows.size());
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            for (int c = 0; c < n; ++c) out[j * n + c] += rows[i * n + c] * T[i * n + j];
    return out;
}

// Whether the rows are (0.99, eta_percent / 100)-LLL-reduced, in exact
// integers. d_i is the Gram determinant of the first i rows and
// lambda_kj = d_{j+1} mu_kj (Cohen, Algorithm 2.6.7, step 2).
bool lll_reduced(int n, const std::vector<int_mp>& b, int eta_percent) {
    std::vector<int_mp> d(static_cast<size_t>(n) + 1), lambda(b.size());
    d[0] = int_mp(1);
    for (int k = 0; k < n; ++k)
        for (int j = 0; j <= k; ++j) {
            int_mp u;
            for (int c = 0; c < n; ++c) u += b[k * n + c] * b[j * n + c];
            for (int i = 0; i < j; ++i)
                u = (u * d[i + 1] - lambda[k * n + i] * lambda[j * n + i]) / d[i];
            (j < k ? lambda[k * n + j] : d[k + 1]) = u;
        }
    for (int k = 0; k < n; ++k) {
        // |mu_kj| <= eta: 100 |lambda_kj| <= eta_percent d_{j+1}
        for (int j = 0; j < k; ++j) {
            const int_mp lhs = lambda[k * n + j] * 100, rhs = d[j + 1] * eta_percent;
            if (lhs > rhs || lhs < -rhs) {
                MESSAGE("|mu(", k, ", ", j, ")| > eta");
                return false;
            }
        }
        // (delta - mu^2) B_{k-1} <= B_k: 99 d_k^2 - 100 lambda^2 <= 100 d_{k+1} d_{k-1}
        if (k == 0) continue;
        const int_mp& lam = lambda[k * n + k - 1];
        if (d[k] * d[k] * 99 - lam * lam * 100 > d[k + 1] * d[k - 1] * 100) {
            MESSAGE("Lovasz fails at ", k);
            return false;
        }
    }
    return true;
}

// T T_inv = I both ways and B T (0.99, eta_percent / 100)-reduced.
void check_reduction(int n, const std::vector<int_mp>& rows, int eta_percent) {
    LllTransform t;
    REQUIRE(lll_transform(n, rows, kLimit, t));
    CHECK(is_identity_product(n, t.T, t.T_inv));
    CHECK(is_identity_product(n, t.T_inv, t.T));
    CHECK(lll_reduced(n, transformed_rows(n, rows, t.T), eta_percent));
}

// n columns, one `thin` times shorter in a skew direction, like scan regions.
std::vector<double> skew_basis(int n, double thin, Rng& rng) {
    std::vector<double> A(static_cast<size_t>(n * n));
    for (double& a : A) a = rng.normal();
    for (int i = 0; i < n; ++i) A[static_cast<size_t>(i * n)] *= thin;
    return A;
}

}  // namespace

TEST_CASE("the LLL reduces random skew bases") {
    Rng rng(29);
    for (int n : {4, 8, 16})
        for (double thin : {1.0, 1e-4, 1e-8}) {
            INFO(n, " ", thin);
            const std::vector<double> rows = lll_rows(n, skew_basis(n, thin, rng).data(), 60);
            check_reduction(n, to_int_mp(rows), 51);
            // the double rows give the same transform as their int_mp copy
            LllTransform a, b;
            REQUIRE(lll_transform(n, rows, kLimit, a));
            REQUIRE(lll_transform(n, to_int_mp(rows), kLimit, b));
            CHECK(a.T == b.T);
        }
}

TEST_CASE("the exact fallback reduces bases too wide for 128 bits") {
    Rng rng(37);
    for (int n : {4, 8, 16})
        for (int bits : {127, 1000}) {
            INFO(n, " ", bits);
            check_reduction(n, to_int_mp(lll_rows(n, skew_basis(n, 1e-3, rng).data(), bits)), 50);
        }
}

TEST_CASE("the LLL reduces the scan's bases") {
    Rng rng(31);
    for (const GateSet& g : {GateSet::sqrt_t(), GateSet::clifford_t()})
        for (double eps : {1e-2, 1e-4, 1e-6, 1e-8}) {
            const SearchRegion region(eps, g);
            for (int dK : {0, 3, 8}) {
                const int K = region.first_scale() + dK;
                const Vec4 v = search_direction(unitary_of(haar_unitary(rng)), 0, g);
                const std::vector<double> A = region.ellipsoid(v, K, embedding(g, K)).to_ball;
                INFO(g.name(), " ", eps, " K ", K);
                check_reduction(g.dim(), to_int_mp(lll_rows(g.dim(), A.data(), kScanBits)), 51);
                check_reduction(g.dim(), to_int_mp(lll_rows(g.dim(), A.data(), kWideBits)), 50);
            }
        }
}

TEST_CASE("entries at 2^120 go to the exact fallback") {
    // b1 - 5 b0 = (2^e, -2^e); the cancelled product <b0, b1> = 10 needs exact_dot.
    for (int e : {119, 121}) {
        INFO(e);
        const int_mp top = shift_left(int_mp(1), e);
        const std::vector<int_mp> rows{int_mp(1), int_mp(1), top + int_mp(5), int_mp(5) - top};
        check_reduction(2, rows, 50);
        LllTransform t;
        REQUIRE(lll_transform(2, rows, kLimit, t));
        CHECK(t.T == std::vector<int64_t>{1, -5, 0, 1});
    }
}

TEST_CASE("linearly dependent rows throw") {
    // the float path meets 0 / 0 at mu(1, 0) and hands over to the exact one
    LllTransform t;
    CHECK_THROWS_AS(lll_transform(2, std::vector<double>{1, 2, 2, 4}, kLimit, t),
                    std::invalid_argument);
}

TEST_CASE("the int_mp inverse covers transforms whose inverse passes 64 bits") {
    // Rows of 140, 94, 48 and 48 random bits, skewed like the z-rotation body
    // at eps = 1e-13; some have T below 2^62 but T^-1 past 2^63.
    std::mt19937_64 rng(5);
    const auto random_bits = [&](int bits) {
        int_mp x;
        for (int b = 0; b < bits; b += 32) x = shift_left(x, 32) + int_mp(int64_t(rng() >> 32));
        return (rng() & 1) != 0 ? -x : x;
    };
    int wide = 0;
    for (int trial = 0; trial < 20; ++trial) {
        std::vector<int_mp> rows(16);
        for (int i = 0; i < 4; ++i)
            for (int c = 0; c < 4; ++c)
                rows[i * 4 + c] = random_bits(c == 0 ? 140 : c == 1 ? 94 : 48);
        LllTransform small;
        LllTransformMp t;
        REQUIRE(lll_transform(4, rows, int64_t{1} << 62, t));
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                int_mp sum;
                for (int l = 0; l < 4; ++l) sum += t.T_inv[i * 4 + l] * t.T[l * 4 + j];
                CHECK(sum == int_mp(i == j ? 1 : 0));
            }
        CHECK(lll_reduced(4, transformed_rows(4, rows, t.T), 51));
        if (lll_transform(4, rows, int64_t{1} << 62, small)) {
            CHECK(small.T == t.T);
            for (size_t k = 0; k < t.T_inv.size(); ++k) CHECK(t.T_inv[k] == int_mp(small.T_inv[k]));
        } else {
            ++wide;
        }
    }
    CHECK(wide > 0);
}

TEST_CASE("an LLL transform beyond the limit is refused") {
    // reducing (1, 0), (2^50 + 1, 1) subtracts 2^50 + 1 times b0
    const std::vector<double> rows{1, 0, 0x1p50 + 1, 1};
    LllTransform t;
    CHECK_FALSE(lll_transform(2, rows, kLimit, t));
    CHECK(lll_transform(2, rows, int64_t{1} << 52, t));
}

TEST_CASE("the product A T is exact where doubles would cancel") {
    // A spans 2^40, A T is O(1). Entries are dyadic, so the exact product is
    // representable and must match bit for bit.
    const int n = 3;
    const std::vector<double> A{0x1p40, 1.0, 0.5, 0.25, 0x1p-10, 3.0, 1.0, 1.0, 1.0};
    const std::vector<int64_t> T{1, 0, 0, -(int64_t{1} << 40), 1, 0, 0, 0, 1};
    std::vector<double> out(static_cast<size_t>(n * n));
    exact_product(n, n, A.data(), T.data(), out.data());
    const std::vector<double> want{0x1p40 - 0x1p40, 1.0, 0.5, 0.25 - 0x1p30, 0x1p-10, 3.0,
                                   1.0 - 0x1p40,    1.0, 1.0};
    CHECK(out == want);
}

TEST_CASE("the product A T is exact where 128-bit sums would overflow") {
    // Row 0 of the golden-122 basis (Clifford+T, scale 29); column 7 has
    // 128-bit terms and partial sums but an O(1) result.
    const int n = 8;
    const std::vector<double> row{0x1.c970c67133e81p+32,
                                  0x1.437594a40fc9p+32,
                                  0x1.ef21427624ef2p-22,
                                  -0x1.437594a40fc9p+32,
                                  0.0,
                                  -0x1.0bf65eaa33dd8p+31,
                                  -0x1.7af4beb229a07p+31,
                                  -0x1.0bf65eaa33dd8p+31};
    const std::vector<int64_t> T{0,     0,     0,     0,     5741,  -2378, -1016705, -1306372,
                                 -2378, 0,     3363,  0,     -5741, 0,     705059,   945027,
                                 3363,  0,     -4756, 0,     2378,  2378,  19601,    -30098,
                                 -2378, 0,     3363,  0,     2378,  -3363, -732779,  -902462,
                                 0,     -3363, 0,     -4756, -2378, 2378,  2378,     2378,
                                 0,     2378,  0,     3363,  3363,  -2378, 1640982,  2359901,
                                 0,     0,     0,     0,     -2378, 985,   -2323077, -3339782,
                                 0,     -2378, 0,     -3363, 0,     985,   1644345,  2363264};
    std::vector<double> out(n);
    exact_product(n, 1, row.data(), T.data(), out.data());
    // exact sums rounded once (Python fractions)
    const std::vector<double> want{0x1.96861ec160731p-10, 0.0,
                                   -0x1.1f74af35d7315p-9, 0.0,
                                   0x1.68b9d2bcd75ccp-4,  0x1.d8cc295e6bae6p-3,
                                   -0x1.21baa1deb9b92p+8, 0x1.993140bdaea06p+8};
    CHECK(out == want);
}

TEST_CASE("the product A T is exact where it exceeds 128 bits") {
    // 2^60 (2^20 + 1) + 1 needs 133 bits at the row's scale 2^52
    const int n = 2;
    const std::vector<double> A{0x1p60, 1.0};
    const std::vector<int64_t> T{(int64_t{1} << 20) + 1, 0, 1, 1};
    std::vector<double> out(n);
    exact_product(n, 1, A.data(), T.data(), out.data());
    CHECK(out == std::vector<double>{0x1p80 + 0x1p60, 1.0});
}

}  // namespace quicksynth::test
