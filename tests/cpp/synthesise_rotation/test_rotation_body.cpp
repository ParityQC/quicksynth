// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The rotation body against two independent enumerations: an exhaustive search
// over the coefficient ball sum c_i^2 <= 2^K (the trace form bounds every body
// point there), and the u entries of the joint walk's hits at the same
// (K, sub, class), which also pins the class convention. Also the enclosing
// ellipse, the level counts pi eps^3 2^l / 12 and pi^3 eps^3 2^l / 768,
// segment corners on the unit circle down to eps = 1e-13, and the reuse of
// one LLL across scales.

#include "quicksynth/synthesise_rotation/rotation_body.hpp"

#include <doctest/doctest.h>

#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"
#include "quicksynth/solve_norm_equation/relative_norm.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <set>
#include <vector>

namespace quicksynth::test {

namespace {

template <size_t N>
using Ring = RingElement<N, __int128>;

constexpr long double kPiL = 3.141592653589793238462643383279502884L;

template <size_t N>
GateSet gateset_of() {
    return N == 4 ? GateSet::clifford_t() : GateSet::sqrt_t();
}

Unitary z_rotation(double theta) {
    return {std::cos(theta / 2), -std::sin(theta / 2), 0, 0, 0, 0,
            std::cos(theta / 2), std::sin(theta / 2)};
}

template <size_t N>
Ring<N> delta_power(int e) {
    Ring<N> out = ring_zeta_power<N, __int128>(0);
    const Ring<N> delta =
        ring_sub(ring_zeta_power<N, __int128>(0), ring_zeta_power<N, __int128>(1));
    for (int i = 0; i < e; ++i) out = ring_mul(out, delta);
    return out;
}

template <size_t N>
bool divides(int e, const Ring<N>& u) {
    return ring_divide(u, delta_power<N>(e)).has_value();
}

// The points of the body, after the exact test.
template <size_t N>
std::set<Ring<N>> body_points(const RotationBody<N>& body, int K) {
    std::set<Ring<N>> out;
    for (const auto& p : body.points(K))
        if (body.contains(p.u, K)) out.insert(p.u);
    return out;
}

// Every point of points(K) once, in delta^sub O_K, with the right flag.
template <size_t N>
void check_points(const RotationBody<N>& body, int K) {
    const auto points = body.points(K);
    std::set<Ring<N>> seen;
    for (const auto& p : points) {
        CHECK(seen.insert(p.u).second);
        CHECK(divides<N>(body.sub(), p.u));
        CHECK(p.new_at_level == !divides<N>(body.sub() + 1, p.u));
    }
}

// Exhaustive: every u of delta^sub O_K in the body at scale K. The outer
// coefficients run over the ball sum c_i^2 <= 2^K; c_0, real at every place,
// then over the interval the segment and the strip |Im(u a*)| <= R eps allow,
// widened by a slack, and the exact test decides.
template <size_t N>
std::set<Ring<N>> exhaustive(const RotationBody<N>& body, const Vec4& a, int K) {
    const long double R = std::sqrt(std::ldexp(1.0L, K)), eps = body.epsilon();
    const long double norm = std::hypot(static_cast<long double>(a[0]), a[1]);
    const long double ar = a[0] / norm, ai = a[1] / norm;
    // sigma_0(zeta^i) a*
    std::array<long double, N> re{}, im{};
    for (size_t i = 0; i < N; ++i) {
        const long double angle = kPiL * static_cast<long double>(i) / N;
        const long double c = std::cos(angle), s = std::sin(angle);
        re[i] = c * ar + s * ai;
        im[i] = s * ar - c * ai;
    }
    const long double floor_re = R * std::sqrt(1 - eps * eps), half = R * eps;
    const long double slack = 1e-6L * R + 1e-6L;
    std::set<Ring<N>> out;
    Ring<N> u{};
    const int64_t total = int64_t{1} << K;
    auto visit = [&](auto&& self, size_t i, int64_t remaining) -> void {
        const auto top = static_cast<int64_t>(std::sqrt(static_cast<long double>(remaining)));
        if (i > 0) {
            for (int64_t c = -top; c <= top; ++c) {
                if (c * c > remaining) continue;
                u[i] = c;
                self(self, i - 1, remaining - c * c);
            }
            u[i] = 0;
            return;
        }
        long double sr = 0, si = 0;
        for (size_t j = 1; j < N; ++j) {
            sr += static_cast<long double>(u[j]) * re[j];
            si += static_cast<long double>(u[j]) * im[j];
        }
        long double lo = -static_cast<long double>(top), hi = static_cast<long double>(top);
        // Re: c0 ar + sr >= floor_re; Im: |si - c0 ai| <= half
        if (std::fabs(ar) > 1e-3L) {
            const long double edge = (floor_re - slack - sr) / ar;
            if (ar > 0)
                lo = std::max(lo, edge);
            else
                hi = std::min(hi, edge);
        }
        if (std::fabs(ai) > 1e-3L) {
            long double e1 = (si - half - slack) / ai, e2 = (si + half + slack) / ai;
            if (e1 > e2) std::swap(e1, e2);
            lo = std::max(lo, e1);
            hi = std::min(hi, e2);
        }
        for (auto c = static_cast<int64_t>(std::ceil(lo));
             c <= static_cast<int64_t>(std::floor(hi)); ++c) {
            if (c * c > remaining) continue;
            const long double x = sr + c * ar, y = si - c * ai;
            if (x < floor_re - slack || std::fabs(y) > half + slack) continue;
            u[0] = c;
            if (divides<N>(body.sub(), u) && body.contains(u, K)) out.insert(u);
        }
        u[0] = 0;
    };
    visit(visit, N - 1, total);
    return out;
}

// The u entries, at scale K, of the joint walk's hits at (K, sub, class).
template <size_t N>
std::set<Ring<N>> walk_entries(Enumerator& enumerator, const Unitary& V, int K, int j, int sub) {
    const int step = 8 / static_cast<int>(N);
    std::set<Ring<N>> out;
    for (const Enumerator::Hit& h : enumerator.enumerate(V, K, j, sub).hits) {
        Z16 u = h.element.matrix[0][0];
        const Z16 sqrt2{0, 0, 1, 0, 0, 0, -1, 0};
        for (int k = h.element.denominator_exponent; k < K; ++k) u = zeta_mul(u, sqrt2);
        Ring<N> r{};
        for (int t = 0; t < kZetaDegree; ++t) {
            if (t % step == 0)
                r[t / step] = u[t];
            else
                CHECK(u[t] == 0);
        }
        out.insert(r);
    }
    return out;
}

// The body points whose xi = 2^K - u u* has a solution t t* = xi.
template <size_t N>
std::set<Ring<N>> solvable(const std::set<Ring<N>>& us, int K, PrimeIdeals<N>& primes) {
    std::set<Ring<N>> out;
    for (const Ring<N>& u : us) {
        RealElement<N / 2, __int128> power{};
        power[0] = static_cast<__int128>(1) << K;
        const NormSolutions<N> s =
            solve_relative_norm<N>(real_sub(power, relative_norm(u)), primes);
        REQUIRE(s.outcome != NormOutcome::undecided);
        if (s.outcome == NormOutcome::solved) out.insert(u);
    }
    return out;
}

// Levels l = L K - sub, L = N / 2, from `first` to `last`.
template <class F>
void for_levels(int N, int first, int last, F&& f) {
    const int L = N / 2;
    for (int l = first; l <= last; ++l) {
        const int K = (l + L - 1) / L;
        f(K, L * K - l);
    }
}

// The first level with one new point expected per class.
int first_useful_level(int N, double eps) {
    const double per_level = N == 4 ? M_PI / 12 : std::pow(M_PI, 3) / 768;
    return static_cast<int>(std::ceil(std::log2(1 / (per_level * eps * eps * eps))));
}

template <size_t N>
void check_exhaustive(double eps, int last_level, int angles) {
    const GateSet g = gateset_of<N>();
    std::mt19937_64 rng(1000 + N);
    std::uniform_real_distribution<double> angle(-M_PI, M_PI);
    for (int trial = 0; trial < angles; ++trial) {
        const Unitary V = z_rotation(angle(rng));
        for (int j = 0; j < GateSet::kDeterminantClasses; ++j) {
            const Vec4 a = search_direction(V, j, g);
            for_levels(static_cast<int>(N), 0, last_level, [&](int K, int sub) {
                INFO("eps ", eps, " trial ", trial, " class ", j, " K ", K, " sub ", sub);
                const RotationBody<N> body(a, sub, eps);
                check_points(body, K);
                CHECK(body_points(body, K) == exhaustive(body, a, K));
            });
        }
    }
}

template <size_t N>
void check_walk(double eps, int first_level, int last_level, int angles) {
    const GateSet g = gateset_of<N>();
    Enumerator enumerator(eps, g, N == 4 ? Price::clifford_t(1) : Price{1, 3}, 1, {});
    PrimeIdeals<N> primes;
    std::mt19937_64 rng(2000 + N);
    std::uniform_real_distribution<double> angle(-M_PI, M_PI);
    int64_t hits = 0;
    for (int trial = 0; trial < angles; ++trial) {
        const Unitary V = z_rotation(angle(rng));
        for (int j = 0; j < GateSet::kDeterminantClasses; ++j) {
            const Vec4 a = search_direction(V, j, g);
            for_levels(static_cast<int>(N), first_level, last_level, [&](int K, int sub) {
                INFO("eps ", eps, " trial ", trial, " class ", j, " K ", K, " sub ", sub);
                const RotationBody<N> body(a, sub, eps);
                const std::set<Ring<N>> walked = walk_entries<N>(enumerator, V, K, j, sub);
                CHECK(solvable<N>(body_points(body, K), K, primes) == walked);
                hits += static_cast<int64_t>(walked.size());
            });
        }
    }
    CHECK(hits > 0);
}

// Mean new body points per (class, level) over the angles, against the area
// of the segment times the conjugate disks over the covolume, halved for
// "new" (N(delta) = 2).
template <size_t N>
void check_level_count(double eps, int l, int angles) {
    const GateSet g = gateset_of<N>();
    const int L = static_cast<int>(N / 2);
    const int K = (l + L - 1) / L, sub = L * K - l;
    const double expected =
        (N == 4 ? M_PI / 12 : std::pow(M_PI, 3) / 768) * std::pow(eps, 3) * std::ldexp(1.0, l);
    std::mt19937_64 rng(3000 + N);
    std::uniform_real_distribution<double> angle(-M_PI, M_PI);
    int64_t count = 0;
    for (int trial = 0; trial < angles; ++trial) {
        const Unitary V = z_rotation(angle(rng));
        for (int j = 0; j < GateSet::kDeterminantClasses; ++j) {
            const RotationBody<N> body(search_direction(V, j, g), sub, eps);
            for (const auto& p : body.points(K)) count += p.new_at_level && body.contains(p.u, K);
        }
    }
    const double mean = static_cast<double>(count) / (2.0 * angles);
    INFO("level ", l, " expected ", expected, " mean ", mean);
    CHECK(std::fabs(mean / expected - 1) < 0.1);
}

// u0 = zeta^m sqrt2^K lies on the boundary circle at every place. Directions
// at angle asin(0.99 eps) from it put u0 next to a corner of the segment, and
// at asin(1.01 eps) just outside. The scales stop below the first useful one
// by `below`, where the ellipsoid's slab at such a direction (see points())
// stays under a few thousand points.
template <size_t N>
void check_corners(double eps, int below) {
    const int L = static_cast<int>(N / 2);
    const Ring<N> sqrt2 =
        N == 4 ? Ring<N>{0, 1, 0, -1}
               : ring_sub(ring_zeta_power<N, __int128>(2), ring_zeta_power<N, __int128>(6));
    const int top = first_useful_level(static_cast<int>(N), eps) / L - below;
    std::mt19937_64 rng(4000 + N);
    for (int K = std::max(1, top - 3); K <= top; ++K) {
        Ring<N> power = ring_zeta_power<N, __int128>(0);
        for (int k = 0; k < K; ++k) power = ring_mul(power, sqrt2);
        for (int trial = 0; trial < 4; ++trial) {
            const int m = static_cast<int>(rng() % (2 * N));
            const int sub = static_cast<int>(rng() % L);
            const Ring<N> u0 = ring_times_zeta(power, m);
            for (double side : {-1.0, 1.0}) {
                INFO("eps ", eps, " K ", K, " m ", m, " sub ", sub, " side ", side);
                const long double base = kPiL * m / N;
                const long double inside = std::asin(0.99L * eps), outside = std::asin(1.01L * eps);
                const long double phi = base + side * inside, psi = base + side * outside;
                const Vec4 a{static_cast<double>(std::cos(phi)), static_cast<double>(std::sin(phi)),
                             0, 0};
                const Vec4 b{static_cast<double>(std::cos(psi)), static_cast<double>(std::sin(psi)),
                             0, 0};
                const RotationBody<N> body(a, sub, eps);
                CHECK(body.contains(u0, K));
                const auto points = body.points(K, 20000);
                CHECK(points.size() <= 20000);
                CHECK(std::any_of(points.begin(), points.end(),
                                  [&](const auto& p) { return p.u == u0; }));
                CHECK_FALSE(RotationBody<N>(b, sub, eps).contains(u0, K));
            }
        }
    }
}

template <size_t N>
void check_reduction_scale(double eps, int last_level) {
    const GateSet g = gateset_of<N>();
    std::mt19937_64 rng(5000 + N);
    std::uniform_real_distribution<double> angle(-M_PI, M_PI);
    for (int trial = 0; trial < 3; ++trial) {
        const Unitary V = z_rotation(angle(rng));
        const Vec4 a = search_direction(V, trial % 2, g);
        for_levels(static_cast<int>(N), 0, last_level, [&](int K, int sub) {
            INFO("K ", K, " sub ", sub);
            const RotationBody<N> shared(a, sub, eps), fresh(a, sub, eps, K);
            std::set<Ring<N>> x, y;
            for (const auto& p : shared.points(K)) x.insert(p.u);
            for (const auto& p : fresh.points(K)) y.insert(p.u);
            CHECK(x == y);
        });
    }
}

}  // namespace

TEST_CASE("the enclosing ellipse covers the circular segment") {
    // In units of R, from the chord: the radial coordinate (cos phi - c) / h,
    // h = 1 - c, with cos phi - c = h - 2 sin^2(phi / 2) free of cancellation.
    for (long double eps : {1e-13L, 1e-10L, 1e-6L, 1e-3L, 0.1L, 0.5L, 0.9L, 0.999L}) {
        const long double c = std::sqrt(1 - eps * eps), h = eps * eps / (1 + c);
        const auto form = [&](long double radial, long double tangential) {
            const long double x = (radial - kSegmentCentre) / kSegmentRadial;
            const long double y = tangential / kSegmentTangential;
            return x * x + y * y;
        };
        long double worst = 0;
        const long double top = std::asin(eps);
        for (int i = 0; i <= 20000; ++i) {
            const long double phi = top * (2.0L * i / 20000 - 1);
            const long double s = std::sin(phi / 2);
            worst = std::max(worst, form(1 - 2 * s * s / h, std::sin(phi) / eps));
            worst = std::max(worst, form(0, 2.0L * i / 20000 - 1));
        }
        INFO("eps ", static_cast<double>(eps));
        CHECK(worst < 0.99975L);
        CHECK(worst > 0.999L);
    }
}

TEST_CASE("the rotation body equals an exhaustive search on Clifford+T") {
    check_exhaustive<4>(1e-2, first_useful_level(4, 1e-2) + 1, 3);
}

TEST_CASE("the rotation body equals an exhaustive search on Clifford+T at 1e-3" *
          doctest::test_suite("slow")) {
    check_exhaustive<4>(1e-3, first_useful_level(4, 1e-3), 1);
}

TEST_CASE("the rotation body equals an exhaustive search on Clifford+sqrt(T)") {
    check_exhaustive<8>(1e-2, 20, 2);
}

TEST_CASE("the rotation body equals an exhaustive search at the first useful sqrt(T) level" *
          doctest::test_suite("slow")) {
    check_exhaustive<8>(1e-2, first_useful_level(8, 1e-2), 2);
}

TEST_CASE("the solvable rotation body points are the joint walk's hits on Clifford+T") {
    check_walk<4>(1e-2, 0, first_useful_level(4, 1e-2) + 2, 8);
    check_walk<4>(1e-3, 0, first_useful_level(4, 1e-3) + 1, 4);
}

TEST_CASE("the solvable rotation body points are the joint walk's hits on Clifford+sqrt(T)") {
    check_walk<8>(1e-2, 0, first_useful_level(8, 1e-2) + 2, 6);
    check_walk<8>(1e-3, 0, first_useful_level(8, 1e-3) + 1, 3);
}

TEST_CASE("the solvable rotation body points are the joint walk's hits at 1e-6" *
          doctest::test_suite("slow")) {
    check_walk<4>(1e-6, first_useful_level(4, 1e-6) - 2, first_useful_level(4, 1e-6) + 1, 2);
    check_walk<8>(1e-6, first_useful_level(8, 1e-6) - 2, first_useful_level(8, 1e-6) + 1, 1);
}

TEST_CASE("the new rotation body points per level match the area") {
    check_level_count<4>(1e-2, 30, 8);
    check_level_count<4>(1e-2, 29, 8);
    check_level_count<8>(1e-2, 33, 4);
    check_level_count<8>(1e-2, 34, 4);
}

TEST_CASE("segment corners on the boundary circle are found down to 1e-13") {
    for (double eps : {1e-3, 1e-6, 1e-10, 1e-13}) {
        check_corners<4>(eps, 10);
        check_corners<8>(eps, 3);
    }
}

TEST_CASE("a point limit cuts the list after one more point") {
    // a = 1 at eps 1e-3: the slab of the exact direction (see points())
    const RotationBody<4> body(Vec4{1, 0, 0, 0}, 0, 1e-3);
    const size_t all = body.points(20).size();
    REQUIRE(all > 100);
    CHECK(body.points(20, 100).size() == 101);
    CHECK(body.points(20, static_cast<int64_t>(all)).size() == all);
}

TEST_CASE("a walk visits the new points of the list in order and stops on request") {
    const RotationBody<4> body(Vec4{1, 0, 0, 0}, 0, 1e-3);
    const auto all = body.points(20);
    std::vector<Ring<4>> fresh;
    for (const auto& p : all)
        if (p.new_at_level) fresh.push_back(p.u);
    REQUIRE(fresh.size() > 10);
    std::vector<Ring<4>> visited;
    const auto collect = [&](const RotationBody<4>::Point& p) {
        CHECK(p.new_at_level);
        visited.push_back(p.u);
        return true;
    };
    CHECK(body.walk(20, kNoPointLimit, false, collect) == static_cast<int64_t>(all.size()));
    CHECK(visited == fresh);
    // stopped at the tenth new point: every point up to it was walked
    int seen = 0;
    const int64_t walked = body.walk(20, kNoPointLimit, false,
                                     [&](const RotationBody<4>::Point&) { return ++seen < 10; });
    CHECK(seen == 10);
    int64_t upto = 0, news = 0;
    while (news < 10) news += all[static_cast<size_t>(upto++)].new_at_level;
    CHECK(walked == upto);
}

TEST_CASE("one LLL serves every scale") {
    check_reduction_scale<4>(1e-2, first_useful_level(4, 1e-2) + 1);
    check_reduction_scale<8>(1e-2, first_useful_level(8, 1e-2) + 1);
}

TEST_CASE("the rotation body checks its arguments") {
    const Vec4 a{1, 0, 0, 0};
    CHECK_THROWS_AS(RotationBody<4>(a, 2, 1e-3), std::invalid_argument);
    CHECK_THROWS_AS(RotationBody<8>(a, -1, 1e-3), std::invalid_argument);
    CHECK_THROWS_AS(RotationBody<4>(a, 0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(RotationBody<4>(a, 0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(RotationBody<4>(Vec4{0, 0, 1, 0}, 0, 1e-3), std::invalid_argument);
    CHECK_THROWS_AS(RotationBody<4>(a, 0, 1e-3, -1), std::invalid_argument);
}

}  // namespace quicksynth::test
