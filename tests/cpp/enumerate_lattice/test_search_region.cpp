// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Finer levels: the sublattice delta^sub of a scale, checked against level 0
// with an independent division by delta = 1 - zeta. The cap's ellipsoid: it
// holds the cap and is never larger than the box's.

#include "fixtures.hpp"

#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>

namespace quicksynth::test {

namespace {

using Ring = std::vector<int64_t>;  // d coefficients of Z[zeta_s], zeta_s^d = -1

Ring ring_of(const Z16& x, const GateSet& g) {
    Ring out(static_cast<size_t>(g.degree()));
    for (int i = 0; i < g.degree(); ++i)
        out[static_cast<size_t>(i)] = x[static_cast<size_t>(g.step * i)];
    return out;
}

Ring times_zeta(const Ring& x, int a) {
    const int d = static_cast<int>(x.size());
    Ring out(x.size(), 0);
    for (int i = 0; i < d; ++i) {
        const int j = (i + a) % (2 * d);
        out[static_cast<size_t>(j % d)] +=
            j < d ? x[static_cast<size_t>(i)] : -x[static_cast<size_t>(i)];
    }
    return out;
}

// sqrt2 = zeta_s^(d/4) - zeta_s^(3d/4)
Ring times_sqrt2(const Ring& x) {
    const int d = static_cast<int>(x.size());
    Ring a = times_zeta(x, d / 4), b = times_zeta(x, 3 * d / 4);
    for (size_t i = 0; i < a.size(); ++i) a[i] -= b[i];
    return a;
}

// Whether delta^r divides x: x = (1 - zeta) q has q_0 = (x_0 - sum_{i>0} x_i) / 2 and
// q_i = x_i + q_{i-1}.
bool divisible(Ring x, int r) {
    for (; r > 0; --r) {
        const int64_t rest = std::accumulate(x.begin() + 1, x.end(), int64_t{0});
        if ((x[0] - rest) % 2 != 0) return false;
        Ring q(x.size());
        q[0] = (x[0] - rest) / 2;
        for (size_t i = 1; i < x.size(); ++i) q[i] = x[i] + q[i - 1];
        x = q;
    }
    return true;
}

auto key_of(const CyclotomicMatrix& e) { return std::make_tuple(e.denominator_exponent, e.matrix); }

}  // namespace

TEST_CASE("a finer level holds the level-0 elements that delta^sub divides") {
    int proper = 0;
    for (const GateSet& g : {GateSet::sqrt_t(), GateSet::clifford_t()}) {
        const int L = levels_per_scale(g);
        Enumerator enumerator(1e-3, g, {1, 3}, 1, {});
        // walks the finest sub first, so its set-ups start from a fresh sub 0
        Enumerator finest_first(1e-3, g, {1, 3}, 1, {});
        for (int t = 0; t < 4; ++t) {
            const Unitary& V = calibration_targets()[static_cast<size_t>(t)];
            for (int K = enumerator.first_scale(); K < enumerator.first_scale() + 3; ++K)
                for (int det_class : {0, 1}) {
                    using Key = decltype(key_of(CyclotomicMatrix{}));
                    std::vector<std::set<Key>> finest(static_cast<size_t>(L));
                    for (int sub = L - 1; sub >= 1; --sub)
                        for (const Enumerator::Hit& h :
                             finest_first.enumerate(V, K, det_class, sub).hits)
                            finest[static_cast<size_t>(sub)].insert(key_of(h.element));
                    const auto all = enumerator.enumerate(V, K, det_class).hits;
                    for (int sub = 1; sub < L; ++sub) {
                        INFO(g.name(), " target ", t, " K ", K, " class ", det_class, " sub ", sub);
                        std::set<Key> want, got;
                        for (const Enumerator::Hit& h : all) {
                            Ring u1 = ring_of(h.element.matrix[0][0], g),
                                 u2 = ring_of(h.element.matrix[1][0], g);
                            for (int s = h.element.denominator_exponent; s < K; ++s) {
                                u1 = times_sqrt2(u1);
                                u2 = times_sqrt2(u2);
                            }
                            if (divisible(u1, sub) && divisible(u2, sub))
                                want.insert(key_of(h.element));
                        }
                        for (const Enumerator::Hit& h :
                             enumerator.enumerate(V, K, det_class, sub).hits)
                            got.insert(key_of(h.element));
                        CHECK(got == want);
                        CHECK(finest[static_cast<size_t>(sub)] == want);
                        proper += !want.empty() && want.size() < all.size();
                    }
                }
        }
    }
    CHECK(proper > 0);
}

TEST_CASE("levels count L per scale and step through the sublattices") {
    CHECK(levels_per_scale(GateSet::sqrt_t()) == 4);
    CHECK(levels_per_scale(GateSet::clifford_t()) == 2);
    CHECK(level_at(12, 4).K == 3);
    CHECK(level_at(12, 4).sub == 0);
    CHECK(level_at(13, 4).K == 4);
    CHECK(level_at(13, 4).sub == 3);
    const SearchRegion region(1e-5, GateSet::sqrt_t());
    CHECK(region.estimate_points(5, 2) == doctest::Approx(region.estimate_points(5) / 16));
    CHECK(level_at(region.first_level(1), 1).K == region.first_scale());
}

TEST_CASE("the cap's ellipsoid holds the cap and is never larger than the box's") {
    // (eps, its semi-axis across over eps): the dome's down to 1e-6, the cylinder's at 3e-7, the
    // box's below
    for (const auto& [eps, across] : std::vector<std::pair<double, double>>{{1e-2, 1.033},
                                                                            {1e-3, 1.033},
                                                                            {1e-5, 1.033},
                                                                            {1e-6, 1.062},
                                                                            {3e-7, 1.155},
                                                                            {1e-7, 2.0},
                                                                            {1e-9, 2.0}}) {
        INFO("eps ", eps);
        const SearchRegion region(eps, GateSet::sqrt_t());
        const CapEllipsoid& cap = region.cap_ellipsoid();
        CHECK(cap.across / eps == doctest::Approx(across).epsilon(1e-3));
        // points of the cap: t = sqrt(1 - r^2) for r from 0 (the tip) to eps (the rim)
        for (const double f : {0.0, 0.25, 0.5, 0.75, 0.999, 1.0}) {
            const double r = f * eps, t = std::sqrt(1.0 - r * r);
            const double along = (t - cap.centre) / cap.along, across_v = r / cap.across;
            CHECK(along * along + across_v * across_v <= 1.0);
        }
        // the box's volume prod(2 h), h floored as in the region
        const double s = std::sqrt(1.0 - eps * eps), h = std::max((1.0 - s) / 2, 2.4e-15);
        CHECK(cap.along * std::pow(cap.across, 3) <= (2 * h) * std::pow(2 * eps, 3));
    }
}

}  // namespace quicksynth::test
