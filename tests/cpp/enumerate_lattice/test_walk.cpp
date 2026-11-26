// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Walks: solutions planted on the boundary of the cap's ellipsoid.
//
// An exact element U is the solution of every target U R(theta) with theta = 2 asin(rho eps),
// planted at rho from 0, the tip of the cap (an exact target), to 1 - 1e-10, its rim, about
// random axes and about z, which keeps a diagonal U's targets z-rotations. Every target the
// acceptance test (chord_distance) takes must come back from the walk at U's scale K, sub 0,
// and from the finest level of scale K + 1, where U lies in every sublattice. An element below
// the first scale is walked at the first scale of its parity instead, at sub 0 and the finest
// level. Both gate sets and both determinant classes, at precisions on both sides of the
// double-double set-up's threshold.

#include "oracles.hpp"

#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace quicksynth::test {

namespace {

// rho: the tip, inside the cap, and up to 1e-10 eps from its rim
const std::vector<double> kFractions{0.0,  1e-6,     1e-3,     0.3,      0.7,
                                     0.99, 1 - 1e-4, 1 - 1e-6, 1 - 1e-8, 1 - 1e-10};

struct Planting {
    std::vector<double> epsilons;
    int elements = 0;     // per gate set, precision and kind (at their scale, below the first)
    int axes = 0;         // random rotation axes per rho > 0, besides z
    double points = 0.0;  // the elements' scales: the first up to the last with this many
                          // estimated points
    uint64_t seed = 0;
};

// Whether the walk at (K, sub) for V in `det_class` returns the element of direction y.
bool returned(Enumerator& enumerator, const Unitary& V, int det_class, int K, int sub,
              const Vector4& y) {
    const std::vector<WalkResult> walks = enumerator.survivors_batch({{V, det_class}}, K, 0, sub);
    const size_t n = static_cast<size_t>(enumerator.gateset().dim());
    for (const CyclotomicMatrix& e : enumerator.matrices(
             walks[0].coefficients.data(), walks[0].coefficients.size() / n, K, det_class))
        if (gap_up_to_sign(matrix_direction(e), y) < 1e-14) return true;
    return false;
}

// A random word over `letters` with `length` letters.
std::string word_of(Rng& rng, const std::string& letters, int length) {
    std::string out;
    for (int i = 0; i < length; ++i)
        out += letters[static_cast<size_t>(rng.below(static_cast<int64_t>(letters.size())))];
    return out;
}

void plant(const Planting& p) {
    Rng rng(p.seed);
    int64_t required = 0;
    std::array<int64_t, 2> classes{};
    for (const GateSet& g : {GateSet::clifford_t(), GateSet::sqrt_t()}) {
        const bool sqrt_t = g.step == 1;
        const std::string letters = sqrt_t ? "HSTQ" : "HST";
        const char magic = sqrt_t ? 'Q' : 'T';
        const int finest = levels_per_scale(g) - 1;
        for (const double eps : p.epsilons) {
            Enumerator enumerator(eps, g, {}, 1, {});
            const int K0 = enumerator.first_scale();
            int K1 = K0;
            while (enumerator.estimate_points(K1 + 1) <= p.points) ++K1;
            for (const bool lifted : {false, true}) {
                if (lifted && K0 == 0) continue;
                int length = lifted ? 2 : 4 * K0, attempts = 0;
                for (int found = 0; found < p.elements;) {
                    REQUIRE(++attempts < 100000);
                    // lifted elements are also drawn without H, so diagonal ones occur
                    const bool diagonal = lifted && found % 2 == 1;
                    const std::string word =
                        word_of(rng, diagonal ? letters.substr(1) : letters, length);
                    const CyclotomicMatrix U = word_matrix(word);
                    const int K = U.denominator_exponent;
                    if (lifted ? K >= K0 : K > K1) {
                        length = std::max(1, length - 1);
                        continue;
                    }
                    if (!lifted && K < K0) {
                        ++length;
                        continue;
                    }
                    ++found;
                    const int det_class =
                        static_cast<int>(std::count(word.begin(), word.end(), magic) % 2);
                    const int walk_scale = K0 + ((K0 - K) & 1);
                    const std::vector<std::pair<int, int>> walks =
                        lifted ? std::vector<std::pair<int, int>>{{walk_scale, 0},
                                                                  {walk_scale, finest}}
                               : std::vector<std::pair<int, int>>{{K, 0}, {K + 1, finest}};
                    const Matrix M = matrix_value(U);
                    const Vector4 y = su2_of(M);
                    std::vector<std::array<double, 3>> axes{{0.0, 0.0, 1.0}};
                    for (int a = 0; a < p.axes; ++a)
                        axes.push_back({rng.normal(), rng.normal(), rng.normal()});
                    for (const double rho : kFractions)
                        for (const auto& axis : axes) {
                            if (rho == 0.0 && axis != axes.front()) continue;
                            const double theta = 2 * std::asin(rho * eps);
                            const Unitary V =
                                unitary_of(mul(M, rotation(theta, axis[0], axis[1], axis[2])));
                            // the walk returns what the acceptance test takes
                            if (chord_distance(matrix_direction(U), su2_vector(V)) > eps) continue;
                            ++required;
                            ++classes[det_class];
                            for (const auto& [K_walk, sub] : walks) {
                                INFO(g.name(), " eps ", eps, " word ", word, " rho ", rho,
                                     " walk K ", K_walk, " sub ", sub);
                                CHECK(returned(enumerator, V, det_class, K_walk, sub, y));
                            }
                        }
                }
            }
        }
    }
    CHECK(required > 0);
    CHECK(classes[0] > 0);
    CHECK(classes[1] > 0);
}

}  // namespace

TEST_CASE("the walk returns solutions planted on the cap's boundary") {
    plant({{1e-3, 1e-6, 3e-7, 1e-7, 1e-8, 1e-9}, 2, 1, 1e4, 1});
}

TEST_CASE("the walk returns many solutions planted on the cap's boundary" *
          doctest::test_suite("slow")) {
    for (const uint64_t seed : {2, 3})
        plant({{1e-3, 1e-4, 1e-5, 1e-6, 5e-7, 3e-7, 2e-7, 1e-7, 5e-8, 3e-8, 1e-8, 3e-9, 1e-9, 3e-10,
                1e-10, 1e-11, 1e-12, 1e-13},
               6,
               3,
               1e5,
               seed});
}

}  // namespace quicksynth::test
