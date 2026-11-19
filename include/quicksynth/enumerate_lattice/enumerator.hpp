// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Enumerates the elements of one denominator scale and determinant class that
// lie within epsilon of a target.
//
// The LLL-reduced search region (search_region.hpp) is walked with Fincke-Pohst (walk.hpp);
// survivors are verified exactly and priced. Per-direction set-ups are cached.
#pragma once

#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/enumerate_lattice/walk.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace quicksynth {

constexpr int64_t kTransformCache = 4096;

class SetupCache;  // enumerator.cpp

class Enumerator {
  public:
    struct Settings {
        // entries per cache (set-ups, LLL transforms)
        int64_t transform_cache = kTransformCache;
    };
    // An element within epsilon and its cost.
    struct Hit {
        Cost cost;
        CyclotomicMatrix element;
    };
    // One walk: one scale and one determinant class.
    struct Result {
        int K = 0;
        int sub = 0;  // sublattice delta^sub (SearchRegion::first_level)
        int det_class = 0;
        int64_t points = 0;  // leaves of the walk
        std::vector<Hit> hits;
        double seconds = 0.0;
    };
    struct Problem {
        Unitary target;
        int det_class = 0;
    };

    // `threads` splits a walk or a batch; survivors do not depend on it.
    Enumerator(double epsilon, GateSet gateset, Price price, int threads, Settings settings);
    ~Enumerator();
    Enumerator(Enumerator&&) noexcept;

    const SearchRegion& region() const { return region_; }
    double epsilon() const { return region_.epsilon(); }
    const GateSet& gateset() const { return region_.gateset(); }
    int threads() const { return threads_; }
    double estimate_points(int K) const { return region_.estimate_points(K); }
    double estimate_points(int K, int sub) const { return region_.estimate_points(K, sub); }
    int first_scale() const { return region_.first_scale(); }

    // Priced elements of class `det_class` at scale K within epsilon of V, in
    // walk order, and the leaf count.
    // `sub` > 0 walks only the sublattice delta^sub at scale K.
    Result enumerate(const Unitary& V, int K, int det_class, int sub = 0);

    // Walks many targets at scale K, each capped at `walk_budget` leaves (0 = no
    // cap). Per walk: flat ring-basis survivor rows and leaf count.
    std::vector<WalkResult> survivors_batch(const std::vector<Problem>& problems, int K,
                                            int64_t walk_budget, int sub = 0);

    // The unpriced elements su2_element(u1, u2, K) Q^det_class (su2.hpp), K the
    // denominator exponent, of the walk's survivors (n power-basis coefficients
    // each, u1's then u2's) with u1 u1* + u2 u2* = 2^K exactly.
    std::vector<CyclotomicMatrix> matrices(const int64_t* walk_survivors, size_t num_survivors,
                                           int denominator_exponent, int det_class) const;

  private:
    SearchRegion region_;
    Pricer pricer_;
    int threads_;
    std::unique_ptr<SetupCache> setups_;
};

}  // namespace quicksynth
