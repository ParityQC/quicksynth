// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/price_candidates/price.hpp"

#include "quicksynth/exact_synthesis/descent.hpp"
#include "quicksynth/data_structures/bloch.hpp"

#include <optional>
#include <stdexcept>
#include <utility>

namespace quicksynth {

namespace {

// The T-count is half the Bloch exponent, which is even exactly on Clifford+T.
Cost by_t_count(const Price& price, int64_t exponent) {
    if (exponent % 2 != 0)
        throw std::invalid_argument("Bloch image is not over Z[1/sqrt2]: not Clifford+T");
    return price.cost(exponent / 2, 0);
}

}  // namespace

Pricer Pricer::for_gateset(const GateSet& gateset, const Price& price) {
    return gateset.step != 1 ? t_count(price.t) : descent(price.t, price.q);
}

Pricer Pricer::descent(int64_t t_cost, int64_t q_cost) {
    return Pricer(false, Price{t_cost, q_cost});
}

Pricer Pricer::t_count(int64_t t_cost) { return Pricer(true, Price{t_cost, 0}); }

// Uses 128-bit digits where they suffice, arbitrary-size digits otherwise.
Cost Pricer::price(const CyclotomicMatrix& e) const {
    if (t_count_) {
        const std::optional<int64_t> exponent =
            try_bloch_exponent_128(e.matrix, e.denominator_exponent);
        return by_t_count(price_,
                          exponent ? *exponent : bloch_exponent(e.matrix, e.denominator_exponent));
    }
    const Descent descent(price_.t, price_.q);
    std::optional<BlochImage128> b = try_bloch_image_128(e.matrix, e.denominator_exponent);
    if (!b) return descent.cost(e.matrix, e.denominator_exponent);
    return descent.cost(std::move(*b));
}

Cost Pricer::price(const CyclotomicMatrixDigits& u) const {
    if (t_count_) return by_t_count(price_, bloch_exponent(u));
    return Descent(price_.t, price_.q).cost(u);
}

}  // namespace quicksynth
