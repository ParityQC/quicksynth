// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The cost of an exact element at a price. On
// Clifford+sqrt(T) it is the weight of the normal form of its Bloch image
// (descent.hpp); on Clifford+T the price is the T-count and nothing else, so it
// is t_cost times the T-count, half the Bloch exponent. Both are read in
// 128-bit digits where they suffice and in digits of any size where they do
// not, which give the same cost.
#pragma once

#include "quicksynth/data_structures/cost.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"

#include <cstdint>

namespace quicksynth {

class Pricer {
  public:
    // By the T-count on Clifford+T, by the normal form on Clifford+sqrt(T).
    static Pricer for_gateset(const GateSet& gateset, const Price& price);
    static Pricer descent(int64_t t_cost, int64_t q_cost);
    // std::invalid_argument for an element that is not Clifford+T.
    static Pricer t_count(int64_t t_cost);
    Cost price(const CyclotomicMatrix& e) const;
    // Same, for an element given by arbitrary-size digits.
    Cost price(const CyclotomicMatrixDigits& u) const;

  private:
    Pricer(bool t_count, const Price& price) : t_count_(t_count), price_(price) {}
    bool t_count_;
    Price price_;
};

}  // namespace quicksynth
