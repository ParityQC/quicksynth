// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Cost and word of a Clifford+sqrt(T) element from the syllable normal form of
// its Bloch image.
//
// The normal form is unique, and for t_cost <= 2 q_cost its weight is the cost,
// so a greedy descent suffices. With lambda = 2 cos(pi/8),
// each step puts the two columns not divisible by lambda first and applies a
// Q syllable of sign s:
//   col0 <- u (col0 + s theta col1),  col1 <- s u (col1 - s theta col0),
//   col2 <- s lambda^3 col2,           u = 2 lambda^2 - 1, theta = lambda^2 - 3.
// It descends iff the common lambda-valuation d of the result is 2 mod 4 (new
// exponent k + 3 - d). Two descents on the same axis and sign form a
// T-syllable, priced min(t_cost, 2 q_cost). Exact because a normal form has
// Bloch exponent 3 n_Q + 2 n_T ([paper: prop:blochexp]; Forest et al., Thm. 4.1).
//
// word() replays the descent (Q for sign +1, X Q X for -1), spells the leftover
// Clifford as the first matching S^l H S^k H S^m H S^n, inverts, merges letters
// (merge_letters) and spells T as Q Q when 2 q_cost <= t_cost. Digits are int_mp:
// a Clifford+T image passes 2^128 near 1e-8. The word is verified to
// reproduce U up to a phase zeta^j.
#pragma once

#include "quicksynth/data_structures/bloch.hpp"
#include "quicksynth/data_structures/cost.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"

#include <cstdint>
#include <string>

namespace quicksynth {

// Methods throw std::invalid_argument if the input is not a unitary over the
// ring or its Bloch image.
class Descent {
  public:
    Descent(int64_t t_cost, int64_t q_cost) : t_(t_cost), q_(q_cost) {}

    // Normal-form weight of a 128-bit image (scan hits, try_bloch_image_128);
    // strips common powers of lambda first.
    Cost cost(BlochImage128 image) const;
    // Same for u / sqrt2^k, via its int_mp image.
    Cost cost(const M2& u, int k) const;
    Cost cost(const CyclotomicMatrixDigits& u) const;

    // Normal-form word of U = u / sqrt2^k.
    std::string word(const M2& u, int k) const;
    std::string word(const CyclotomicMatrixDigits& u) const;

  private:
    int64_t t_, q_;
};

}  // namespace quicksynth
