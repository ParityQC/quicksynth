// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The relative norm equation t t* = xi for xi in O_K+, K = Q(zeta_2N), N = 4
// (Clifford+T) or 8 (Clifford+sqrt(T)) [Ross, Selinger, arXiv:1403.2975,
// section 6; Kliuchnikov, Bocharov, Roetteler, Yard, arXiv:1510.03888].
//
// Both rings have class number 1, the units of relative norm 1 are the roots
// of unity (Hasse unit index 1) and every totally positive unit is a square
// [Washington, GTM 83, Theorem 4.12 and Corollary 4.13 and section 11.1;
// Weber]. So a solution exists iff xi is 0, or totally positive with even
// valuation at every K+ prime inert in K/K+. The solutions up to roots of unity
// are then
//   t = eta delta^v2 prod_inert pi^(v/2) prod_split pi^a pi*^(v-a), 0 <= a <= v,
// one per choice of the exponents a, with delta = 1 - zeta above 2 and eta the
// unit square root that makes t t* = xi exact.
//
// N_{K+/Q}(xi) is factored in 64 bits (integer_factor.hpp). Valuations come
// from exact trial division by the prime generators, which works for every
// residue degree. Odd primes p = -1 mod 2N are the inert ones; an odd power of
// p in the norm rejects xi before any prime is generated.
#pragma once

#include "quicksynth/data_structures/number_field.hpp"
#include "quicksynth/solve_norm_equation/prime_ideals.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace quicksynth {

#pragma GCC visibility push(hidden)

enum class NormOutcome {
    solved,
    none,       // no t in O_K has t t* = xi
    undecided,  // past what the solver decides; see NormUndecided
};

enum class NormUndecided {
    no,
    norm_too_large,  // N_{K+/Q}(xi) without its power of 2 is 2^63 or more
    factor_gave_up,  // Pollard-Brent left a composite cofactor
    no_generator,    // no generator found for a prime above p
    digit_overflow,  // a solution or an intermediate passed 127 bits
    unit_root,       // the unit square root is beyond double-double rounding
};

// Choices returned per xi; realistic xi have at most a few dozen, and a cut
// list is reported as capped.
inline constexpr int kChoiceCap = 64;

template <size_t N>
struct NormSolutions {
    NormOutcome outcome = NormOutcome::none;
    NormUndecided undecided = NormUndecided::no;
    // Number of choices, prod (v + 1) over the split pairs, saturating at
    // INT64_MAX; halved, rounded up, with transpose_halving.
    int64_t count = 0;
    // t with t t* = xi exactly, one per choice, the first min(count, cap) in
    // mixed-radix order of the exponents a, the first split pair fastest.
    std::vector<RingElement<N, __int128>> choices;
    bool capped = false;  // count > cap
};

// The solutions of t t* = xi.
//
// transpose_halving keeps the first half of the choices: choice c and choice
// count - 1 - c (every a -> v - a) give t and zeta^m t* for some m. A candidate
// [[u, -t* w], [t, u* w]] and its transpose, whose lower-left entry is -t* w,
// have the same distance to a diagonal target and the same cost, since H, S,
// T and sqrt(T) are their own transposes. The halving is exact only for a
// caller that tries every phase t zeta^i of each choice it keeps (the cost
// test in test_relative_norm.cpp prices both halves).
template <size_t N>
NormSolutions<N> solve_relative_norm(const RealElement<N / 2, __int128>& xi, PrimeIdeals<N>& primes,
                                     int cap = kChoiceCap, bool transpose_halving = false);

#pragma GCC visibility pop

}  // namespace quicksynth
