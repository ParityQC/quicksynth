// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Polynomials over F_p of degree at most 8 and the factorisation of the
// cyclotomic polynomials Phi_8 = x^4 + 1 and Phi_16 = x^8 + 1 mod p. Each
// irreducible factor h of Phi_{2n} mod p gives one prime (p, h(zeta)) of
// Z[zeta] above p (Kummer-Dedekind), all of residue degree f, the order of p
// mod 2n [Washington, GTM 83, Theorem 2.13]. The factors come in closed form;
// the tests check them against the general splitting (fp_poly_oracle.hpp).
#pragma once

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace quicksynth {

#pragma GCC visibility push(hidden)

// A polynomial over F_p, coefficients c[0] .. c[degree] reduced mod p and
// c[degree] nonzero; degree -1 is the zero polynomial. The prime p is passed
// to each operation.
struct FpPoly {
    static constexpr int kMaxDegree = 8;
    std::array<uint64_t, kMaxDegree + 1> c{};
    int degree = -1;
    bool operator==(const FpPoly& o) const { return degree == o.degree && c == o.c; }
    bool operator!=(const FpPoly& o) const { return !(*this == o); }
    bool operator<(const FpPoly& o) const;  // by degree, then coefficients from the top
};

// (q, r) with a = q b + r and deg r < deg b; p prime, std::domain_error for b = 0.
std::pair<FpPoly, FpPoly> poly_divmod(const FpPoly& a, const FpPoly& b, uint64_t p);

// The residue degree of p in Q(zeta_{2n}): the least f with p^f = 1 mod 2n;
// p odd, n a power of two (std::invalid_argument otherwise).
int cyclotomic_residue_degree(uint64_t p, int n);

// The irreducible monic factors, ascending, of x^n + 1 mod an odd prime p,
// n = 4 or 8 (std::invalid_argument otherwise; primality is not checked). All
// have degree cyclotomic_residue_degree(p, n). They are binomials x^k - w for
// p = 1 mod 4 and trinomials for p = 3 mod 4, with coefficients from roots of
// unity and square roots mod p in closed form per class of p mod 16 [Meyn,
// Finite Fields Appl. 2 (1996)], which costs a few modular powers instead of
// the O(log p) polynomial powers of a random splitting.
std::vector<FpPoly> cyclotomic_factors(int n, uint64_t p);

#pragma GCC visibility pop

}  // namespace quicksynth
