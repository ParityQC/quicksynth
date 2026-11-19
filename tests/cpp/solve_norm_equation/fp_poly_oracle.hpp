// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// F_p polynomial arithmetic and the general factorisation of a squarefree
// polynomial: distinct-degree, then Cantor-Zassenhaus equal-degree splitting
// [Cohen, GTM 138, section 3.4] with a fixed-seed generator, so that the
// factors are the same on every run. The oracle for cyclotomic_factors
// (prime_field.hpp).
#pragma once

#include "quicksynth/solve_norm_equation/prime_field.hpp"

#include <cstdint>
#include <initializer_list>
#include <utility>
#include <vector>

namespace quicksynth::test {

// The polynomial with the given coefficients from c[0] up, reduced mod p;
// std::length_error past degree 8.
FpPoly fp_poly(std::initializer_list<uint64_t> coefficients, uint64_t p);

FpPoly poly_add(const FpPoly& a, const FpPoly& b, uint64_t p);
FpPoly poly_sub(const FpPoly& a, const FpPoly& b, uint64_t p);
// a b; std::length_error if the product passes degree 8.
FpPoly poly_mul(const FpPoly& a, const FpPoly& b, uint64_t p);
// a b mod m and a^e mod m, deg m >= 1.
FpPoly poly_mulmod(const FpPoly& a, const FpPoly& b, const FpPoly& m, uint64_t p);
FpPoly poly_powmod(const FpPoly& a, uint64_t e, const FpPoly& m, uint64_t p);
// The monic gcd; zero for a = b = 0.
FpPoly poly_gcd(const FpPoly& a, const FpPoly& b, uint64_t p);
FpPoly poly_monic(const FpPoly& a, uint64_t p);

// For a monic squarefree f: the pairs (d, product of the irreducible factors
// of f of degree d), by increasing d, nonempty products only.
std::vector<std::pair<int, FpPoly>> distinct_degree_factors(const FpPoly& f, uint64_t p);

// The irreducible monic factors, ascending, of a monic squarefree f whose
// irreducible factors all have degree d; p odd prime. Throws
// std::runtime_error if the splitting does not converge, which only an input
// that breaks these conditions causes.
std::vector<FpPoly> equal_degree_factors(const FpPoly& f, int d, uint64_t p);

}  // namespace quicksynth::test
