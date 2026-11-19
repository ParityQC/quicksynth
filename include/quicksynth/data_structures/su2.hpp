// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// SU(2) in double precision. Operation order is fixed and part of the result:
// lattice reduction reads search_direction, recorded distances chord_distance.
#pragma once

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"

#include <array>

namespace quicksynth {

// A 2x2 complex matrix, row-major, as (re, im) pairs.
using Unitary = std::array<double, 8>;
using Vec4 = std::array<double, 4>;  // SU(2) vector (Re a, Im a, Re b, Im b)

// (Re, Im) of zeta16^j, from square roots alone.
std::array<double, 2> zeta_power(int j);

// 2^(n/2): a power of two, times sqrt(2) or sqrt(1/2) for odd n.
double sqrt2_power(int n);

// Global phase and determinant classes. With zeta_n the gate set's root of unity
// (zeta16 on Clifford+sqrt(T), zeta8 on Clifford+T) and Q = diag(1, zeta_n),
// every exact answer is
//
//     U = zeta_n^p S Q^j:
//
// a ring phase zeta_n^p, a global phase and so irrelevant; an exact S in SU(2),
// [[u1, -u2*], [u2, u1*]] / sqrt2^k, which is what the lattice search finds, from
// its first column (su2_element(), cyclotomic_matrix.hpp); and j in {0, 1}, the
// determinant class. det U = zeta_n^(2p + j): a ring phase changes the
// determinant by even powers only, so Q^j supplies the odd one. Hence U is within
// epsilon of V up to phase exactly when S is within epsilon of V Q^-j up to
// sign: the search runs once per class, around search_direction(V, j), and
// multiplies Q^j back onto what it finds (times_q_power(), cyclotomic_matrix.hpp).
// Unitaries are compared up to phase through their SU(2) vectors V / sqrt(det V),
// which are fixed up to sign (su2_vector(), chord_distance()).

// The SU(2) vector the lattice search for determinant class j aims at: that of
// V Q^-j.
Vec4 search_direction(const Unitary& V, int det_class, const GateSet& gateset);

// V / sqrt(det V) (principal root) as (Re a, Im a, Re b, Im b), (a, b) its first
// column: V up to phase, fixed up to sign.
Vec4 su2_vector(const Unitary& V);

// Diamond distance of the PU(2) classes, from the shorter chord; accurate at
// all distances.
double chord_distance(const Vec4& a, const Vec4& b);

// An element as a complex matrix, and its SU(2) vector. The matrix also from
// the exponent and the 32 coefficients as doubles (row-major, each over 1, zeta,
// ..., zeta^7), for an element whose coefficients pass 64 bits.
Unitary matrix_unitary(const CyclotomicMatrix& e);
Unitary matrix_unitary(int k, const std::array<double, 32>& coefficients);
Vec4 matrix_direction(const CyclotomicMatrix& e);

// Conjugate transpose and product.
Unitary dagger(const Unitary& M);
Unitary matmul(const Unitary& A, const Unitary& B);

// Normalised Haar measure of an epsilon-ball in PU(2): two caps of angular
// radius asin(epsilon) on S^3.
double ball_measure(double epsilon);

}  // namespace quicksynth
