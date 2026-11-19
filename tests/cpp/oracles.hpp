// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Independent long-double reference code for the C++ tests. The seeded sampler
// is not numpy's; numpy draws come from calibration_targets.json.
#pragma once

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace quicksynth::test {

using Complex = std::complex<long double>;
using Matrix = std::array<Complex, 4>;  // row-major
using Vector4 = std::array<long double, 4>;

Matrix matrix_of(const Unitary& t);
Unitary unitary_of(const Matrix& m);
Matrix mul(const Matrix& a, const Matrix& b);
Matrix dagger_of(const Matrix& a);
Matrix scaled(const Matrix& a, Complex s);
Matrix transpose_of(const Matrix& a);
Complex trace_of(const Matrix& a);

// M / sqrt2^k.
Matrix matrix_value(const CyclotomicMatrix& e);
// (Re a, Im a, Re b, Im b) of V / sqrt(det V); defined up to sign.
Vector4 su2_of(const Matrix& V);
// [[a, -conj b], [b, conj a]] from (Re a, Im a, Re b, Im b).
Matrix su2_matrix(const Vector4& v);

// PU(2) diamond distance sqrt(s - s^2 / 4), s = min(|a - b|^2, |a + b|^2).
double chord_diamond_distance(const Matrix& u, const Matrix& v);
// sqrt(1 - |tr(u v^dagger)|^2 / 4); inaccurate below about 1e-8.
double trace_diamond_distance(const Matrix& u, const Matrix& v);
// max |U^dagger U - I|.
double unitarity_error(const Matrix& U);

// |a - b| <= rel * max(|a|, |b|) or <= abs.
bool close(double a, double b, double rel = 1e-10, double abs = 1e-13);

// min(max |a - b|, max |a + b|) for 4-vectors defined up to sign.
template <class A, class B>
long double gap_up_to_sign(const A& a, const B& b) {
    long double minus = 0, plus = 0;
    for (int i = 0; i < 4; ++i) {
        minus = std::max(minus, std::abs(static_cast<long double>(a[i]) - b[i]));
        plus = std::max(plus, std::abs(static_cast<long double>(a[i]) + b[i]));
    }
    return std::min(minus, plus);
}

// exp(-i theta/2 n.sigma), n normalised; axis is 'x', 'y' or 'z'.
Matrix rotation(double theta, double nx, double ny, double nz);
Matrix rotation(double theta, char axis);
extern const Matrix kH, kS, kX, kY, kZ, kI;
// The word's matrix, its letters multiplied in long double.
Matrix letters_multiplied(const std::string& word);

// (n_T, n_sqrtT); inverses count as the gate.
std::pair<int64_t, int64_t> magic_counts(const std::string& word);

// Canonical Bloch image; equal iff same coset U * Clifford modulo phase.
using CosetKey = std::pair<int64_t, std::vector<std::vector<std::string>>>;
CosetKey coset_key(const CyclotomicMatrix& e);
CosetKey coset_key(const CyclotomicMatrixDigits& e);

// zeta^j * e, using zeta^8 = -1.
CyclotomicMatrix times_zeta(const CyclotomicMatrix& e, int j);

class Rng {
  public:
    explicit Rng(uint64_t seed) : engine_(seed) {}
    double uniform();  // [0, 1), 53 bits
    double normal();   // Box-Muller
    int64_t below(int64_t n);
    std::mt19937_64& engine() { return engine_; }

  private:
    std::mt19937_64 engine_;
};

// Haar U(2) element; uniform unit 4-vector.
Matrix haar_unitary(Rng& rng);
Vector4 unit_vector(Rng& rng);
// Random word over `letters`, length < max_length.
std::string random_word(Rng& rng, const std::string& letters, int max_length);

}  // namespace quicksynth::test
