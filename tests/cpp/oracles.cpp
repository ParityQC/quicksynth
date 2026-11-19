// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "oracles.hpp"

#include "quicksynth/data_structures/bloch.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quicksynth::test {

namespace {

const long double kPi = std::acos(-1.0L);

Vector4 canonical_sign(Vector4 v) {
    for (long double x : v) {
        if (x == 0) continue;
        if (x < 0)
            for (long double& y : v) y = -y;
        break;
    }
    return v;
}

CosetKey key_of(const BlochDigits& b) {
    std::vector<std::vector<std::string>> columns;
    for (int c = 0; c < 3; ++c) {
        std::vector<std::string> column;
        for (int r = 0; r < 3; ++r)
            for (int p = 0; p < 4; ++p) column.push_back(at(b, r, c, p));
        const auto pivot = std::find_if(column.begin(), column.end(),
                                        [](const std::string& d) { return d != "0"; });
        if (pivot != column.end() && (*pivot)[0] == '-')
            for (std::string& d : column)
                if (d != "0") d = d[0] == '-' ? d.substr(1) : "-" + d;
        columns.push_back(column);
    }
    std::sort(columns.begin(), columns.end());
    return {b.exponent, columns};
}

}  // namespace

const Matrix kI{Complex(1), Complex(0), Complex(0), Complex(1)};
const Matrix kH{Complex(1 / std::sqrt(2.0L)), Complex(1 / std::sqrt(2.0L)),
                Complex(1 / std::sqrt(2.0L)), Complex(-1 / std::sqrt(2.0L))};
const Matrix kS{Complex(1), Complex(0), Complex(0), Complex(0, 1)};
const Matrix kX{Complex(0), Complex(1), Complex(1), Complex(0)};
const Matrix kY{Complex(0), Complex(0, -1), Complex(0, 1), Complex(0)};
const Matrix kZ{Complex(1), Complex(0), Complex(0), Complex(-1)};

Matrix matrix_of(const Unitary& t) {
    Matrix m;
    for (int i = 0; i < 4; ++i) m[i] = Complex(t[2 * i], t[2 * i + 1]);
    return m;
}

Unitary unitary_of(const Matrix& m) {
    Unitary t;
    for (int i = 0; i < 4; ++i) {
        t[2 * i] = static_cast<double>(m[i].real());
        t[2 * i + 1] = static_cast<double>(m[i].imag());
    }
    return t;
}

Matrix mul(const Matrix& a, const Matrix& b) {
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3], a[2] * b[0] + a[3] * b[2],
            a[2] * b[1] + a[3] * b[3]};
}

Matrix dagger_of(const Matrix& a) {
    return {std::conj(a[0]), std::conj(a[2]), std::conj(a[1]), std::conj(a[3])};
}

Matrix scaled(const Matrix& a, Complex s) { return {a[0] * s, a[1] * s, a[2] * s, a[3] * s}; }

Matrix transpose_of(const Matrix& a) { return {a[0], a[2], a[1], a[3]}; }

Complex trace_of(const Matrix& a) { return a[0] + a[3]; }

Matrix matrix_value(const CyclotomicMatrix& e) {
    Matrix m;
    const long double scale = std::pow(std::sqrt(2.0L), -e.denominator_exponent);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            Complex z = 0;
            for (int t = 0; t < kZetaDegree; ++t)
                z += static_cast<long double>(e.matrix[i][j][t]) *
                     Complex(std::cos(kPi * t / 8), std::sin(kPi * t / 8));
            m[2 * i + j] = z * scale;
        }
    return m;
}

Vector4 su2_of(const Matrix& V) {
    const Complex root = std::sqrt(V[0] * V[3] - V[1] * V[2]);
    const Complex a = V[0] / root, b = V[2] / root;
    return {a.real(), a.imag(), b.real(), b.imag()};
}

Matrix su2_matrix(const Vector4& v) {
    const Complex a(v[0], v[1]), b(v[2], v[3]);
    return {a, -std::conj(b), b, std::conj(a)};
}

double chord_diamond_distance(const Matrix& u, const Matrix& v) {
    const Vector4 a = su2_of(u), b = su2_of(v);
    long double minus = 0, plus = 0;
    for (int i = 0; i < 4; ++i) {
        minus += (a[i] - b[i]) * (a[i] - b[i]);
        plus += (a[i] + b[i]) * (a[i] + b[i]);
    }
    const long double s = std::min(minus, plus);
    return static_cast<double>(std::sqrt(std::max(s - s * s / 4, 0.0L)));
}

double trace_diamond_distance(const Matrix& u, const Matrix& v) {
    const long double overlap = std::min(std::abs(trace_of(mul(u, dagger_of(v)))) / 2, 1.0L);
    return static_cast<double>(std::sqrt(std::max(0.0L, 1 - overlap * overlap)));
}

double unitarity_error(const Matrix& U) {
    const Matrix p = mul(dagger_of(U), U);
    long double worst = 0;
    for (int i = 0; i < 4; ++i) worst = std::max(worst, std::abs(p[i] - kI[i]));
    return static_cast<double>(worst);
}

bool close(double a, double b, double rel, double abs) {
    return std::abs(a - b) <= std::max(rel * std::max(std::abs(a), std::abs(b)), abs);
}

Matrix rotation(double theta, double nx, double ny, double nz) {
    const long double norm =
        std::sqrt(static_cast<long double>(nx) * nx + static_cast<long double>(ny) * ny +
                  static_cast<long double>(nz) * nz);
    const long double c = std::cos(theta / 2.0L), s = std::sin(theta / 2.0L);
    const Complex x(nx / norm), y(ny / norm), z(nz / norm), i(0, 1);
    Matrix out;
    for (int e = 0; e < 4; ++e) out[e] = c * kI[e] - i * s * (x * kX[e] + y * kY[e] + z * kZ[e]);
    return out;
}

Matrix rotation(double theta, char axis) {
    return rotation(theta, axis == 'x' ? 1 : 0, axis == 'y' ? 1 : 0, axis == 'z' ? 1 : 0);
}

Matrix letters_multiplied(const std::string& word) {
    const auto phase = [](int j) {  // diag(1, zeta^j), zeta = exp(i pi / 8)
        return Matrix{Complex(1), Complex(0), Complex(0),
                      std::polar(1.0L, std::acos(-1.0L) * j / 8)};
    };
    Matrix m = kI;
    for (char c : word) {
        switch (c) {
            case 'H': m = mul(m, kH); break;
            case 'X': m = mul(m, kX); break;
            case 'Y': m = mul(m, kY); break;
            case 'Z': m = mul(m, kZ); break;
            case 'S': m = mul(m, phase(4)); break;
            case 's': m = mul(m, phase(-4)); break;
            case 'T': m = mul(m, phase(2)); break;
            case 't': m = mul(m, phase(-2)); break;
            case 'Q': m = mul(m, phase(1)); break;
            case 'q': m = mul(m, phase(-1)); break;
            default: throw std::invalid_argument(std::string("no letter ") + c);
        }
    }
    return m;
}

std::pair<int64_t, int64_t> magic_counts(const std::string& word) {
    const auto n = [&word](const char* letters) {
        return static_cast<int64_t>(std::count_if(word.begin(), word.end(), [letters](char c) {
            return std::string(letters).find(c) != std::string::npos;
        }));
    };
    return {n("Tt"), n("Qq")};
}

CosetKey coset_key(const CyclotomicMatrix& e) {
    return key_of(bloch_digits(e.matrix, e.denominator_exponent));
}

CosetKey coset_key(const CyclotomicMatrixDigits& e) { return key_of(bloch_digits(e)); }

CyclotomicMatrix times_zeta(const CyclotomicMatrix& e, int j) {
    CyclotomicMatrix out{M2{}, e.denominator_exponent};
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) out.matrix[a][b] = quicksynth::times_zeta(e.matrix[a][b], j);
    return out;
}

double Rng::uniform() { return static_cast<double>(engine_() >> 11) * 0x1.0p-53; }

double Rng::normal() {
    const double u = 1.0 - uniform(), v = uniform();
    return std::sqrt(-2.0 * std::log(u)) * std::cos(2.0 * static_cast<double>(kPi) * v);
}

int64_t Rng::below(int64_t n) { return static_cast<int64_t>(engine_() % static_cast<uint64_t>(n)); }

Matrix haar_unitary(Rng& rng) {
    std::array<Complex, 4> z;
    for (Complex& x : z) x = Complex(rng.normal(), rng.normal());
    // Gram-Schmidt on columns; positive real diag(R) makes Q Haar-distributed.
    Complex c0[2] = {z[0], z[2]}, c1[2] = {z[1], z[3]};
    const long double n0 = std::sqrt(std::norm(c0[0]) + std::norm(c0[1]));
    for (Complex& x : c0) x /= n0;
    const Complex proj = std::conj(c0[0]) * c1[0] + std::conj(c0[1]) * c1[1];
    for (int i = 0; i < 2; ++i) c1[i] -= proj * c0[i];
    const long double n1 = std::sqrt(std::norm(c1[0]) + std::norm(c1[1]));
    for (Complex& x : c1) x /= n1;
    return {c0[0], c1[0], c0[1], c1[1]};
}

Vector4 unit_vector(Rng& rng) {
    Vector4 v;
    long double norm = 0;
    for (long double& x : v) {
        x = rng.normal();
        norm += x * x;
    }
    for (long double& x : v) x /= std::sqrt(norm);
    return canonical_sign(v);
}

std::string random_word(Rng& rng, const std::string& letters, int max_length) {
    std::string word;
    const int64_t n = rng.below(max_length);
    for (int64_t i = 0; i < n; ++i)
        word += letters[static_cast<size_t>(rng.below(static_cast<int64_t>(letters.size())))];
    return word;
}

}  // namespace quicksynth::test
