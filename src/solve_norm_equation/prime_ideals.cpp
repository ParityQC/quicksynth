// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/solve_norm_equation/prime_ideals.hpp"

#include "quicksynth/enumerate_lattice/lll.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace quicksynth {

namespace {

// Bound on the LLL transform's entries; with basis entries below 2^63 the
// reduced vectors are exact in __int128.
constexpr int64_t kTransformLimit = int64_t{1} << 40;
// Generators have entries near p^(f/N) < 2^32; a reduced vector past this
// bound is not short, and below it the Gram matrix is exact in __int128.
constexpr __int128 kShortEntry = __int128{1} << 40;
// Division steps before Euclid's algorithm gives up; each one lowers the norm.
constexpr int kEuclidSteps = 4096;

template <size_t N>
using Ring128 = RingElement<N, __int128>;
template <size_t N>
using RingMp = RingElement<N, int_mp>;

int_mp power(uint64_t p, int f) {
    int_mp out(1);
    for (int i = 0; i < f; ++i) out = out * int_mp(static_cast<__int128>(p));
    return out;
}

void check_prime_input(uint64_t p, const FpPoly& h, size_t n) {
    if ((p & 1) == 0 || p >= (uint64_t{1} << 63))
        throw std::invalid_argument("prime_ideals: p is an odd prime below 2^63");
    if (h.degree < 1 || static_cast<size_t>(h.degree) >= n || h.c[h.degree] != 1)
        throw std::invalid_argument("prime_ideals: h is a monic proper factor of x^N + 1");
}

// The Hermite basis of P = (p, h(zeta)): p zeta^i for i < f and
// zeta^i - (zeta^i mod h) for f <= i < N, residues centred mod p. A short
// vector of P has small coordinates in it, so the LLL transform stays small
// [Cohen, GTM 138, sections 2.4.3 and 4.8.2].
template <size_t N>
std::vector<int_mp> ideal_basis(uint64_t p, const FpPoly& h) {
    const int n = static_cast<int>(N), f = h.degree;
    std::vector<int_mp> rows(N * N);
    for (int i = 0; i < f; ++i)
        rows[static_cast<size_t>(i * n + i)] = int_mp(static_cast<__int128>(p));
    for (int i = f; i < n; ++i) {
        FpPoly monomial;
        monomial.c[static_cast<size_t>(i)] = 1;
        monomial.degree = i;
        const FpPoly residue = poly_divmod(monomial, h, p).second;
        rows[static_cast<size_t>(i * n + i)] = int_mp(1);
        for (int k = 0; k <= residue.degree; ++k) {
            const __int128 c = static_cast<__int128>(residue.c[static_cast<size_t>(k)]);
            rows[static_cast<size_t>(i * n + k)] = int_mp(2 * c > p ? p - c : -c);
        }
    }
    return rows;
}

// Norm of x estimated from its embeddings; 0 when x is 0.
template <size_t N>
double approximate_norm(const Ring128<N>& x) {
    double out = 1;
    for (int k = 1; k < static_cast<int>(N); k += 2) out *= std::norm(embed(x, k));
    return out;
}

// Fincke-Pohst enumeration [Cohen, GTM 138, Algorithm 2.7.5] of the integer
// combinations of a reduced basis with squared length at most `bound`.
template <size_t N>
class ShortVectors {
  public:
    explicit ShortVectors(const std::array<Ring128<N>, N>& basis) : basis_(basis) {
        // Gram-Schmidt from the Gram matrix: r_ij = G_ij - sum_k<j mu_jk r_ik.
        std::array<std::array<double, N>, N> r{};
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j <= i; ++j) {
                __int128 g = 0;
                for (size_t k = 0; k < N; ++k) g += basis_[i][k] * basis_[j][k];
                double s = static_cast<double>(g);
                for (size_t k = 0; k < j; ++k) s -= mu_[j][k] * r[i][k];
                r[i][j] = s;
                if (j < i) mu_[i][j] = s / r[j][j];
            }
        for (size_t i = 0; i < N; ++i) b_[i] = r[i][i];
    }

    // Calls visit(v) for each nonzero lattice vector in the ball; false if
    // the node budget ran out first.
    template <class Visit>
    bool run(double bound, int64_t nodes, Visit&& visit) {
        nodes_ = nodes;
        x_.fill(0);
        return search(static_cast<int>(N) - 1, bound, visit);
    }

  private:
    template <class Visit>
    bool search(int i, double remaining, Visit& visit) {
        double centre = 0;
        for (size_t j = static_cast<size_t>(i) + 1; j < N; ++j) centre -= mu_[j][i] * x_[j];
        const double half = std::sqrt(std::max(remaining, 0.0) / b_[i]);
        const int64_t lo = static_cast<int64_t>(std::ceil(centre - half));
        const int64_t hi = static_cast<int64_t>(std::floor(centre + half));
        for (int64_t xi = lo; xi <= hi; ++xi) {
            if (--nodes_ < 0) return false;
            x_[i] = xi;
            const double d = static_cast<double>(xi) - centre;
            const double left = remaining - b_[i] * d * d;
            if (i > 0) {
                if (!search(i - 1, left, visit)) return false;
            } else {
                Ring128<N> v{};
                bool nonzero = false;
                for (size_t j = 0; j < N; ++j) {
                    if (x_[j] == 0) continue;
                    nonzero = true;
                    for (size_t k = 0; k < N; ++k) v[k] += x_[j] * basis_[j][k];
                }
                if (nonzero) visit(v);
            }
        }
        x_[i] = 0;
        return true;
    }

    std::array<Ring128<N>, N> basis_;
    std::array<std::array<double, N>, N> mu_{};
    std::array<double, N> b_{};
    std::array<int64_t, N> x_{};
    int64_t nodes_ = 0;
};

template <size_t N>
__int128 length2(const Ring128<N>& v) {
    __int128 s = 0;
    for (const __int128& c : v) s += c * c;
    return s;
}

template <size_t N>
bool is_zero(const RingMp<N>& x) {
    for (const int_mp& c : x)
        if (!c.zero()) return false;
    return true;
}

// round(a / n) for n > 0, halves up.
int_mp nearest_quotient(const int_mp& a, const int_mp& n) { return floor_div(a * 2 + n, n * 2); }

// a - q b with N(a - q b) < N(b) for q near a / b, or nullopt.
template <size_t N>
std::optional<RingMp<N>> euclid_remainder(const RingMp<N>& a, const RingMp<N>& b,
                                          const int_mp& norm_b) {
    // a / b = a c / N(b), c the product of b's other conjugates
    RingMp<N> c = ring_zeta_power<N, int_mp>(0);
    for (int j = 3; j < 2 * static_cast<int>(N); j += 2) c = ring_mul(c, galois(b, j));
    const RingMp<N> numerator = ring_mul(a, c);
    RingMp<N> q;
    for (size_t i = 0; i < N; ++i) q[i] = nearest_quotient(numerator[i], norm_b);
    RingMp<N> best = ring_sub(a, ring_mul(q, b));
    int_mp best_norm = ring_norm(best);
    // The rounded quotient can miss the Euclidean minimum; its unit neighbours
    // usually reach it.
    for (int j = 0; j < 2 * static_cast<int>(N) && !(best_norm < norm_b); ++j) {
        const RingMp<N> r = ring_sub(a, ring_mul(ring_add(q, ring_zeta_power<N, int_mp>(j)), b));
        const int_mp norm_r = ring_norm(r);
        if (norm_r < best_norm) {
            best = r;
            best_norm = norm_r;
        }
    }
    if (!(best_norm < norm_b)) return std::nullopt;
    return best;
}

}  // namespace

template <size_t N>
std::optional<Ring128<N>> short_generator(uint64_t p, const FpPoly& h) {
    check_prime_input(p, h, N);
    const int n = static_cast<int>(N), f = h.degree;
    const std::vector<int_mp> rows = ideal_basis<N>(p, h);
    LllTransform t;
    if (!lll_transform(n, rows, kTransformLimit, t)) return std::nullopt;
    // Reduced vector j = sum_i T_ij row_i, exact: |row| < 2^63, |T| < 2^40.
    std::array<Ring128<N>, N> basis{};
    for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k) {
            __int128 s = 0;
            for (int i = 0; i < n; ++i) {
                const int_mp& r = rows[static_cast<size_t>(i * n + k)];
                if (!r.zero()) s += r.to_int128() * t.T[static_cast<size_t>(i * n + j)];
            }
            if (s >= kShortEntry || s <= -kShortEntry) return std::nullopt;
            basis[j][k] = s;
        }
    const int_mp target = power(p, f);
    const double target_d = target.to_double(0);
    const double least = std::pow(static_cast<double>(p), 2.0 * f / n);
    ShortVectors<N> scan(basis);
    for (const double radius : kGeneratorRadii) {
        std::optional<Ring128<N>> best;
        __int128 best_length = 0;
        const auto visit = [&](const Ring128<N>& v) {
            const double estimate = approximate_norm<N>(v);
            // N(v) is a positive multiple of p^f; the estimate only screens.
            if (!(estimate > 0.5 * target_d && estimate < 1.5 * target_d)) return;
            const __int128 length = length2<N>(v);
            if (best && (length > best_length || (length == best_length && !(v < *best)))) return;
            if (ring_norm(v) != target) return;
            best = v;
            best_length = length;
        };
        const bool complete = scan.run(radius * least * (1 + 1e-9), kGeneratorScanNodes, visit);
        if (best) return best;
        if (!complete) return std::nullopt;
    }
    return std::nullopt;
}

template <size_t N>
std::optional<Ring128<N>> euclid_generator(uint64_t p, const FpPoly& h) {
    check_prime_input(p, h, N);
    RingMp<N> a{}, b{};
    a[0] = int_mp(static_cast<__int128>(p));
    for (int k = 0; k <= h.degree; ++k)
        b[static_cast<size_t>(k)] = int_mp(static_cast<__int128>(h.c[k]));
    int_mp norm_b = ring_norm(b);
    for (int step = 0; step < kEuclidSteps && !is_zero<N>(b); ++step) {
        std::optional<RingMp<N>> r = euclid_remainder<N>(a, b, norm_b);
        if (!r) return std::nullopt;
        a = std::move(b);
        b = std::move(*r);
        norm_b = ring_norm(b);
    }
    if (!is_zero<N>(b) || ring_norm(a) != power(p, h.degree)) return std::nullopt;
    Ring128<N> out;
    for (size_t i = 0; i < N; ++i) {
        if (a[i].bits() > 126) return std::nullopt;
        out[i] = a[i].to_int128();
    }
    return out;
}

template <size_t N>
std::optional<Ring128<N>> prime_generator(uint64_t p, const FpPoly& h) {
    if (std::optional<Ring128<N>> g = short_generator<N>(p, h)) return g;
    return euclid_generator<N>(p, h);
}

template <size_t N>
PrimeDecomposition<N> decompose(uint64_t p, const Ring128<N>& pi) {
    const int m = 2 * static_cast<int>(N);
    const int residue = static_cast<int>(p % static_cast<uint64_t>(m));
    PrimeDecomposition<N> out;
    out.f = cyclotomic_residue_degree(p, static_cast<int>(N));
    bool inert = false;
    for (int q = residue, k = 0; k < out.f; q = q * residue % m, ++k) inert |= q == m - 1;
    std::array<bool, 2 * N> covered{};
    const auto cover = [&](int j) {
        for (int q = j, k = 0; k < out.f; q = q * residue % m, ++k)
            covered[static_cast<size_t>(q)] = true;
    };
    for (int j = 1; j < m; j += 2) {
        if (covered[static_cast<size_t>(j)]) continue;
        cover(j);
        if (!inert) cover(m - j);  // sigma_-j P = (sigma_j P)*, the pair's other prime
        PrimeAbove<N> prime;
        prime.pi = galois(pi, j);
        prime.inert = inert;
        prime.norm = relative_norm(prime.pi);
        out.primes.push_back(prime);
    }
    return out;
}

template <size_t N>
std::optional<PrimeDecomposition<N>> PrimeIdeals<N>::above(uint64_t p) {
    if ((p & 1) == 0) throw std::invalid_argument("PrimeIdeals::above: p is odd");
    if (PrimeDecomposition<N>* hit = cache_.find(p)) return *hit;
    const std::optional<Ring128<N>> pi =
        prime_generator<N>(p, cyclotomic_factors(static_cast<int>(N), p).front());
    if (!pi) return std::nullopt;
    return cache_.put(p, decompose<N>(p, *pi));
}

#define QUICKSYNTH_PRIME_IDEALS(N)                                                   \
    template std::optional<Ring128<N>> short_generator<N>(uint64_t, const FpPoly&);  \
    template std::optional<Ring128<N>> euclid_generator<N>(uint64_t, const FpPoly&); \
    template std::optional<Ring128<N>> prime_generator<N>(uint64_t, const FpPoly&);  \
    template PrimeDecomposition<N> decompose<N>(uint64_t, const Ring128<N>&);        \
    template class PrimeIdeals<N>;

QUICKSYNTH_PRIME_IDEALS(4)
QUICKSYNTH_PRIME_IDEALS(8)

#undef QUICKSYNTH_PRIME_IDEALS

}  // namespace quicksynth
