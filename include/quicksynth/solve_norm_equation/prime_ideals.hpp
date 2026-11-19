// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Generators of the primes of O_K = Z[zeta_2N] (N = 4, 8) above an odd rational
// prime p. Both rings have class number 1 [Washington, GTM 83, Theorem 11.1 and
// the table of section 11.1], so every prime P = (p, h(zeta)), h an irreducible
// factor of x^N + 1 mod p, is principal. A generator is a lattice vector of P
// of norm exactly p^f; the power basis is orthogonal for the trace form, so the
// shortest such vectors are short in coefficient space, where an exact LLL and
// a Fincke-Pohst scan find them [Cohen, GTM 138, Algorithms 2.6.7 and 2.7.5].
// Euclid's algorithm on p and h(zeta) is the fallback; both rings are
// norm-Euclidean [Lenstra, J. London Math. Soc. (2) 10 (1975)].
//
// The other primes above p are Galois images sigma_j(P), one per coset of the
// decomposition group <p> in (Z/2N)^* [Washington, GTM 83, Theorem 2.13].
// Complex conjugation is sigma_-1, so P is its own conjugate, and the prime of
// K+ below it inert in K/K+, exactly when -1 is a power of p mod 2N.
#pragma once

#include "quicksynth/data_structures/number_field.hpp"
#include "quicksynth/solve_norm_equation/prime_field.hpp"
#include "quicksynth/support/lru.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace quicksynth {

#pragma GCC visibility push(hidden)

// The Fincke-Pohst scan's radii, as multiples of p^(2f/N), the least squared
// coefficient norm of an element of norm p^f (AM-GM on the trace form).
inline constexpr double kGeneratorRadii[] = {2, 4, 8, 16};
// Lattice nodes the scan visits per radius before it gives up on that radius.
inline constexpr int64_t kGeneratorScanNodes = int64_t{1} << 16;

// Generators of P = (p, h(zeta)) for an odd prime p < 2^63 and a monic
// irreducible factor h of x^N + 1 mod p, or nullopt if the method finds none;
// std::invalid_argument for even or larger p or an h of degree 0 or N. Each
// result lies in P and has norm exactly p^f, f = deg h.
//
// short_generator: the shortest norm-p^f vector of the first scan radius that
// has one, ties broken by the coefficients. euclid_generator: gcd(p, h(zeta))
// by nearest-integer division, which may leave a large unit factor.
// prime_generator: short_generator, else euclid_generator.
template <size_t N>
std::optional<RingElement<N, __int128>> short_generator(uint64_t p, const FpPoly& h);
template <size_t N>
std::optional<RingElement<N, __int128>> euclid_generator(uint64_t p, const FpPoly& h);
template <size_t N>
std::optional<RingElement<N, __int128>> prime_generator(uint64_t p, const FpPoly& h);

// A prime of K+ above p, given by a prime P of K above it.
template <size_t N>
struct PrimeAbove {
    RingElement<N, __int128> pi;  // generates P
    // P = P*: pi* / pi is a root of unity. Otherwise pi* generates P* != P.
    bool inert = false;
    // pi pi*: generates the K+ prime below P if it splits, its square if inert.
    RealElement<N / 2, __int128> norm;
};

template <size_t N>
struct PrimeDecomposition {
    int f = 0;  // residue degree of p in K; N(pi) = p^f
    // One per K+ prime above p: N / (2 f) split pairs, or N / f inert primes.
    std::vector<PrimeAbove<N>> primes;
};

// The decomposition from a generator of one prime above p.
template <size_t N>
PrimeDecomposition<N> decompose(uint64_t p, const RingElement<N, __int128>& pi);

// Prime decompositions by p, cached across calls (the small primes recur from
// one norm equation to the next). Not thread-safe.
inline constexpr int64_t kPrimeCacheCapacity = 256;

template <size_t N>
class PrimeIdeals {
  public:
    explicit PrimeIdeals(int64_t capacity = kPrimeCacheCapacity) : cache_(capacity) {}
    // The primes above an odd prime p; std::invalid_argument for even p. Nullopt
    // if no generator was found, which is not cached.
    std::optional<PrimeDecomposition<N>> above(uint64_t p);

  private:
    Lru<uint64_t, PrimeDecomposition<N>> cache_;
};

#pragma GCC visibility pop

}  // namespace quicksynth
