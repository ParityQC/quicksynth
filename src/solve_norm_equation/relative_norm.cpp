// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/solve_norm_equation/relative_norm.hpp"

#include "quicksynth/solve_norm_equation/integer_factor.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <optional>
#include <stdexcept>
#include <vector>

namespace quicksynth {

namespace {

template <size_t N>
using Ring = RingElement<N, __int128>;
template <size_t N>
using Real = RealElement<N / 2, __int128>;

[[noreturn]] void broken(const char* what) {
    throw std::logic_error(std::string("solve_relative_norm: ") + what);
}

template <size_t N>
NormSolutions<N> decided(NormOutcome outcome) {
    NormSolutions<N> out;
    out.outcome = outcome;
    return out;
}

template <size_t N>
NormSolutions<N> undecided(NormUndecided why) {
    NormSolutions<N> out;
    out.outcome = NormOutcome::undecided;
    out.undecided = why;
    return out;
}

template <size_t N>
Ring<N> one() {
    return ring_zeta_power<N, __int128>(0);
}

// x^e by repeated multiplication; e is a small valuation.
template <size_t N>
Ring<N> power(const Ring<N>& x, int e) {
    Ring<N> out = one<N>();
    for (int i = 0; i < e; ++i) out = ring_mul(out, x);
    return out;
}

// x divided by d as often as d divides it; returns the count.
template <size_t N>
int strip(Real<N>& x, const Real<N>& d, int most) {
    int v = 0;
    while (v < most) {
        std::optional<Real<N>> q = real_divide(x, d);
        if (!q) break;
        x = *q;
        ++v;
    }
    return v;
}

// A split pair (P, P*) at valuation v.
template <size_t N>
struct SplitPair {
    std::vector<Ring<N>> pi_powers, conjugate_powers;  // index 0 .. v
    int v = 0;
};

int64_t saturating_product(int64_t a, int64_t b) {
    int64_t out;
    if (__builtin_mul_overflow(a, b, &out)) return std::numeric_limits<int64_t>::max();
    return out;
}

template <size_t N>
NormSolutions<N> solve(const Real<N>& xi, PrimeIdeals<N>& primes, int cap, bool halving) {
    constexpr int m = 2 * static_cast<int>(N);
    if (xi == Real<N>{}) {
        NormSolutions<N> out = decided<N>(NormOutcome::solved);
        out.count = 1;
        out.choices.push_back(Ring<N>{});
        return out;
    }
    if (!totally_positive(xi)) return decided<N>(NormOutcome::none);

    // The prime above 2 is ramified, delta delta* of norm 2, so its
    // valuation is the power of 2 in the norm.
    int_mp n = real_norm(xi);
    int v2 = 0;
    while (n.even()) {
        n = half_of(n);
        ++v2;
    }
    if (n.bits() > 63) return undecided<N>(NormUndecided::norm_too_large);
    const std::optional<std::vector<PrimePower>> factors =
        factor_u64(static_cast<uint64_t>(n.to_int64()));
    if (!factors) return undecided<N>(NormUndecided::factor_gave_up);
    for (const PrimePower& q : *factors)
        if (q.e % 2 == 1 && q.p % m == static_cast<uint64_t>(m - 1))
            return decided<N>(NormOutcome::none);

    const Ring<N> delta = ring_sub(one<N>(), ring_zeta_power<N, __int128>(1));
    Real<N> rest = xi;
    if (strip<N>(rest, relative_norm(delta), v2) != v2) broken("the 2-part does not divide xi");
    Ring<N> fixed = power<N>(delta, v2);
    std::vector<SplitPair<N>> pairs;
    for (const PrimePower& q : *factors) {
        const std::optional<PrimeDecomposition<N>> d = primes.above(q.p);
        if (!d) return undecided<N>(NormUndecided::no_generator);
        // Each division by pi pi* removes p^f from the norm: one K+ prime if
        // split, the square of one if inert.
        int left = q.e;
        for (const PrimeAbove<N>& prime : d->primes) {
            if (left == 0) break;
            const int w = strip<N>(rest, prime.norm, left / d->f);
            left -= w * d->f;
            if (w == 0) continue;
            if (prime.inert) {
                fixed = ring_mul(fixed, power<N>(prime.pi, w));
            } else {
                SplitPair<N> pair;
                pair.v = w;
                pair.pi_powers.push_back(one<N>());
                pair.conjugate_powers.push_back(one<N>());
                const Ring<N> conj = conjugate(prime.pi);
                for (int a = 1; a <= w; ++a) {
                    pair.pi_powers.push_back(ring_mul(pair.pi_powers.back(), prime.pi));
                    pair.conjugate_powers.push_back(ring_mul(pair.conjugate_powers.back(), conj));
                }
                pairs.push_back(std::move(pair));
            }
        }
        // An inert K+ prime at odd valuation leaves p^(f/2) behind.
        if (left != 0) {
            if (d->primes.front().inert) return decided<N>(NormOutcome::none);
            broken("the split primes above p do not account for its power in the norm");
        }
    }
    if (real_norm(rest) != int_mp(1)) broken("xi is not a unit after its primes are removed");

    NormSolutions<N> out = decided<N>(NormOutcome::solved);
    out.count = 1;
    for (const SplitPair<N>& pair : pairs) out.count = saturating_product(out.count, pair.v + 1);
    if (halving) out.count = out.count / 2 + out.count % 2;
    out.capped = out.count > cap;
    const int64_t wanted = out.capped ? cap : out.count;
    for (int64_t c = 0; c < wanted; ++c) {
        Ring<N> t0 = fixed;
        int64_t digits = c;
        for (const SplitPair<N>& pair : pairs) {
            const int a = static_cast<int>(digits % (pair.v + 1));
            digits /= pair.v + 1;
            t0 = ring_mul(t0, ring_mul(pair.pi_powers[a], pair.conjugate_powers[pair.v - a]));
        }
        const std::optional<Real<N>> unit = real_divide(xi, relative_norm(t0));
        if (!unit) broken("t0 t0* does not divide xi");
        const std::optional<Real<N>> eta = unit_sqrt(*unit);
        if (!eta) return undecided<N>(NormUndecided::unit_root);
        Ring<N> t = ring_mul(from_real(*eta), t0);
        if (relative_norm(t) != xi) broken("a solution fails t t* = xi");
        out.choices.push_back(t);
    }
    return out;
}

}  // namespace

template <size_t N>
NormSolutions<N> solve_relative_norm(const RealElement<N / 2, __int128>& xi, PrimeIdeals<N>& primes,
                                     int cap, bool transpose_halving) {
    try {
        return solve<N>(xi, primes, cap, transpose_halving);
    } catch (const std::overflow_error&) {
        // __int128 digits throw past 127 bits; the solutions are not exact.
        return undecided<N>(NormUndecided::digit_overflow);
    }
}

template NormSolutions<4> solve_relative_norm<4>(const RealElement<2, __int128>&, PrimeIdeals<4>&,
                                                 int, bool);
template NormSolutions<8> solve_relative_norm<8>(const RealElement<4, __int128>&, PrimeIdeals<8>&,
                                                 int, bool);

}  // namespace quicksynth
