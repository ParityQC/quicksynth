// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/enumerate_lattice/lll.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace quicksynth {

namespace {

using Int128 = __int128;
using UInt128 = unsigned __int128;

// L2's tolerances, midway between (0.51, 0.99) and (0.5, 1). A test passed
// under rounding still leaves the basis (0.99, 0.51)-reduced.
constexpr double kEta = (0.51 + 0.5) / 2;
constexpr double kDelta = (0.99 + 1.0) / 2;
// A double Gram product below this fraction of the norms has cancelled and is
// recomputed exactly.
constexpr double kCancelled = 0x1p-26;
// Bound on entries and size-reduction steps of the 128-bit basis; exact_dot
// needs entries below it.
constexpr double kEntryBound = 0x1p120;
// Size-reduction passes per vector and total steps before the float LLL gives up.
constexpr int kMaxPasses = 32;
constexpr int64_t kMaxSteps = 200000;

// x -= q y; false on overflow.
bool sub_mul(Int128& x, Int128 q, Int128 y) {
    Int128 p;
    return !__builtin_mul_overflow(q, y, &p) && !__builtin_sub_overflow(x, p, &x);
}

// x += q y; false on overflow.
bool add_mul(Int128& x, Int128 q, Int128 y) {
    Int128 p;
    return !__builtin_mul_overflow(q, y, &p) && !__builtin_add_overflow(x, p, &x);
}

bool sub_mul(int_mp& x, const int_mp& q, const int_mp& y) {
    x -= q * y;
    return true;
}

bool add_mul(int_mp& x, const int_mp& q, const int_mp& y) {
    x += q * y;
    return true;
}

// x into out if |x| < bound.
bool narrow(Int128 x, int64_t bound, int64_t& out) {
    if (!(x < bound && x > -bound)) return false;
    out = static_cast<int64_t>(x);
    return true;
}

bool narrow(const int_mp& x, int64_t bound, int64_t& out) {
    if (!x.fits_int64()) return false;
    out = x.to_int64();
    return out < bound && out > -bound;
}

// Basis rows b = U b_in with U and U^-1, n x n row-major. Each step is an
// elementary row operation E on b and U; U^-1 takes E^-1 on the right.
template <class Int>
class TrackedBasis {
  public:
    TrackedBasis(int n, std::vector<Int> rows)
        : n_(n), b_(std::move(rows)), u_(b_.size()), u_inv_(b_.size()) {
        for (int i = 0; i < n; ++i) u(i, i) = u_inv(i, i) = Int(1);
    }

    int n() const { return n_; }
    const Int& b(int i, int j) const { return b_[cell(i, j)]; }
    const Int* row(int i) const { return b_.data() + cell(i, 0); }

    // b_k -= q b_j; U^-1 adds q times its column k to column j. False on
    // overflow, which leaves the basis invalid.
    bool subtract(int k, int j, const Int& q) {
        for (int c = 0; c < n_; ++c)
            if (!sub_mul(b_[cell(k, c)], q, b(j, c)) || !sub_mul(u(k, c), q, u(j, c)) ||
                !add_mul(u_inv(c, j), q, u_inv(c, k)))
                return false;
        return true;
    }

    // b_{k-1} <-> b_k, its own inverse.
    void exchange(int k) {
        for (int c = 0; c < n_; ++c) {
            std::swap(b_[cell(k - 1, c)], b_[cell(k, c)]);
            std::swap(u(k - 1, c), u(k, c));
            std::swap(u_inv(c, k - 1), u_inv(c, k));
        }
    }

    // T = U^t and T_inv = (U^-1)^t, since the reduced rows are U b_in and the
    // reduced columns of B are B U^t. False if an entry of T is not below `limit`
    // or one of T_inv not below 2^63 - 1.
    bool transform(int64_t limit, LllTransform& out) {
        out.T.assign(b_.size(), 0);
        out.T_inv.assign(b_.size(), 0);
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < n_; ++j)
                if (!narrow(u(j, i), limit, out.T[cell(i, j)]) ||
                    !narrow(u_inv(j, i), std::numeric_limits<int64_t>::max(),
                            out.T_inv[cell(i, j)]))
                    return false;
        return true;
    }

    // The same with T_inv of any size.
    bool transform(int64_t limit, LllTransformMp& out) {
        out.T.assign(b_.size(), 0);
        out.T_inv.assign(b_.size(), int_mp());
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < n_; ++j) {
                if (!narrow(u(j, i), limit, out.T[cell(i, j)])) return false;
                out.T_inv[cell(i, j)] = int_mp(u_inv(j, i));
            }
        return true;
    }

  private:
    size_t cell(int i, int j) const { return static_cast<size_t>(i) * n_ + j; }
    Int& u(int i, int j) { return u_[cell(i, j)]; }
    Int& u_inv(int i, int j) { return u_inv_[cell(i, j)]; }

    int n_;
    std::vector<Int> b_, u_, u_inv_;
};

// <a, b> summed exactly in 256 bits, then rounded to a double within two
// roundings. NaN if an entry is not below 2^120 in absolute value.
double exact_dot(const Int128* a, const Int128* b, int n) {
    constexpr Int128 kBound = Int128{1} << 120;
    constexpr UInt128 kLow64 = ~uint64_t{0};
    UInt128 lo = 0;  // the sum is hi 2^128 + lo
    Int128 hi = 0;
    for (int c = 0; c < n; ++c) {
        if (a[c] >= kBound || a[c] <= -kBound || b[c] >= kBound || b[c] <= -kBound)
            return std::nan("");
        const UInt128 x = a[c] < 0 ? -static_cast<UInt128>(a[c]) : static_cast<UInt128>(a[c]);
        const UInt128 y = b[c] < 0 ? -static_cast<UInt128>(b[c]) : static_cast<UInt128>(b[c]);
        // x y = high 2^128 + low from 64-bit halves; each cross product is below 2^120
        const UInt128 cross = (x >> 64) * (y & kLow64) + (x & kLow64) * (y >> 64);
        UInt128 low = (x & kLow64) * (y & kLow64);
        UInt128 high = (x >> 64) * (y >> 64) + (cross >> 64);
        low += cross << 64;
        if (low < cross << 64) ++high;
        if ((a[c] < 0) != (b[c] < 0)) {
            hi -= static_cast<Int128>(high) + (lo < low ? 1 : 0);
            lo -= low;
        } else {
            lo += low;
            hi += static_cast<Int128>(high) + (lo < low ? 1 : 0);
        }
    }
    // the magnitude, so that its two halves add without cancelling
    const bool negative = hi < 0;
    if (negative) {
        lo = ~lo + 1;
        hi = ~hi + (lo == 0 ? 1 : 0);
    }
    const double magnitude = std::ldexp(static_cast<double>(hi), 128) + static_cast<double>(lo);
    return negative ? -magnitude : magnitude;
}

// Schnorr-Euchner floating-point LLL on the 128-bit basis. r_kj = <b_k, b_j*>
// and mu_kj are doubles, recomputed from the Gram products as in L2's Cholesky
// step; products that cancel are taken exactly. Size reduction is lazy, pass by
// pass, until every |mu| <= eta. run() returns false, for the exact fallback,
// on a non-finite value, a pass that does not shrink the largest |mu|,
// overflow or the step limit.
class FloatLll {
  public:
    explicit FloatLll(TrackedBasis<Int128>& basis)
        : basis_(basis),
          n_(basis.n()),
          approx_(static_cast<size_t>(n_) * n_),
          norm_sq_(n_),
          r_(static_cast<size_t>(n_) * n_),
          mu_(static_cast<size_t>(n_) * n_),
          step_(n_) {}

    bool run() {
        for (int i = 0; i < n_; ++i) refresh(i);
        r(0, 0) = norm_sq_[0];
        int k = 1;
        int64_t steps = 0;
        while (k < n_) {
            if (++steps > kMaxSteps) return false;
            if (!size_reduce(k)) return false;
            // The Lovasz test uses s = |b_k|^2 - sum_{j < k-1} mu_kj r_kj, b_k*
            // after the exchange, not r_kk = s - mu r. If b_k* is much shorter
            // than b_k-1*, r_kk has cancelled entirely and the exchange is due.
            double s = norm_sq_[k];
            for (int j = 0; j < k - 1; ++j) s -= mu(k, j) * r(k, j);
            if (!std::isfinite(s)) return false;
            if (kDelta * r(k - 1, k - 1) > s) {
                exchange(k);
                if (--k == 0) {
                    r(0, 0) = norm_sq_[0];
                    k = 1;
                }
                continue;
            }
            r(k, k) = s - mu(k, k - 1) * r(k, k - 1);
            ++k;
        }
        return true;
    }

  private:
    double& r(int i, int j) { return r_[static_cast<size_t>(i) * n_ + j]; }
    double& mu(int i, int j) { return mu_[static_cast<size_t>(i) * n_ + j]; }
    double* approx(int i) { return approx_.data() + static_cast<size_t>(i) * n_; }

    // The double copy of row i and its squared norm.
    void refresh(int i) {
        double* a = approx(i);
        double sq = 0.0;
        for (int c = 0; c < n_; ++c) {
            a[c] = static_cast<double>(basis_.b(i, c));
            sq += a[c] * a[c];
        }
        norm_sq_[i] = sq;
    }

    double gram(int i, int j) {
        const double* a = approx(i);
        const double* b = approx(j);
        double g = 0.0;
        for (int c = 0; c < n_; ++c) g += a[c] * b[c];
        if (std::fabs(g) >= kCancelled * std::sqrt(norm_sq_[i] * norm_sq_[j])) return g;
        return exact_dot(basis_.row(i), basis_.row(j), n_);
    }

    // r_kj and mu_kj for j < k; returns the largest |mu_kj|, NaN if any is NaN.
    double cholesky_row(int k) {
        double largest = 0.0;
        for (int j = 0; j < k; ++j) {
            double x = gram(k, j);
            for (int i = 0; i < j; ++i) x -= mu(j, i) * r(k, i);
            r(k, j) = x;
            mu(k, j) = x / r(j, j);
            const double m = std::fabs(mu(k, j));
            if (std::isnan(m)) return m;
            largest = std::max(largest, m);
        }
        return largest;
    }

    bool size_reduce(int k) {
        double previous = HUGE_VAL;
        for (int pass = 0;; ++pass) {
            const double largest = cholesky_row(k);
            if (!std::isfinite(largest)) return false;
            if (largest <= kEta) return true;
            if (pass == kMaxPasses || largest >= previous) return false;
            previous = largest;
            // round from the last coefficient down, carrying each step into the lower mu
            for (int j = k - 1; j >= 0; --j) {
                const double x = std::nearbyint(mu(k, j));
                step_[j] = x;
                if (x == 0.0) continue;
                for (int i = 0; i < j; ++i) mu(k, i) -= x * mu(j, i);
                mu(k, j) -= x;
            }
            for (int j = k - 1; j >= 0; --j) {
                if (step_[j] == 0.0) continue;
                if (!(std::fabs(step_[j]) < kEntryBound)) return false;
                if (!basis_.subtract(k, j, static_cast<Int128>(step_[j]))) return false;
            }
            refresh(k);
        }
    }

    void exchange(int k) {
        basis_.exchange(k);
        std::swap_ranges(approx(k - 1), approx(k - 1) + n_, approx(k));
        std::swap(norm_sq_[k - 1], norm_sq_[k]);
    }

    TrackedBasis<Int128>& basis_;
    int n_;
    std::vector<double> approx_, norm_sq_, r_, mu_, step_;
};

// Cohen's integral LLL (Algorithm 2.6.7), delta = 99/100. With d_0 = 1 and d_i
// the Gram determinant of b_0..b_{i-1}, B_i = d_{i+1} / d_i and
// lambda_kj = d_{j+1} mu_kj are integers, so every step is exact.
class ExactLll {
  public:
    explicit ExactLll(TrackedBasis<int_mp>& basis)
        : basis_(basis),
          n_(basis.n()),
          lambda_(static_cast<size_t>(n_) * n_),
          d_(static_cast<size_t>(n_) + 1) {}

    void run() {
        d_[0] = int_mp(1);
        gram_schmidt(0);
        int k = 1, k_max = 0;
        while (k < n_) {
            if (k > k_max) {
                k_max = k;
                gram_schmidt(k);
            }
            reduce(k, k - 1);
            if (lovasz_fails(k)) {
                exchange(k, k_max);
                k = std::max(1, k - 1);
                continue;
            }
            for (int l = k - 2; l >= 0; --l) reduce(k, l);
            ++k;
        }
    }

  private:
    int_mp& lambda(int k, int j) { return lambda_[static_cast<size_t>(k) * n_ + j]; }

    int_mp dot(int i, int j) const {
        int_mp s;
        for (int c = 0; c < n_; ++c) s += basis_.b(i, c) * basis_.b(j, c);
        return s;
    }

    // lambda_kj for j < k and d_{k+1} (Cohen, step 2).
    void gram_schmidt(int k) {
        for (int j = 0; j <= k; ++j) {
            int_mp u = dot(k, j);
            for (int i = 0; i < j; ++i) u = (u * d_[i + 1] - lambda(k, i) * lambda(j, i)) / d_[i];
            (j < k ? lambda(k, j) : d_[k + 1]) = std::move(u);
        }
        if (d_[k + 1].zero())
            throw std::invalid_argument("lll: the basis vectors are linearly dependent");
    }

    // RED(k, l): b_k -= q b_l with q the integer nearest mu_kl, if |mu_kl| > 1/2.
    void reduce(int k, int l) {
        int_mp& lam = lambda(k, l);
        const int_mp& d = d_[l + 1];
        const int_mp twice = lam + lam;
        if (!(twice > d) && !(twice < -d)) return;
        const int_mp q = floor_div(twice + d, d + d);
        basis_.subtract(k, l, q);
        lam -= q * d;
        for (int i = 0; i < l; ++i) lambda(k, i) -= q * lambda(l, i);
    }

    // B_k < (delta - mu^2) B_{k-1}: 100 d_{k+1} d_{k-1} < 99 d_k^2 - 100 lambda^2.
    bool lovasz_fails(int k) {
        const int_mp& lam = lambda(k, k - 1);
        return d_[k + 1] * d_[k - 1] * 100 < d_[k] * d_[k] * 99 - lam * lam * 100;
    }

    // SWAP(k): exchanges b_{k-1} and b_k and updates d_k and the lambdas (Cohen, step 4).
    void exchange(int k, int k_max) {
        basis_.exchange(k);
        for (int j = 0; j < k - 1; ++j) std::swap(lambda(k, j), lambda(k - 1, j));
        const int_mp& lam = lambda(k, k - 1);
        // the new d_k
        int_mp b = (d_[k - 1] * d_[k + 1] + lam * lam) / d_[k];
        for (int i = k + 1; i <= k_max; ++i) {
            int_mp& a = lambda(i, k);
            int_mp& c = lambda(i, k - 1);
            const int_mp t = a;
            a = (d_[k + 1] * c - lam * t) / d_[k];
            c = (b * t + lam * a) / d_[k + 1];
        }
        d_[k] = std::move(b);
    }

    TrackedBasis<int_mp>& basis_;
    int n_;
    std::vector<int_mp> lambda_, d_;
};

// The rows in 128-bit integers; false if an entry is not below 2^120 in
// absolute value.
bool to_int128(const std::vector<double>& rows, std::vector<Int128>& out) {
    out.resize(rows.size());
    for (size_t i = 0; i < rows.size(); ++i) {
        if (!(std::fabs(rows[i]) < kEntryBound)) return false;
        out[i] = static_cast<Int128>(rows[i]);
    }
    return true;
}

bool to_int128(const std::vector<int_mp>& rows, std::vector<Int128>& out) {
    out.resize(rows.size());
    for (size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].bits() > 120) return false;
        out[i] = rows[i].to_int128();
    }
    return true;
}

std::vector<int_mp> to_int_mp(const std::vector<double>& rows) {
    std::vector<int_mp> out(rows.size());
    for (size_t i = 0; i < rows.size(); ++i) out[i] = int_mp::from_double(rows[i]);
    return out;
}

std::vector<int_mp> to_int_mp(const std::vector<int_mp>& rows) { return rows; }

template <class Rows, class Out>
bool reduce(int n, const Rows& rows, int64_t limit, Out& out) {
    std::vector<Int128> wide;
    if (to_int128(rows, wide)) {
        TrackedBasis<Int128> basis(n, std::move(wide));
        if (FloatLll(basis).run()) return basis.transform(limit, out);
    }
    TrackedBasis<int_mp> basis(n, to_int_mp(rows));
    ExactLll(basis).run();
    return basis.transform(limit, out);
}

}  // namespace

std::vector<double> lll_rows(int n, const double* A, int precision_bits) {
    std::vector<double> rows(static_cast<size_t>(n) * n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            rows[static_cast<size_t>(j) * n + i] =
                std::nearbyint(std::ldexp(A[static_cast<size_t>(i) * n + j], precision_bits));
    return rows;
}

bool lll_transform(int n, const std::vector<double>& rows, int64_t limit, LllTransform& out) {
    return reduce(n, rows, limit, out);
}

bool lll_transform(int n, const std::vector<int_mp>& rows, int64_t limit, LllTransform& out) {
    return reduce(n, rows, limit, out);
}

bool lll_transform(int n, const std::vector<int_mp>& rows, int64_t limit, LllTransformMp& out) {
    return reduce(n, rows, limit, out);
}

namespace {

constexpr int kWideProductBits = 64;  // fractional bits of A in wide_product

// A T over int_mp, with A rounded to `bits` fractional bits.
void wide_product(int n, int rows, const double* A, const int64_t* T, int bits, double* out) {
    std::vector<int_mp> a(static_cast<size_t>(n));
    for (int i = 0; i < rows; ++i) {
        for (int k = 0; k < n; ++k)
            a[static_cast<size_t>(k)] = int_mp::from_double(
                std::nearbyint(std::ldexp(A[static_cast<size_t>(i) * n + k], bits)));
        for (int j = 0; j < n; ++j) {
            int_mp s;
            for (int k = 0; k < n; ++k)
                s += a[static_cast<size_t>(k)] * T[static_cast<size_t>(k) * n + j];
            // one rounding: s 2^-bits to the nearest double
            out[static_cast<size_t>(i) * n + j] = s.to_double(-bits);
        }
    }
}

// Bit width of x; ceil(log2 n) is bit_width(n - 1).
int bit_width(uint64_t x) { return x == 0 ? 0 : 64 - __builtin_clzll(x); }

// A T in 128-bit integers, each row of A scaled by a power of two to convert exactly.
// A row whose terms or partial sums could reach 2^127 goes to wide_product at the same
// scale, which gives the same rounding. False if a row spans more than 60 binades.
bool product_128(int n, int rows, const double* A, const int64_t* T, double* out) {
    uint64_t t_max = 0;
    for (size_t k = 0; k < static_cast<size_t>(n) * n; ++k)
        t_max = std::max(t_max, T[k] < 0 ? uint64_t{0} - static_cast<uint64_t>(T[k])
                                         : static_cast<uint64_t>(T[k]));
    // partial sums of a row are below 2^(a bits + sum_bits)
    const int sum_bits = bit_width(t_max) + bit_width(static_cast<uint64_t>(n) - 1);
    std::vector<__int128> a(n);
    for (int i = 0; i < rows; ++i) {
        const double* row = A + static_cast<size_t>(i) * n;
        double* row_out = out + static_cast<size_t>(i) * n;
        int lo = 0, hi = 0;
        bool any = false;
        for (int j = 0; j < n; ++j) {
            if (row[j] == 0.0) continue;
            int e;
            std::frexp(row[j], &e);
            if (!any) {
                lo = hi = e;
                any = true;
            }
            lo = std::min(lo, e);
            hi = std::max(hi, e);
        }
        if (!any) {
            for (int j = 0; j < n; ++j) row_out[j] = 0.0;
            continue;
        }
        // every entry exact at 2**-shift, the largest then below 2**(53 + range)
        const int shift = 53 - lo;
        if (hi - lo > 60) return false;
        if (53 + hi - lo + sum_bits > 127) {
            wide_product(n, 1, row, T, shift, row_out);
            continue;
        }
        const double scale = std::ldexp(1.0, shift);
        for (int j = 0; j < n; ++j) a[j] = static_cast<__int128>(row[j] * scale);
        for (int j = 0; j < n; ++j) {
            __int128 s = 0;
            for (int k = 0; k < n; ++k) s += a[k] * static_cast<__int128>(T[k * n + j]);
            row_out[j] = static_cast<double>(s) / scale;
        }
    }
    return true;
}

}  // namespace

void exact_product(int n, int rows, const double* A, const int64_t* T, double* out) {
    if (!product_128(n, rows, A, T, out)) wide_product(n, rows, A, T, kWideProductBits, out);
}

// QR decomposition of A to turn the corresponding elipsoid equations into the right form. With
// A= Q R, since Q does not change lengths, for || A (c -c0) ||^2 <= 1 we can write
// || R (c -c0) ||^2 <= 1, with R upper triangular.
QrFactors qr_positive(int m, int n, const double* A) {
    if (m > n) throw std::invalid_argument("qr_positive: more rows than columns");
    const size_t rows = static_cast<size_t>(m), cols = static_cast<size_t>(n);
    QrFactors f;
    // R starts as A and is reduced in place; Q collects the reflections,
    // H_0 H_1 ... H_{m-1}.
    f.R.reserve(cols * cols);  // room for the zero rows
    f.R.assign(A, A + rows * cols);
    f.Q.assign(rows * rows, 0.0);
    auto R = [&](int i, int j) -> double& { return f.R[i * cols + j]; };
    auto Q = [&](int i, int j) -> double& { return f.Q[i * rows + j]; };
    for (int i = 0; i < m; ++i) Q(i, i) = 1.0;
    std::vector<double> u(rows);
    for (int j = 0; j < m; ++j) {
        // The reflection I - 2 u u^t / |u|^2 taking column j from row j down to
        // (alpha, 0, ..., 0). alpha has the sign opposite R(j, j), so that
        // u_j = R(j, j) - alpha does not cancel.
        double norm = 0.0;
        for (int i = j; i < m; ++i) norm += R(i, j) * R(i, j);
        norm = std::sqrt(norm);
        if (norm == 0.0) continue;
        const double alpha = R(j, j) > 0.0 ? -norm : norm;
        for (int i = j; i < m; ++i) u[i] = R(i, j);
        u[j] -= alpha;
        double unorm = 0.0;
        for (int i = j; i < m; ++i) unorm += u[i] * u[i];
        if (unorm == 0.0) continue;
        const double scale = 2.0 / unorm;
        for (int c = j; c < n; ++c) {
            double dot = 0.0;
            for (int i = j; i < m; ++i) dot += u[i] * R(i, c);
            dot *= scale;
            for (int i = j; i < m; ++i) R(i, c) -= dot * u[i];
        }
        for (int i = j + 1; i < m; ++i) R(i, j) = 0.0;
        R(j, j) = alpha;
        // Q <- Q (I - scale u u^t)
        for (int r = 0; r < m; ++r) {
            double dot = 0.0;
            for (int i = j; i < m; ++i) dot += Q(r, i) * u[i];
            dot *= scale;
            for (int i = j; i < m; ++i) Q(r, i) -= dot * u[i];
        }
    }
    // Flip each row of R whose diagonal is negative, and its column of Q.
    for (int r = 0; r < m; ++r) {
        const double sign = R(r, r) < 0.0 ? -1.0 : 1.0;
        for (int c = 0; c < n; ++c) R(r, c) = c < r ? 0.0 : sign * R(r, c);
        for (int i = 0; i < m; ++i) Q(i, r) = sign * Q(i, r);
    }
    f.R.resize(cols * cols, 0.0);  // the zero rows m to n-1
    return f;
}

}  // namespace quicksynth
