// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/synthesise_rotation/rotation_body.hpp"

#include "quicksynth/enumerate_lattice/lll.hpp"
#include "quicksynth/data_structures/int_mp.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace quicksynth {

namespace {

// Fractional bits of the integer LLL rows at reduction scale 0, one more per
// two scales above it. The conjugate rows are O(1), so rounding moves the form
// by about 2^-48 relative; the radial row, about 2 log2(1/eps) bits larger,
// stays within the 128-bit LLL down to eps = 1e-10 and goes to the exact
// int_mp LLL below that.
constexpr int kRowBits = 48;
constexpr int64_t kTransformLimit = int64_t{1} << 62;
// The bound is inflated by these on top of kBodyInflation, so that the
// double-precision walk cannot drop a point of the body. The double-double
// rows and centre place the radial coordinate to about 2^-106 / eps^2 of the
// radial semi-axis (1e-5 at eps = 1e-13), which kBodyInflation covers.
constexpr double kRelativeMargin = 1e-6;
constexpr double kAbsoluteMargin = 1e-9;
// Walk coordinates stay well inside the exactly representable doubles.
constexpr double kCoordinateLimit = 0x1p52;

[[noreturn]] void overflow(const char* what) {
    throw std::overflow_error(std::string("RotationBody: ") + what);
}

int64_t to_int64(const int_mp& x) {
    if (!x.fits_int64()) overflow("a lattice coordinate passes 64 bits");
    return x.to_int64();
}

__int128 to_int128(const int_mp& x) {
    if (x.bits() > 126) overflow("a lattice coordinate passes 126 bits");
    return x.to_int128();
}

ComplexDD times(const ComplexDD& a, const ComplexDD& b) {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}

// 2^(K/2).
DD root_two_power(int K) {
    const int whole = K >= 0 ? K / 2 : -((1 - K) / 2);
    return dd_ldexp(K - 2 * whole == 1 ? sqrt2_dd() : DD{1, 0}, whole);
}

// Every z in Z^N with |T (z - f)|^2 <= bound, T upper triangular with positive
// diagonal [Fincke, Pohst, Math. Comp. 44 (1985)], last coordinate first,
// until emit(z, step) returns false; step is true when z is the previous z
// plus e_0. Coordinates 1 .. N - 1 run outward from the centre of their
// interval [Schnorr, Euchner, Math. Programming 66 (1994)], so that a walk cut
// short has listed the central points, those of the body, first; coordinate 0
// runs upward.
template <size_t N, class Emit>
void enumerate_ellipsoid(const std::array<double, N * N>& T, const std::array<double, N>& f,
                         double bound, Emit&& emit) {
    std::array<int64_t, N> z{};
    std::array<double, N + 1> partial{};  // partial[i]: the terms of coordinates i ..
    bool going = true;
    auto visit = [&](auto&& self, size_t i) -> void {
        double shift = 0;
        for (size_t j = i + 1; j < N; ++j)
            shift += T[i * N + j] * (static_cast<double>(z[j]) - f[j]);
        const double diagonal = T[i * N + i];
        const double rest = bound - partial[i + 1];
        if (rest < 0) return;
        const double centre = f[i] - shift / diagonal;
        const double width = std::sqrt(rest) / diagonal;
        if (!(std::fabs(centre) + width < kCoordinateLimit)) overflow("the walk leaves 52 bits");
        const auto first = static_cast<int64_t>(std::ceil(centre - width));
        const auto last = static_cast<int64_t>(std::floor(centre + width));
        bool step = false;
        const auto at = [&](int64_t v) {
            const double term = diagonal * (static_cast<double>(v) - f[i]) + shift;
            const double total = partial[i + 1] + term * term;
            if (total > bound) {
                step = false;
                return;
            }
            z[i] = v;
            partial[i] = total;
            if (i == 0) {
                going = emit(z, step);
                step = true;
            } else {
                self(self, i - 1);
            }
        };
        if (i == 0) {
            for (int64_t v = first; going && v <= last; ++v) at(v);
            return;
        }
        if (first > last) return;
        // centre, centre + 1, centre - 1, centre + 2, ...
        const int64_t middle = std::clamp(static_cast<int64_t>(std::llround(centre)), first, last);
        for (int64_t k = 0; going && (middle + k <= last || middle - k >= first); ++k) {
            if (middle + k <= last) at(middle + k);
            if (k > 0 && going && middle - k >= first) at(middle - k);
        }
    };
    visit(visit, N - 1);
}

}  // namespace

template <size_t N>
RotationBody<N>::RotationBody(const Vec4& direction, int sub, double epsilon, int reduction_scale)
    : sub_(sub), epsilon_(epsilon), reduction_scale_(reduction_scale) {
    constexpr int P = static_cast<int>(N / 2);
    if (sub < 0 || sub >= P) throw std::invalid_argument("RotationBody: sub outside [0, N / 2)");
    if (!(epsilon > 0 && epsilon < 1))
        throw std::invalid_argument("RotationBody: epsilon outside (0, 1)");
    if (reduction_scale < 0) throw std::invalid_argument("RotationBody: negative reduction scale");
    const DD ar{direction[0], 0}, ai{direction[1], 0};
    const DD length = dd_sqrt(ar * ar + ai * ai);
    if (!(length.hi > 0)) throw std::invalid_argument("RotationBody: zero direction");
    a_ = {ar / length, ai / length};
    const DD eps{epsilon, 0};
    cos_squared_ = DD{1, 0} - two_prod(epsilon, epsilon);
    const DD cosine = dd_sqrt(cos_squared_);
    // h / R = 1 - cos without cancellation.
    const DD depth = two_prod(epsilon, epsilon) / (cosine + 1.0);

    // delta^sub zeta^i
    using Ring128 = RingElement<N, __int128>;
    Ring128 delta_power = ring_zeta_power<N, __int128>(0);
    const Ring128 delta =
        ring_sub(ring_zeta_power<N, __int128>(0), ring_zeta_power<N, __int128>(1));
    for (int s = 0; s < sub; ++s) delta_power = ring_mul(delta_power, delta);
    std::array<Ring128, N> b;
    for (size_t i = 0; i < N; ++i) {
        b[i] = ring_times_zeta(delta_power, static_cast<int>(i));
        for (size_t k = 0; k < N; ++k) basis_[k][i] = static_cast<int64_t>(b[i][k]);
    }

    // The whitened basis W / R at reduction_scale: the identity place rotated by
    // a*, radial over kSegmentRadial h, tangential over kSegmentTangential eps;
    // the conjugate places as they are. Row k, column i, rounded to integers.
    const DD radial = depth * kSegmentRadial, tangential = eps * kSegmentTangential;
    const DD scale = root_two_power(reduction_scale);
    const int bits = kRowBits + (reduction_scale + 1) / 2;
    std::vector<int_mp> rows(N * N);
    for (size_t i = 0; i < N; ++i) {
        std::array<DD, N> w;
        const ComplexDD s = embed_dd(b[i], 1);
        w[0] = (s.re * a_.re + s.im * a_.im) / radial;
        w[1] = (s.im * a_.re - s.re * a_.im) / tangential;
        for (int p = 1; p < P; ++p) {
            const ComplexDD c = embed_dd(b[i], 2 * p + 1);
            w[2 * p] = c.re;
            w[2 * p + 1] = c.im;
        }
        for (size_t k = 0; k < N; ++k)
            rows[i * N + k] = nearest_integer(dd_ldexp(w[k] / scale, bits));
    }
    LllTransformMp t;
    if (!lll_transform(static_cast<int>(N), rows, kTransformLimit, t))
        overflow("the LLL transform passes 62 bits");
    transform_ = std::move(t.T);
    inverse_ = std::move(t.T_inv);
    // W T exactly, then rounded once.
    std::array<double, N * N> reduced;
    for (size_t k = 0; k < N; ++k)
        for (size_t j = 0; j < N; ++j) {
            int_mp sum;
            for (size_t i = 0; i < N; ++i) sum += rows[i * N + k] * transform_[i * N + j];
            reduced[k * N + j] = sum.to_double(-bits);
        }
    const std::vector<double> R =
        qr_positive(static_cast<int>(N), static_cast<int>(N), reduced.data()).R;
    std::copy(R.begin(), R.end(), triangle_.begin());

    // The centre u* has sigma_0(u*) = R (cos + kSegmentCentre h / R) a and
    // sigma_p(u*) = 0 elsewhere. With w* = u* / delta^sub = sum x_i zeta^i the
    // inverse Fourier transform over the N embeddings gives
    // x_i = (2 / N) Re(sigma_0(w*) zeta^-i).
    ComplexDD c{cosine + depth * kSegmentCentre, DD{}};
    c = times(c, a_);
    const int q = 16 / static_cast<int>(N);       // zeta = exp(i pi q / 16)
    const DD gap = DD{2, 0} - cos_pi16(q) * 2.0;  // |1 - zeta|^2
    const ComplexDD delta_inverse{(DD{1, 0} - cos_pi16(q)) / gap, sin_pi16(q) / gap};
    for (int s = 0; s < sub; ++s) c = times(c, delta_inverse);
    for (size_t i = 0; i < N; ++i) {
        const int m = static_cast<int>(i) * q;
        centre_[i] = (c.re * cos_pi16(m) + c.im * sin_pi16(m)) * (2.0 / static_cast<double>(N));
    }
}

template <size_t N>
std::vector<typename RotationBody<N>::Point> RotationBody<N>::points(int K, int64_t limit) const {
    std::vector<Point> out;
    walk(K, limit, true, [&](const Point& p) {
        out.push_back(p);
        return true;
    });
    return out;
}

template <size_t N>
int64_t RotationBody<N>::walk(int K, int64_t limit, bool all,
                              const std::function<bool(const Point&)>& visit) const {
    constexpr double P = static_cast<double>(N / 2);
    // x = n0 + T y, centre x* = n0 + r, y = m + z, centre T^-1 r = m + f.
    const DD R = root_two_power(K);
    std::array<int64_t, N> n0;
    std::array<__int128, N> m;
    std::array<DD, N> r;
    for (size_t i = 0; i < N; ++i) {
        const DD x = centre_[i] * R;
        n0[i] = to_int64(nearest_integer(x));
        r[i] = x - dd(n0[i]);
    }
    std::array<double, N> f;
    for (size_t j = 0; j < N; ++j) {
        DD y{};
        for (size_t i = 0; i < N; ++i) y = y + dd(inverse_[j * N + i]) * r[i];
        m[j] = to_int128(nearest_integer(y));
        f[j] = to_double(y - dd(m[j]));
    }
    const double bound = std::ldexp(kBodyInflation * P * (1 + kRelativeMargin) + kAbsoluteMargin,
                                    K - reduction_scale_);
    // A step z -> z + e_0 adds column 0 of T to x and its image to u.
    std::array<__int128, N> x_step, u_step;
    for (size_t i = 0; i < N; ++i) x_step[i] = transform_[i * N];
    for (size_t k = 0; k < N; ++k) {
        u_step[k] = 0;
        for (size_t i = 0; i < N; ++i)
            u_step[k] = checked_add(u_step[k], checked_mul(basis_[k][i], x_step[i]));
    }
    std::array<__int128, N> x;
    Point p;
    int64_t walked = 0;
    enumerate_ellipsoid<N>(triangle_, f, bound, [&](const std::array<int64_t, N>& z, bool step) {
        if (step) {
            for (size_t i = 0; i < N; ++i) x[i] = checked_add(x[i], x_step[i]);
            for (size_t k = 0; k < N; ++k) p.u[k] = checked_add(p.u[k], u_step[k]);
        } else {
            for (size_t i = 0; i < N; ++i) {
                x[i] = n0[i];
                for (size_t j = 0; j < N; ++j)
                    x[i] = checked_add(x[i],
                                       checked_mul(transform_[i * N + j], checked_add(m[j], z[j])));
            }
            for (size_t k = 0; k < N; ++k) {
                __int128 sum = 0;
                for (size_t i = 0; i < N; ++i)
                    if (basis_[k][i] != 0) sum = checked_add(sum, checked_mul(basis_[k][i], x[i]));
                p.u[k] = sum;
            }
        }
        __int128 parity = 0;
        for (size_t i = 0; i < N; ++i) parity ^= x[i] & 1;
        p.new_at_level = parity != 0;
        ++walked;
        const bool more = !(all || p.new_at_level) || visit(p);
        return more && walked <= limit;
    });
    return walked;
}

template <size_t N>
bool RotationBody<N>::contains(const Ring& u, int K) const {
    if (K < 0 || K > 125) throw std::invalid_argument("RotationBody: scale outside [0, 125]");
    RealElement<N / 2, __int128> power{};
    power[0] = static_cast<__int128>(1) << K;
    const RealElement<N / 2, __int128> xi = real_sub(power, relative_norm(u));
    if (xi != RealElement<N / 2, __int128>{} && !totally_positive(xi)) return false;
    const ComplexDD s = embed_dd(u, 1);
    const DD x = s.re * a_.re + s.im * a_.im;
    return x.hi > 0 && x * x >= dd_ldexp(cos_squared_, K);
}

template class RotationBody<4>;
template class RotationBody<8>;

}  // namespace quicksynth
