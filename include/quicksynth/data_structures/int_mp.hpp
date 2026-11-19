// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Signed integers of any size (int_mp) and exact elements M / sqrt2^k over
// Z[zeta16] in them. Values past 2^128 occur in the Bloch-image digits of
// Clifford+T words near 1e-8 (Descent::word).
#pragma once

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

namespace quicksynth {

// Hidden: library-internal, so calls bypass the symbol table.
#pragma GCC visibility push(hidden)

// A value that fits int64 is stored in place, any other as sign and magnitude
// in 64-bit limbs on the heap (least significant first, the top limb nonzero).
// Operations on in-place values check for overflow and take no allocation.
class int_mp {
  public:
    int_mp() noexcept {}
    // From any built-in signed integer of at most 64 bits: int64_t is long long
    // on macOS and long on Linux, so one overload for long would leave
    // int_mp(int64_t) ambiguous with the __int128 one on macOS.
    template <class I, std::enable_if_t<std::is_integral_v<I> && std::is_signed_v<I> &&
                                            sizeof(I) <= sizeof(int64_t),
                                        int> = 0>
    explicit int_mp(I x) : small_(static_cast<int64_t>(x)) {}
    explicit int_mp(__int128 x) {
        if (x >= std::numeric_limits<int64_t>::min() && x <= std::numeric_limits<int64_t>::max())
            small_ = static_cast<int64_t>(x);
        else
            set_wide(x);
    }
    int_mp(const int_mp& o) : small_(o.small_) {
        if (o.wide_ != nullptr) copy_wide(o);
    }
    int_mp(int_mp&& o) noexcept : small_(o.small_), wide_(o.wide_) { o.wide_ = nullptr; }
    int_mp& operator=(const int_mp& o) {
        if (o.wide_ != nullptr) {
            if (this != &o) copy_wide(o);
        } else {
            release();
            small_ = o.small_;
        }
        return *this;
    }
    int_mp& operator=(int_mp&& o) noexcept {
        std::swap(small_, o.small_);
        std::swap(wide_, o.wide_);
        return *this;
    }
    ~int_mp() { release(); }

    // x truncated toward zero; x finite.
    static int_mp from_double(double x);
    // Parses an optional '-' and digits of base 10 or 16 into out; white space
    // before the sign and between digits is skipped. False, with out
    // unchanged, on any other character or no digit.
    static bool parse(const char* text, int base, int_mp& out);

    bool zero() const { return wide_ == nullptr && small_ == 0; }
    bool even() const {
        return ((wide_ != nullptr ? wide_->limbs()[0] : static_cast<uint64_t>(small_)) & 1) == 0;
    }
    int sign() const {
        if (wide_ != nullptr) return wide_->negative ? -1 : 1;
        return (small_ > 0) - (small_ < 0);
    }
    int bits() const;  // bit length of |x|; 0 for zero
    bool fits_int64() const { return wide_ == nullptr; }
    int64_t to_int64() const { return small_; }  // requires fits_int64()
    __int128 to_int128() const;                  // requires bits() <= 127
    // Base 10 or 16, lowercase, with '-' for a negative value.
    std::string str(int base) const;
    // x 2^exponent rounded to 53 bits, ties to even, then to a double, which
    // rounds again only below the normal range.
    double to_double(int exponent) const;
    // m and e with x truncated to 53 bits = m 2^e, 0.5 <= |m| < 1; 0 and
    // e = 0 for zero.
    double truncated_frexp(int& exponent) const;

    bool operator==(const int_mp& o) const {
        if (wide_ == nullptr || o.wide_ == nullptr)
            return wide_ == nullptr && o.wide_ == nullptr && small_ == o.small_;
        return equal_wide(o);
    }
    bool operator!=(const int_mp& o) const { return !(*this == o); }
    friend bool operator==(const int_mp& a, int64_t s);
    friend int_mp operator+(const int_mp& a, const int_mp& b);
    friend int_mp operator-(const int_mp& a, const int_mp& b);
    friend int_mp operator-(const int_mp& a);
    friend int_mp operator*(const int_mp& a, const int_mp& b);
    friend int_mp operator*(const int_mp& a, int64_t s);
    friend int_mp& operator+=(int_mp& a, const int_mp& b);
    friend int_mp& operator-=(int_mp& a, const int_mp& b);
    friend int_mp half_of(const int_mp& x);
    friend int_mp shift_left(const int_mp& x, int bits);
    friend int compare(const int_mp& a, const int_mp& b);
    friend int_mp operator/(const int_mp& a, const int_mp& b);
    friend int_mp floor_div(const int_mp& a, const int_mp& b);

  private:
    // Heap header; the limbs follow it.
    struct alignas(8) Wide {
        uint32_t size;  // limbs in use
        uint32_t capacity;
        bool negative;
        uint64_t* limbs() { return reinterpret_cast<uint64_t*>(this + 1); }
        const uint64_t* limbs() const { return reinterpret_cast<const uint64_t*>(this + 1); }
    };
    struct Digits;  // a magnitude view of either form (int_mp.cpp)

    void release() {
        if (wide_ != nullptr) {
            free_wide(wide_);
            wide_ = nullptr;
        }
    }
    void set_wide(__int128 x);
    void copy_wide(const int_mp& o);
    bool equal_wide(const int_mp& o) const;
    static Wide* allocate(uint32_t capacity);
    static void free_wide(Wide* w) noexcept;
    // Takes ownership of w and returns its value, in place if it fits int64.
    static int_mp adopt(Wide* w);
    static int_mp sum(const int_mp& a, const int_mp& b, bool subtract);
    static int_mp product(const int_mp& a, const int_mp& b);
    static int_mp quotient(const int_mp& a, const int_mp& b, bool floor);
    static int compare_wide(const int_mp& a, const int_mp& b);
    int_mp negated() const;
    int_mp halved() const;

    int64_t small_ = 0;     // the value if wide_ is null
    Wide* wide_ = nullptr;  // owned; only for values outside int64
};

inline int_mp operator+(const int_mp& a, const int_mp& b) {
    int64_t r;
    if (a.wide_ == nullptr && b.wide_ == nullptr && !__builtin_add_overflow(a.small_, b.small_, &r))
        return int_mp(r);
    return int_mp::sum(a, b, false);
}

inline int_mp operator-(const int_mp& a, const int_mp& b) {
    int64_t r;
    if (a.wide_ == nullptr && b.wide_ == nullptr && !__builtin_sub_overflow(a.small_, b.small_, &r))
        return int_mp(r);
    return int_mp::sum(a, b, true);
}

inline int_mp operator-(const int_mp& a) {
    if (a.wide_ == nullptr && a.small_ != std::numeric_limits<int64_t>::min())
        return int_mp(-a.small_);
    return a.negated();
}

inline int_mp operator*(const int_mp& a, const int_mp& b) {
    int64_t r;
    if (a.wide_ == nullptr && b.wide_ == nullptr) {
        if (!__builtin_mul_overflow(a.small_, b.small_, &r)) return int_mp(r);
        return int_mp(static_cast<__int128>(a.small_) * b.small_);
    }
    return int_mp::product(a, b);
}

inline int_mp operator*(const int_mp& a, int64_t s) {
    int64_t r;
    if (a.wide_ == nullptr) {
        if (!__builtin_mul_overflow(a.small_, s, &r)) return int_mp(r);
        return int_mp(static_cast<__int128>(a.small_) * s);
    }
    return int_mp::product(a, int_mp(s));
}

inline int_mp operator*(int64_t s, const int_mp& a) { return a * s; }

inline int_mp& operator+=(int_mp& a, const int_mp& b) {
    int64_t r;
    if (a.wide_ == nullptr && b.wide_ == nullptr && !__builtin_add_overflow(a.small_, b.small_, &r))
        a.small_ = r;
    else
        a = int_mp::sum(a, b, false);
    return a;
}

inline int_mp& operator-=(int_mp& a, const int_mp& b) {
    int64_t r;
    if (a.wide_ == nullptr && b.wide_ == nullptr && !__builtin_sub_overflow(a.small_, b.small_, &r))
        a.small_ = r;
    else
        a = int_mp::sum(a, b, true);
    return a;
}

inline bool operator==(const int_mp& a, int64_t s) { return a.wide_ == nullptr && a.small_ == s; }

// -1, 0 or 1 as a < b, a == b or a > b.
inline int compare(const int_mp& a, const int_mp& b) {
    if (a.wide_ == nullptr && b.wide_ == nullptr)
        return (a.small_ > b.small_) - (a.small_ < b.small_);
    return int_mp::compare_wide(a, b);
}
inline bool operator<(const int_mp& a, const int_mp& b) { return compare(a, b) < 0; }
inline bool operator>(const int_mp& a, const int_mp& b) { return compare(a, b) > 0; }

// a / b truncated toward zero; b nonzero.
inline int_mp operator/(const int_mp& a, const int_mp& b) {
    if (a.wide_ == nullptr && b.wide_ == nullptr &&
        !(a.small_ == std::numeric_limits<int64_t>::min() && b.small_ == -1))
        return int_mp(a.small_ / b.small_);
    return int_mp::quotient(a, b, false);
}

// floor(a / b); b nonzero.
inline int_mp floor_div(const int_mp& a, const int_mp& b) { return int_mp::quotient(a, b, true); }

// Parity and halving for the Z[lambda] templates (data_structures/bloch.hpp).
inline bool is_even(const int_mp& x) { return x.even(); }
// x / 2 truncated toward zero.
inline int_mp half_of(const int_mp& x) {
    return x.wide_ == nullptr ? int_mp(x.small_ / 2) : x.halved();
}

// x 2^bits, bits >= 0.
int_mp shift_left(const int_mp& x, int bits);

using Z16Mp = Z16Of<int_mp>;

// M / sqrt2^k with M over Z[zeta16].
struct CyclotomicMatrixMp {
    M2Of<int_mp> matrix;           // M
    int denominator_exponent = 0;  // k
};

// Multiplies out a word over H S X Y Z T Q s t q, reducing after each letter.
// Returns false on an unknown letter.
bool word_matrix_mp(const std::string& word, CyclotomicMatrixMp& out);

// Throws unknown_letter() on an unknown letter.
CyclotomicMatrixMp word_matrix_mp(const std::string& word);

// Widens u / sqrt2^k to int_mp and reduces it.
CyclotomicMatrixMp widen(const M2& u, int k);

// Parses hex digits; std::invalid_argument on a non-hex coefficient.
CyclotomicMatrixMp from_digits(const CyclotomicMatrixDigits& u);

// The j in [0, 16) with a = zeta^j b, if any.
std::optional<int> phase_between(const CyclotomicMatrixMp& a, const CyclotomicMatrixMp& b);

bool equal_up_to_phase(const CyclotomicMatrixMp& a, const CyclotomicMatrixMp& b);

#pragma GCC visibility pop

}  // namespace quicksynth
