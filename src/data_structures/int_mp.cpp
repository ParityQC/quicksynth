// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/int_mp.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <new>
#include <vector>

#if defined(__SANITIZE_ADDRESS__)
#define QUICKSYNTH_ASAN 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define QUICKSYNTH_ASAN 1
#endif
#endif
#ifdef QUICKSYNTH_ASAN
// ASan's interface (sanitizer/asan_interface.h), declared here so no header is needed.
extern "C" void __asan_poison_memory_region(void const volatile* addr, size_t size);
extern "C" void __asan_unpoison_memory_region(void const volatile* addr, size_t size);
#endif

namespace quicksynth {

namespace {

using u128 = unsigned __int128;

constexpr uint64_t kTen19 = 10'000'000'000'000'000'000u;  // largest power of ten in a limb

int bit_width(uint64_t x) { return x == 0 ? 0 : 64 - __builtin_clzll(x); }

uint64_t magnitude(int64_t x) {
    return x < 0 ? 0 - static_cast<uint64_t>(x) : static_cast<uint64_t>(x);
}

// Magnitude comparison: -1, 0 or 1.
int compare_magnitudes(const uint64_t* x, uint32_t nx, const uint64_t* y, uint32_t ny) {
    if (nx != ny) return nx < ny ? -1 : 1;
    for (uint32_t i = nx; i-- > 0;)
        if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
    return 0;
}

// d = x + y; returns the size of d.
uint32_t add_magnitudes(const uint64_t* x, uint32_t nx, const uint64_t* y, uint32_t ny,
                        uint64_t* d) {
    if (nx < ny) {
        std::swap(x, y);
        std::swap(nx, ny);
    }
    uint64_t carry = 0;
    uint32_t i = 0;
    for (; i < ny; ++i) {
        const u128 s = static_cast<u128>(x[i]) + y[i] + carry;
        d[i] = static_cast<uint64_t>(s);
        carry = static_cast<uint64_t>(s >> 64);
    }
    for (; i < nx; ++i) {
        const uint64_t s = x[i] + carry;
        carry = s < carry;
        d[i] = s;
    }
    d[nx] = carry;
    return nx + static_cast<uint32_t>(carry);
}

// d = x - y for |x| >= |y|; returns the size of d.
uint32_t subtract_magnitudes(const uint64_t* x, uint32_t nx, const uint64_t* y, uint32_t ny,
                             uint64_t* d) {
    uint64_t borrow = 0;
    uint32_t i = 0;
    for (; i < ny; ++i) {
        const uint64_t a = x[i], b = y[i];
        d[i] = a - b - borrow;
        borrow = (a < b) || (a - b < borrow);
    }
    for (; i < nx; ++i) {
        const uint64_t a = x[i];
        d[i] = a - borrow;
        borrow = a < borrow;
    }
    while (nx > 0 && d[nx - 1] == 0) --nx;
    return nx;
}

// q = x / y truncated, for nx >= ny >= 1 and y[ny - 1] nonzero; q has nx - ny + 1
// limbs. Returns whether the remainder is nonzero. Knuth, TAOCP vol. 2, 4.3.1 D.
bool divide_magnitudes(const uint64_t* x, uint32_t nx, const uint64_t* y, uint32_t ny,
                       uint64_t* q) {
    if (ny == 1) {
        u128 rest = 0;
        for (uint32_t i = nx; i-- > 0;) {
            const u128 cur = (rest << 64) | x[i];
            q[i] = static_cast<uint64_t>(cur / y[0]);
            rest = cur % y[0];
        }
        return rest != 0;
    }
    // normalised copies: the top bit of v set, u one limb longer
    const int s = __builtin_clzll(y[ny - 1]);
    std::vector<uint64_t> v(ny), u(nx + 1);
    for (uint32_t i = ny; i-- > 0;)
        v[i] = (y[i] << s) | (s != 0 && i > 0 ? y[i - 1] >> (64 - s) : 0);
    u[nx] = s != 0 ? x[nx - 1] >> (64 - s) : 0;
    for (uint32_t i = nx; i-- > 0;)
        u[i] = (x[i] << s) | (s != 0 && i > 0 ? x[i - 1] >> (64 - s) : 0);
    for (uint32_t j = nx - ny + 1; j-- > 0;) {
        // estimate from the top two limbs, at most one too large afterwards
        const u128 top = (static_cast<u128>(u[j + ny]) << 64) | u[j + ny - 1];
        u128 qhat = top / v[ny - 1];
        u128 rhat = top % v[ny - 1];
        while ((qhat >> 64) != 0 || qhat * v[ny - 2] > ((rhat << 64) | u[j + ny - 2])) {
            --qhat;
            rhat += v[ny - 1];
            if ((rhat >> 64) != 0) break;
        }
        // u[j..j+ny] -= qhat v
        uint64_t carry = 0, borrow = 0;
        for (uint32_t i = 0; i < ny; ++i) {
            const u128 p = qhat * v[i] + carry;
            carry = static_cast<uint64_t>(p >> 64);
            const uint64_t a = u[i + j], b = static_cast<uint64_t>(p);
            u[i + j] = a - b - borrow;
            borrow = (a < b) || (a - b < borrow);
        }
        const uint64_t a = u[j + ny];
        u[j + ny] = a - carry - borrow;
        if ((a < carry) || (a - carry < borrow)) {
            // qhat was one too large: add v back
            --qhat;
            uint64_t c = 0;
            for (uint32_t i = 0; i < ny; ++i) {
                const u128 t = static_cast<u128>(u[i + j]) + v[i] + c;
                u[i + j] = static_cast<uint64_t>(t);
                c = static_cast<uint64_t>(t >> 64);
            }
            u[j + ny] += c;
        }
        q[j] = static_cast<uint64_t>(qhat);
    }
    return std::any_of(u.begin(), u.begin() + ny, [](uint64_t limb) { return limb != 0; });
}

// The 64 bits of x from bit `from` up.
uint64_t bits_from(const uint64_t* x, uint32_t n, int from) {
    const uint32_t q = static_cast<uint32_t>(from / 64);
    const int r = from % 64;
    uint64_t out = q < n ? x[q] >> r : 0;
    if (r != 0 && q + 1 < n) out |= x[q + 1] << (64 - r);
    return out;
}

// True if any of the bits of x below bit `below` is set.
bool any_bit_below(const uint64_t* x, int below) {
    const uint32_t q = static_cast<uint32_t>(below / 64);
    for (uint32_t i = 0; i < q; ++i)
        if (x[i] != 0) return true;
    const int r = below % 64;
    return r != 0 && (x[q] & ((uint64_t{1} << r) - 1)) != 0;
}

bool is_space(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r';
}

// Value of a digit of base <= 16, or 16 if c is none.
int digit_value(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return 16;
}

}  // namespace

// Sign and magnitude limbs of an int_mp in either form. May point into itself, so
// it is not copied, and it does not outlive x.
struct int_mp::Digits {
    explicit Digits(const int_mp& x) {
        if (x.wide_ != nullptr) {
            d = x.wide_->limbs();
            n = x.wide_->size;
            negative = x.wide_->negative;
        } else {
            one = magnitude(x.small_);
            d = &one;
            n = one != 0 ? 1 : 0;
            negative = x.small_ < 0;
        }
    }
    Digits(const Digits&) = delete;
    Digits& operator=(const Digits&) = delete;

    int bits() const { return n == 0 ? 0 : 64 * static_cast<int>(n - 1) + bit_width(d[n - 1]); }

    const uint64_t* d;
    uint32_t n;
    bool negative;
    uint64_t one;
};

// --- storage ----------------------------------------------------------------------

namespace {

// Freed blocks of 4 << c limbs (c < kClasses), kept per thread for reuse.
// The exact layer allocates and frees blocks of a few sizes at a high rate,
// which is slow through malloc. Trivially destructible, so it stays valid
// while other thread-local objects are destroyed. At most kKept blocks per
// class, 518 KiB per thread. A thread whose first allocation or free happens
// after its thread_local destructors ran (in a pthread key destructor, or
// in a static destructor on the main thread) never returns its kept blocks.
constexpr int kClasses = 8;
constexpr int kKept = 64;  // blocks per class
struct FreeBlocks {
    void* blocks[kClasses][kKept];
    int count[kClasses];
    bool armed;   // the reaper below is registered
    bool closed;  // the reaper ran; blocks go back to malloc
};
thread_local FreeBlocks free_blocks;

// Bytes of a block of class c: the 16-byte header and its limbs.
constexpr size_t block_bytes(int c) { return 16 + (size_t{4} << c) * sizeof(uint64_t); }

// Under ASan a kept block is poisoned, so that a use of a freed int_mp is still reported.
#ifdef QUICKSYNTH_ASAN
void poison(void* block, int c) { __asan_poison_memory_region(block, block_bytes(c)); }
void unpoison(void* block, int c) { __asan_unpoison_memory_region(block, block_bytes(c)); }
#else
void poison(void*, int) {}
void unpoison(void*, int) {}
#endif

// Removes a kept block of class c; count[c] > 0.
void* take_kept(int c) {
    void* block = free_blocks.blocks[c][--free_blocks.count[c]];
    unpoison(block, c);
    return block;
}

// Returns the kept blocks to malloc when the thread ends.
struct FreeBlocksReaper {
    bool used = false;
    ~FreeBlocksReaper() {
        free_blocks.closed = true;
        for (int c = 0; c < kClasses; ++c)
            while (free_blocks.count[c] > 0) ::operator delete(take_kept(c));
    }
};
thread_local FreeBlocksReaper free_blocks_reaper;

// Registers this thread's reaper.
void arm() {
    free_blocks.armed = true;
    free_blocks_reaper.used = true;  // the first use registers its destructor
}

// A new block from malloc. Arms the reaper early, so that it also runs before
// pthread key destructors and static destructors that free this thread's int_mp values.
// Out of line so that allocate stays inlined.
[[gnu::noinline]] void* fresh_block(size_t bytes) {
    if (!free_blocks.armed) arm();
    return ::operator new(bytes);
}

// The class of a block of at least `limbs` limbs, kClasses if none.
int size_class(uint32_t limbs) {
    int c = 0;
    while (c < kClasses && (uint32_t{4} << c) < limbs) ++c;
    return c;
}

}  // namespace

int_mp::Wide* int_mp::allocate(uint32_t capacity) {
    static_assert(sizeof(Wide) == 16, "block_bytes assumes a 16-byte header");
    const int c = size_class(capacity);
    if (c < kClasses) capacity = uint32_t{4} << c;
    void* memory = c < kClasses && free_blocks.count[c] > 0
                       ? take_kept(c)
                       : fresh_block(sizeof(Wide) + capacity * sizeof(uint64_t));
    Wide* w = new (memory) Wide;
    w->size = 0;
    w->capacity = capacity;
    w->negative = false;
    return w;
}

void int_mp::free_wide(Wide* w) noexcept {
    const int c = size_class(w->capacity);
    if (c < kClasses && w->capacity == (uint32_t{4} << c) && !free_blocks.closed &&
        free_blocks.count[c] < kKept) {
        if (!free_blocks.armed) arm();
        free_blocks.blocks[c][free_blocks.count[c]++] = w;
        poison(w, c);
        return;
    }
    ::operator delete(w);
}

int_mp int_mp::adopt(Wide* w) {
    while (w->size > 0 && w->limbs()[w->size - 1] == 0) --w->size;
    const uint64_t top = uint64_t{1} << 63;
    int_mp out;
    if (w->size == 0 ||
        (w->size == 1 && (w->limbs()[0] < top || (w->negative && w->limbs()[0] == top)))) {
        const uint64_t m = w->size == 0 ? 0 : w->limbs()[0];
        out.small_ = static_cast<int64_t>(w->negative ? 0 - m : m);
        free_wide(w);
    } else {
        out.wide_ = w;
    }
    return out;
}

void int_mp::set_wide(__int128 x) {
    const u128 m = x < 0 ? 0 - static_cast<u128>(x) : static_cast<u128>(x);
    Wide* w = allocate(2);
    w->limbs()[0] = static_cast<uint64_t>(m);
    w->limbs()[1] = static_cast<uint64_t>(m >> 64);
    w->size = 2;
    w->negative = x < 0;
    *this = adopt(w);
}

void int_mp::copy_wide(const int_mp& o) {
    if (wide_ == nullptr || wide_->capacity < o.wide_->size) {
        release();
        wide_ = allocate(o.wide_->size);
    }
    std::memcpy(wide_->limbs(), o.wide_->limbs(), o.wide_->size * sizeof(uint64_t));
    wide_->size = o.wide_->size;
    wide_->negative = o.wide_->negative;
    small_ = 0;
}

bool int_mp::equal_wide(const int_mp& o) const {
    return wide_->negative == o.wide_->negative &&
           compare_magnitudes(wide_->limbs(), wide_->size, o.wide_->limbs(), o.wide_->size) == 0;
}

// --- arithmetic -------------------------------------------------------------------

int_mp int_mp::sum(const int_mp& a, const int_mp& b, bool subtract) {
    const Digits x(a), y(b);
    const bool y_negative = y.negative != subtract;
    Wide* w = allocate(std::max(x.n, y.n) + 1);
    uint64_t* d = w->limbs();
    if (x.negative == y_negative) {
        w->size = add_magnitudes(x.d, x.n, y.d, y.n, d);
        w->negative = x.negative;
    } else if (compare_magnitudes(x.d, x.n, y.d, y.n) >= 0) {
        w->size = subtract_magnitudes(x.d, x.n, y.d, y.n, d);
        w->negative = x.negative;
    } else {
        w->size = subtract_magnitudes(y.d, y.n, x.d, x.n, d);
        w->negative = y_negative;
    }
    return adopt(w);
}

int_mp int_mp::product(const int_mp& a, const int_mp& b) {
    const Digits x(a), y(b);
    if (x.n == 0 || y.n == 0) return int_mp();
    // the longer operand in x, so that a scalar factor takes the one-limb loop
    const Digits& u = x.n >= y.n ? x : y;
    const Digits& v = x.n >= y.n ? y : x;
    Wide* w = allocate(u.n + v.n);
    uint64_t* d = w->limbs();
    if (v.n == 1) {
        const uint64_t m = v.d[0];
        uint64_t carry = 0;
        for (uint32_t i = 0; i < u.n; ++i) {
            const u128 t = static_cast<u128>(u.d[i]) * m + carry;
            d[i] = static_cast<uint64_t>(t);
            carry = static_cast<uint64_t>(t >> 64);
        }
        d[u.n] = carry;
    } else {
        std::fill(d, d + v.n, uint64_t{0});
        for (uint32_t i = 0; i < u.n; ++i) {
            uint64_t carry = 0;
            for (uint32_t j = 0; j < v.n; ++j) {
                const u128 t = static_cast<u128>(u.d[i]) * v.d[j] + d[i + j] + carry;
                d[i + j] = static_cast<uint64_t>(t);
                carry = static_cast<uint64_t>(t >> 64);
            }
            d[i + v.n] = carry;
        }
    }
    w->size = u.n + v.n;
    w->negative = x.negative != y.negative;
    return adopt(w);
}

int_mp int_mp::quotient(const int_mp& a, const int_mp& b, bool floor) {
    const Digits x(a), y(b);
    const bool negative = x.negative != y.negative;
    if (x.n < y.n) {
        // |a| < |b|: 0, or -1 when flooring a negative quotient
        return int_mp(floor && negative && x.n != 0 ? -1 : 0);
    }
    Wide* w = allocate(x.n - y.n + 2);
    uint64_t* d = w->limbs();
    const bool rest = divide_magnitudes(x.d, x.n, y.d, y.n, d);
    w->size = x.n - y.n + 1;
    d[w->size] = 0;
    if (floor && negative && rest) {
        // floor of a negative quotient is one further from zero
        uint32_t i = 0;
        while (++d[i] == 0) ++i;
        w->size = std::max(w->size, i + 1);
    }
    w->negative = negative;
    return adopt(w);
}

int int_mp::compare_wide(const int_mp& a, const int_mp& b) {
    const int sa = a.sign(), sb = b.sign();
    if (sa != sb) return sa < sb ? -1 : 1;
    const Digits x(a), y(b);
    const int m = compare_magnitudes(x.d, x.n, y.d, y.n);
    return sa < 0 ? -m : m;
}

int_mp int_mp::negated() const {
    const Digits x(*this);
    Wide* w = allocate(x.n);
    std::memcpy(w->limbs(), x.d, x.n * sizeof(uint64_t));
    w->size = x.n;
    w->negative = !x.negative;
    return adopt(w);
}

int_mp int_mp::halved() const {
    const Digits x(*this);
    Wide* w = allocate(x.n);
    uint64_t* d = w->limbs();
    for (uint32_t i = 0; i + 1 < x.n; ++i) d[i] = (x.d[i] >> 1) | (x.d[i + 1] << 63);
    d[x.n - 1] = x.d[x.n - 1] >> 1;
    w->size = x.n;
    w->negative = x.negative;
    return adopt(w);
}

int_mp shift_left(const int_mp& x, int bits) {
    const int_mp::Digits s(x);
    if (s.n == 0 || bits == 0) return x;
    const uint32_t q = static_cast<uint32_t>(bits / 64);
    const int r = bits % 64;
    int_mp::Wide* w = int_mp::allocate(s.n + q + 1);
    uint64_t* d = w->limbs();
    std::fill(d, d + q, uint64_t{0});
    uint64_t carry = 0;
    for (uint32_t i = 0; i < s.n; ++i) {
        d[q + i] = (s.d[i] << r) | carry;
        carry = r == 0 ? 0 : s.d[i] >> (64 - r);
    }
    d[s.n + q] = carry;
    w->size = s.n + q + 1;
    w->negative = s.negative;
    return int_mp::adopt(w);
}

// --- conversions ------------------------------------------------------------------

int int_mp::bits() const { return Digits(*this).bits(); }

__int128 int_mp::to_int128() const {
    const Digits x(*this);
    u128 m = 0;
    for (uint32_t i = x.n; i-- > 0;) m = (m << 64) | x.d[i];
    return x.negative ? -static_cast<__int128>(m) : static_cast<__int128>(m);
}

int_mp int_mp::from_double(double x) {
    int exponent;
    const double fraction = std::frexp(std::fabs(x), &exponent);  // |x| = fraction 2^exponent
    if (exponent <= 0) return int_mp();                           // |x| < 1, or 0
    // 53-bit integer m with |x| = m 2^(exponent - 53)
    int64_t m = static_cast<int64_t>(std::ldexp(fraction, 53));
    if (x < 0) m = -m;
    if (exponent < 53) return int_mp(m / (int64_t{1} << (53 - exponent)));
    return shift_left(int_mp(m), exponent - 53);
}

bool int_mp::parse(const char* text, int base, int_mp& out) {
    const char* p = text;
    while (is_space(*p)) ++p;
    const bool negative = *p == '-';
    if (negative) ++p;
    if (digit_value(*p) >= base) return false;
    // magnitude limbs; digits are gathered in chunks that fit a limb, 19 decimal or 15 hex
    std::vector<uint64_t> limbs;
    const int chunk_digits = base == 16 ? 15 : 19;
    uint64_t chunk = 0, scale = 1;
    int digits = 0;
    auto flush = [&] {
        uint64_t carry = chunk;
        for (uint64_t& limb : limbs) {
            const u128 t = static_cast<u128>(limb) * scale + carry;
            limb = static_cast<uint64_t>(t);
            carry = static_cast<uint64_t>(t >> 64);
        }
        if (carry != 0) limbs.push_back(carry);
        chunk = 0;
        scale = 1;
        digits = 0;
    };
    for (; *p != '\0'; ++p) {
        if (is_space(*p)) continue;
        const int v = digit_value(*p);
        if (v >= base) return false;
        chunk = chunk * static_cast<uint64_t>(base) + static_cast<uint64_t>(v);
        scale *= static_cast<uint64_t>(base);
        if (++digits == chunk_digits) flush();
    }
    if (digits != 0) flush();
    Wide* w = allocate(static_cast<uint32_t>(limbs.size()));
    if (!limbs.empty()) std::memcpy(w->limbs(), limbs.data(), limbs.size() * sizeof(uint64_t));
    w->size = static_cast<uint32_t>(limbs.size());
    w->negative = negative;
    out = adopt(w);
    return true;
}

std::string int_mp::str(int base) const {
    const Digits x(*this);
    if (x.n == 0) return "0";
    std::string out;
    if (base == 16) {
        static const char kHex[] = "0123456789abcdef";
        for (uint32_t i = x.n; i-- > 0;)
            for (int shift = 60; shift >= 0; shift -= 4) out += kHex[(x.d[i] >> shift) & 15];
    } else {
        // remainders by 10^19, least significant first, each 19 digits wide
        std::vector<uint64_t> q(x.d, x.d + x.n);
        while (!q.empty()) {
            u128 rest = 0;
            for (size_t i = q.size(); i-- > 0;) {
                const u128 cur = (rest << 64) | q[i];
                q[i] = static_cast<uint64_t>(cur / kTen19);
                rest = cur % kTen19;
            }
            while (!q.empty() && q.back() == 0) q.pop_back();
            uint64_t chunk = static_cast<uint64_t>(rest);
            for (int k = 0; k < 19; ++k) {
                out += static_cast<char>('0' + chunk % 10);
                chunk /= 10;
            }
        }
        std::reverse(out.begin(), out.end());
    }
    out.erase(0, out.find_first_not_of('0'));
    if (x.negative) out.insert(out.begin(), '-');
    return out;
}

double int_mp::to_double(int exponent) const {
    const Digits x(*this);
    if (x.n == 0) return 0.0;
    const int n = x.bits();
    uint64_t m = x.d[0];
    int drop = 0;
    if (n > 53) {
        drop = n - 53;
        m = bits_from(x.d, x.n, drop);
        const bool half = (bits_from(x.d, x.n, drop - 1) & 1) != 0;
        if (half && ((m & 1) != 0 || any_bit_below(x.d, drop - 1))) ++m;  // m <= 2^53, exact
    }
    const double out = std::ldexp(static_cast<double>(m), exponent + drop);
    return x.negative ? -out : out;
}

double int_mp::truncated_frexp(int& exponent) const {
    const Digits x(*this);
    exponent = x.bits();
    if (x.n == 0) return 0.0;
    const uint64_t m = exponent > 53 ? bits_from(x.d, x.n, exponent - 53) : x.d[0];
    const double out = std::ldexp(static_cast<double>(m), -std::min(exponent, 53));
    return x.negative ? -out : out;
}

}  // namespace quicksynth
