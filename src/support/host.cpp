// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/support/host.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>

namespace quicksynth {

namespace {

std::function<void()>& check() {
    static std::function<void()> f;
    return f;
}

std::function<void(const std::string&, bool)>& printer() {
    static std::function<void(const std::string&, bool)> f;
    return f;
}

}  // namespace

void set_interrupt_check(std::function<void()> f) { check() = std::move(f); }

void set_printer(std::function<void(const std::string&, bool)> f) { printer() = std::move(f); }

void poll_interrupt() {
    if (check()) check()();
}

void print_line(const std::string& line, bool flush) {
    if (printer()) return printer()(line, flush);
    std::printf("%s\n", line.c_str());
    if (flush) std::fflush(stdout);
}

std::string float_repr(double x) {
    if (x == 0) return std::signbit(x) ? "-0.0" : "0.0";
    char buf[40];
    for (int p = 1; p <= 17; ++p) {
        std::snprintf(buf, sizeof buf, "%.*e", p - 1, x);
        if (std::strtod(buf, nullptr) == x) break;
    }
    std::string mantissa(buf, std::strchr(buf, 'e'));
    const int exp10 = std::atoi(std::strchr(buf, 'e') + 1);
    const std::string sign = mantissa[0] == '-' ? "-" : "";
    std::string digits;
    for (char ch : mantissa)
        if (ch >= '0' && ch <= '9') digits += ch;
    while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
    const int point = exp10 + 1;  // digits before the decimal point
    const int n = static_cast<int>(digits.size());
    if (point > -4 && point <= 16) {
        if (point <= 0) return sign + "0." + std::string(-point, '0') + digits;
        if (point >= n) return sign + digits + std::string(point - n, '0') + ".0";
        return sign + digits.substr(0, point) + "." + digits.substr(point);
    }
    char e[16];
    std::snprintf(e, sizeof e, "e%c%02d", exp10 < 0 ? '-' : '+', std::abs(exp10));
    return sign + digits.substr(0, 1) + (n > 1 ? "." + digits.substr(1) : "") + e;
}

}  // namespace quicksynth
