// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Host-language hooks, set once before any call.
//
// The interrupt check runs between long steps (QuickSynth stages, cost-ball and
// prefix-collection levels, prefix radii) and throws to stop the call; nothing
// from a stopped call is cached. Unset: no check, lines go to stdout.
#pragma once

#include <functional>
#include <string>

namespace quicksynth {

void set_interrupt_check(std::function<void()> check);
void set_printer(std::function<void(const std::string& line, bool flush)> printer);

void poll_interrupt();
// Prints one progress line (`line` has no newline).
void print_line(const std::string& line, bool flush);

// Python's repr of a finite float: shortest round-trip digits, fixed-point
// when the decimal point is within 16 digits.
std::string float_repr(double x);

}  // namespace quicksynth
