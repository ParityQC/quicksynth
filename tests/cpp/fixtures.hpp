// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Readers for the JSON fixtures in tests/fixtures. Integers past 64 bits are
// decimal strings; infinities are bare Infinity / -Infinity.
#pragma once

#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <nlohmann/json.hpp>

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace quicksynth::test {

using json = nlohmann::json;

// Parsed repository file; cached.
const json& repository_json(const std::string& relative);
// tests/fixtures/golden_<group>.json
const json& golden(const std::string& group);

// Value as double; accepts Infinity, -Infinity and NaN.
double number(const json& value);
// Integer of any size as decimal digits.
std::string decimal(const json& value);
// Integer, or nullopt for null.
std::optional<int64_t> optional_of(const json& value);

// Target stored as {"re": [[a, b], [c, d]], "im": ...}.
Unitary decode_target(const json& target);
std::vector<Unitary> decode_targets(const json& targets);
// Haar targets of tests/fixtures/calibration_targets.json.
const std::vector<Unitary>& calibration_targets();

// Gate set for a ring name, "sqrtT" or "T".
GateSet ring_gateset(const std::string& ring);

// CyclotomicMatrix <-> row [k, 32 coefficients] (int64).
CyclotomicMatrix matrix_of_row(const json& row);
json row_of(const CyclotomicMatrix& element);
// Arbitrary-size element <-> row of decimal strings.
CyclotomicMatrixDigits digits_of_row(const json& row);
std::vector<std::string> decimal_row(const CyclotomicMatrixDigits& element);
CyclotomicMatrixDigits digits_of(const CyclotomicMatrix& element);
// The element unreduced, as 2^bits M / sqrt2^(k + 2 bits).
CyclotomicMatrixDigits widened(const CyclotomicMatrix& element, int bits);

}  // namespace quicksynth::test
