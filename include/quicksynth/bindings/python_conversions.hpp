// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Conversions shared by the quicksynth._core binders to the python API.
#pragma once

#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/data_structures/bloch.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <array>
#include <complex>
#include <cstdint>
#include <string>
#include <utility>

namespace quicksynth::bindings {

namespace nb = nanobind;

using Matrix =
    nb::ndarray<const std::complex<double>, nb::shape<2, 2>, nb::c_contig, nb::device::cpu>;
// (k, 32 hex coefficients), arbitrary precision.
using DigitsArg = std::pair<int, std::array<std::string, 32>>;

Unitary unitary_of(const Matrix& m);
CyclotomicMatrixDigits digits_of(const DigitsArg& element);
// (36 Bloch image digits as Python ints, least lambda exponent).
nb::tuple bloch_tuple(const BlochDigits& b);

// ValueError unless `name` is fixed or model.
StoppingRule stopping_rule_named(const std::string& name);

// Module binders, one per library part.
void bind_exact(nb::module_& mod);
void bind_synthesis(nb::module_& mod);

}  // namespace quicksynth::bindings
