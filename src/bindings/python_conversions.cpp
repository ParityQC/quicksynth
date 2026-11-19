// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/bindings/python_conversions.hpp"

namespace quicksynth::bindings {

Unitary unitary_of(const Matrix& m) {
    const std::complex<double>* d = m.data();
    return {d[0].real(), d[0].imag(), d[1].real(), d[1].imag(),
            d[2].real(), d[2].imag(), d[3].real(), d[3].imag()};
}

CyclotomicMatrixDigits digits_of(const DigitsArg& element) {
    return {element.first, element.second};
}

nb::tuple bloch_tuple(const BlochDigits& b) {
    nb::list digits;
    for (const std::string& d : b.m)
        digits.append(nb::steal(PyLong_FromString(d.c_str(), nullptr, 10)));
    return nb::make_tuple(digits, b.exponent);
}

StoppingRule stopping_rule_named(const std::string& name) {
    if (name == "fixed") return StoppingRule::fixed;
    if (name == "model") return StoppingRule::model;
    throw nb::value_error(("unknown stopping rule " + name).c_str());
}

}  // namespace quicksynth::bindings
