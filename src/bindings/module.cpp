// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// nanobind module quicksynth._core: the exact layer (exact_bindings.cpp) and the
// one-call synthesis (synthesis_bindings.cpp). Built by scikit-build-core
// (pyproject.toml).
#include "quicksynth/bindings/python_conversions.hpp"

#include "quicksynth/support/host.hpp"

#include <nanobind/stl/string.h>

#include <string>

namespace nb = nanobind;

namespace {

// Installs the host.hpp hooks: KeyboardInterrupt stops long calls between
// steps; progress goes to Python's print.
void bind_host() {
    quicksynth::set_interrupt_check([] {
        nb::gil_scoped_acquire acquire;
        if (PyErr_CheckSignals() != 0) throw nb::python_error();
    });
    quicksynth::set_printer([](const std::string& line, bool flush) {
        nb::gil_scoped_acquire acquire;
        nb::module_::import_("builtins").attr("print")(line, nb::arg("flush") = flush);
    });
}

}  // namespace

NB_MODULE(_core, mod) {
    mod.doc() = "Native core of quicksynth: the exact layer and the one-call synthesis.";
    quicksynth::bindings::bind_exact(mod);
    quicksynth::bindings::bind_synthesis(mod);
    bind_host();
}
