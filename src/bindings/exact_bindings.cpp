// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Bindings for the exact layer (cyclotomic_matrix.hpp, bloch.hpp, price.hpp, descent.hpp,
// normal_form.hpp).
#include "quicksynth/bindings/python_conversions.hpp"

#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/exact_synthesis/descent.hpp"
#include "quicksynth/exact_synthesis/normal_form.hpp"

#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include <tuple>

namespace quicksynth::bindings {

namespace {

Pricer pricer(bool t_count, int64_t t_cost, int64_t q_cost) {
    return t_count ? Pricer::t_count(t_cost) : Pricer::descent(t_cost, q_cost);
}

std::tuple<int64_t, int64_t, int64_t> triple(const Cost& c) { return {c.value, c.n_t, c.n_q}; }

}  // namespace

void bind_exact(nb::module_& mod) {
    mod.def(
        "word_matrix",
        [](const std::string& word) {
            const CyclotomicMatrixDigits e = word_matrix_digits(word);
            nb::list row;
            row.append(nb::int_(e.denominator_exponent));
            for (const std::string& d : e.hex)
                row.append(nb::steal(PyLong_FromString(d.c_str(), nullptr, 16)));
            return row;
        },
        nb::arg("word"),
        "The exact element of a word over H S X Y Z T Q s t q, reduced, as a list (k, then the "
        "32 coefficients as Python ints of any size); ValueError on another letter.");
    mod.def(
        "bloch_image",
        [](const DigitsArg& element) { return bloch_tuple(bloch_digits(digits_of(element))); },
        nb::arg("element"),
        "(36 integers, exponent) of an element given as (k, 32 hexadecimal coefficients): its "
        "Bloch image over Z[lambda] (row, column, power of lambda) and its least lambda "
        "exponent.");
    mod.def(
        "price",
        [](const DigitsArg& element, int64_t t_cost, int64_t q_cost, bool t_count) {
            return triple(pricer(t_count, t_cost, q_cost).price(digits_of(element)));
        },
        nb::arg("element"), nb::arg("t_cost"), nb::arg("q_cost"), nb::arg("t_count") = false,
        "(cost, n_T, n_Q) of an element given as (k, 32 hexadecimal coefficients): the weight "
        "of its normal form at (t_cost, q_cost), or with t_count t_cost times its T-count; "
        "ValueError where the descent finds no unitary over the ring, or with t_count where "
        "the element is not Clifford+T.");
    mod.def(
        "matrix_word",
        [](const DigitsArg& element, int64_t t_cost, int64_t q_cost) {
            return Descent(t_cost, q_cost).word(digits_of(element));
        },
        nb::arg("element"), nb::arg("t_cost"), nb::arg("q_cost"),
        "The word of the normal form at (t_cost, q_cost) of an element given as (k, 32 "
        "hexadecimal coefficients), checked against it up to phase; ValueError where the "
        "descent finds no unitary over the ring.");
    mod.def("normal_form", &normal_form, nb::arg("word"),
            "The normal-form spelling of a word; ValueError on an unknown letter.");
}

}  // namespace quicksynth::bindings
