#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Exact layer through the Python package: README example, errors, spelling,
elements past 64 bits and past the double range, and long words against
Python ints.

The layer itself is tested in C++ (tests/cpp/data_structures,
price_candidates, exact_synthesis)."""

from __future__ import annotations

import json
import random
import re

import mpmath
import numpy as np
import pytest

from conftest import ROOT
from oracles import exact_bloch_image, exact_word_matrix, magic_counts
from quicksynth import (
    CyclotomicMatrix,
    bloch_image,
    matrix_word,
    normal_form,
    price,
    word_unitary,
)


def _matrix_at_50_digits(element: CyclotomicMatrix) -> mpmath.matrix:
    """``M / sqrt2**k`` from the exact coefficients, at 50 digits."""
    with mpmath.workdps(50):
        zeta = mpmath.expjpi(mpmath.mpf(1) / 8)
        scale = mpmath.sqrt(2) ** -element.sqrt2_exponent
        c = element.coefficients
        return mpmath.matrix(
            [
                [
                    scale * sum(c[16 * i + 8 * j + t] * zeta**t for t in range(8))
                    for j in range(2)
                ]
                for i in range(2)
            ]
        )


def _as_numpy(W: mpmath.matrix) -> np.ndarray:
    return np.array([[complex(W[i, j]) for j in range(2)] for i in range(2)])


def test_readme_exact_layer_example() -> None:
    """The README's exact-layer example runs as written."""
    blocks = re.findall(
        r"```python\n(.*?)```", (ROOT / "README.md").read_text(encoding="utf-8"), re.S
    )
    (example,) = [b for b in blocks if "matrix_word" in b]
    exec(example, {})


def test_invalid_input_raises() -> None:
    """Unknown letters and non-unitary elements raise ValueError.

    The non-unitary element has a Bloch image; the descent rejects it.
    """
    for spell in (word_unitary, normal_form):
        with pytest.raises(ValueError, match="unknown letter 'R'"):
            spell("HTR")
    good = word_unitary("HTSHQ")
    c = good.coefficients
    bad = CyclotomicMatrix(good.sqrt2_exponent, c[:8] * 2 + c[16:])  # [[a, a], [c, d]]
    bloch_image(bad)
    wide = CyclotomicMatrix(
        bad.sqrt2_exponent + 140, tuple(x << 70 for x in bad.coefficients)
    )
    for call in (lambda e: matrix_word(e, 1, 3), price):
        for element in (bad, wide):
            with pytest.raises(ValueError, match="descent"):
                call(element)


def test_normal_form_keeps_the_gate_and_the_magic_gates() -> None:
    rng = random.Random(4)
    for _ in range(200):
        w = "".join(rng.choice("HSXYZTQsqt") for _ in range(rng.randint(0, 50)))
        s = normal_form(w)
        assert bloch_image(word_unitary(s)) == bloch_image(word_unitary(w))
        assert magic_counts(s) == magic_counts(w)
        # at most two Clifford letters before each magic gate
        runs = [r for r in "".join("|" if c in "TtQq" else c for c in s).split("|")]
        assert all(len(r) <= 2 for r in runs[:-1])


def test_an_element_past_64_bits_crosses_as_digits() -> None:
    """An unreduced element past 64 bits keeps its Bloch image, price and word."""
    good = word_unitary("HTSHQ")
    wide = CyclotomicMatrix(
        good.sqrt2_exponent + 140, tuple(x << 70 for x in good.coefficients)
    )
    assert bloch_image(wide) == bloch_image(good)
    assert price(wide) == price(good)
    assert matrix_word(wide, 1, 3) == matrix_word(good, 1, 3)


def test_a_t_count_is_refused_off_clifford_t() -> None:
    """The T-count price rejects elements containing sqrt(T)."""
    assert price(word_unitary("HTHT"), 2, 3, t_count=True) == (4, 2, 0)
    with pytest.raises(ValueError, match="not Clifford\\+T"):
        price(word_unitary("HQ"), 1, 3, t_count=True)


def test_to_numpy_of_the_reference_elements() -> None:
    """Reference element matrices match to 4e-16 (as test_distance.cpp).

    The reference used CPython complex arithmetic; the package sums in the
    core's order.
    """
    path = ROOT / "tests" / "fixtures" / "golden_words.json"
    for row in json.loads(path.read_text(encoding="utf-8"))["elements"]:
        element = CyclotomicMatrix(row["element"][0], tuple(row["element"][1:]))
        want = np.array(row["matrix"][0::2]) + 1j * np.array(row["matrix"][1::2])
        got = element.to_numpy()
        assert got.shape == (2, 2)
        assert np.abs(got.ravel() - want).max() <= 4e-16, row["word"]


def test_elements_past_64_bits_to_numpy() -> None:
    """The matrices of the reference's elements past 64 bits, which it holds
    without matrices, against a 50-digit evaluation."""
    path = ROOT / "tests" / "fixtures" / "golden_words.json"
    rows = json.loads(path.read_text(encoding="utf-8"))["past_64_bits"]["rows"]
    for row in rows:
        element = CyclotomicMatrix(row["element"][0], tuple(row["element"][1:]))
        assert element == word_unitary(row["word"])
        want = _as_numpy(_matrix_at_50_digits(element))
        assert np.abs(element.to_numpy() - want).max() <= 1e-15, row["word"]


def test_a_word_past_the_double_range() -> None:
    """A word whose coefficients pass 2^1024, where a coefficient converted to
    a double would overflow: its matrix agrees with a 50-digit evaluation of its
    exact element."""
    rng = random.Random(16)
    word = "".join(rng.choice(("HT", "HQ", "SHT", "HQX")) for _ in range(3400))
    element = word_unitary(word)
    assert element.sqrt2_exponent > 2100
    assert max(abs(c) for c in element.coefficients).bit_length() > 1024
    W = _matrix_at_50_digits(element)
    assert np.abs(element.to_numpy() - _as_numpy(W)).max() <= 1e-15


def _check_against_python_ints(word: str) -> None:
    """CyclotomicMatrix, Bloch image, price and word of ``word`` against oracles.py.

    The Bloch image identifies an element up to phase, so the word is checked
    through it.
    """
    element = word_unitary(word)
    assert element == exact_word_matrix(word)
    image = exact_bloch_image(element)
    assert bloch_image(element) == image
    for t_cost, q_cost in ((1, 3), (4, 1)):
        spelt = matrix_word(element, t_cost, q_cost)
        assert exact_bloch_image(exact_word_matrix(spelt)) == image
        cost, n_t, n_q = price(element, t_cost, q_cost)
        assert (n_t, n_q) == magic_counts(spelt)
        assert cost == t_cost * n_t + q_cost * n_q


def test_random_words_against_python_ints() -> None:
    """Random words over every letter, from one syllable to digits of a few
    thousand bits."""
    rng = random.Random(25)
    syllables = ("HT", "HQ", "SHT", "HQX", "HtY", "HqZ", "Hs", "HTQ")
    for n in (1, 2, 5, 20, 60, 200, 600, 1500):
        _check_against_python_ints("".join(rng.choice(syllables) for _ in range(n)))


def test_the_word_past_the_double_range_against_python_ints() -> None:
    rng = random.Random(16)
    _check_against_python_ints(
        "".join(rng.choice(("HT", "HQ", "SHT", "HQX")) for _ in range(3400))
    )
