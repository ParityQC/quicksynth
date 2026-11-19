#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Tests for the oracles.py SU(2) vector and chord distance.

The library versions are tested in tests/cpp/data_structures/test_su2.cpp."""

from __future__ import annotations

import math

import numpy as np

from oracles import chord_diamond_distance, diamond_distance, su2_vector


def test_su2_vector_is_phase_invariant() -> None:
    rng = np.random.default_rng(0)
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, _ = np.linalg.qr(z)
    v1, v2 = su2_vector(q), su2_vector(np.exp(0.7j) * q)
    assert np.allclose(v1, v2) or np.allclose(v1, -v2)


def test_chord_distance_is_accurate_where_the_trace_formula_is_not() -> None:
    for theta in (0.3, 1e-3, 1e-6, 1e-9):
        rz = np.diag([np.exp(-0.5j * theta), np.exp(0.5j * theta)])
        exact = abs(math.sin(theta / 2))
        chord = chord_diamond_distance(rz, np.eye(2))
        assert abs(chord - exact) <= 1e-12 * exact
        if theta >= 1e-3:  # trace formula still accurate here
            assert abs(diamond_distance(rz, np.eye(2)) - exact) <= 1e-6 * exact
    assert diamond_distance(np.eye(2), np.eye(2)) == 0.0
    assert chord_diamond_distance(np.eye(2), -np.eye(2)) == 0.0  # phase invariant
