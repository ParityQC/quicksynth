#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Paper counting (tests/fixtures/counting.py): level recurrence, Poisson mean
versus union bound, and epsilon-ball measure."""

from __future__ import annotations

import math

from counting import (
    EXACT_LEVELS,
    RECURRENCE,
    ball_measure,
    mean_cost_lower_bound,
    poisson_mean_cost,
)

from conftest import fixture_json


def test_exact_levels_satisfy_the_weighted_recurrence() -> None:
    # ball levels are tested in test_cost_ball.cpp
    a = EXACT_LEVELS
    for n in range(4, len(a)):
        assert a[n] == sum(c * a[n - i] for i, c in enumerate(RECURRENCE, start=1))


def test_poisson_mean_sits_a_constant_above_the_union_bound() -> None:
    gaps = [poisson_mean_cost(e) - mean_cost_lower_bound(e) for e in (1e-3, 1e-5, 1e-7)]
    assert all(0.3 < g < 0.6 for g in gaps), gaps
    # slope 3 / log2(rho), rho the dominant root of x^3 - 2x^2 - 4
    slope = (poisson_mean_cost(1e-8) - poisson_mean_cost(1e-4)) / (4 * math.log2(10))
    assert abs(slope - 2.181) < 0.01, slope


def test_ball_measure_equals_the_reference() -> None:
    for eps, mu in fixture_json("golden_window.json")["ball_measure"]:
        assert math.isclose(ball_measure(eps), mu, rel_tol=1e-14, abs_tol=1e-13)
