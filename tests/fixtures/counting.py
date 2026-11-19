#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""The paper's counting: the volume lower bound on the mean cost of approximating
Haar-random unitaries, and the Poisson means it is compared with.

If at most ``N(C)`` elements of Clifford+sqrt(T) (cosets times the 24
Cliffords) have cost at most ``C``, a Haar-random target lies within
``epsilon`` of one of them with probability at most ``N(C) * mu(epsilon)``,
where ``mu`` is the measure of an ``epsilon``-ball in PU(2). Hence

    E[cost] >= sum_C max(0, 1 - N(C) mu(epsilon)).

Everything is at the catalysis price (T = 1, sqrt(T) = 3), except
:func:`clifford_t_mean_t_count`, which counts T gates on Clifford+T. ``N(C)``
uses the enumerated level sizes to cost 20 and :data:`RECURRENCE` beyond. The
recurrence counts normal-form tuples, which bound the level sizes from above,
so the bound is unconditional.

The library's copy of ``ball_measure`` is in C++, in
``include/quicksynth/data_structures/su2.hpp``.
The two Poisson means are the expected cost that sizes CertSynth's ball;
CertSynth sums them in C++ as numpy sums them here, to the bit.
"""

from __future__ import annotations

import math

import numpy as np

EXACT_LEVELS = [
    1,
    3,
    6,
    18,
    48,
    120,
    312,
    816,
    2112,
    5472,
    14208,
    36864,
    95616,
    248064,
    643584,
    1669632,
    4331520,
    11237376,
    29153280,
    75632640,
    196214784,
]
"""Cosets of the syllable graph at each cost from 0 to 20, enumerated."""

CLIFFORD_ORDER = 24


def ball_measure(epsilon: float) -> float:
    """Normalised Haar measure of ``{U : d_diamond(U, V) < epsilon}`` in PU(2).

    The ball lifts to two caps of angular radius ``arcsin(epsilon)`` on ``S^3``,
    each of measure ``(2a - sin 2a) / (2 pi)``.
    """
    a = math.asin(epsilon)
    if a < 1e-3:  # series, to avoid cancellation
        x = 2 * a
        return (x**3 / 6 - x**5 / 120 + x**7 / 5040) / math.pi
    return (2 * a - math.sin(2 * a)) / math.pi


RECURRENCE = (2, 0, 4)
"""``a_n = 2 a_(n-1) + 4 a_(n-3)`` for ``n >= 4`` (the paper's Corollary
``cor:recurrence``): each vertex of the normal-form tree has two children per
syllable class, T at weight 1 and Q, Q^-1 at weight 3.

The normal-form tuples satisfy it by construction and bound the level sizes
from above, with equality by the tree of clusters; :data:`EXACT_LEVELS`
satisfies it from cost 4 to 20. Its dominant root is 2.5943
(``x**3 - 2 x**2 - 4 = 0``), so the bound's slope is ``3 / log2(2.5943) = 2.181``
per bit.
"""


def level_sizes(max_cost: int) -> list[float]:
    """Level sizes to ``max_cost``: :data:`EXACT_LEVELS`, then the recurrence."""
    sizes = list(map(float, EXACT_LEVELS))
    while len(sizes) <= max_cost:
        sizes.append(sum(c * sizes[-i] for i, c in enumerate(RECURRENCE, start=1)))
    return sizes


def mean_cost_lower_bound(epsilon: float, max_cost: int = 200) -> float:
    """``sum_C max(0, 1 - N(C) mu(epsilon))`` over ``C <= max_cost``."""
    mu = ball_measure(epsilon)
    cumulative = np.cumsum(level_sizes(max_cost)) * CLIFFORD_ORDER
    covered = np.minimum(1.0, cumulative * mu)
    return float(np.sum(1.0 - covered))


def poisson_mean_cost(epsilon: float, max_cost: int = 200) -> float:
    """Mean cost if the ball elements were independent Haar-random points.

    Then the number of elements within ``epsilon`` of a target is Poisson and
    ``Pr[cost > C] = exp(-N(C) mu(epsilon))``. The paper's matching upper bound
    (Proposition ``prop:poisson``) is conditional on a second-moment hypothesis
    on that number, which the Poisson case satisfies. The result lies a
    constant, about 0.45 T states, above :func:`mean_cost_lower_bound`.
    """
    mu = ball_measure(epsilon)
    cumulative = np.cumsum(level_sizes(max_cost)) * CLIFFORD_ORDER
    return float(np.sum(np.exp(-cumulative * mu)))


def clifford_t_mean_t_count(epsilon: float) -> float:
    """:func:`poisson_mean_cost` on Clifford+T, by T-count: the cosets of the
    Clifford+T syllable graph number 1, 3, 6, 12, 24, ... at T-count 0, 1, 2,
    ..., and the sum runs to T-count 120."""
    mu = ball_measure(epsilon)
    levels = np.array([1.0] + [3.0 * 2.0**c for c in range(120)])
    cumulative = np.cumsum(levels) * CLIFFORD_ORDER
    return float(np.sum(np.exp(-cumulative * mu)))
