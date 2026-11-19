#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""QuickSynth settings."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, get_args

GateSetName = Literal["sqrt_t", "clifford_t"]
StoppingRule = Literal["model", "fixed"]


@dataclass(frozen=True)
class SynthesisConfig:
    """Settings of :class:`quicksynth.QuickSynth`; frozen.

    QuickSynth minimises ``t_cost * n_T + q_cost * n_sqrtT`` in magic states
    (Weiden et al.: 1 and 3 under catalysis).
    """

    gateset: GateSetName = "sqrt_t"
    """``"clifford_t"`` restricts words to Clifford+T; cost is ``t_cost`` per T."""
    t_cost: int = 1
    """Price of a ``T`` gate; positive integer."""
    q_cost: int = 3
    """Price of a ``sqrt(T)`` gate; positive integer. Not used on Clifford+T."""
    extra_scales: int | None = None
    """Fixed stopping rule on whole scales: scan scales past the first answer
    (default: rate rule). No effect with ``fine_levels``."""
    point_budget: float | None = None
    """Scan entry price in estimated lattice points per determinant class; larger
    scales are skipped."""
    workers: int = 1
    """Lattice-walk threads; the walk is split over its top coordinates."""
    time_limit: float | None = None
    """Seconds after which no new stage starts (then :class:`TimeoutError`). Stage
    0 counts against it and stops at it; the first stage starts unless stage 0
    left a word within epsilon, which is then the answer, and a running stage
    always finishes."""
    verbose: bool = False
    """Print a line per scan walk and per stage to ``sys.stdout``."""
    stopping_rule: StoppingRule = "fixed"
    """How far the first stage searches past its first answer.

    - ``"fixed"``: with ``fine_levels`` (the default on Clifford+sqrt(T)) the
      search walks the level of best expected value past the first answer. On
      whole scales the scan continues while a fixed fraction of the best cost
      saved is worth the time, at most ``extra_scales`` past the first hit;
      below ``POLICY.split_below`` the split walks a fixed number of composite
      scales. Later stages always use it.
    - ``"model"``: a counting model of the gates within ``epsilon`` decides
      (``search_scales/stopping_rule.hpp``). The scan continues while the next scale's
      expected saving is worth ``rate`` times its time, up to the integer scale
      bound; the split picks the best-value prefix scale per composite scale.
      Clifford-symmetric targets and prices in the hundreds (model too large)
      use the fixed stopping rule."""
    fast: bool = False
    """Rate rule at ``POLICY.fast_rate`` instead of ``POLICY.rate`` magic states a
    second, in the scan and, under the fixed stopping rule, the split: on whole
    scales no walk past the first answer whose predicted time exceeds about 1 s at
    cost 50. Faster near a split scale's point budget, at a slightly higher cost.
    With ``fine_levels`` LevelSearch's rate is ``POLICY.fine_fast_rate``."""
    fine_levels: bool | None = None
    """Clifford+sqrt(T): search the four finer levels of each denominator scale
    (sublattices delta^sub, delta = 1 - zeta16) and, past the first answer, walk
    the scan or split level whose expected saving most exceeds the rate times
    its predicted time, ``POLICY.fine_rate`` (``POLICY.fine_fast_rate`` with
    ``fast``). ``None`` is on; ``False`` walks whole scales only, at
    ``POLICY.rate`` (``POLICY.fast_rate`` with ``fast``). No effect on
    Clifford+T. Only stages under the fixed stopping rule use finer levels: with
    ``stopping_rule="model"`` the first stage walks whole scales."""
    rotation_path: bool = True
    """A target within ``rotation_tolerance * epsilon`` of a Clifford-equivalent
    z-rotation is first searched on the z-rotation path (stage 0,
    ``include/quicksynth/synthesise_rotation/``): its diagonal entries level by
    level, with the norm equation solved exactly. On Clifford+T the answer to an
    exact one is certified T-optimal; on Clifford+sqrt(T) the search stops by
    ``rotation_solvable``. No answer falls through to the stages; a
    Clifford+sqrt(T) word found before a cap stopped the search is returned if
    the stages find nothing cheaper."""
    rotation_tolerance: float = 1e-2
    """Largest off-axis part, in units of ``epsilon``, that stage 0 takes; it
    searches at ``epsilon`` less that part. Exact targets, whose off-axis part
    is rounding, are taken at any value. A target just beyond it goes to the
    stages, which are far slower on near-diagonal targets at small epsilon."""
    rotation_solvable: int | None = None
    """Clifford+sqrt(T) count rule of stage 0: whole levels are searched until
    this many diagonal entries had a solvable norm equation (default
    ``POLICY.rotation_solvable``, ``POLICY.rotation_solvable_fast`` with
    ``fast``). More costs proportionally more time for a slightly lower mean
    cost."""

    def __post_init__(self) -> None:
        if self.t_cost < 1 or self.q_cost < 1:
            raise ValueError("gate prices must be positive integers")
        if self.gateset not in get_args(GateSetName):
            raise ValueError(f"unknown gateset {self.gateset!r}")
        if self.stopping_rule not in get_args(StoppingRule):
            raise ValueError(f"unknown stopping rule {self.stopping_rule!r}")
        if self.time_limit is not None and not self.time_limit > 0:
            raise ValueError("time_limit must be positive")
        if self.workers < 1:
            raise ValueError("workers must be at least 1")
        if not 0 <= self.rotation_tolerance < 1:
            raise ValueError("rotation_tolerance must lie in [0, 1)")
        if self.rotation_solvable is not None and self.rotation_solvable < 1:
            raise ValueError("rotation_solvable must be at least 1")
        object.__setattr__(self, "verbose", bool(self.verbose))
        object.__setattr__(self, "fast", bool(self.fast))
        if self.fine_levels is not None:
            object.__setattr__(self, "fine_levels", bool(self.fine_levels))
        object.__setattr__(self, "rotation_path", bool(self.rotation_path))
