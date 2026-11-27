#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Approximate synthesis into priced gate words and exact synthesis algorithm.

    >>> import numpy as np, quicksynth as qs
    >>> r = qs.QuickSynth().synthesise(qs.rotation(np.pi / 4), 1e-3)
    >>> r.word, r.cost, r.certified   # R_z(pi/4) is T up to phase
    ('tS', 1, False)

Results are verified (distance and price) and deterministic. The searches run
in the native core (``include/quicksynth/synthesise.hpp``).

An exact element is ``M / sqrt2**k`` with ``M`` over ``Z[zeta16]``,
``zeta = exp(i pi / 8)``, and ``k`` minimal. :func:`normal_form` gives the normal
form

    A_1 M_1  A_2 M_2  ...  A_n M_n  C

with ``M_i`` in ``T Q t q``, ``A_i`` in ``"" H SH`` (left cosets of <S, X>) and
``C`` a shortest Clifford word; cost is preserved, global phase is not. The
exact layer wraps ``data_structures/``, ``price_candidates/`` and
``exact_synthesis/``.
"""

from __future__ import annotations

import dataclasses
import math
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any, NamedTuple

import numpy as np

import quicksynth._core as _core
from quicksynth.config import SynthesisConfig

GATE_QASM = {
    "H": "h",
    "S": "s",
    "s": "sdg",
    "X": "x",
    "Y": "y",
    "Z": "z",
    "T": "t",
    "t": "tdg",
    "Q": "p(pi/8)",  # sqrt(T) = diag(1, e^{i pi/8})
    "q": "p(-pi/8)",
}

UNITARITY_TOLERANCE = 1e-6
"""Largest accepted ``||U^dagger U - I||``; ``epsilon`` is measured from the
nearest unitary."""

MIN_EPSILON = 1e-13
"""Smallest accepted ``epsilon``. A double target is uncertain at about 1e-16,
and the search runs at ``epsilon - POLICY.epsilon_margin``."""


@dataclass(frozen=True)
class _Policy:
    epsilon_margin: float
    split_below: float
    rate: float
    fast_rate: float
    fine_fast_rate: float
    fine_rate: float
    prefix_radius: int
    prefix_count: int
    rotation_solvable: int
    rotation_solvable_fast: int
    scan_throughput: Callable[[str], float]


POLICY = _Policy(
    epsilon_margin=_core.EPSILON_MARGIN,
    split_below=_core.SPLIT_BELOW,
    rate=_core.RATE,
    fast_rate=_core.FAST_RATE,
    fine_fast_rate=_core.FINE_FAST_RATE,
    fine_rate=_core.FINE_RATE,
    prefix_radius=_core.PREFIX_RADIUS,
    prefix_count=_core.PREFIX_COUNT,
    rotation_solvable=_core.ROTATION_SOLVABLE,
    rotation_solvable_fast=_core.ROTATION_SOLVABLE_FAST,
    scan_throughput=_core.reference_throughput,
)
"""Read-only core policy constants (``policy.hpp``, ``search_scales/``,
``synthesise_rotation/``).

- ``epsilon_margin``: searches run at ``epsilon - epsilon_margin``.
- ``split_below``, ``rate``: split threshold, rate-rule factor.
- ``fast_rate``: rate-rule factor of ``fast``.
- ``fine_rate``, ``fine_fast_rate``: LevelSearch's rate on Clifford+sqrt(T), by
  default and under ``fast``.
- ``prefix_radius``, ``prefix_count``: split prefix ball.
- ``rotation_solvable``, ``rotation_solvable_fast``: stage 0's Clifford+sqrt(T)
  count rule, by default and with ``fast``.
- ``scan_throughput(gateset)``: single-thread walk rate the fixed stopping rule uses.
"""


# --- the exact layer ------------------------------------------------------------------

_SQRT2 = math.sqrt(2.0)
_COS1 = math.sqrt(2.0 + _SQRT2) / 2  # cos(pi/8)
_SIN1 = math.sqrt(2.0 - _SQRT2) / 2  # sin(pi/8)
_ZETA_POWERS = (
    (1.0, 0.0),
    (_COS1, _SIN1),
    (_SQRT2 / 2, _SQRT2 / 2),
    (_SIN1, _COS1),
    (0.0, 1.0),
    (-_SIN1, _COS1),
    (-_SQRT2 / 2, _SQRT2 / 2),
    (-_COS1, _SIN1),
)
"""(Re, Im) of ``zeta**t`` for t < 8, as ``zeta_power`` in
``data_structures/su2.cpp``."""


def _sqrt2_power(n: int) -> float:
    """``2**(n/2)``, as ``sqrt2_power`` in ``data_structures/su2.cpp``."""
    if n < 0:
        return math.ldexp(1.0 if n % 2 == 0 else math.sqrt(0.5), -(-n // 2))
    return math.ldexp(1.0 if n % 2 == 0 else _SQRT2, n // 2)


def _as_double(c: int, shift: int) -> float:
    """``c * 2**-shift`` converted as the core does: rounded within 64 bits,
    truncated to 53 bits beyond (``int_mp::truncated_frexp``)."""
    if -(2**63) <= c < 2**63:
        return math.ldexp(float(c), -shift)
    drop = abs(c).bit_length() - 53
    x = math.ldexp(float(abs(c) >> drop), drop - shift)
    return -x if c < 0 else x


class CyclotomicMatrix(NamedTuple):
    """An exact element ``M / sqrt2**k`` with ``k`` least."""

    sqrt2_exponent: int
    coefficients: tuple[int, ...]
    """32 coefficients of ``M``, row-major, each entry over ``zeta**0..7``."""

    def to_numpy(self) -> np.ndarray:
        """Returns ``M / sqrt2**k`` as doubles, bit-identical to the core
        (``data_structures/su2.cpp``, ``data_structures/distance.cpp``)."""
        # coefficients above 2**1000 are scaled by 2**-shift to stay finite
        widest = max(abs(c).bit_length() for c in self.coefficients)
        shift = max(widest - 1000, 0)
        x = [_as_double(c, shift) for c in self.coefficients]
        scale = _sqrt2_power(2 * shift - self.sqrt2_exponent)
        entries = []
        for ij in range(4):
            # core's summation order, for identical bits
            re = im = 0.0
            for t in range(8):
                re += x[8 * ij + t] * _ZETA_POWERS[t][0]
                im += x[8 * ij + t] * _ZETA_POWERS[t][1]
            entries.append(complex(re * scale, im * scale))
        return np.array(entries, dtype=np.complex128).reshape(2, 2)


def word_unitary(word: str) -> CyclotomicMatrix:
    """Returns the exact element of a word; :class:`ValueError` on an unknown letter."""
    k, *coefficients = _core.word_matrix(word)
    return CyclotomicMatrix(k, tuple(coefficients))


def _digits(element: CyclotomicMatrix) -> tuple[int, list[str]]:
    """``k`` and the 32 coefficients in hexadecimal."""
    return element.sqrt2_exponent, [format(c, "x") for c in element.coefficients]


def normal_form(word: str) -> str:
    """Returns the shortest normal-form spelling of ``word``, verified exactly."""
    return _core.normal_form(word)


def matrix_word(element: CyclotomicMatrix, t_cost: int, q_cost: int) -> str:
    """Returns a normal-form word of cost ``price(element)[0]``, verified up to phase.

    ``T`` is spelt ``Q Q`` when ``2 * q_cost <= t_cost``. Cliffords are spelt
    out, so ``T^dagger`` becomes ``X T X``."""
    return _core.matrix_word(_digits(element), t_cost, q_cost)


def price(
    element: CyclotomicMatrix, t_cost: int = 1, q_cost: int = 3, t_count: bool = False
) -> tuple[int, int, int]:
    """Returns ``(cost, n_T, n_Q)`` of the element's normal form.

    With ``t_count``, cost is ``t_cost`` times the T-count; :class:`ValueError`
    if the element is not Clifford+T."""
    return _core.price(_digits(element), t_cost, q_cost, t_count)


def bloch_image(element: CyclotomicMatrix) -> tuple[tuple[int, ...], int]:
    """Returns the Bloch image over ``Z[lambda]``, ``lambda = 2 cos(pi / 8)``: 36
    coefficients (row, column, power) and the minimal lambda-exponent."""
    digits, exponent = _core.bloch_image(_digits(element))
    return tuple(digits), int(exponent)


# --- the result -----------------------------------------------------------------------


@dataclass(frozen=True)
class SynthesisResult:
    """Verified result of :meth:`QuickSynth.synthesise`."""

    word: str
    """Letters ``H S X Y Z T Q`` and inverses ``s t q`` as a matrix product:
    ``H Q S`` is ``H @ Q @ S``, so the rightmost acts first. At most two
    Clifford letters between magic gates (see :func:`normal_form`)."""
    cost: int
    """``t_cost * n_t + q_cost * n_q``."""
    n_t: int
    """Number of ``T`` and ``t`` letters."""
    n_q: int
    """Number of ``Q`` and ``q`` letters."""
    distance: float
    """Diamond distance to the target up to phase, from the exact element.
    Error <= ~5e-16, so only answers within ``epsilon - 5e-16`` are accepted."""
    epsilon: float
    """Requested tolerance."""
    certified: bool
    """True if the search proved no cheaper word exists (on Clifford+T: the scan
    reached the scale bound, or stage 0 its level bound on an exactly
    Clifford-diagonal target); otherwise ``cost`` is an upper bound."""
    method: str
    """``"quicksynth"``."""
    seconds: float
    """Wall-clock seconds of the call."""
    t_cost: int
    """``T`` price used."""
    q_cost: int
    """``sqrt(T)`` price used."""
    gateset: str = "sqrt_t"
    """``"sqrt_t"`` or ``"clifford_t"``."""
    stages: int = 1
    """Escalation stages used (1: default settings; 0: the z-rotation path of
    ``rotation_path`` answered)."""
    log: list[dict[str, Any]] = field(default_factory=list, compare=False, repr=False)
    """Walk log for time-model calibration; not compared or serialised. One dict
    per scan walk (scale, determinant class) and split composite scale: ``stage``,
    ``search``, ``scale``, ``parity`` or ``prefix_scale``, ``seconds``,
    ``estimated_points``, ``points`` (leaves), ``best_cost`` (None before any
    answer). The z-rotation path (stage 0, ``search`` ``"rotation"``) logs one
    per level and class, with ``sub``, ``points`` (body points) and the counts
    ``candidates``, ``solvable``, ``undecided``, ``capped`` and ``pricings``."""
    search: dict[str, Any] | None = field(default=None, compare=False, repr=False)
    """Answering-stage summary; not compared or serialised. Keys: ``scales``,
    ``first_hit_scale``, ``last_full_scale``, ``lower_bound`` (from fully walked
    scales), ``estimated_points``, ``scan_cost``, ``scan_certified``,
    ``split_ran``, ``source`` (``"scan"``, ``"split"``, ``"rotation"`` or
    ``"clifford_t"``: on Clifford+sqrt(T), a cheaper pure-T word from the
    Clifford+T scan)."""

    @property
    def gates(self) -> list[str]:
        """Gates in time order; the first acts first."""
        return list(reversed(self.word))

    @property
    def unitary(self) -> np.ndarray:
        """The word's matrix, from its exact element."""
        return word_unitary(self.word).to_numpy()

    def to_qasm(self, qubit: str = "q[0]") -> str:
        """OpenQASM 3 with ``sqrt(T)`` as the phase gate ``p(pi/8)``."""
        lines = ["OPENQASM 3.0;", 'include "stdgates.inc";', "qubit[1] q;"]
        lines += [f"{GATE_QASM[g]} {qubit};" for g in self.gates]
        return "\n".join(lines) + "\n"

    def as_dict(self) -> dict[str, Any]:
        """Returns the compared fields and :attr:`gates`, for JSON."""
        d = {
            f.name: getattr(self, f.name) for f in dataclasses.fields(self) if f.compare
        }
        return {"word": d.pop("word"), "gates": self.gates, **d}


# --- the targets ----------------------------------------------------------------------


def rotation(theta: float, axis: str = "z") -> np.ndarray:
    """``exp(-i theta sigma_axis / 2)``."""
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    if axis == "z":
        return np.array([[c - 1j * s, 0], [0, c + 1j * s]])
    if axis == "x":
        return np.array([[c, -1j * s], [-1j * s, c]])
    if axis == "y":
        return np.array([[c, -s], [s, c]], dtype=complex)
    raise ValueError(f"axis must be x, y or z, not {axis!r}")


def nearest_unitary(target: Any) -> np.ndarray:
    """The target as a 2x2 complex matrix, replaced by its nearest unitary.

    Raises :class:`ValueError` unless it is a finite 2x2 matrix within
    :data:`UNITARITY_TOLERANCE` of unitary.

    It stays in numpy: the target is replaced by numpy's SVD polar factor, whose
    last bits depend on the LAPACK build, and every search, distance and
    final check downstream starts from those bits.
    """
    try:
        U = np.asarray(target, dtype=complex)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"target is not a numeric 2x2 matrix: {exc}") from exc
    if U.shape != (2, 2):
        raise ValueError(f"target must be a 2x2 unitary, not shape {U.shape}")
    if not np.all(np.isfinite(U)):
        raise ValueError("target has non-finite entries")
    defect = float(np.linalg.norm(U.conj().T @ U - np.eye(2)))
    if defect > UNITARITY_TOLERANCE:
        raise ValueError(
            f"target is not unitary: ||U^dagger U - I|| = {defect:.1e} exceeds "
            f"{UNITARITY_TOLERANCE:g}"
        )
    W, _, Vh = np.linalg.svd(U)
    return W @ Vh


def _epsilon(epsilon: Any) -> float:
    try:
        epsilon = float(epsilon)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"epsilon must be a number: {exc}") from exc
    if not 0 < epsilon < 1:
        raise ValueError("epsilon must lie in (0, 1)")
    if epsilon < MIN_EPSILON:
        raise ValueError(
            f"epsilon below {MIN_EPSILON:g} is beneath what a double-precision "
            "target defines"
        )
    return epsilon


def _config(
    config: SynthesisConfig | None, overrides: dict[str, Any]
) -> SynthesisConfig:
    cfg = config or SynthesisConfig()
    return dataclasses.replace(cfg, **overrides) if overrides else cfg


# --- the synthesiser ------------------------------------------------------------------


class QuickSynth:
    """Finds a cheap word fast at any precision.

    Runs the lattice scan (whole denominator scales), then the prefix split
    (deeper scales for the same work); the cheaper answer wins. Stages escalate
    prefixes and budgets until one answers. ``overrides`` are
    :class:`SynthesisConfig` fields. :meth:`synthesise` raises
    :class:`TimeoutError` at ``time_limit`` and :class:`RuntimeError` when no
    stage answers."""

    method = "quicksynth"

    def __init__(self, config: SynthesisConfig | None = None, **overrides: Any):
        self.config = cfg = _config(config, overrides)
        self._core = _core.QuickSynth(
            cfg.gateset,
            cfg.t_cost,
            cfg.q_cost,
            extra_scales=None if cfg.extra_scales is None else int(cfg.extra_scales),
            point_budget=None if cfg.point_budget is None else float(cfg.point_budget),
            workers=cfg.workers,
            time_limit=None if cfg.time_limit is None else float(cfg.time_limit),
            time_limit_text=""
            if cfg.time_limit is None
            else format(cfg.time_limit, ""),
            verbose=cfg.verbose,
            stopping_rule=cfg.stopping_rule,
            fast=cfg.fast,
            rate=cfg.rate,
            fine_levels=cfg.fine_levels,
            rotation_path=cfg.rotation_path,
            rotation_tolerance=float(cfg.rotation_tolerance),
            rotation_solvable=None
            if cfg.rotation_solvable is None
            else int(cfg.rotation_solvable),
        )

    def synthesise(self, target: Any, epsilon: float) -> SynthesisResult:
        U = nearest_unitary(target)
        epsilon = _epsilon(epsilon)
        start = time.time()
        r = self._core.synthesise(U, epsilon)
        cfg = self.config
        return SynthesisResult(
            word=r.word,
            cost=r.cost,
            n_t=r.n_t,
            n_q=r.n_q,
            distance=r.distance,
            epsilon=epsilon,
            certified=r.certified,
            method=self.method,
            seconds=time.time() - start,
            t_cost=cfg.t_cost,
            q_cost=cfg.q_cost,
            gateset=cfg.gateset,
            stages=r.stages,
            log=r.log,
            search=r.search,
        )

    def prepare(self, epsilon: float) -> None:
        """Builds every stage's prefix ball and prefix sets in memory, so a timed
        call at ``epsilon`` pays for neither."""
        self._core.prepare(_epsilon(epsilon))

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.config!r})"


def synthesise(
    target: Any,
    epsilon: float,
    config: SynthesisConfig | None = None,
    **overrides: Any,
) -> SynthesisResult:
    """Returns ``QuickSynth(config, **overrides).synthesise(target, epsilon)``.

    ``target``: 2x2 array-like within 1e-6 of unitary. ``overrides``:
    :class:`SynthesisConfig` fields, e.g. ``synthesise(U, 1e-3, q_cost=2)``.
    """
    return QuickSynth(config, **overrides).synthesise(target, epsilon)


def synthesise_rotation(
    theta: float,
    epsilon: float,
    axis: str = "z",
    config: SynthesisConfig | None = None,
    **overrides: Any,
) -> SynthesisResult:
    """:func:`synthesise` for a rotation ``exp(-i theta sigma_axis / 2)``."""
    return synthesise(rotation(theta, axis), epsilon, config, **overrides)
