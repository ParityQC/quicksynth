#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Model stopping rule through the one-call API: price check, fallback, default, split
log and worker independence. Unit tests: tests/cpp/search_scales/test_stopping_rule.cpp.
"""

from __future__ import annotations

import numpy as np
import pytest

import quicksynth as sq
from oracles import haar_random_unitary
from quicksynth.cli import build_parser

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1]).astype(complex)


def _rotation(theta: float, n: tuple[float, float, float]) -> np.ndarray:
    """exp(-i theta/2 n.sigma), n normalised."""
    x, y, z = np.asarray(n, dtype=float) / np.linalg.norm(n)
    return np.cos(theta / 2) * np.eye(2) - 1j * np.sin(theta / 2) * (
        x * X + y * Y + z * Z
    )


def _walks(r: sq.SynthesisResult, points: bool = False) -> list[tuple]:
    """Search, scales and optionally points of each logged walk."""
    return [
        (e["search"], e["scale"], e.get("prefix_scale"))
        + ((e["points"],) if points else ())
        for e in r.log
    ]


# --- the price ------------------------------------------------------------------------


def test_the_model_stopping_rule_refuses_a_t_price_above_two_q() -> None:
    with pytest.raises(ValueError, match="t_cost"):
        sq.synthesise(np.eye(2), 1e-3, t_cost=3, q_cost=1, stopping_rule="model")


# --- the fallback ---------------------------------------------------------------------


def test_a_fallback_target_synthesises_as_under_the_fixed_stopping_rule() -> None:
    V = _rotation(0.7, (0.0, 0.6, 0.5))
    runs = [sq.synthesise(V, 1e-4, stopping_rule=w) for w in ("fixed", "model")]
    assert runs[0].word == runs[1].word
    assert _walks(runs[0], points=True) == _walks(runs[1], points=True)


def test_a_price_without_a_model_synthesises_as_under_the_fixed_stopping_rule(
    calibration_targets: list[np.ndarray],
) -> None:
    # model tables would need over 1e10 costs per scale at this price
    V = calibration_targets[2]
    runs = [
        sq.synthesise(V, 1e-4, q_cost=10**9, stopping_rule=w)
        for w in ("fixed", "model")
    ]
    assert runs[0].word == runs[1].word
    assert _walks(runs[0], points=True) == _walks(runs[1], points=True)


# --- the stopping rule in the searches ------------------------------------------------


def test_fixed_is_the_default_and_unknown_stopping_rules_raise() -> None:
    assert sq.SynthesisConfig().stopping_rule == "fixed"
    assert sq.QuickSynth().config.stopping_rule == "fixed"
    assert build_parser().parse_args(["1", "1e-3"]).stopping_rule == "fixed"
    with pytest.raises(ValueError, match="stopping rule"):
        sq.SynthesisConfig(stopping_rule="adaptive")  # type: ignore[arg-type]


def test_combined_model_stopping_rule_logs_the_split(
    calibration_targets: list[np.ndarray],
) -> None:
    V = calibration_targets[1]
    r = sq.synthesise(V, 1e-4, stopping_rule="model")
    fixed = sq.synthesise(V, 1e-4, stopping_rule="fixed")
    assert r.distance <= 1e-4 and r.cost <= fixed.cost
    searches = {e["search"] for e in r.log}
    assert "scan" in searches and all(e["stage"] == 1 for e in r.log)
    for e in r.log:
        if e["search"] == "split":
            assert e["prefix_scale"] in (2, 3) and e["scale"] > max(
                x["scale"] for x in r.log if x["search"] == "scan"
            )
    # log is not serialised
    assert "log" not in r.as_dict()


@pytest.mark.slow
@pytest.mark.parametrize("eps", [1e-4, 1e-6])
def test_model_stopping_rule_is_deterministic(
    eps: float, calibration_targets: list[np.ndarray]
) -> None:
    for V in calibration_targets[:3]:
        runs = [
            sq.synthesise(V, eps, stopping_rule="model", workers=w) for w in (1, 1, 4)
        ]
        assert len({r.word for r in runs}) == 1
        assert _walks(runs[0]) == _walks(runs[1]) == _walks(runs[2])
        assert runs[0].distance <= eps


@pytest.mark.slow
def test_model_stopping_rule_on_clifford_t_below_1e_7_answers() -> None:
    V = haar_random_unitary(np.random.default_rng(2))
    r = sq.synthesise(V, 1e-8, gateset="clifford_t", stopping_rule="model")
    assert r.distance <= 1e-8 and r.n_q == 0
    assert any(e["search"] == "split" for e in r.log)
