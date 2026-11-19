#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""One-call API against the frozen runs in ``golden_synthesise.json``."""

from __future__ import annotations

import dataclasses
import platform
import sys

import numpy as np
import pytest

import quicksynth as qs
from conftest import fixture_json
from oracles import decode_target
from quicksynth import bloch_image, word_unitary


def _synthesiser(run: dict) -> qs.QuickSynth:
    """Returns QuickSynth with the run's settings, the fixed stopping rule, whole
    scales (``fine_levels=False``) and no z-rotation path, which the recorded
    runs predate.

    The recorded ``method`` and ``deterministic`` are not settings any more.
    """
    config = {
        k: v for k, v in run["config"].items() if k not in ("method", "deterministic")
    }
    pinned = {"stopping_rule": "fixed", "fine_levels": False, "rotation_path": False}
    return qs.QuickSynth(**{**config, **pinned})


def _fields(r: qs.SynthesisResult) -> dict:
    """Result fields as recorded: all but seconds, log and search."""
    out = dataclasses.asdict(r)
    for field in ("seconds", "log", "search"):
        out.pop(field)
    return out


def _recorded(result: dict) -> dict:
    """Recorded result with method ``"scan"`` renamed to ``"quicksynth"``."""
    assert result["method"] == "scan"
    return {**result, "method": "quicksynth"}


DISTANCE_TOLERANCE = 1e-15
"""A recorded distance is mpmath's at 40 digits, and the core's the chord in
double precision; both start from the bits of the target's nearest unitary,
which numpy's SVD computes with the host's LAPACK. They part by a few 1e-16,
so the distance is compared to this absolute tolerance, the rest exactly."""


def _assert_matches(got: dict, want: dict) -> None:
    got, want = dict(got), dict(want)
    if "distance" in got or "distance" in want:
        assert got.pop("distance") == pytest.approx(
            want.pop("distance"), abs=DISTANCE_TOLERANCE
        )
    assert got == want


STRICT = sys.platform == "linux" and platform.machine() == "x86_64"
"""Whether every run must return its recorded word (Linux x86-64 only).

Elsewhere LAPACK changes the last bits of numpy's SVD polar factor; on macOS
x86-64 run 107 (a symmetric target) returns the transposed word.
"""


def _transpose_tie(U: np.ndarray, got: dict, want: dict) -> bool:
    """Whether ``got`` is ``want`` transposed, for a symmetric target off ``STRICT``.

    Both words are equally far from the target, so SVD rounding picks one. Every
    letter is symmetric up to sign, so the reversed word is the transpose.
    """
    return (
        not STRICT
        and "word" in got
        and np.array_equal(U, U.T)
        and bloch_image(word_unitary(got["word"][::-1]))
        == bloch_image(word_unitary(want["word"]))
    )


def _synthesise_cases() -> list:
    return [
        pytest.param(
            r,
            id=f"{r['config'].get('gateset', 'sqrt_t')}-{r['config']['method']}-"
            f"{r['epsilon']:.0e}-{i}",
            marks=[pytest.mark.slow] if r.get("slow") else [],
        )
        for i, r in enumerate(fixture_json("golden_synthesise.json")["runs"])
    ]


@pytest.mark.parametrize("run", _synthesise_cases())
def test_synthesise_end_to_end(run: dict) -> None:
    """Every recorded run reproduces all result fields.

    The word must match on ``STRICT`` hosts; distances within
    ``DISTANCE_TOLERANCE``.
    """
    U = decode_target(run["target"])
    synthesiser = _synthesiser(run)
    r = synthesiser.synthesise(U, run["epsilon"])
    got = _fields(r)
    want = _recorded(run["result"])
    if _transpose_tie(U, got, want):
        got["word"] = want["word"]
    _assert_matches(got, want)
