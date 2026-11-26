#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Deployment tests: any 2x2 unitary in, a checked word out.

Covers near-Clifford targets, Clifford+T below 1e-4 and at 1e-9, concurrent cold
starts, interrupts and load independence.
"""

from __future__ import annotations

import dataclasses
import json
import os
import signal
import subprocess
import sys
import textwrap
import threading
import time
from fractions import Fraction

import numpy as np
import pytest
from fingerprint import H, near_clifford

import quicksynth as sq
from conftest import fixture_json
from oracles import chord_diamond_distance, haar_random_unitary
from test_golden import _fields

# ------------------------------------------------------------- the inputs


def test_nearly_unitary_input_is_rounded_and_bad_input_rejected() -> None:
    U = haar_random_unitary(np.random.default_rng(1))
    r = sq.QuickSynth().synthesise(U.astype(np.complex64), 1e-3)
    assert r.distance <= 1e-3
    for bad in (np.eye(3), [[np.nan, 0], [0, 1]], [[1, 1], [0, 1]], "x"):
        with pytest.raises(ValueError):
            sq.synthesise(bad, 1e-3)
    for eps in (0.0, 1.0, -1e-3, float("nan"), 1e-14, 1e-15):
        with pytest.raises(ValueError):
            sq.synthesise(U, eps)


def test_answer_at_exactly_epsilon_is_never_outside() -> None:
    # the identity is at distance epsilon, up to rounding
    eps = 1e-3
    U = near_clifford(np.eye(2), 1.0, eps)
    r = sq.QuickSynth().synthesise(U, eps)
    assert r.distance <= eps


# --------------------------------------------------------- the hard targets


@pytest.mark.slow
@pytest.mark.parametrize("f", [1.01, 3.0])
def test_a_few_epsilon_from_a_clifford_is_answered(f: float) -> None:
    eps = 1e-6
    for C in (np.eye(2), H):
        r = sq.synthesise(near_clifford(C, f, eps), eps)
        assert r.distance <= eps and r.cost > 0


def test_clifford_t_below_1e_4_is_answered_and_certified() -> None:
    U = haar_random_unitary(np.random.default_rng(2))
    for eps in (1e-4, 1e-6):
        r = sq.synthesise(U, eps, gateset="clifford_t")
        assert r.distance <= eps and r.n_q == 0 and r.certified


def _u3(theta: float, phi: float, lam: float) -> np.ndarray:
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array(
        [
            [c, -np.exp(1j * lam) * s],
            [np.exp(1j * phi) * s, np.exp(1j * (phi + lam)) * c],
        ]
    )


@pytest.mark.slow
@pytest.mark.parametrize(
    ("target", "n_t"),
    [
        (haar_random_unitary(np.random.default_rng(104)), 81),
        # no answer before the split moved to prefix scale 6 below 1e-8
        (_u3(1.1, 0.7, -2.3), 83),
    ],
    ids=["haar104", "u3"],
)
def test_clifford_t_at_1e_9_is_answered_in_stage_1(
    target: np.ndarray, n_t: int
) -> None:
    r = sq.QuickSynth(gateset="clifford_t").synthesise(target, 1e-9)
    assert (r.n_t, r.n_q, r.stages) == (n_t, 0, 1)
    assert chord_diamond_distance(r.unitary, target) <= 1e-9
    assert {e["prefix_scale"] for e in r.log if e["search"] == "split"} == {6}


def test_time_limit_stops_escalation() -> None:
    # stage 1 cannot answer; no time for stage 2
    U = near_clifford(np.eye(2), 1.01, 1e-6)
    with pytest.raises(TimeoutError):
        sq.synthesise(U, 1e-6, time_limit=1e-9)


# -------------------------------------------------------- the same answer


def test_same_target_same_word_under_load() -> None:
    rng = np.random.default_rng(5)
    jobs = [haar_random_unitary(rng) for _ in range(6)]
    first = [sq.synthesise(U, 1e-6).word for U in jobs]
    out: list[str | None] = [None] * len(jobs)

    def run(i: int) -> None:
        out[i] = sq.synthesise(jobs[i], 1e-6).word

    threads = [threading.Thread(target=run, args=(i,)) for i in range(len(jobs))]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert out == first


def test_settings_follow_the_requested_epsilon_not_the_margin(capsys):
    """The search margin does not change progress output at 1e-7.

    Stage settings and prefix radii are pinned in the C++ tests.
    """
    U = haar_random_unitary(np.random.default_rng(8))
    sq.synthesise(U, 1e-7, verbose=True, stopping_rule="fixed")
    assert "\nstage 1/4 {}: cost " in capsys.readouterr().out


def test_errors_keep_their_order_and_text() -> None:
    """Errors: configuration before target before epsilon.

    The time limit keeps the caller's spelling; prices >= 2^31 are rejected.
    """
    with pytest.raises(ValueError, match="gate prices must be positive integers"):
        sq.synthesise(np.eye(3), 2.0, q_cost=0)
    with pytest.raises(TypeError):
        sq.synthesise(np.eye(3), 2.0, bogus=1)
    with pytest.raises(ValueError, match="target must be a 2x2 unitary"):
        sq.synthesise(np.eye(3), 2.0)
    U = near_clifford(np.eye(2), 1.01, 1e-6)
    for limit, text in ((1e-9, "1e-09"), (Fraction(1, 10**9), "1/1000000000")):
        with pytest.raises(
            TimeoutError,
            match=r"no word within epsilon after 1 of 4 stages and \d+ s "
            rf"\(time_limit {text} s\)",
        ):
            sq.synthesise(U, 1e-6, time_limit=limit)
    V = haar_random_unitary(np.random.default_rng(9))
    with pytest.raises(ValueError, match="t_cost exceeds two sqrt"):
        sq.synthesise(V, 1e-3, t_cost=7)
    with pytest.raises(ValueError, match="a price must be below 2\\^31"):
        sq.synthesise(V, 1e-3, q_cost=2**62)
    r = sq.synthesise(V, 1e-3, t_cost=2, q_cost=1)
    assert r.n_t == 0 and r.distance <= 1e-3
    # on Clifford+T, T stays T at any price; the descent at (6, 3) would give Q Q
    one = sq.synthesise(V, 1e-3, t_cost=1, gateset="clifford_t")
    r = sq.synthesise(V, 1e-3, t_cost=6, gateset="clifford_t")
    assert (r.n_q, r.cost, r.n_t) == (0, 6 * r.n_t, one.n_t) and r.distance <= 1e-3
    assert sq.SynthesisResult(
        "", 0, 0, 0, 0.0, 1e-3, True, "quicksynth", 0.0, 1, 3
    ).to_qasm() == ('OPENQASM 3.0;\ninclude "stdgates.inc";\nqubit[1] q;\n')


def test_verbose_lines_go_through_sys_stdout(capsys) -> None:
    """Scan and stage lines keep the format of the former Python loops."""
    U = haar_random_unitary(np.random.default_rng(3))
    r = sq.synthesise(U, 1e-2, verbose=True)
    lines = capsys.readouterr().out.splitlines()
    assert lines[-1] == "stage 1/4 {}: cost 12"
    scan = [e for e in r.log if e["search"] == "scan"]
    assert len(lines) == len(scan) + 1
    for line, e in zip(lines, scan, strict=False):
        assert line.startswith(f"    K={e['scale']} parity {e['parity']}: ")
        assert f": {e['points']} points, hits [" in line
    r = sq.synthesise(U, 1e-8, gateset="clifford_t", verbose=True)
    assert capsys.readouterr().out.splitlines()[-1] == "stage 1/3 {}: cost 72"
    # whole scales, which the former loops walked
    sq.synthesise(U, 1e-3, verbose=True, fine_levels=False)
    lines = capsys.readouterr().out.splitlines()
    assert lines[4].startswith(
        "    K=8 parity 0: 78 points, "
        "hits [16, 16, 20, 24, 24, 24, 28, 29]...  est 5408  ["
    )


def test_the_answer_ignores_the_global_phase() -> None:
    """A global phase leaves the word unchanged; the distance moves by rounding only."""
    U = haar_random_unitary(np.random.default_rng(4))
    for eps in (1e-2, 1e-5):
        a = sq.QuickSynth().synthesise(U, eps)
        b = sq.QuickSynth().synthesise(np.exp(0.7j) * U, eps)
        assert dataclasses.replace(a, seconds=0, distance=0) == dataclasses.replace(
            b, seconds=0, distance=0
        )
        assert abs(a.distance - b.distance) < 1e-15


# ------------------------------------------------------------- interrupts


def _interrupted(
    script: str, trigger: str, delay: float, *args: str
) -> tuple[list[str], dict]:
    """Runs ``script`` and sends SIGINT ``delay`` s after a ``trigger`` line.

    Returns:
        The printed lines and the JSON of the last one.
    """
    child = subprocess.Popen(
        [sys.executable, "-c", textwrap.dedent(script), *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        # line-by-line output
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
    )
    assert child.stdout is not None
    # loose timeout: a missed signal fails the test instead of hanging it
    timer = threading.Timer(600, child.kill)
    timer.start()
    try:
        lines = []
        for line in child.stdout:
            lines.append(line.rstrip("\n"))
            if line.startswith(trigger):
                time.sleep(delay)
                child.send_signal(signal.SIGINT)
                break
        lines += child.stdout.read().splitlines()
        child.wait()
    finally:
        timer.cancel()
        child.kill()
        child.stdout.close()
    assert child.returncode == 0, "\n".join(lines)
    return lines, json.loads(lines[-1])


@pytest.mark.slow
def test_an_interrupted_synthesis_leaves_nothing_that_changes_an_answer() -> None:
    """SIGINT after an empty stage 1 stops the call before stage 3.

    The same process then reproduces the uninterrupted answer and the frozen
    end-to-end result. The target is a z-rotation, which stage 0 would answer,
    so the z-rotation path is off.
    """
    eps = 1e-6
    U = near_clifford(np.eye(2), 1.01, eps)
    reference = _fields(sq.synthesise(U, eps, rotation_path=False))
    assert reference["stages"] == 3
    frozen = next(
        r
        for r in fixture_json("golden_synthesise.json")["runs"]
        if r["epsilon"] == eps and r["config"].get("gateset") == "sqrt_t"
    )
    script = """
        import dataclasses, json, math, signal, sys
        import numpy as np
        import quicksynth as qs

        signal.signal(signal.SIGINT, signal.default_int_handler)
        frozen = json.loads(sys.argv[1])

        def fields(r):
            out = dataclasses.asdict(r)
            for field in ("seconds", "log", "search"):
                out.pop(field)
            return out

        # near_clifford(np.eye(2), 1.01, 1e-6)
        U = np.eye(2) @ qs.rotation(2 * math.asin(1.01 * 1e-6), "z")
        try:
            qs.synthesise(U, 1e-6, verbose=True, rotation_path=False)
            stopped = False
        except KeyboardInterrupt:
            stopped = True
        print("stopped")
        again = fields(qs.synthesise(U, 1e-6, verbose=True, rotation_path=False))
        t = frozen["target"]
        V = np.array(t["re"]) + 1j * np.array(t["im"])
        recorded = ("method", "deterministic")
        config = {k: v for k, v in frozen["config"].items() if k not in recorded}
        # the frozen run was recorded on whole scales under the fixed stopping rule
        config = {
            **config,
            "stopping_rule": "fixed",
            "fine_levels": False,
            "rotation_path": False,
        }
        run = fields(qs.synthesise(V, frozen["epsilon"], **config))
        print(json.dumps({"stopped": stopped, "again": again, "frozen": run}))
    """
    lines, got = _interrupted(script, "stage 1/4", 0.0, json.dumps(frozen))
    # scan lines precede each stage line
    before = [
        line for line in lines[: lines.index("stopped")] if line.startswith("stage")
    ]
    assert got["stopped"] and before[0].startswith("stage 1/4")
    assert not any(line.startswith("stage 3/4") for line in before)
    assert got["again"] == reference
    # the frozen run's distance is mpmath's at 40 digits, the core's the chord
    # in double precision, both from the host's LAPACK bits (test_golden.py)
    want = {**frozen["result"], "method": "quicksynth"}
    assert got["frozen"].pop("distance") == pytest.approx(
        want.pop("distance"), abs=1e-15
    )
    assert got["frozen"] == want
