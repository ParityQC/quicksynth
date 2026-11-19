#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Public API: QuickSynth, the one-call functions and results."""

from __future__ import annotations

import dataclasses
import json
import math

import numpy as np
import pytest

import quicksynth as sq
from oracles import (
    chord_diamond_distance,
    diamond_distance,
    haar_random_unitary,
    magic_counts,
)
from quicksynth import api


@pytest.fixture
def quicksynth() -> sq.QuickSynth:
    return sq.QuickSynth()


def test_rotation_by_pi_over_4_is_one_t_gate(quicksynth) -> None:
    r = quicksynth.synthesise(sq.rotation(math.pi / 4), 1e-3)
    assert (r.cost, r.n_t, r.n_q, r.certified, r.method) == (
        1,
        1,
        0,
        False,
        "quicksynth",
    )
    assert r.distance < 1e-9
    # T^dagger S is T up to phase
    assert r.word == "tS" and r.gates == ["S", "t"]
    assert "tdg q[0];" in r.to_qasm()


def test_rotation_by_pi_over_8_is_one_sqrt_t_gate(quicksynth) -> None:
    r = quicksynth.synthesise(sq.rotation(math.pi / 8), 1e-3)
    assert (r.cost, r.n_t, r.n_q) == (3, 0, 1)
    assert "p(pi/8) q[0];" in r.to_qasm()
    cheap = sq.QuickSynth(quicksynth.config, q_cost=1)
    assert cheap.synthesise(sq.rotation(math.pi / 8), 1e-3).cost == 1


def test_scan_method_gives_a_correct_word_too(quicksynth) -> None:
    r = sq.synthesise_rotation(math.pi / 4, 1e-3, config=quicksynth.config)
    assert (r.cost, r.certified, r.method) == (1, False, "quicksynth")
    assert r.distance < 1e-6  # chord formula float floor ~1e-8


def test_synthesise_is_quicksynths_one_call() -> None:
    """synthesise equals QuickSynth.synthesise; the config has no method."""
    U = haar_random_unitary(np.random.default_rng(3))
    quick = sq.QuickSynth().synthesise(U, 1e-2)
    assert quick.method == "quicksynth" and quick.search is not None
    again = sq.synthesise(U, 1e-2)
    assert dataclasses.replace(again, seconds=0) == dataclasses.replace(
        quick, seconds=0
    )
    assert "method" not in {f.name for f in dataclasses.fields(sq.SynthesisConfig)}
    with pytest.raises(TypeError, match="method"):
        sq.SynthesisConfig(method="optimal")  # type: ignore[call-arg]


def test_random_target_within_epsilon(quicksynth) -> None:
    rng = np.random.default_rng(7)
    target = haar_random_unitary(rng)
    for eps in (0.2, 0.05):
        r = quicksynth.synthesise(target, eps)
        assert r.distance <= eps
        assert r.method == "quicksynth"
        n_t, n_q = magic_counts(r.word)
        assert r.cost == r.n_t + 3 * r.n_q == n_t + 3 * n_q
        assert diamond_distance(r.unitary, target) == pytest.approx(r.distance)
        assert r.as_dict()["cost"] == r.cost


def test_config_validation() -> None:
    with pytest.raises(ValueError):
        sq.SynthesisConfig(q_cost=0)
    with pytest.raises(ValueError):
        sq.synthesise(np.eye(3), 1e-2)
    with pytest.raises(ValueError):
        sq.synthesise(np.array([[1, 1], [0, 1]]), 1e-2)
    with pytest.raises(ValueError):
        sq.synthesise(np.eye(2), 2.0)
    with pytest.raises(TypeError, match="cache_dir"):
        sq.SynthesisConfig(cache_dir="balls")  # type: ignore[call-arg]
    with pytest.raises(ValueError, match="unknown gate set x"):
        api._core.QuickSynth(gateset="x", t_cost=1, q_cost=3)


def test_settings_reach_the_core_as_bool() -> None:
    """verbose reaches the core as bool; a very large price still answers."""
    U = haar_random_unitary(np.random.default_rng(5))
    for eps in (1e-3, 1e-4):
        word = sq.QuickSynth().synthesise(U, eps).word
        for kw in (
            {"verbose": 1},
            {"verbose": 0},
            {"verbose": np.bool_(True)},
        ):
            assert sq.QuickSynth(**kw).synthesise(U, eps).word == word, kw
    default = sq.synthesise(U, 1e-4, q_cost=10**9)
    assert default.method == "quicksynth"


def test_clifford_t_gateset_gives_t_count_optimal_words() -> None:
    quicksynth = sq.QuickSynth(gateset="clifford_t")
    r = quicksynth.synthesise(sq.rotation(math.pi / 4), 1e-3)
    assert (r.cost, r.n_t, r.n_q, r.certified, r.gateset) == (
        1,
        1,
        0,
        True,
        "clifford_t",
    )
    target = haar_random_unitary(np.random.default_rng(9))
    r = quicksynth.synthesise(target, 0.1)
    assert r.n_q == 0 and not set("Qq") & set(r.word) and r.cost == r.n_t
    assert r.certified and r.distance <= 0.1


def test_the_policy_is_read_from_the_core() -> None:
    """POLICY mirrors the core's constants; stage settings are pinned in C++."""
    P = api.POLICY
    assert (P.epsilon_margin, P.split_below, P.rate, P.fast_rate) == (
        1e-15,
        1e-5,
        0.3,
        3.0,
    )
    assert (P.fine_rate, P.fine_fast_rate) == (0.2, 6.0)
    assert (P.prefix_radius, P.prefix_count) == (6, 600)
    assert (P.rotation_solvable, P.rotation_solvable_fast) == (256, 40)
    assert (P.scan_throughput("sqrt_t"), P.scan_throughput("clifford_t")) == (
        8.3e6,
        2.1e7,
    )
    assert (api.MIN_EPSILON, api.UNITARITY_TOLERANCE) == (1e-13, 1e-6)
    with pytest.raises(AttributeError):
        P.rate = 1.0  # type: ignore[misc]


def _walks(r: sq.SynthesisResult) -> list[tuple]:
    return [
        (e["stage"], e["search"], e["scale"], e.get("parity"), e.get("prefix_scale"))
        for e in r.log
    ]


@pytest.mark.parametrize(
    ("gateset", "eps"),
    [("sqrt_t", 2e-7), ("clifford_t", 1e-6)],
    ids=["sqrt_t", "clifford_t"],
)
def test_fast_walks_a_prefix_of_the_default_walks(gateset: str, eps: float) -> None:
    """fast on whole scales stops the same search sooner: its walks begin those
    of the default on whole scales, its cost is no lower. At 2e-7 on
    Clifford+sqrt(T) the default's last composite scale (about 8e6 points) is
    past the fast rate."""
    default = sq.QuickSynth(gateset=gateset, fine_levels=False)
    fast = sq.QuickSynth(gateset=gateset, fast=True, fine_levels=False)
    fewer = 0
    for seed in (11, 12, 13):
        U = haar_random_unitary(np.random.default_rng(seed))
        d, f = default.synthesise(U, eps), fast.synthesise(U, eps)
        assert f.distance <= eps and f.cost >= d.cost
        assert _walks(d)[: len(_walks(f))] == _walks(f)
        fewer += len(_walks(f)) < len(_walks(d))
    if gateset == "sqrt_t":
        assert fewer > 0


@pytest.mark.parametrize("gateset", ["sqrt_t", "clifford_t"])
def test_fast_walks_the_finer_levels(gateset: str) -> None:
    """fast keeps the finer levels on: on Clifford+sqrt(T) the log holds walks of
    sublattices (sub > 0), as with fine_levels=True, and every answer is within
    epsilon; fine_levels=False walks whole scales; Clifford+T is unchanged.
    """
    assert sq.SynthesisConfig().fine_levels is None
    assert sq.SynthesisConfig(fine_levels=1).fine_levels is True  # type: ignore[arg-type]
    assert sq.SynthesisConfig(fine_levels=0).fine_levels is False  # type: ignore[arg-type]
    fast = sq.QuickSynth(gateset=gateset, fast=True)
    fine = sq.QuickSynth(gateset=gateset, fast=True, fine_levels=True)
    whole = sq.QuickSynth(gateset=gateset, fast=True, fine_levels=False)
    for seed in (21, 22):
        U = haar_random_unitary(np.random.default_rng(seed))
        r, f, w = (q.synthesise(U, 1e-6) for q in (fast, fine, whole))
        assert r.distance <= 1e-6 and w.distance <= 1e-6
        untimed = [
            [{k: v for k, v in e.items() if k != "seconds"} for e in x.log]
            for x in (r, f)
        ]
        assert r.word == f.word and untimed[0] == untimed[1]
        assert any(e.get("sub", 0) > 0 for e in r.log) == (gateset == "sqrt_t")
        assert not any(e.get("sub", 0) > 0 for e in w.log)
        if gateset == "clifford_t":
            assert r.word == w.word


def test_fine_levels_without_fast() -> None:
    """The default walks finer levels, as fine_levels=True does, with LevelSearch
    at POLICY.fine_rate; fine_levels=False walks whole scales at POLICY.rate."""
    U = haar_random_unitary(np.random.default_rng(23))
    d = sq.QuickSynth().synthesise(U, 1e-6)
    r = sq.QuickSynth(fine_levels=True).synthesise(U, 1e-6)
    w = sq.QuickSynth(fine_levels=False).synthesise(U, 1e-6)
    assert d.distance <= 1e-6 and w.distance <= 1e-6
    assert any(e.get("sub", 0) > 0 for e in d.log)
    assert not any(e.get("sub", 0) > 0 for e in w.log)
    assert d.word == r.word
    core = api._core.QuickSynth
    target = np.ascontiguousarray(api.nearest_unitary(U), dtype=np.complex128)
    P = api.POLICY
    for answer, explicit in (
        (d, core("sqrt_t", 1, 3, fine_levels=True, rate=P.fine_rate)),
        (w, core("sqrt_t", 1, 3, fine_levels=False, rate=P.rate)),
    ):
        assert explicit.synthesise(target, 1e-6).word == answer.word


@pytest.mark.parametrize("stopping_rule", ["fixed", "model"])
def test_quicksynth_prepares_ahead_of_a_call(stopping_rule: str) -> None:
    """prepare builds the balls of a call at epsilon and checks epsilon as a call does.

    tests/cpp/test_synthesise.cpp checks which balls prepare builds.
    """
    quick = sq.QuickSynth(stopping_rule=stopping_rule)
    quick.prepare(1e-7)
    r = quick.synthesise(haar_random_unitary(np.random.default_rng(8)), 1e-7)
    assert r.distance <= 1e-7
    with pytest.raises(ValueError, match="epsilon"):
        quick.prepare(0.0)


def test_result_fields_keep_their_types() -> None:
    """All fields have builtin types and serialise to JSON."""
    U = haar_random_unitary(np.random.default_rng(6))
    types = {
        "word": str,
        "cost": int,
        "n_t": int,
        "n_q": int,
        "distance": float,
        "epsilon": float,
        "certified": bool,
        "method": str,
        "seconds": float,
        "t_cost": int,
        "q_cost": int,
        "gateset": str,
        "stages": int,
    }
    for synthesiser in (
        sq.QuickSynth(),
        sq.QuickSynth(gateset="clifford_t"),
    ):
        r = synthesiser.synthesise(U, 1e-2)
        assert {f: type(getattr(r, f)) for f in types} == types
        assert json.loads(json.dumps(r.as_dict()))["certified"] in (True, False)
        assert list(r.as_dict()) == ["word", "gates", *list(types)[1:]]
    fields = [f.name for f in dataclasses.fields(sq.SynthesisResult)]
    assert fields == [*types, "log", "search"]


# --- the z-rotation path (stage 0) ----------------------------------------------------

_H = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
_S = np.diag([1, 1j])


def _comparable(r: sq.SynthesisResult) -> tuple:
    """All fields but the timings."""
    fields = {k: v for k, v in dataclasses.asdict(r).items() if k != "seconds"}
    fields["log"] = [{k: v for k, v in e.items() if k != "seconds"} for e in r.log]
    return tuple(sorted(fields.items(), key=lambda kv: kv[0]))


@pytest.mark.parametrize("axis", ["z", "x", "y"])
def test_clifford_t_rotations_take_the_rotation_path(axis: str) -> None:
    """A rotation about a coordinate axis is a z-rotation up to Cliffords: stage
    0 answers it certified, at no more than the joint walk's T-count and at the
    same one when the joint walk certifies its answer."""
    for theta, eps in ((0.7, 1e-4), (2.1, 1e-6)):
        U = sq.rotation(theta, axis)
        r = sq.synthesise(U, eps, gateset="clifford_t")
        assert (r.stages, r.certified, r.n_q) == (0, True, 0)
        assert magic_counts(r.word) == (r.n_t, 0) and r.cost == r.n_t
        assert r.distance <= eps and chord_diamond_distance(r.unitary, U) <= eps
        assert r.search is not None and r.search["source"] == "rotation"
        assert r.log and all(e["search"] == "rotation" for e in r.log)
        assert all(e["stage"] == 0 for e in r.log)
        counts = ("candidates", "solvable", "undecided", "capped", "pricings")
        assert set(counts) <= set(r.log[0])
        assert sum(e["solvable"] for e in r.log) > 0
        if eps < 1e-4:
            continue  # the joint walk takes seconds on z-rotations below 1e-4
        joint = sq.synthesise(U, eps, gateset="clifford_t", rotation_path=False)
        assert joint.stages >= 1 and r.cost <= joint.cost
        assert joint.cost == r.cost or not joint.certified


def test_exact_clifford_t_rotations_cost_their_t_gates() -> None:
    quick = sq.QuickSynth(gateset="clifford_t")
    for k in range(-4, 5):
        r = quick.synthesise(sq.rotation(k * math.pi / 4), 1e-6)
        assert (r.cost, r.stages, r.certified) == (k % 2, 0, True)
        assert r.distance < 1e-9


def test_near_diagonal_targets_take_the_path_up_to_the_tolerance() -> None:
    """An off-axis tilt within rotation_tolerance * epsilon takes stage 0,
    uncertified; a larger one falls through to the stages unchanged."""
    eps = 1e-4
    C = _H @ _S
    for f, path in ((1e-3, True), (1e-1, False)):
        U = C @ sq.rotation(0.7) @ sq.rotation(2 * math.asin(f * eps), "x")
        r = sq.synthesise(U, eps, gateset="clifford_t")
        assert r.distance <= eps and chord_diamond_distance(r.unitary, U) <= eps
        assert (r.stages == 0) == path
        if path:
            assert not r.certified
        else:
            off = sq.synthesise(U, eps, gateset="clifford_t", rotation_path=False)
            assert _comparable(r) == _comparable(off)


@pytest.mark.parametrize("axis", ["z", "x"])
def test_sqrt_t_rotations_take_the_rotation_path(axis: str) -> None:
    """On Clifford+sqrt(T) stage 0 answers a rotation about a coordinate axis,
    uncertified, after rotation_solvable solvable diagonal entries (fewer with
    fast), its cost the word's."""
    for theta, eps in ((0.7, 1e-4), (2.1, 1e-7)):
        U = sq.rotation(theta, axis)
        default = sq.synthesise(U, eps)
        fast = sq.synthesise(U, eps, fast=True)
        chosen = sq.synthesise(U, eps, rotation_solvable=8)
        budgets = (api.POLICY.rotation_solvable, api.POLICY.rotation_solvable_fast, 8)
        for r, budget in zip((default, fast, chosen), budgets, strict=True):
            assert (r.stages, r.certified) == (0, False)
            assert r.search is not None and r.search["source"] == "rotation"
            assert all(e["stage"] == 0 for e in r.log)
            assert sum(e["solvable"] for e in r.log) >= budget
            assert magic_counts(r.word) == (r.n_t, r.n_q)
            assert r.cost == r.n_t + 3 * r.n_q
            assert r.distance <= eps and chord_diamond_distance(r.unitary, U) <= eps
        assert default.cost <= fast.cost <= chosen.cost


@pytest.mark.parametrize("gateset", ["clifford_t", "sqrt_t"])
def test_the_rotation_path_leaves_other_targets_unchanged(gateset: str) -> None:
    """Targets off the path give the answer, log and summary of
    rotation_path=False."""
    rng = np.random.default_rng(31)
    targets = [haar_random_unitary(rng) for _ in range(3)]
    for U in targets:
        for eps in (1e-3, 1e-4):
            on = sq.synthesise(U, eps, gateset=gateset)
            off = sq.synthesise(U, eps, gateset=gateset, rotation_path=False)
            assert on.stages >= 1
            assert _comparable(on) == _comparable(off)


def test_rotation_settings_are_checked() -> None:
    for tolerance in (-0.1, 1.0):
        with pytest.raises(ValueError, match="rotation_tolerance"):
            sq.SynthesisConfig(rotation_tolerance=tolerance)
    with pytest.raises(ValueError, match="rotation_solvable"):
        sq.SynthesisConfig(rotation_solvable=0)
    assert sq.SynthesisConfig(rotation_path=0).rotation_path is False  # type: ignore[arg-type]
    # no tolerance: only targets diagonal up to rounding take the path
    U = _H @ sq.rotation(0.7) @ sq.rotation(2 * math.asin(1e-7), "x")
    assert (
        sq.synthesise(U, 1e-4, gateset="clifford_t", rotation_tolerance=0).stages >= 1
    )
