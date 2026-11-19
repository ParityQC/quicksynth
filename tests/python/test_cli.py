#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Tests for the quicksynth command."""

from __future__ import annotations

import argparse
import json
import math

import numpy as np
import pytest

import quicksynth as sq
from oracles import haar_random_unitary
from quicksynth.cli import _config, build_parser, main, parse_angle, parse_unitary


def test_angle_parsing() -> None:
    assert parse_angle("0.5") == 0.5
    assert parse_angle("pi/8") == pytest.approx(math.pi / 8)
    assert parse_angle("3*pi/16") == pytest.approx(3 * math.pi / 16)
    with pytest.raises(argparse.ArgumentTypeError):
        parse_angle("__import__('os')")
    with pytest.raises(argparse.ArgumentTypeError):
        parse_angle("9**9**9**9")


def test_unitary_parsing() -> None:
    u = parse_unitary('[["0.6+0.8j", 0], [0, [0, 1]]]')
    assert np.allclose(u, [[0.6 + 0.8j, 0], [0, 1j]])


def test_rotation_json_output(capsys) -> None:
    # on Clifford+T the scan certifies the T
    code = main(["pi/4", "1e-3", "--json", "--gateset", "clifford_t"])
    assert code == 0
    out = json.loads(capsys.readouterr().out)
    assert out["cost"] == 1 and out["n_t"] == 1 and out["certified"] is True


def test_unitary_qasm_output(capsys) -> None:
    main(
        [
            "--unitary",
            '[[1, 0], [0, "0.7071067811865476+0.7071067811865476j"]]',
            "1e-3",
            "--qasm",
        ]
    )
    out = capsys.readouterr().out
    assert out.startswith("OPENQASM 3.0;") and "tdg q[0];" in out


def test_fine_levels_flag(capsys) -> None:
    """--fine-levels and --no-fine-levels set fine_levels; else it is None (on).
    The flag reaches the search: a non-rotation target (seed 22, where the two
    settings give different words) matches the API word for word."""
    parse = build_parser().parse_args
    assert _config(parse(["pi/7", "1e-5", "--fast"])).fine_levels is None
    assert _config(parse(["pi/7", "1e-5", "--fine-levels"])).fine_levels is True
    assert _config(parse(["pi/7", "1e-5", "--no-fine-levels"])).fine_levels is False
    U = haar_random_unitary(np.random.default_rng(22))
    text = json.dumps([[[z.real, z.imag] for z in row] for row in U])
    words = []
    for flag, fine_levels in (([], None), (["--no-fine-levels"], False)):
        assert main(["--unitary", text, "1e-6", "--fast", *flag, "--json"]) == 0
        out = json.loads(capsys.readouterr().out)
        r = sq.QuickSynth(fast=True, fine_levels=fine_levels).synthesise(U, 1e-6)
        assert out["word"] == r.word and out["distance"] <= 1e-6
        words.append(r.word)
    assert words[0] != words[1]


def test_the_method_is_quicksynth(capsys) -> None:
    """The CLI runs QuickSynth and rejects other methods' and removed options."""
    assert main(["pi/4", "1e-3", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["method"] == "quicksynth"
    for argv in (
        ["pi/4", "1e-3", "--method", "optimal"],
        ["pi/4", "1e-3", "--max-radius", "6"],
        ["1e-3", "--prepare"],
        ["pi/4", "1e-3", "--cache-dir", "balls"],
    ):
        with pytest.raises(SystemExit):
            main(argv)
    assert build_parser().prog == "quicksynth"


def test_rotation_path_flag(capsys) -> None:
    """Clifford+T z-rotations take stage 0 unless --no-rotation-path is given."""
    for flag, stages in (([], 0), (["--no-rotation-path"], 1)):
        argv = ["0.7", "1e-4", "--json", "--gateset", "clifford_t", *flag]
        assert main(argv) == 0
        out = json.loads(capsys.readouterr().out)
        assert (out["stages"], out["certified"]) == (stages, True)
    args = build_parser().parse_args(["0.7", "1e-4", "--rotation-tolerance", "0.05"])
    assert _config(args).rotation_tolerance == 0.05 and _config(args).rotation_path
    assert _config(args).rotation_solvable is None
    args = build_parser().parse_args(["0.7", "1e-4", "--rotation-solvable", "64"])
    assert _config(args).rotation_solvable == 64
