#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Command line ``quicksynth``: Clifford+T and Clifford+sqrt(T) synthesis.

quicksynth THETA EPSILON              # R_z(THETA), THETA may be "pi/8"
quicksynth --unitary '[["0.6+0.8j","0"],["0","1"]]' 1e-4
quicksynth pi/8 1e-3 --price 1 2 --json
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import sys
from pathlib import Path
from typing import get_args

import numpy as np

from quicksynth import __version__
from quicksynth.api import POLICY, QuickSynth, rotation
from quicksynth.config import GateSetName, StoppingRule, SynthesisConfig

_DESCRIPTION = "QuickSynth: a cheap word fast, over Clifford+sqrt(T) or Clifford+T."


def parse_angle(text: str) -> float:
    """Parses a float or an expression in ``pi``, ``sqrt``, ``e`` (``3*pi/16``)."""
    try:
        return float(text)
    except ValueError:
        pass
    allowed = {"pi": math.pi, "e": math.e, "sqrt": math.sqrt}
    if not set(text) <= set("0123456789.+-*/() pisqrte") or "**" in text:
        # no powers: "9**9**9" would hang eval
        raise argparse.ArgumentTypeError(f"cannot parse angle {text!r}")
    try:
        return float(eval(text, {"__builtins__": {}}, allowed))
    except Exception as exc:  # pragma: no cover
        raise argparse.ArgumentTypeError(f"cannot parse angle {text!r}: {exc}") from exc


def parse_unitary(text: str) -> np.ndarray:
    """Parses a JSON 2x2 matrix of numbers, ``"a+bj"`` strings or ``[re, im]``."""
    raw = json.loads(text)

    def entry(x):
        if isinstance(x, str):
            return complex(x.replace(" ", "").replace("i", "j"))
        if isinstance(x, list | tuple):
            return complex(x[0], x[1])
        return complex(x)

    return np.array([[entry(x) for x in row] for row in raw], dtype=complex)


def build_parser() -> argparse.ArgumentParser:
    """The command line of ``quicksynth``."""
    p = argparse.ArgumentParser(
        prog="quicksynth",
        description=_DESCRIPTION,
        epilog=__doc__.split("\n\n")[1] if __doc__ else None,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("theta", nargs="?", type=parse_angle, help="angle of R_z(theta)")
    p.add_argument("epsilon", type=float, help="tolerance in diamond distance")
    p.add_argument("--axis", choices=("x", "y", "z"), default="z")
    p.add_argument("--unitary", type=parse_unitary, help="2x2 unitary as JSON")
    p.add_argument("--unitary-file", type=Path, help=".npy file with a 2x2 unitary")
    p.add_argument(
        "--price",
        nargs=2,
        type=int,
        default=(1, 3),
        metavar=("T", "SQRT_T"),
        help="magic-state price of T and sqrt(T) (default 1 3)",
    )
    p.add_argument(
        "--gateset",
        choices=get_args(GateSetName),
        default="sqrt_t",
        help="clifford_t restricts to Clifford+T (cost = T-count)",
    )
    p.add_argument(
        "--extra-scales",
        type=int,
        default=None,
        help="override the scan's stopping",
    )
    p.add_argument(
        "--point-budget", type=float, default=None, help="override the entry price"
    )
    p.add_argument(
        "--time-limit",
        type=float,
        default=None,
        help="seconds after which no further escalation stage starts; stage 0 "
        "counts against it",
    )
    p.add_argument(
        "--stopping-rule",
        choices=get_args(StoppingRule),
        default="fixed",
        help="how far the search goes past its first answer (default fixed)",
    )
    p.add_argument(
        "--fast",
        action="store_true",
        help="rate rule at POLICY.fast_rate (finer levels: POLICY.fine_fast_rate): "
        "stop past the first answer sooner",
    )
    p.add_argument(
        "--fine-levels",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Clifford+sqrt(T): four finer levels per scale, by saving per second; "
        "fixed stopping rule only (default: on; --no-fine-levels: whole scales)",
    )
    p.add_argument(
        "--no-rotation-path",
        dest="rotation_path",
        action="store_false",
        help="no stage 0: search z-rotations up to Cliffords with the stages only",
    )
    p.add_argument(
        "--rotation-tolerance",
        type=float,
        default=1e-2,
        help="largest off-axis part, in units of epsilon, that stage 0 takes "
        "(default 0.01)",
    )
    p.add_argument(
        "--rotation-solvable",
        type=int,
        default=None,
        metavar="N",
        help="Clifford+sqrt(T): stage 0 stops after whole levels with N diagonal "
        "entries of solvable norm equation (default "
        f"{POLICY.rotation_solvable}, {POLICY.rotation_solvable_fast} with --fast)",
    )
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--json", action="store_true", help="machine-readable output")
    p.add_argument("--qasm", action="store_true", help="print OpenQASM 3")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--version", action="version", version=f"quicksynth {__version__}")
    return p


def _config(args: argparse.Namespace) -> SynthesisConfig:
    # settings other than the price map to same-named flags
    rest = {
        f.name: getattr(args, f.name)
        for f in dataclasses.fields(SynthesisConfig)
        if f.name not in ("t_cost", "q_cost")
    }
    return SynthesisConfig(t_cost=args.price[0], q_cost=args.price[1], **rest)


def main(argv: list[str] | None = None) -> int:
    """The ``quicksynth`` command."""
    parser = build_parser()
    args = parser.parse_args(argv)
    synthesiser = QuickSynth(_config(args))
    if args.unitary is not None:
        target = args.unitary
    elif args.unitary_file is not None:
        target = np.load(args.unitary_file)
    elif args.theta is not None:
        target = rotation(args.theta, args.axis)
    else:
        parser.error("give THETA, --unitary or --unitary-file")
    result = synthesiser.synthesise(target, args.epsilon)
    if args.json:
        print(json.dumps(result.as_dict()))
    elif args.qasm:
        sys.stdout.write(result.to_qasm())
    else:
        status = "certified optimal" if result.certified else "upper bound"
        print(f"word:  {result.word}")
        print(
            f"cost:  {result.cost} T states  "
            f"(n_T = {result.n_t}, n_sqrtT = {result.n_q}; "
            f"T = {result.t_cost}, sqrt(T) = {result.q_cost})  {status}"
        )
        print(f"error: {result.distance:.3e}  (epsilon {args.epsilon:g})")
        print(f"time:  {result.seconds:.2f} s")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
