#!/usr/bin/env python3

#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Records and checks QuickSynth's answers on a fixed set of calls
(``tests/fixtures/fingerprint.json``).

    python tools/fingerprint.py check [--fast] [-n PROCESSES] [--case ID ...]
    python tools/fingerprint.py record [-n PROCESSES] [--records internal/data]

Refactors that must not change answers have to pass ``check``. A case is one
call: target, epsilon and settings (gate set, price and stopping rule always given).
Targets are stored as the core's input unitary in hexadecimal floats, so the
check does not depend on the host LAPACK. Compared: every result field
(distance bit-exact) and log field except time; exceptions by type and message.
``check`` prints each difference and exits 1 on any.

``record`` takes cases from ``tests/fixtures/golden_synthesise.json``, the paper
records under ``internal/data`` (C1, campaign_post/C1t, window evaluation) and
its own targets; only ``record`` reads the records. Cases marked ``fast`` run
in CI (tests/python/test_fingerprint.py).
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import subprocess
import sys
from collections.abc import Iterator, Sequence
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

import quicksynth as qs
from quicksynth.api import nearest_unitary

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "fingerprint.json"
RECORDS = ROOT / "internal" / "data"

LOG_FIELDS = (
    "stage",
    "scale",
    "parity",
    "prefix_scale",
    "estimated_points",
    "points",
    "best_cost",
)
"""Stored log-entry fields, in order. ``parity`` is None for split entries and
``prefix_scale`` for scan entries, which encodes ``search``; ``seconds`` is
omitted. ``estimated_points`` round-trips exactly as a JSON number."""

ORIGINS = {
    "golden": "golden-<run>-<stopping_rule>: tests/fixtures/golden_synthesise.json",
    "c1": "c1-<t><q>-<epsilon>-<target>: internal/data/campaign/C1.jsonl, "
    "targets_2027.json",
    "c1t": "c1t-<epsilon>-<target>: internal/data/campaign_post/C1t.jsonl, "
    "targets_2027.json",
    "window": "window-<t><q> or ct-<epsilon>-<target>: the model-rule rows of "
    "internal/data/window_eval/rotations.jsonl.gz, targets_2031.json",
    "deploy": "tests/python/test_deploy.py's near_clifford targets",
    "symmetric": "targets with a Clifford symmetry",
    "deep": "deep-<word>-<gateset>-<epsilon>: the word moved by epsilon/3",
    "whole": "whole-<c1 case>: the first C1 case of each setting with fine_levels "
    "false, whole scales",
}
"""Case sources, keyed by the first word of the case id."""

FAST_CASES = 200
"""Size of the fast subset."""

Case = dict[str, Any]


# --- the values ---------------------------------------------------------------


def _hex(x: float | None) -> str | None:
    return None if x is None else float(x).hex()


def encode_matrix(U: np.ndarray) -> list[str]:
    """``U`` row-major as (re, im) pairs of hexadecimal floats."""
    entries = np.asarray(U, dtype=np.complex128).reshape(4)
    return [float(v).hex() for z in entries for v in (z.real, z.imag)]


def decode_matrix(values: Sequence[str]) -> np.ndarray:
    x = [float.fromhex(v) for v in values]
    return np.array([complex(x[i], x[i + 1]) for i in range(0, 8, 2)]).reshape(2, 2)


def _search(search: dict[str, Any] | None) -> dict[str, Any] | None:
    if search is None:
        return None
    out = dict(search)
    out["lower_bound"] = _hex(search["lower_bound"])
    out["estimated_points"] = _hex(search["estimated_points"])
    return out


def _log_entry(e: dict[str, Any]) -> list[Any]:
    return [e.get(field) for field in LOG_FIELDS]


def _config(settings: dict[str, Any]) -> dict[str, Any]:
    # the fixture stores the stopping rule under its recorded key "window"
    return {"stopping_rule" if k == "window" else k: v for k, v in settings.items()}


def run_case(case: Case) -> dict[str, Any]:
    """Runs a case; returns ``{"result": ...}`` or ``{"raises": ...}``."""
    core = qs.QuickSynth(**_config(case["settings"]))._core
    U = np.ascontiguousarray(decode_matrix(case["target"]), dtype=np.complex128)
    try:
        r = core.synthesise(U, float.fromhex(case["epsilon"]))
    except Exception as exc:
        return {"raises": {"type": type(exc).__name__, "message": str(exc)}}
    return {
        "result": {
            "word": r.word,
            "cost": r.cost,
            "n_t": r.n_t,
            "n_q": r.n_q,
            "distance": _hex(r.distance),
            "certified": r.certified,
            "stages": r.stages,
            "search": _search(r.search),
            "log": [_log_entry(e) for e in r.log],
        }
    }


def outcome(case: Case) -> dict[str, Any]:
    """Recorded outcome of a case, in :func:`run_case` form."""
    return {k: case[k] for k in ("result", "raises") if k in case}


def diff(want: Any, got: Any, path: str = "") -> list[str]:
    """Lists differing fields as ``path: want -> got``."""
    if isinstance(want, dict) and isinstance(got, dict):
        out = []
        for key in sorted(want.keys() | got.keys()):
            sub = f"{path}.{key}" if path else key
            if key not in got:
                out.append(f"{sub}: {want[key]!r} -> missing")
            elif key not in want:
                out.append(f"{sub}: missing -> {got[key]!r}")
            else:
                out += diff(want[key], got[key], sub)
        return out
    if isinstance(want, list) and isinstance(got, list):
        if path.endswith(".log") or path == "log":
            return _log_diff(want, got, path)
        if len(want) != len(got) or any(
            not isinstance(w, list | dict) for w in want + got
        ):
            return [] if want == got else [f"{path}: {want!r} -> {got!r}"]
        out = []
        for i, (w, g) in enumerate(zip(want, got, strict=True)):
            out += diff(w, g, f"{path}[{i}]")
        return out
    if want == got and type(want) is type(got):
        return []
    return [f"{path}: {want!r} -> {got!r}"]


def _log_diff(want: list, got: list, path: str) -> list[str]:
    out = []
    for i in range(max(len(want), len(got))):
        w = want[i] if i < len(want) else None
        g = got[i] if i < len(got) else None
        if w is None or g is None:
            out.append(f"{path}[{i}]: {w!r} -> {g!r}")
            continue
        for field, a, b in zip(LOG_FIELDS, w, g, strict=True):
            if a != b or type(a) is not type(b):
                out.append(f"{path}[{i}].{field}: {a!r} -> {b!r}")
    return out


# --- the build ----------------------------------------------------------------


def _commit() -> str:
    def git(*args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()

    head = git("rev-parse", "HEAD")
    dirty = git("status", "--porcelain", "--untracked-files=no", "--", "src", "python")
    return head + ("-dirty" if dirty else "")


# --- the cases ----------------------------------------------------------------


def _settings(
    gateset: str, t_cost: int, q_cost: int, stopping_rule: str, workers: int = 1
) -> dict[str, Any]:
    s: dict[str, Any] = {
        "gateset": gateset,
        "t_cost": t_cost,
        "q_cost": q_cost,
        "window": stopping_rule,  # recorded key; see _config
    }
    if workers != 1:
        s["workers"] = workers
    return s


def _case(
    case_id: str,
    target: np.ndarray,
    epsilon: float,
    settings: dict[str, Any],
    expect: dict[str, Any] | None = None,
) -> Case:
    """Builds a case. ``expect`` (recorded word, cost, stages) is checked by
    ``record`` and not stored."""
    case: Case = {
        "id": case_id,
        "target": encode_matrix(nearest_unitary(target)),
        "epsilon": _hex(epsilon),
        "settings": settings,
    }
    if expect is not None:
        case["expect"] = expect
    return case


def _decode_target(t: dict[str, Any]) -> np.ndarray:
    return np.array(t["re"]) + 1j * np.array(t["im"])


def golden_cases() -> list[Case]:
    """QuickSynth runs of golden_synthesise.json (``"scan"``; ``"auto"`` below
    1e-3) under the recorded fixed stopping rule and, except test_golden.py's three
    slow runs, the model stopping rule."""
    runs = json.loads((ROOT / "tests/fixtures/golden_synthesise.json").read_text())
    out = []
    for i, run in enumerate(runs["runs"]):
        config = run["config"]
        quick = config["method"] == "scan" or (
            config["method"] == "auto" and run["epsilon"] < 1e-3
        )
        if not quick or "raises" in run:
            continue
        for stopping_rule in ("fixed",) if run.get("slow") else ("fixed", "model"):
            settings = _settings(
                config.get("gateset", "sqrt_t"),
                config.get("t_cost", 1),
                config.get("q_cost", 3),
                stopping_rule,
                config.get("workers", 1),
            )
            expect = None
            if stopping_rule == "fixed":
                r = run["result"]
                expect = {"word": r["word"], "cost": r["cost"], "stages": r["stages"]}
            out.append(
                _case(
                    f"golden-{i:03d}-{stopping_rule}",
                    _decode_target(run["target"]),
                    run["epsilon"],
                    settings,
                    expect,
                )
            )
    return out


CAMPAIGN_STAGE2 = {("clifford+t", 5)}
"""(gate set, target) of included campaign stage-2 calls: one of C1t's 13
(Clifford+T, 1e-8, ~35 CPU-s each), since the window evaluation's 17 cover that
setting. C1's one stage-2 call (target 78, 173 CPU-s) is excluded."""

CAMPAIGN_PER_SETTING = 7
"""Stage-1 calls per campaign setting, from those recorded at <= 2 CPU-s. None
qualify at (1,1) and (1,3) at 1e-8 (~18 s each); the quickest (1,3) call is
taken instead."""


def campaign_cases(records: Path) -> list[Case]:
    """Campaign C1 and its Clifford+T re-run C1t under the recorded fixed stopping rule:
    the ``CAMPAIGN_STAGE2`` calls and ``CAMPAIGN_PER_SETTING`` stage-1 calls per
    setting (24 settings)."""
    targets = json.loads((records / "campaign/targets_2027.json").read_text())
    rows = []
    for path in ("campaign/C1.jsonl", "campaign_post/C1t.jsonl"):
        rows += [
            (path, json.loads(line))
            for line in (records / path).read_text().splitlines()
        ]
    rng = np.random.default_rng(2026)
    settings_of: dict[tuple, list] = {}
    for path, r in rows:
        key = (r["gate_set"], tuple(r["price"]), r["eps"])
        settings_of.setdefault(key, []).append((path, r))
    chosen = []
    for key in sorted(settings_of):
        group = settings_of[key]
        chosen += [
            (path, r)
            for path, r in group
            if r["stages"] > 1 and (r["gate_set"], r["target"]) in CAMPAIGN_STAGE2
        ]
        quick = [(p, r) for p, r in group if r["stages"] == 1 and r["cpu_seconds"] <= 2]
        if quick:
            n = min(CAMPAIGN_PER_SETTING, len(quick))
            picks = rng.choice(len(quick), size=n, replace=False)
            chosen += [quick[i] for i in sorted(picks.tolist())]
        elif key[1] == (1, 3):
            chosen.append(min(group, key=lambda pr: pr[1]["cpu_seconds"]))
    out = []
    for _, r in chosen:
        clifford_t = r["gate_set"] == "clifford+t"
        t, q = (1, 3) if clifford_t else r["price"]
        name = "c1t" if clifford_t else f"c1-{r['price'][0]}{r['price'][1]}"
        out.append(
            _case(
                f"{name}-{r['eps']:.0e}-{r['target']:02d}",
                _decode_target(targets[r["target"]]),
                r["eps"],
                _settings("clifford_t" if clifford_t else "sqrt_t", t, q, "fixed"),
                {"word": r["word"], "cost": r["cost"], "stages": r["stages"]},
            )
        )
    return out


WINDOW_EVAL_PER_SETTING = 3
"""Stage-1 model-rule rows per window_eval setting, drawn from those at <= 1
CPU-s; all stage-2 rows are taken too."""


def window_eval_cases(records: Path) -> list[Case]:
    """Window-evaluation model-rule rows: all 18 stage-2 rows and
    ``WINDOW_EVAL_PER_SETTING`` per setting (24 settings)."""
    targets = json.loads((records / "window_eval/targets_2031.json").read_text())
    with gzip.open(records / "window_eval/rotations.jsonl.gz", "rt") as f:
        rows = [json.loads(line) for line in f]
    rows = [r for r in rows if r["window"] == "model"]
    rng = np.random.default_rng(2031)
    groups: dict[tuple, list] = {}
    for r in rows:
        groups.setdefault(
            (r["gateset"], r["t_cost"], r["q_cost"], r["eps"]), []
        ).append(r)
    chosen = []
    for key in sorted(groups):
        group = groups[key]
        chosen += [r for r in group if r["stages"] > 1]
        quick = [r for r in group if r["stages"] == 1 and r["cpu"] <= 1]
        n = min(WINDOW_EVAL_PER_SETTING, len(quick))
        picks = rng.choice(len(quick), size=n, replace=False)
        chosen += [quick[i] for i in sorted(picks.tolist())]
    out = []
    for r in chosen:
        name = "ct" if r["gateset"] == "clifford_t" else f"{r['t_cost']}{r['q_cost']}"
        out.append(
            _case(
                f"window-{name}-{r['eps']:.0e}-{r['target']:02d}",
                _decode_target(targets[r["target"]]),
                r["eps"],
                _settings(r["gateset"], r["t_cost"], r["q_cost"], "model"),
                {"word": r["word"], "cost": r["cost"], "stages": r["stages"]},
            )
        )
    return out


H = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
S = np.diag([1, 1j])


def near_clifford(C: np.ndarray, f: float, eps: float) -> np.ndarray:
    """``C`` times a z-rotation that puts the target ``f * eps`` from ``C``."""
    return C @ qs.rotation(2 * math.asin(f * eps), "z")


def own_cases() -> list[Case]:
    """Own targets, model stopping rule: near-Clifford (test_deploy.py; stages 3 and 2,
    and one exactly epsilon away), Clifford-symmetric, and exact words moved by
    epsilon/3 at 1e-10 to 1e-13."""
    out = []
    model = "model"
    for name, C, f, eps in (
        ("identity", np.eye(2), 1.01, 1e-6),
        ("h", H, 1.01, 1e-6),
        ("identity-at-epsilon", np.eye(2), 1.0, 1e-3),
    ):
        out.append(
            _case(
                f"deploy-{name}",
                near_clifford(C, f, eps),
                eps,
                _settings("sqrt_t", 1, 3, model),
            )
        )
    haar = np.random.default_rng(7)
    z = haar.normal(size=(2, 2)) + 1j * haar.normal(size=(2, 2))
    W, _ = np.linalg.qr(z)
    symmetric = {
        "rz0.1": qs.rotation(0.1, "z"),
        "rz1": qs.rotation(1.0, "z"),
        "rzpi7": qs.rotation(math.pi / 7, "z"),
        "rx0.7": qs.rotation(0.7, "x"),
        "ry1.3": qs.rotation(1.3, "y"),
        "srx0.4": S @ qs.rotation(0.4, "x") @ S.conj().T,
        "wwt": W @ W.T,
    }
    for name, U in symmetric.items():
        for gateset in ("sqrt_t", "clifford_t"):
            for eps in (1e-3, 1e-5):
                out.append(
                    _case(
                        f"symmetric-{name}-{gateset}-{eps:.0e}",
                        U,
                        eps,
                        _settings(gateset, 1, 3, model),
                    )
                )
    for gateset in ("sqrt_t", "clifford_t"):
        out.append(
            _case(
                f"symmetric-rz0.1-{gateset}-1e-07",
                symmetric["rz0.1"],
                1e-7,
                _settings(gateset, 1, 3, model),
            )
        )
    for t, q in ((1, 1), (2, 3)):
        out.append(
            _case(
                f"symmetric-rz1-sqrt_t-{t}{q}-1e-04",
                symmetric["rz1"],
                1e-4,
                _settings("sqrt_t", t, q, model),
            )
        )
    words = {"HTSHTHT": True, "HQSHTQHSQ": False, "SHTHTHTSHTHTHTSHT": True}
    for (word, clifford_t), eps, axis in zip(
        [(w, c) for w, c in words.items() for _ in range(2)],
        (1e-10, 1e-13, 1e-11, 1e-12, 1e-10, 1e-13),
        "xzzxxz",
        strict=True,
    ):
        E = qs.word_unitary(word).to_numpy()
        for gateset in ("sqrt_t", "clifford_t") if clifford_t else ("sqrt_t",):
            out.append(
                _case(
                    f"deep-{word}-{gateset}-{eps:.0e}",
                    E @ qs.rotation(eps / 3, axis),
                    eps,
                    _settings(gateset, 1, 3, model),
                )
            )
    return out


def whole_scale_cases(campaign: list[Case]) -> list[Case]:
    """The first C1 case of each setting again with ``fine_levels=False``: whole
    scales, the default before the finer levels, which the C1 records predate."""
    out = []
    seen = set()
    for case in campaign:
        setting = case["id"].rsplit("-", 1)[0]
        if not case["id"].startswith("c1-") or setting in seen:
            continue
        seen.add(setting)
        settings = {**case["settings"], "fine_levels": False}
        out.append({**case, "id": f"whole-{case['id']}", "settings": settings})
    return out


def all_cases(records: Path) -> list[Case]:
    campaign = campaign_cases(records)
    cases = (
        golden_cases()
        + campaign
        + window_eval_cases(records)
        + own_cases()
        + whole_scale_cases(campaign)
    )
    ids = [c["id"] for c in cases]
    duplicates = {i for i in ids if ids.count(i) > 1}
    if duplicates:
        raise RuntimeError(f"duplicate case ids: {sorted(duplicates)}")
    return cases


def _points(case: Case) -> int:
    """Leaves walked by a case; ranks the fast subset."""
    if "result" not in case:
        return 1 << 62
    return sum(entry[LOG_FIELDS.index("points")] for entry in case["result"]["log"])


def mark_fast(cases: list[Case]) -> None:
    """Marks the ``FAST_CASES`` stage-1, single-thread cases with fewest leaves.

    Half the subset is split equally among kinds (first word of the id), so
    every kind is represented."""
    eligible = [
        c
        for c in cases
        if "result" in c
        and c["result"]["stages"] == 1
        and c["settings"].get("workers", 1) == 1
    ]
    eligible.sort(key=lambda c: (_points(c), c["id"]))
    groups: dict[str, list[Case]] = {}
    for c in eligible:
        groups.setdefault(c["id"].split("-")[0], []).append(c)
    fast: list[Case] = []
    share = FAST_CASES // (2 * len(groups))
    for group in groups.values():
        fast += group[:share]
    chosen = {c["id"] for c in fast}
    for c in eligible:
        if len(chosen) >= FAST_CASES:
            break
        chosen.add(c["id"])
    for c in cases:
        c["fast"] = c["id"] in chosen


# --- running ------------------------------------------------------------------


def run_all(cases: Sequence[Case], processes: int) -> Iterator[dict[str, Any]]:
    """Yields each case's outcome in order, over ``processes`` processes."""
    if processes <= 1:
        yield from map(run_case, cases)
        return
    with ProcessPoolExecutor(processes) as pool:
        yield from pool.map(run_case, cases, chunksize=1)


def _dump(fixture: dict[str, Any]) -> str:
    """Serialises the fixture, one case per line."""
    lines = [json.dumps(c, separators=(",", ":")) for c in fixture["cases"]]
    meta = json.dumps(fixture["meta"], indent=2)
    return '{\n"meta": ' + meta + ',\n"cases": [\n' + ",\n".join(lines) + "\n]\n}\n"


def record(args: argparse.Namespace) -> int:
    cases = all_cases(Path(args.records))
    print(f"recording {len(cases)} cases on {args.processes} processes", flush=True)
    mismatches = 0
    for case, got in zip(cases, run_all(cases, args.processes), strict=True):
        case.update(got)
        expect = case.pop("expect", None)
        if expect is None:
            continue
        if "raises" in got:
            print(f"{case['id']}: raised {got['raises']}, recorded {expect}")
            mismatches += 1
            continue
        mine = {k: got["result"][k] for k in expect}
        if mine != expect:
            print(f"{case['id']}: answers {mine}, recorded {expect}")
            mismatches += 1
    mark_fast(cases)
    order = (
        "id",
        "fast",
        "target",
        "epsilon",
        "settings",
        "result",
        "raises",
    )
    cases = [{k: c[k] for k in order if k in c} for c in cases]
    meta: dict[str, Any] = {
        "commit": _commit(),
        "cases": len(cases),
        "fast_cases": sum(c["fast"] for c in cases),
        "log_fields": list(LOG_FIELDS),
        "origins": ORIGINS,
        "recorded_with": "python tools/fingerprint.py record",
    }
    Path(args.out).write_text(_dump({"meta": meta, "cases": cases}))
    print(
        f"wrote {args.out}: {len(cases)} cases, {meta['fast_cases']} fast; "
        f"{mismatches} differ from their records"
    )
    return 0


def load(path: Path = FIXTURE) -> dict[str, Any]:
    return json.loads(path.read_text())


def check(args: argparse.Namespace) -> int:
    fixture = load(Path(args.fixture))
    cases = fixture["cases"]
    if args.case:
        wanted = set(args.case)
        cases = [c for c in cases if c["id"] in wanted]
        missing = wanted - {c["id"] for c in cases}
        if missing:
            print(f"no such case: {sorted(missing)}")
            return 2
    elif args.fast:
        cases = [c for c in cases if c["fast"]]
    print(
        f"fingerprint of {fixture['meta']['commit']}; checking {len(cases)} cases "
        f"on {args.processes} processes",
        flush=True,
    )
    differing = 0
    for case, got in zip(cases, run_all(cases, args.processes), strict=True):
        lines = diff(outcome(case), got)
        if lines:
            differing += 1
            print(f"{case['id']}:")
            for line in lines:
                print(f"  {line}")
            sys.stdout.flush()
    print(f"{differing} of {len(cases)} cases differ")
    return 1 if differing else 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check", help="compare this build with the fixture")
    p.add_argument("--fixture", default=str(FIXTURE))
    p.add_argument("--fast", action="store_true", help="the fast subset only")
    p.add_argument("--case", nargs="+", help="these case ids only")
    p.add_argument("-n", "--processes", type=int, default=1)
    p.set_defaults(func=check)
    p = sub.add_parser("record", help="record the fixture from this build")
    p.add_argument("--records", default=str(RECORDS))
    p.add_argument("--out", default=str(FIXTURE))
    p.add_argument("-n", "--processes", type=int, default=1)
    p.set_defaults(func=record)
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
