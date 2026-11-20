#!/usr/bin/env python3

#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Writes the data of the interactive figures in docs/how-it-works/mechanism.md.

    python tools/showcase_data.py [--out docs/_static/showcase]

``points.bin``: every circuit up to cost 9 on Clifford+T and on Clifford+sqrt(T)
at the price (1,3), found cheapest first and priced by the library, as the
states it prepares from the six stabiliser states. ``mechanism.json``: the
circuit counts per cost, the scan of one Haar-random target scale by scale, a
search the split improves, and one target at seven prices.

The scan's times are measured: run it on an otherwise idle machine, one
thread. The data of docs/how-it-works/numbers.md (``results.json``) comes from
the paper's figure script instead.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import unitary_group

import quicksynth as qs

ROOT = Path(__file__).resolve().parents[1]
MAX_COST = 9
SCALE = 32767  # int16 coordinates of the unit vectors
SCAN_SEED, SPLIT_SEED = 11, 9
PRICES = [(2, 1), (1, 1), (2, 3), (1, 2), (2, 5), (1, 3), (1, 4)]
PRICE_EPSILON = 1e-4
WHOLE_SCALES = {"fine_levels": False, "stopping_rule": "fixed"}

H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
S = np.diag([1, 1j])
T = np.diag([1, np.exp(1j * np.pi / 4)])
Q = np.diag([1, np.exp(1j * np.pi / 8)])


def su2(m: np.ndarray) -> np.ndarray:
    m = np.asarray(m, dtype=complex)
    return m / np.sqrt(np.linalg.det(m))


def quaternion(m: np.ndarray) -> np.ndarray:
    """[[a, -b*], [b, a*]] in SU(2) as (Re a, Im a, Re b, Im b)."""
    return np.stack(
        [m[..., 0, 0].real, m[..., 0, 0].imag, m[..., 1, 0].real, m[..., 1, 0].imag], -1
    )


def cliffords() -> list[tuple[np.ndarray, str]]:
    """The 24 single-qubit Cliffords modulo phase, as words over H and S."""

    def key(m: np.ndarray) -> tuple[float, ...]:
        q = quaternion(su2(m))
        q = q if q[np.argmax(np.abs(q) > 1e-9)] > 0 else -q
        return tuple(np.round(q, 6))

    found = {key(np.eye(2)): (su2(np.eye(2)), "")}
    frontier = [(np.eye(2), "")]
    while frontier:
        nxt = []
        for m, w in frontier:
            for g, letter in ((H, "H"), (S, "S")):
                k = key(m @ g)
                if k not in found:
                    found[k] = (su2(m @ g), w + letter)
                    nxt.append((m @ g, w + letter))
        frontier = nxt
    assert len(found) == 24
    return list(found.values())


CLIFFORDS = cliffords()
CLIFFORD_MATRICES = np.array([m for m, _ in CLIFFORDS])
FUNCTIONAL = np.array(
    [0.7071, 0.3141, 0.5772, 0.2718]
)  # generic: picks one coset member


def coset_keys(ms: np.ndarray) -> np.ndarray:
    """A key per coset m C (C Clifford, up to phase): its member maximising
    |q . FUNCTIONAL|."""
    q = quaternion(np.einsum("nij,cjk->ncik", ms, CLIFFORD_MATRICES))
    f = q @ FUNCTIONAL
    q = q * np.sign(f)[..., None]
    rep = q[np.arange(len(ms)), np.argmax(np.abs(f), axis=1)]
    return np.round(rep * 1e6).astype(np.int64)


def circuits(
    t_cost: int, q_cost: int, sqrt_t: bool
) -> list[tuple[np.ndarray, str, int]]:
    """Every circuit up to MAX_COST modulo a final Clifford, cheapest first, as
    (matrix, word, cost)."""
    gates = [(su2(T), "T", t_cost)] + ([(su2(Q), "Q", q_cost)] if sqrt_t else [])
    seen = {tuple(coset_keys(su2(np.eye(2))[None])[0])}
    layers: dict[int, list[tuple[np.ndarray, str]]] = {0: [(su2(np.eye(2)), "")]}
    out = []
    for cost in range(MAX_COST + 1):
        layer = layers.pop(cost, [])
        out += [(m, w, cost) for m, w in layer]
        if not layer:
            continue
        ms = np.array([m for m, _ in layer])
        for g, letter, price in gates:
            if cost + price > MAX_COST:
                continue
            for cm, cw in CLIFFORDS:
                prod = ms @ cm @ g
                for i, k in enumerate(coset_keys(prod)):
                    if tuple(k) not in seen:
                        seen.add(tuple(k))
                        layers.setdefault(cost + price, []).append(
                            (prod[i], layer[i][1] + cw + letter)
                        )
    return out


def sphere_points(
    gate_set: str, t_cost: int, q_cost: int
) -> tuple[bytes, dict[str, Any]]:
    """The states the circuits prepare from the stabiliser states, each once, at
    its least cost."""
    sqrt_t = gate_set == "clifford_sqrt_t"
    found = circuits(t_cost, q_cost, sqrt_t)
    n_q = []
    for _, word, cost in found:  # every cost from the library's exact pricer
        price, _, nq = qs.price(qs.word_unitary(word), t_cost, q_cost, not sqrt_t)
        assert price == cost, (word, cost, price)
        n_q.append(nq)
    stabiliser = (
        np.array([[1, 0], [0, 1], [1, 1], [1, -1], [1, 1j], [1, -1j]])
        / np.array([1, 1, *[np.sqrt(2)] * 4])[:, None]
    )
    psi = np.einsum("nij,sj->nsi", np.array([m for m, _, _ in found]), stabiliser)
    a, b = psi[..., 0], psi[..., 1]
    xyz = np.stack(
        [2 * (a.conj() * b).real, 2 * (a.conj() * b).imag, abs(a) ** 2 - abs(b) ** 2],
        -1,
    ).reshape(-1, 3)
    cost = np.repeat([c for _, _, c in found], 6)
    nq = np.repeat(n_q, 6)
    order = np.argsort(cost, kind="stable")
    _, first = np.unique(
        np.round(xyz[order] * 1e5).astype(np.int64), axis=0, return_index=True
    )
    keep = order[first]
    keep = keep[np.argsort(cost[keep], kind="stable")]
    rec = np.zeros(
        len(keep),
        dtype=[("x", "<i2"), ("y", "<i2"), ("z", "<i2"), ("cost", "u1"), ("nq", "u1")],
    )
    v = xyz[keep] / np.linalg.norm(xyz[keep], axis=1, keepdims=True)
    v = np.round(v * SCALE).astype(np.int16)
    rec["x"], rec["y"], rec["z"], rec["cost"], rec["nq"] = (
        v[:, 0],
        v[:, 1],
        v[:, 2],
        cost[keep],
        nq[keep],
    )
    per_cost = np.bincount([c for _, _, c in found], minlength=MAX_COST + 1)
    meta = {
        "gate_set": gate_set,
        "price": [t_cost, q_cost],
        "label": "Clifford+T, T-count"
        if not sqrt_t
        else f"Clifford+√T, catalysis ({t_cost},{q_cost})",
        "count": len(rec),
        "max_cost": MAX_COST,
        "circuits": per_cost.tolist(),
    }
    return rec.tobytes(), meta


def scale_bound(gate_set: str, price: list[int], k: int) -> float:
    """No circuit of scale k or deeper costs less (search_scales/scale_bound.hpp)."""
    t, q = price
    return (
        t * (2 * k - 3) if gate_set == "clifford_t" else min(t / 2, q / 3) * (4 * k - 7)
    )


def walks(result: qs.SynthesisResult) -> list[dict[str, Any]]:
    """The first stage's walks in order; a scale's determinant classes merged."""
    out: list[dict[str, Any]] = []
    for e in result.log:
        if e["stage"] != 1:
            continue
        key = (e["search"], e["scale"], e.get("prefix_scale"))
        if not out or out[-1]["key"] != key:
            out.append(
                {
                    "key": key,
                    "search": e["search"],
                    "scale": e["scale"],
                    "prefix_scale": e.get("prefix_scale"),
                    "points": 0,
                    "estimated": 0.0,
                    "seconds": 0.0,
                    "best": None,
                }
            )
        w = out[-1]
        w["points"] += e["points"]
        w["estimated"] += e["estimated_points"]
        w["seconds"] += e["seconds"]
        if e["best_cost"] is not None:
            w["best"] = (
                e["best_cost"] if w["best"] is None else min(w["best"], e["best_cost"])
            )
    for w in out:
        del w["key"]
    return out


def scan_runs() -> list[dict[str, Any]]:
    target = unitary_group.rvs(2, random_state=SCAN_SEED)
    runs = []
    for gate_set, name in [("clifford_t", "clifford_t"), ("sqrt_t", "clifford_sqrt_t")]:
        for eps in [1e-3, 1e-5]:
            r = qs.synthesise(target, eps, gateset=gate_set, **WHOLE_SCALES)
            sc = walks(r)
            assert all(w["search"] == "scan" for w in sc)
            last, prev = sc[-1], sc[-2]
            growth = last["estimated"] / prev["estimated"]
            price = [r.t_cost, r.q_cost]
            nxt = {
                "scale": last["scale"] + 1,
                "points": round(last["points"] * growth),
                "seconds": last["seconds"] * growth,
            }
            runs.append(
                {
                    "gate_set": name,
                    "epsilon": eps,
                    "price": price,
                    "cost": r.cost,
                    "n_t": r.n_t,
                    "n_q": r.n_q,
                    "word": r.word,
                    "certified": r.certified,
                    "seconds": round(r.seconds, 4),
                    "scales": [
                        {
                            "scale": w["scale"],
                            "points": w["points"],
                            "estimated": round(w["estimated"], 3),
                            "seconds": round(w["seconds"], 6),
                            "best": w["best"],
                            "bound": scale_bound(name, price, w["scale"]),
                        }
                        for w in sc
                    ],
                    "next": nxt | {"bound": scale_bound(name, price, nxt["scale"])},
                }
            )
    return runs


def split_runs() -> list[dict[str, Any]]:
    target = unitary_group.rvs(2, random_state=SPLIT_SEED)
    runs = []
    for eps in [1e-7, 1e-8]:
        r = qs.synthesise(target, eps, gateset="sqrt_t", **WHOLE_SCALES)
        ws = walks(r)
        scan = [w for w in ws if w["search"] == "scan"]
        growth = scan[-1]["estimated"] / scan[-2]["estimated"]
        out = []
        for w in ws:
            v = {k: w[k] for k in ("search", "scale", "prefix_scale", "points", "best")}
            if w["search"] == "split":  # what the scan alone would walk at this scale
                v["scan_points"] = round(
                    scan[-1]["points"] * growth ** (w["scale"] - scan[-1]["scale"])
                )
            out.append(v)
        runs.append(
            {
                "epsilon": eps,
                "cost": r.cost,
                "n_t": r.n_t,
                "n_q": r.n_q,
                "word": r.word,
                "source": r.search["source"],
                "price": [r.t_cost, r.q_cost],
                "walks": out,
            }
        )
    return runs


def price_rows() -> list[dict[str, Any]]:
    target = unitary_group.rvs(2, random_state=SCAN_SEED)
    rows = []
    for t, q in PRICES:
        r = qs.synthesise(target, PRICE_EPSILON, t_cost=t, q_cost=q)
        rows.append(
            {
                "price": [t, q],
                "cost": r.cost,
                "n_t": r.n_t,
                "n_q": r.n_q,
                "word": r.word,
            }
        )
    r = qs.synthesise(target, PRICE_EPSILON, gateset="clifford_t")
    rows.append(
        {"price": [1, None], "cost": r.cost, "n_t": r.n_t, "n_q": 0, "word": r.word}
    )
    for a in rows[
        :-1
    ]:  # the figure shows each answer as the cheapest of them all at its price
        t, q = a["price"]
        assert all(b["n_t"] * t + b["n_q"] * q >= a["cost"] for b in rows), a
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--out", type=Path, default=ROOT / "docs" / "_static" / "showcase"
    )
    out = parser.parse_args().out
    commit = subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    blobs, sets, offset = [], [], 0
    for gate_set, t, q in [("clifford_t", 1, 0), ("clifford_sqrt_t", 1, 3)]:
        blob, meta = sphere_points(gate_set, t, q)
        sets.append(meta | {"offset": offset})
        offset += meta["count"]
        blobs.append(blob)
    (out / "points.bin").write_bytes(b"".join(blobs))
    data = {
        "source": f"tools/showcase_data.py at quicksynth {commit}: Haar-random "
        f"targets (scipy unitary_group, seed {SCAN_SEED}; seed {SPLIT_SEED} for the "
        "split); scan and split on whole scales (fine_levels=False), fixed stopping "
        "rule, one thread; prices with the default settings",
        "sphere": {"scale": SCALE, "record_bytes": 8, "sets": sets},
        "scan": scan_runs(),
        "split": {"seed": SPLIT_SEED, "runs": split_runs()},
        "prices": {"seed": SCAN_SEED, "epsilon": PRICE_EPSILON, "rows": price_rows()},
    }
    (out / "mechanism.json").write_text(json.dumps(data, separators=(",", ":")) + "\n")
    print(f"wrote {out / 'points.bin'} ({offset} points) and {out / 'mechanism.json'}")


if __name__ == "__main__":
    main()
