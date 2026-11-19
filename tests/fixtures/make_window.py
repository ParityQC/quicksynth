#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Record the model window's reference values as a golden fixture.

    uv run --no-sync python tests/fixtures/make_window.py REFERENCE_DIR

``REFERENCE_DIR`` holds the scripts of the window study: ``census.py`` and
``census.json`` (the closed-form census N(K, c) and its exact table),
``window_rule.py`` (the model of one setting: ``Census``, ``Setting``,
``ball_measure``, the dynamic programme ``Setting.dp``) and ``verify_extra.py``
(the one-step test ``myopic_Q``). It writes ``tests/fixtures/golden_window.json``
with

* ``census``: N(K, c) for K <= 20 and c <= 40 at (1, 1), (1, 3) and (2, 3) on
  Clifford+sqrt(T) and on Clifford+T, exact integers from ``census.json``;
* ``bounds``: the least cost of any gate of scale >= K (``Census.cmin_from``) at
  nine prices, K = 3..31 (Clifford+sqrt(T)) and 3..43 (Clifford+T);
* ``ball_measure`` at 1e-2..1e-9;
* ``settings``: for five settings, the seconds of each scale at the reference
  throughput (None over the point budget), and on a grid of states (K, b) the
  one-step value S_1 - rate sec, the dynamic programme's Q, the savings and the
  chances of something cheaper (next group and every deeper one), and the
  split's saving at the study's phi.

tests/cpp/search_scales/test_stopping_rule.cpp compares the core with it.
Nothing here is timed. The reference counts the groups to scale 34 (46 on
Clifford+T) and the core 16 (28) scales past the first; the grid keeps to costs
where the two agree.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

CENSUS_KEYS = {
    "clifford+sqrt(t)|(1, 1)": (1, 1, 1),
    "clifford+sqrt(t)|(1, 3)": (1, 1, 3),
    "clifford+sqrt(t)|(2, 3)": (1, 2, 3),
    "clifford+t|(1, 1)": (2, 1, 1),
}
CENSUS_K, CENSUS_C = 20, 40
BOUND_PRICES = [(2, 1), (1, 1), (2, 3), (1, 2), (2, 5), (1, 3), (1, 4), (1, 5)]
SETTINGS = [
    ("Q", (1, 1), 1e-3),
    ("Q", (1, 3), 1e-4),
    ("Q", (2, 3), 1e-5),
    ("Q", (1, 3), 1e-6),
    ("T", (1, 1), 1e-4),
]
RATE = 0.3
STATES = 6  # K = s .. s + STATES - 1
SPLIT_PHI = {
    (1, 1): [[0.4838, 0.5162, 0.0], [0.1721, 0.5973, 0.0127]],
    (1, 3): [[0.4349, 0.5310, 0.0], [0.2637, 0.6184, 0.0126]],
    (2, 3): [[0.4813, 0.4475, 0.0], [0.0, 0.1969, 0.0179]],
}
"""split_rule.json's phi_k[j], k = 2 and 3, as the core rounds them."""


def main() -> None:
    ref = Path(sys.argv[1]).resolve()
    sys.path.insert(0, str(ref))
    import numpy as np
    import verify_extra
    import window_rule as wr

    out: dict = {"census": [], "bounds": [], "ball_measure": [], "settings": []}
    table = json.loads((ref / "census.json").read_text())
    for key, (step, t, q) in CENSUS_KEYS.items():
        c = table[key]
        out["census"].append(
            {
                "step": step,
                "t": t,
                "q": q,
                "N": [row[: CENSUS_C + 1] for row in c["N"][: CENSUS_K + 1]],
            }
        )
    for gs, prices in (("Q", BOUND_PRICES), ("T", [(1, 1)])):
        for p in prices:
            cen = wr.Census(gs, p)
            out["bounds"].append(
                {
                    "step": 1 if gs == "Q" else 2,
                    "t": p[0],
                    "q": p[1],
                    "K": list(range(3, wr.KMAX[gs] - 2)),
                    "c_min": [int(cen.cmin_from[K]) for K in range(3, wr.KMAX[gs] - 2)],
                }
            )
    for e in range(2, 10):
        out["ball_measure"].append([10.0**-e, wr.ball_measure(10.0**-e)])
    for gs, p, eps in SETTINGS:
        st = wr.Setting(wr.Census(gs, p), eps)
        my = verify_extra.myopic_Q(st, RATE, wr.POINT_BUDGET, False)
        dp, _ = st.dp(RATE, budget=wr.POINT_BUDGET)
        groups = 14 if gs == "Q" else 20
        seconds = [
            float(st.sec[i]) if st.est[i] <= wr.POINT_BUDGET else None
            for i in range(groups)
        ]
        # costs below the least of any gate past the core's last group
        top = int(st.cen.cmin_from[st.s + (16 if gs == "Q" else 28) + 1])
        costs = sorted({b for b in range(0, top, 3)} | {top - 1})
        grid = []
        for i in range(STATES):
            for b in costs:
                sv, pc = st.saving(i, b)
                sd, pd = st.saving(i, b, deeper=True)
                grid.append(
                    {
                        "K": st.s + i,
                        "b": b,
                        "one_step": float(my[i, b]),
                        "dp": float(dp[i, b]),
                        "saving": sv,
                        "p_cheaper": pc,
                        "saving_deeper": sd,
                        "p_deeper": pd,
                    }
                )
        split = []
        if gs == "Q" and p in SPLIT_PHI:
            for L in (st.s, st.s + 1, st.s + 2):
                for Kc in (L + 1, L + 2, L + 3):
                    for k, phi in zip((2, 3), SPLIT_PHI[p], strict=True):
                        for b in costs[::2]:
                            lam = np.zeros(st.nc + 1)
                            for j in range(3):
                                g = Kc - j - st.s
                                if Kc - j > L and 0 <= g < st.nG:
                                    lam += phi[j] * st.cum[g]
                            bb = int(min(b, st.nc))
                            S = float(np.sum(1 - np.exp(-lam[1 : bb + 1])))
                            split.append(
                                {
                                    "Kc": Kc,
                                    "k": k,
                                    "phi": phi,
                                    "b": b,
                                    "L": L,
                                    "saving": S,
                                }
                            )
        out["settings"].append(
            {
                "step": 1 if gs == "Q" else 2,
                "t": p[0],
                "q": p[1],
                "epsilon": eps,
                "start": st.s,
                "mu": st.mu,
                "rate": RATE,
                "seconds": seconds,
                "grid": grid,
                "split": split,
            }
        )
    path = HERE / "golden_window.json"
    path.write_text(json.dumps(out, separators=(",", ":")) + "\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
