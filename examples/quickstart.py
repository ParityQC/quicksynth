#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Synthesise a few targets with QuickSynth.

    python examples/quickstart.py

The first call at a price builds a small cost ball in memory; later calls in
the same process reuse it.
"""

from __future__ import annotations

import numpy as np

import quicksynth as qs


def show(label: str, r: qs.SynthesisResult) -> None:
    kind = "certified optimal" if r.certified else "upper bound"
    print(f"{label}: cost {r.cost} ({kind}), T {r.n_t}, sqrt(T) {r.n_q}")
    print(f"  distance {r.distance:.3e} <= {r.epsilon:.0e}, {r.seconds:.2f} s")
    print(f"  {r.word}")


def main() -> None:
    # R_z(0.7) at the catalysis price (T = 1, sqrt(T) = 3), to 1e-3 and to 1e-7,
    # and on Clifford+T, where the cost is the T-count
    show("R_z(0.7), 1e-3", qs.QuickSynth().synthesise(qs.rotation(0.7), 1e-3))
    show("R_z(0.7), 1e-7", qs.synthesise_rotation(0.7, 1e-7))
    r = qs.synthesise_rotation(0.7, 1e-7, gateset="clifford_t")
    show("R_z(0.7), 1e-7, Clifford+T", r)
    # a general unitary at the gate count (T = sqrt(T) = 1)
    rng = np.random.default_rng(7)
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    U, _ = np.linalg.qr(z)
    show("random U, 1e-5, gate count", qs.synthesise(U, 1e-5, q_cost=1))


if __name__ == "__main__":
    main()
