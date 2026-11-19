#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Recompute the counts of golden_relative_norm.json with PARI/GP.

    python tests/fixtures/check_relative_norm.py [FIXTURE]

Needs cypari2 (not a dependency of the package). For each case, xi over
zeta^0..zeta^(d-1) in Q(zeta_2d), PARI's ``idealfactor`` gives the prime ideals
of xi; complex conjugation, zeta -> zeta^(2d-1), pairs them. The number of
solutions of t t* = xi up to roots of unity is the product of (e + 1) over the
conjugate pairs, and 0 if a self-conjugate prime above an odd prime has odd
exponent. Prints the number of cases that differ from the fixture.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import cypari2

FIXTURE = Path(__file__).with_name("golden_relative_norm.json")


def count(pari: cypari2.Pari, field: object, degree: int, xi: list[int]) -> int:
    element = pari(" + ".join(f"({c})*x^{i}" for i, c in enumerate(xi)))
    conjugation = pari(f"lift(Mod(x^{2 * degree - 1}, polcyclo({2 * degree})))")
    factors = pari.idealfactor(field, element)
    primes, exponents = list(factors[0]), [int(e) for e in factors[1]]
    forms = [pari.idealhnf(field, p) for p in primes]
    paired: set[int] = set()
    n = 1
    for i, prime in enumerate(primes):
        if i in paired:
            continue
        image = pari.idealhnf(field, pari.nfgaloisapply(field, conjugation, prime))
        j = forms.index(image)
        paired.update((i, j))
        if j != i:
            n *= exponents[i] + 1
        elif int(prime[0]) != 2 and exponents[i] % 2:
            return 0
    return n


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else FIXTURE
    cases = json.loads(path.read_text())["cases"]
    pari = cypari2.Pari()
    fields = {d: pari.nfinit(pari.polcyclo(2 * d)) for d in {c[0] for c in cases}}
    wrong = [c for c in cases if count(pari, fields[c[0]], c[0], c[3]) != c[2]]
    for case in wrong[:10]:
        print("differs:", case)
    version = ".".join(str(v) for v in pari.version())
    print(f"{len(cases)} cases, {len(wrong)} differ (PARI/GP {version})")


if __name__ == "__main__":
    main()
