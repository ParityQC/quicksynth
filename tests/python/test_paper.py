#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Independent checks of the statements behind the paper's lower bound.

Exact arithmetic wherever the statement is exact:

* eps-ball measure in PU(2) against Haar sampling and quadrature over S^3;
* growth recurrence from alternating-word counts against the ball levels;
* dominant roots (6, 2.5943, 2) are simple and strictly dominant;
* orders and intersections of the groups in ``G / phases = C *_K D``;
* existence half of the normal form (proof of Corollary cost): value kept,
  cost never raised, cost equals the descent's;
* the bound's asymptotic slope 3 / log2 rho.
"""

from __future__ import annotations

import heapq
import math
import random

import numpy as np
import pytest
from counting import EXACT_LEVELS, ball_measure

from quicksynth import price, word_unitary

# --------------------------------------------------------------- ball measure


def _distance_to_identity(u: np.ndarray) -> np.ndarray:
    """d(U, I) = sqrt(1 - a0^2) for SU(2) vectors u = (a0, a1, a2, a3)."""
    return np.sqrt(np.maximum(0.0, 1.0 - u[:, 0] ** 2))


@pytest.mark.parametrize("eps", [0.2, 0.4, 0.7])
def test_ball_measure_matches_haar_sampling(eps: float) -> None:
    rng = np.random.default_rng(7)
    u = rng.standard_normal((400_000, 4))
    u /= np.linalg.norm(u, axis=1, keepdims=True)  # Haar on SU(2) = uniform on S^3
    p = float(np.mean(_distance_to_identity(u) < eps))
    se = math.sqrt(p * (1 - p) / len(u))
    assert abs(p - ball_measure(eps)) < 5 * se


@pytest.mark.parametrize("eps", [1e-3, 0.05, 0.3, 0.9])
def test_ball_measure_matches_quadrature(eps: float) -> None:
    from scipy.integrate import quad

    # a0 on S^3 has density (2/pi) sqrt(1-t^2); d < eps iff |a0| > sqrt(1 - eps^2)
    c = math.sqrt(1 - eps * eps)
    p, _ = quad(lambda t: (2 / math.pi) * math.sqrt(1 - t * t), c, 1, epsabs=1e-15)
    assert 2 * p == pytest.approx(ball_measure(eps), rel=1e-9, abs=1e-15)


# ----------------------------------------------------------------- the count


def _word_counts(t: int, q: int, upto: int) -> list[int]:
    """Counts alternating words a1 b1 ... an bn (one per coset) by cost.

    a1 in {1, H, SH}, a_i in {H, SH}, b_i in {sqrtT, T, sqrtT^-1} at q, t, q.
    """
    prices = [q, t, q]
    counts = [0] * (upto + 1)
    counts[0] = 1  # identity coset
    # f[c]: words of cost c ending in a syllable (n >= 1)
    f = [0] * (upto + 1)
    for c in range(1, upto + 1):
        total = 0
        for p in prices:
            if c - p == 0:
                total += 3  # first syllable: three choices of a1
            elif c - p > 0:
                total += 2 * f[c - p]  # later syllables: two choices
        f[c] = total
        counts[c] = total
    return counts


def test_word_count_is_the_enumerated_ball() -> None:
    counts = _word_counts(1, 3, len(EXACT_LEVELS) - 1)
    assert counts == list(EXACT_LEVELS)


def test_word_count_satisfies_the_recurrence() -> None:
    counts = _word_counts(1, 3, 60)
    for n in range(4, 61):  # n above every syllable price
        assert counts[n] == 2 * counts[n - 1] + 4 * counts[n - 3]


@pytest.mark.parametrize("t, q, rho", [(1, 1, 6.0), (1, 3, 2.5943), (2, 3, 2.0)])
def test_dominant_root_is_simple_and_strict(t: int, q: int, rho: float) -> None:
    coef = [0.0] * (q + 1)  # x^q - 2 x^(q-t) - 4
    coef[0] = 1.0
    coef[t] -= 2.0
    coef[q] -= 4.0
    roots = np.roots(coef)
    top = max(roots, key=abs)
    assert abs(top.imag) < 1e-12 and top.real == pytest.approx(rho, abs=5e-5)
    others = [r for r in roots if abs(r - top) > 1e-9]
    assert len(others) == len(roots) - 1  # simple
    assert all(abs(r) < abs(top) - 0.3 for r in others)  # strictly dominant


def test_bound_has_the_stated_slope() -> None:
    """Sum_C max(0, 1 - N(C) mu(eps)) - (3/log2 rho) log2(1/eps) stays bounded."""
    rho = max(abs(r) for r in np.roots([1, -2, 0, -4]))
    counts = _word_counts(1, 3, 400)
    cum = np.cumsum([24 * a for a in counts], dtype=float)
    gaps = []
    for k in range(3, 31):
        eps = 10.0**-k
        mu = ball_measure(eps)
        s = sum(max(0.0, 1 - n * mu) for n in cum)
        gaps.append(s - 3 / math.log2(rho) * math.log2(1 / eps))
    assert max(gaps) - min(gaps) < 2.0  # O(1) remainder, no slope


# --------------------------------------------------- the groups modulo phase
#
# Group elements are words; the product is concatenation.

INVERSE = {"H": "H", "X": "X", "Y": "Y", "Z": "Z", "S": "s", "T": "t", "Q": "q"}
INVERSE.update({v: k for k, v in INVERSE.items()})


def inverse(word: str) -> str:
    return "".join(INVERSE[ch] for ch in reversed(word))


def phase_key(word: str) -> tuple:
    """Key of ``word``'s exact element modulo phase: its least phase multiple.

    Multiplying by zeta shifts coefficients up and negates the wrapped one.
    """
    e = word_unitary(word)
    entries = [e.coefficients[8 * i : 8 * i + 8] for i in range(4)]
    keys = []
    for _ in range(16):
        keys.append((e.sqrt2_exponent, tuple(entries)))
        entries = [(-x[7], *x[:7]) for x in entries]
    return min(keys)


def group_closure(generators: list[str]) -> dict[tuple, str]:
    """Returns one word per element of the generated group modulo phase."""
    seen = {phase_key(""): ""}
    frontier = [""]
    while frontier:
        new = []
        for w in frontier:
            for g in generators:
                k = phase_key(w + g)
                if k not in seen:
                    seen[k] = w + g
                    new.append(w + g)
        frontier = new
    return seen


def named_groups() -> dict[str, dict[tuple, str]]:
    """``D = <Q, X>``, ``K = <S, X>``, ``KT = <S, X, T>``, ``C = <H, S>``."""
    return {
        "D": group_closure(["Q", "X"]),
        "K": group_closure(["S", "X"]),
        "KT": group_closure(["S", "X", "T"]),
        "C": group_closure(["H", "S"]),
    }


GROUPS = named_groups()
D, K, C = GROUPS["D"], GROUPS["K"], GROUPS["C"]


def test_named_groups_have_the_orders_of_the_amalgam() -> None:
    D, K, KT, C = (set(GROUPS[g]) for g in ("D", "K", "KT", "C"))
    assert (len(D), len(K), len(KT), len(C)) == (32, 8, 16, 24)
    assert D & C == K and K < KT < D


# ----------------------------------------- the existence half of the normal form

A_FIRST = {"1": "", "H": "H", "SH": "SH"}
B = {"Q": "Q", "T": "T", "Qi": "q"}  # non-trivial D/K representatives
PRICE = {"T": 1, "Q": 3}
P_CLASS = {"Q": 3, "T": 1, "Qi": 3}  # p(b) at the catalysis price


def _split(x: str, reps: dict[str, str]) -> tuple[str, str]:
    """Splits ``x = r * k``, ``r`` in ``reps``, ``k`` in K; returns (name, K word)."""
    for name, r in reps.items():
        k = phase_key(inverse(r) + x)
        if k in K:
            return name, K[k]
    raise AssertionError("no representative")


def _reduce_to_alternating(word: str) -> tuple[list[tuple[str, str]], str, int]:
    """Reduces ``word`` to a1 b1 ... an bn kappa as in the proof of Corollary cost.

    Returns:
        (syllables, kappa, word cost).
    """
    blocks: list[list] = []  # [kind, value, cost]
    for ch in word:
        kind = "D" if ch in PRICE else "C"
        if blocks and blocks[-1][0] == kind:
            blocks[-1][1] += ch
            blocks[-1][2] += PRICE.get(ch, 0)
        else:
            blocks.append([kind, ch, PRICE.get(ch, 0)])
    changed = True
    while changed:
        changed = False
        for i, (kind, val, _cost) in enumerate(blocks):
            key = phase_key(val)
            if kind == "D" and key in K:  # magic block equal to a Clifford
                blocks[i] = ["C", val, 0]
            elif (
                kind == "C"
                and 0 < i < len(blocks) - 1
                and key in K
                and blocks[i - 1][0] == "D"
                and blocks[i + 1][0] == "D"
            ):  # merge D K D
                left, right = blocks[i - 1], blocks[i + 1]
                blocks[i - 1 : i + 2] = [
                    ["D", left[1] + val + right[1], left[2] + right[2]]
                ]
            else:
                continue
            # coalesce equal-kind neighbours and rescan
            merged: list[list] = []
            for blk in blocks:
                if merged and merged[-1][0] == blk[0]:
                    merged[-1][1] += blk[1]
                    merged[-1][2] += blk[2]
                else:
                    merged.append(list(blk))
            blocks = merged
            changed = True
            break
    if not blocks or blocks[0][0] == "D":
        blocks.insert(0, ["C", "", 0])
    if blocks[-1][0] == "D":
        blocks.append(["C", "", 0])
    k = ""
    syllables = []
    for i in range(0, len(blocks) - 1, 2):
        c, d = blocks[i][1], blocks[i + 1][1]
        reps = A_FIRST if i == 0 else {n: A_FIRST[n] for n in ("H", "SH")}
        a, k = _split(k + c, reps)
        b, k = _split(k + d, B)
        syllables.append((a, b))
    kappa = k + blocks[-1][1]
    assert phase_key(kappa) in C
    return syllables, kappa, sum(b[2] for b in blocks if b[0] == "D")


@pytest.mark.parametrize("seed", range(4))
def test_reduction_keeps_the_value_and_never_raises_the_cost(seed: int) -> None:
    rng = random.Random(seed)
    for _ in range(60):
        word = "".join(rng.choice("HSXTQQ") for _ in range(rng.randint(5, 40)))
        syllables, kappa, word_cost = _reduce_to_alternating(word)
        # same value
        alternating = "".join(A_FIRST[a] + B[b] for a, b in syllables) + kappa
        assert phase_key(alternating) == phase_key(word)
        # cost not raised
        nf_cost = sum(P_CLASS[b] for _, b in syllables)
        assert nf_cost <= word_cost
        # cost equals the descent's
        assert price(word_unitary(word), 1, 3)[0] == nf_cost


# ------------------------------------------------- the class prices p(e)


def _class_prices(k: int, prices: dict[int, float]) -> dict[int, float]:
    """Returns the class prices p(e) by shortest path on Z / 2^(k-1).

    Each letter R_j^{+-1} adds +-2^(k-j) at price prices[j].
    """
    mod = 2 ** (k - 1)
    best = {0: 0.0}
    frontier = [(0.0, 0)]

    while frontier:
        c, r = heapq.heappop(frontier)
        if c > best.get(r, math.inf):
            continue
        for j, cj in prices.items():
            for sign in (1, -1):
                r2 = (r + sign * 2 ** (k - j)) % mod
                if c + cj < best.get(r2, math.inf):
                    best[r2] = c + cj
                    heapq.heappush(frontier, (c + cj, r2))
    return best


@pytest.mark.parametrize("t, q", [(1, 1), (1, 3), (2, 3)])
def test_clifford_sqrt_t_class_prices(t: int, q: int) -> None:
    p = _class_prices(3, {2: t, 3: q})
    assert p[1] == p[3] == q and p[2] == min(t, 2 * q)


def test_class_of_a_monomial_word_is_its_signed_sum() -> None:
    """At k = 3 a word in sqrtT^{+-1}, T, S, X lies in Q^s K, s its signed sum.

    Q -> +-1, T -> +-2, S -> 4; X flips the sign of what follows.
    """
    rng = random.Random(11)
    letters = {"Q": 1, "Qi": -1, "T": 2, "S": 4, "X": 0}
    spelt = {"Q": "Q", "Qi": "q", "T": "T", "S": "S", "X": "X"}
    reps = {e: "Q" * e for e in range(4)}
    for _ in range(200):
        word = [rng.choice(list(letters)) for _ in range(rng.randint(1, 25))]
        s, sign = 0, 1
        for ch in word:
            if ch == "X":
                sign = -sign
            else:
                s += sign * letters[ch]
        M = "".join(spelt[ch] for ch in word)
        assert phase_key(inverse(reps[s % 4]) + M) in K


def test_explicit_bound_constants_and_gap() -> None:
    """Checks the explicit bound: a_n <= B rho^n, B = 3/rho, A = 72/(rho-1).

    log_rho(1/(A mu)) - rho/(rho-1) is 1.1 to 1.25 below the covering sum.
    """
    rho = max(np.roots([1, -2, 0, -4]).real)
    counts = _word_counts(1, 3, 400)
    B = max(counts[n] / rho**n for n in range(4))
    assert math.isclose(B, 3 / rho)
    assert all(counts[n] <= B * rho**n * (1 + 1e-12) for n in range(len(counts)))
    A = 24 * B * rho / (rho - 1)
    assert math.isclose(A, 72 / (rho - 1)) and round(A, 1) == 45.2
    cum = np.cumsum([24 * a for a in counts], dtype=float)
    for k in range(3, 9):
        mu = ball_measure(10.0**-k)
        explicit = math.log(1 / (A * mu), rho) - rho / (rho - 1)
        covering = sum(max(0.0, 1 - n * mu) for n in cum)
        assert 1.1 <= covering - explicit <= 1.25
