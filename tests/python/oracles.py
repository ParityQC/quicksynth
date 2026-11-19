#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Reference helpers for the tests: diamond distance, magic counts, Clifford
coset keys, targets, and exact elements and Bloch images in Python ints.

For ``SU(2)`` vectors ``a, b`` the distance is ``sqrt(1 - (a . b)**2)``, which
cancels below about ``1e-8``. The chord form is accurate to machine precision
at all distances.
"""

from __future__ import annotations

import math

import numpy as np

from quicksynth import CyclotomicMatrix, bloch_image


def su2_vector(V: np.ndarray) -> np.ndarray:
    """``(Re a, Im a, Re b, Im b)`` of ``V / sqrt(det V)``; defined up to sign."""
    V = V / np.sqrt(np.linalg.det(V))
    return np.array([V[0, 0].real, V[0, 0].imag, V[1, 0].real, V[1, 0].imag])


def chord_diamond_distance(u: np.ndarray, v: np.ndarray) -> float:
    """``PU(2)`` diamond distance ``sqrt(s - s**2 / 4)``.

    ``s = min(|a - b|**2, |a + b|**2)``, the shorter chord.
    """
    a, b = su2_vector(u), su2_vector(v)
    s = min(float(np.sum((a - b) ** 2)), float(np.sum((a + b) ** 2)))
    return math.sqrt(max(s - s * s / 4, 0.0))


def diamond_distance(u: np.ndarray, v: np.ndarray) -> float:
    """Trace form ``sqrt(1 - |tr(u v^dagger)|^2 / 4)``; inaccurate below ``1e-8``.

    Independent check of :func:`chord_diamond_distance` at coarse precisions.
    """
    overlap = abs(np.trace(u @ v.conj().T)) / 2
    return math.sqrt(max(0.0, 1.0 - min(overlap, 1.0) ** 2))


def haar_random_unitary(rng: np.random.Generator) -> np.ndarray:
    """Returns a Haar-random element of ``U(2)``."""
    z = (rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))) / math.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def decode_target(t: dict[str, list]) -> np.ndarray:
    """Decodes a stored ``{"re", "im"}`` target."""
    return np.array(t["re"]) + 1j * np.array(t["im"])


def magic_counts(word: str) -> tuple[int, int]:
    """``(n_T, n_sqrtT)``; inverses count as the gate."""
    return sum(word.count(c) for c in "Tt"), sum(word.count(c) for c in "Qq")


def coset_key(element: CyclotomicMatrix) -> tuple:
    """Exact key of ``element``'s Clifford coset modulo phase.

    The Bloch image with canonical column signs and order.
    """
    digits, exponent = bloch_image(element)
    columns = []
    for c in range(3):
        column = [digits[12 * r + 4 * c : 12 * r + 4 * c + 4] for r in range(3)]
        pivot = next((x for entry in column for x in entry if x != 0), 1)
        if pivot < 0:
            column = [tuple(-x for x in entry) for entry in column]
        columns.append(tuple(column))
    return (exponent, tuple(sorted(columns)))


# --- exact elements and Bloch images in Python ints ------------------------------

Zeta = list[int]
"""Element of Z[zeta16] as the coefficients of zeta**0..7, zeta**8 = -1."""


def _zeta_mul(a: Zeta, b: Zeta) -> Zeta:
    out = [0] * 8
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j < 8:
                    out[i + j] += x * y
                else:
                    out[i + j - 8] -= x * y
    return out


def _zeta_add(a: Zeta, b: Zeta, b_scale: int = 1) -> Zeta:
    """``a + b_scale * b``."""
    return [x + b_scale * y for x, y in zip(a, b, strict=True)]


def _monomial(sign: int, power: int) -> Zeta:
    """``sign * zeta**power``."""
    power %= 16
    out = [0] * 8
    out[power % 8] = sign if power < 8 else -sign
    return out


# Each letter's nonzero entries, (row, column): (sign, power) for sign * zeta**power.
_LETTERS: dict[str, dict[tuple[int, int], tuple[int, int]]] = {
    "H": {(0, 0): (1, 0), (0, 1): (1, 0), (1, 0): (1, 0), (1, 1): (-1, 0)},
    "X": {(0, 1): (1, 0), (1, 0): (1, 0)},
    "Y": {(0, 1): (-1, 4), (1, 0): (1, 4)},
    "Z": {(0, 0): (1, 0), (1, 1): (-1, 0)},
    "S": {(0, 0): (1, 0), (1, 1): (1, 4)},
    "s": {(0, 0): (1, 0), (1, 1): (1, 12)},
    "T": {(0, 0): (1, 0), (1, 1): (1, 2)},
    "t": {(0, 0): (1, 0), (1, 1): (1, 14)},
    "Q": {(0, 0): (1, 0), (1, 1): (1, 1)},
    "q": {(0, 0): (1, 0), (1, 1): (1, 15)},
}


def exact_word_matrix(word: str) -> CyclotomicMatrix:
    """Exact element of ``word`` in Python ints.

    Multiplies the letters' matrices (H without its 1/sqrt2) and divides by
    sqrt2 only at the end, unlike the core, which reduces after each letter.
    The reduced element is unique.
    """
    m = [[_monomial(1, 0), [0] * 8], [[0] * 8, _monomial(1, 0)]]
    k = 0
    for letter in word:
        product = [[[0] * 8 for _ in range(2)] for _ in range(2)]
        for (r, j), (sign, power) in _LETTERS[letter].items():
            factor = _monomial(sign, power)
            for i in range(2):
                term = _zeta_mul(m[i][r], factor)
                product[i][j] = _zeta_add(product[i][j], term)
        m = product
        k += letter == "H"
    sqrt2 = [0, 0, 1, 0, 0, 0, -1, 0]  # zeta**2 - zeta**6
    while k > 0:
        doubled = [[_zeta_mul(entry, sqrt2) for entry in row] for row in m]
        if any(c % 2 for row in doubled for entry in row for c in entry):
            break
        m = [[[c // 2 for c in entry] for entry in row] for row in doubled]
        k -= 1
    return CyclotomicMatrix(k, tuple(c for row in m for entry in row for c in entry))


Lambda = list[int]
"""Element of Z[lambda], lambda = zeta + zeta**-1, as the coefficients of
lambda**0..3, lambda**4 = 4 lambda**2 - 2."""


def _lambda_mul(a: Lambda, b: Lambda) -> Lambda:
    full = [0] * 7
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            full[i + j] += x * y
    for p in (6, 5, 4):  # lambda**p = 4 lambda**(p-2) - 2 lambda**(p-4)
        full[p - 2] += 4 * full[p]
        full[p - 4] -= 2 * full[p]
    return full[:4]


def exact_bloch_image(element: CyclotomicMatrix) -> tuple[tuple[int, ...], int]:
    """:func:`quicksynth.bloch_image` in Python ints: the 36 digits of
    ``N u**(k+1)`` (row, column, power of lambda) and the least exponent.

    ``N[i][j] = tr(sigma_i M sigma_j M^dagger)``, common powers of lambda
    stripped, ``u = 2 lambda**2 - 1``.
    """
    c = element.coefficients
    m = [
        [list(c[16 * i + 8 * j : 16 * i + 8 * j + 8]) for j in range(2)]
        for i in range(2)
    ]

    def conj(x: Zeta) -> Zeta:
        return [x[0]] + [-x[8 - t] for t in range(1, 8)]

    def matmul(a: list[list[Zeta]], b: list[list[Zeta]]) -> list[list[Zeta]]:
        return [
            [
                _zeta_add(_zeta_mul(a[i][0], b[0][j]), _zeta_mul(a[i][1], b[1][j]))
                for j in range(2)
            ]
            for i in range(2)
        ]

    zero, one = [0] * 8, _monomial(1, 0)
    pauli = (
        [[zero, one], [one, zero]],
        [[zero, _monomial(-1, 4)], [_monomial(1, 4), zero]],
        [[one, zero], [zero, _monomial(-1, 0)]],
    )
    dagger = [[conj(m[j][i]) for j in range(2)] for i in range(2)]
    lambda_powers = [[1, 0, 0, 0, 0, 0, 0, 0]]
    for _ in range(3):
        lambda_powers.append(_zeta_mul(lambda_powers[-1], [0, 1, 0, 0, 0, 0, 0, -1]))
    cells: list[Lambda] = []
    for i in range(3):
        for j in range(3):
            x = matmul(matmul(matmul(pauli[i], m), pauli[j]), dagger)
            trace = _zeta_add(x[0][0], x[1][1])
            # peel lambda**3 .. lambda**0 off the trace
            digits = [0] * 4
            for p in (3, 2, 1, 0):
                digits[p] = trace[p]
                trace = _zeta_add(trace, lambda_powers[p], -digits[p])
            assert trace == [0] * 8, "a trace that is not real"
            cells.append(digits)
    # 1 / lambda = (4 lambda - lambda**3) / 2
    stripped = 0
    while any(any(cell) for cell in cells):
        doubled = [_lambda_mul(cell, [0, 4, 0, -1]) for cell in cells]
        if any(d % 2 for cell in doubled for d in cell):
            break
        cells = [[d // 2 for d in cell] for cell in doubled]
        stripped += 1
    k = element.sqrt2_exponent
    for _ in range(k + 1):
        cells = [_lambda_mul(cell, [-1, 0, 2, 0]) for cell in cells]
    return tuple(d for cell in cells for d in cell), 4 * (k + 1) - stripped
