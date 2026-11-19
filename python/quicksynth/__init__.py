#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Cost-optimal single-qubit synthesis over Clifford+sqrt(T) and Clifford+T.

:class:`QuickSynth` finds a cheap word at any precision (command line:
``quicksynth``); :func:`synthesise` and :func:`synthesise_rotation` are one-call
forms. The exact layer (:class:`CyclotomicMatrix`, :func:`word_unitary`,
:func:`price`, :func:`matrix_word`, :func:`normal_form`, :func:`bloch_image`) maps
between words, exact elements, prices and normal forms.
Everything runs in the native core ``quicksynth._core``.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from quicksynth.api import (
    CyclotomicMatrix,
    QuickSynth,
    SynthesisResult,
    bloch_image,
    matrix_word,
    nearest_unitary,
    normal_form,
    price,
    rotation,
    synthesise,
    synthesise_rotation,
    word_unitary,
)
from quicksynth.config import SynthesisConfig

try:
    __version__ = version("quicksynth")
except PackageNotFoundError:  # source tree, not installed
    __version__ = "0.0.0+unknown"

__all__ = [
    "CyclotomicMatrix",
    "QuickSynth",
    "SynthesisConfig",
    "SynthesisResult",
    "__version__",
    "bloch_image",
    "matrix_word",
    "nearest_unitary",
    "normal_form",
    "price",
    "rotation",
    "synthesise",
    "synthesise_rotation",
    "word_unitary",
]
