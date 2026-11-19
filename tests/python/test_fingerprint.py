#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""QuickSynth against ``tests/fixtures/fingerprint.json`` (``tools/fingerprint.py``).

Every recorded field must match. Non-fast cases are marked ``slow``. Linux
x86-64 only, where the fixture was recorded.
"""

from __future__ import annotations

import platform
import sys

import fingerprint
import pytest

pytestmark = pytest.mark.skipif(
    not (sys.platform == "linux" and platform.machine() == "x86_64"),
    reason="the fingerprint was recorded on Linux x86-64",
)


def _cases() -> list:
    return [
        pytest.param(c, id=c["id"], marks=[] if c["fast"] else [pytest.mark.slow])
        for c in fingerprint.load()["cases"]
    ]


@pytest.mark.parametrize("case", _cases())
def test_quicksynth_answers_as_fingerprinted(case: dict) -> None:
    got = fingerprint.run_case(case)
    lines = fingerprint.diff(fingerprint.outcome(case), got)
    assert not lines, f"{case['id']}:\n" + "\n".join(lines)
