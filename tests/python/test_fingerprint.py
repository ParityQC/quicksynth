#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""QuickSynth against ``tests/fixtures/fingerprint.json`` (``tools/fingerprint.py``).

Every recorded field must match, except a log entry's ``points`` (the walk's
leaves), which may move by two or 5%: on the GitHub runners they differ from the
recording host by up to 3% (cause not yet traced), while every answer, distance
and estimate stays bit-exact. ``tools/fingerprint.py check`` compares exactly.
Non-fast cases are marked ``slow``. Linux x86-64 only, where the fixture was
recorded.
"""

from __future__ import annotations

import platform
import re
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


POINTS = re.compile(r"\.log\[\d+\]\.points: (\d+) -> (\d+)$")


def _beyond_tolerance(line: str) -> bool:
    """Whether a difference is more than a leaf count moved by two or 5%."""
    match = POINTS.search(line)
    if match is None:
        return True
    want, got = int(match[1]), int(match[2])
    return abs(got - want) > max(2, want // 20)


@pytest.mark.parametrize("case", _cases())
def test_quicksynth_answers_as_fingerprinted(case: dict) -> None:
    got = fingerprint.run_case(case)
    lines = [
        line
        for line in fingerprint.diff(fingerprint.outcome(case), got)
        if _beyond_tolerance(line)
    ]
    assert not lines, f"{case['id']}:\n" + "\n".join(lines)
