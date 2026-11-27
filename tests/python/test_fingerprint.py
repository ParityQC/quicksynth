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
recorded; the test of ``update`` runs on any host.
"""

from __future__ import annotations

import platform
import re
import sys

import fingerprint
import pytest

recording_host = pytest.mark.skipif(
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


@recording_host
@pytest.mark.parametrize("case", _cases())
def test_quicksynth_answers_as_fingerprinted(case: dict) -> None:
    got = fingerprint.run_case(case)
    lines = [
        line
        for line in fingerprint.diff(fingerprint.outcome(case), got)
        if _beyond_tolerance(line)
    ]
    assert not lines, f"{case['id']}:\n" + "\n".join(lines)


def test_update_rewrites_only_the_cases_that_differ(tmp_path) -> None:
    """update re-runs a fixture's cases in place: a case whose recorded outcome
    no longer holds gets this build's, the others and meta stay byte for byte.
    Any host: both cases carry this host's outcomes."""
    fixture = fingerprint.load()
    cases = [c for c in fixture["cases"] if c["fast"] and "result" in c][:2]
    for case in cases:
        case["result"] = fingerprint.run_case(case)["result"]
    fixture["cases"] = cases
    expected = fingerprint._dump(fixture)
    cases[0]["result"]["cost"] += 1
    path = tmp_path / "fingerprint.json"
    path.write_text(fingerprint._dump(fixture))
    assert fingerprint.main(["update", "--fixture", str(path)]) == 0
    assert path.read_text() == expected
    assert fingerprint.main(["update", "--fixture", str(path), "--case", "none"]) == 2
