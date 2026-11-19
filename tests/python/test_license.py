#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""License headers: every source file opens with the MIT banner, and the old
proprietary header ("All rights ...") appears nowhere."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from conftest import ROOT

RULE = "===----------------------------------------------------------------------===//"
BANNER = [
    RULE,
    "",
    "Part of QuickSynth, under the MIT License.",
    "See <repo-root>/LICENSE.txt for license information.",
    "Copyright (c) 2026 Parity Quantum Computing GmbH",
    "SPDX-License-Identifier: MIT",
    "",
    RULE,
]
# the old header's phrase, assembled so that this file does not carry it
RESERVED = " ".join(("All", "rights", "reserved")).encode()
# files that describe the old header rather than carry it
MENTION_RESERVED = {"tools/release/export_public.py", "docs/releasing.md"}


def tracked_files() -> list[str]:
    """The tracked files; the test is skipped outside a git checkout (sdist)."""
    try:
        out = subprocess.run(
            ["git", "-C", str(ROOT), "ls-files", "-z"],
            capture_output=True,
            check=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout")
    files = [path for path in out.decode().split("\0") if path]
    if not files:
        pytest.skip("not a git checkout")
    return files


def banner(path: str) -> list[str] | None:
    """The banner in the comment syntax of a source file, None for other files."""
    if path.endswith((".cpp", ".hpp")):
        return ["// " + line if line else "//" for line in BANNER]
    if path.endswith((".py", ".sh")) or Path(path).name == "CMakeLists.txt":
        return ["#  " + line if line else "#" for line in BANNER]
    return None


def test_every_source_file_opens_with_the_license_banner() -> None:
    missing = []
    for path in tracked_files():
        expected = banner(path)
        if expected is None or not (ROOT / path).is_file():
            continue
        lines = (ROOT / path).read_text(encoding="utf-8").split("\n")
        if lines[0].startswith("#!"):
            # shebang, blank line, banner
            if lines[1:2] != [""]:
                missing.append(path)
                continue
            lines = lines[2:]
        n = len(expected)
        if (
            lines[:n] != expected
            or lines[n : n + 1] != [""]
            or lines[n + 1 : n + 2] == [""]
        ):
            missing.append(path)
    assert not missing, (
        "files without the license banner (followed by one blank line): "
        f"{', '.join(missing)}"
    )


def test_no_file_reserves_all_rights() -> None:
    reserved = [
        path
        for path in tracked_files()
        if path not in MENTION_RESERVED
        and (ROOT / path).is_file()
        and RESERVED in (ROOT / path).read_bytes()
    ]
    assert not reserved, f'"{RESERVED.decode()}" in {", ".join(reserved)}'
