#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Shared test helpers: repository root, fixture reader and calibration targets."""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from oracles import decode_target

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures"


@cache
def fixture_json(name: str) -> Any:
    """Parsed ``tests/fixtures/<name>``; cached, do not mutate."""
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture(scope="session")
def calibration_targets() -> list[np.ndarray]:
    """Targets of ``tests/fixtures/calibration_targets.json``."""
    return [decode_target(t) for t in fixture_json("calibration_targets.json")]
