#!/usr/bin/env bash

#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

# Runs once after the dev container is created.
set -euo pipefail

# The uv cache volume keeps the image's owner; the container user's UID may
# have been changed to the host user's since.
cache="${HOME}/.cache/uv"
if [ -d "${cache}" ] && [ ! -w "${cache}" ]; then
    sudo chown -R "$(id -u):$(id -g)" "${cache}"
fi

echo "==> Syncing Python environment into ${UV_PROJECT_ENVIRONMENT:-.venv}"
uv sync --all-extras

echo "==> Installing pre-commit hooks"
if [ -e .git ]; then
    uv run pre-commit install
else
    echo "    (not a git repository yet - skipping)"
fi

echo "==> Toolchain"
uv run python --version
uv --version
cmake --version | head -n 1
ninja --version | sed 's/^/ninja /'
g++ --version | head -n 1

echo "==> Smoke test"
uv run python -c "import quicksynth as s; print('quicksynth', s.__version__, 'from', s.__file__)"

echo "==> Ready. Try: make demo"
