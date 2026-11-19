<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# Contributing

This guide covers the development setup, the checks a change must pass and the
continuous integration. [docs/development.md](docs/development.md) describes the
build, the tests and their fixtures in detail, and
[docs/architecture.md](docs/architecture.md) describes how the library fits
together.

## Setup

Development needs Python 3.12 or later, a C++17 compiler and
[uv](https://docs.astral.sh/uv/) 0.12 or later.

```bash
make sync                   # uv sync --all-extras: editable install and dev tools
uv run pre-commit install   # the hooks of .pre-commit-config.yaml
```

The Python sources are imported in place. After an edit under `include/`,
`src/` or to `CMakeLists.txt`, run `make cpp` to rebuild the C++ core. On
Windows the core builds with clang-cl only
([docs/development.md](docs/development.md#prerequisites)).

## Checks

```bash
make check             # lint, typecheck and test-fast: run this before a pull request
make lint              # ruff format --check and ruff check
make format            # ruff format and ruff check --fix
make typecheck         # mypy
make test-fast         # the C++ and Python tests without the slow ones
make test              # both suites, slow tests included
make test-sanitizers   # the C++ tests without the slow ones, under ASan and UBSan
```

The pre-commit hooks run ruff, clang-format 14.0.6 and the standard file
checks; `uv run pre-commit run --all-files` runs them on the whole tree.

A change to a golden fixture follows the rules in
[docs/development.md](docs/development.md#golden-fixtures).

## Continuous integration

`.github/workflows/ci.yml` runs on pull requests, on pushes to `main`, weekly
and on demand:

| Job | What it runs |
| --- | --- |
| Lint | `uv lock --check`, `make lint`, `make typecheck`, the remaining pre-commit hooks and actionlint |
| Test | `make test-fast` on Ubuntu (Python 3.12 to 3.14) and on macOS arm64 and x86_64 (3.12 and 3.14), and the same ctest and pytest commands on Windows x64 (3.12 and 3.14); the weekly run tests Ubuntu only |
| Full tests and coverage | `make test` on Ubuntu, slow tests included; not on pull requests |
| Sanitizers | `make test-sanitizers` on Ubuntu |

`.github/workflows/wheels.yml` builds the wheels and the sdist
([docs/development.md](docs/development.md#wheels)). It runs on tags `v*`, on
pushes to `main`, weekly, on demand, and on pull requests that change
`pyproject.toml`, `CMakeLists.txt` or the workflow. Its publish job uploads a
tagged release to PyPI; it is switched on for the first release
([docs/releasing.md](docs/releasing.md)).

The workflows pin each action to a release commit and uv to the version that
wrote `uv.lock`. Dependabot proposes updates of the actions and of `uv.lock`
once a week.
