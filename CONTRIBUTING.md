<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# Contributing

[docs/development.md](docs/development.md) describes the build and the tests
in depth, and [docs/architecture.md](docs/architecture.md) the design.

## Setup

You need Python 3.12 or later, a C++17 compiler (clang-cl on Windows), CMake
3.18 or later and [uv](https://docs.astral.sh/uv/) 0.12 or later
([prerequisites](docs/development.md#prerequisites)). The
[dev container](.devcontainer/README.md) comes with all of them.

```bash
make sync                   # editable install with the dev tools
uv run pre-commit install   # install the hooks of .pre-commit-config.yaml
```

After editing `include/`, `src/` or `CMakeLists.txt`, run `make cpp`. Python
edits need no rebuild.

## Checks

Run `make check`, which mirrors the Lint and Test jobs of CI, before you open a
pull request.

| Command | What it runs |
| --- | --- |
| `make check` | `lint`, `typecheck` and `test-fast` |
| `make test-fast` | the C++ and Python tests without the slow ones |
| `make test` | both test suites, slow tests included |
| `make test-sanitizers` | the C++ tests without the slow ones, under ASan and UBSan |
| `make lint` | `ruff format --check` and `ruff check` |
| `make format` | `ruff format` and `ruff check --fix` |
| `make typecheck` | mypy |
| `make cpp` | rebuilds the C++ core after a C++ edit |

The pre-commit hooks run `ruff check --fix`, `ruff format`, clang-format on the
C++ files and the standard file checks.

```bash
uv run pre-commit run --all-files                # every hook on the whole tree
uv run pre-commit run clang-format --all-files   # format every C++ file
```

The hook pins clang-format 14.0.6, the dev container's version, because other
versions lay out the same `.clang-format` differently. The end-of-file fixer
skips `tests/fixtures/calibration_targets.json`, which is kept byte for byte as
it was written, without a final newline.

## Golden fixtures

The tests compare the library with results of earlier versions, stored in
`tests/fixtures/`. A change that should not alter any answer must reproduce them
unchanged. Only a change meant to alter an answer may re-record, and only the
runs or cases it changes. Its commit states which changed and why. The full
rules are in [docs/development.md](docs/development.md#golden-fixtures).

## Continuous integration

`.github/workflows/ci.yml` runs on pull requests, on pushes to `main`, weekly
and on demand.

| Job | What it runs |
| --- | --- |
| Lint | `uv lock --check`, `make lint`, `make typecheck`, the other pre-commit hooks and actionlint |
| Test | `make test-fast` on Ubuntu (Python 3.12, 3.13 and 3.14) and on macOS arm64 and x86_64 (3.12 and 3.14); the equivalent ctest and pytest commands on Windows x64 (3.12 and 3.14) |
| Full tests and coverage | `make test` on Ubuntu, slow tests included; not on pull requests |
| C++ tests under ASan and UBSan | `make test-sanitizers` on Ubuntu |

The weekly run tests on Ubuntu only.

`.github/workflows/wheels.yml` builds the wheels and the sdist
([Wheels](docs/development.md#wheels)). It runs on `v*` tags, on pushes to
`main`, weekly, on demand, and on pull requests that change `pyproject.toml`,
`CMakeLists.txt` or the workflow itself.

The workflows pin each action to a release commit and uv to the version that
wrote `uv.lock`. Dependabot proposes updates of the actions and of `uv.lock`
once a week.
