<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# Development

Run every command on this page from the repository root.

## Prerequisites

| For | Needed |
| --- | --- |
| the package | Python ≥ 3.12, a C++17 compiler, CMake ≥ 3.18 |
| development | [uv](https://docs.astral.sh/uv/) 0.12 or later, which reads `uv.lock` (lockfile revision 3) without rewriting it |

The core uses only the C++ standard library. nanobind, which the build
fetches, binds it to Python. The dev container ships the C++ toolchain and
uv 0.12.17 ([.devcontainer/README.md](../.devcontainer/README.md)).

On Windows, build with clang-cl from Visual Studio's "C++ Clang tools for
Windows" component. MSVC's own compiler lacks `__int128`, so CMake stops under
it. `pip install .` selects the ClangCL toolset. With another generator, set
`CMAKE_CXX_COMPILER=clang-cl`.

## Building

```bash
make sync    # uv sync --all-extras
make cpp     # uv sync --all-extras --reinstall-package quicksynth
```

`make sync` installs the `dev` dependency group, the extras `plot` (matplotlib)
and `cpp` (scikit-build-core, nanobind, CMake, Ninja), and the package as an
editable install. Edits under `python/quicksynth/` take effect at once.

scikit-build-core compiles `quicksynth._core` in Release, in
`build/<wheel tag>/`, and does not rebuild it on import. After an edit under
`include/` or `src/`, or to `CMakeLists.txt`, run `make cpp`.

uv installs into `$UV_PROJECT_ENVIRONMENT` if it is set and into `.venv/`
otherwise. Give each checkout its own environment, because checkouts that
share one also share one `_core`. The dev container does this for you. To see
which module is in use:

```bash
uv run python -c "import quicksynth._core as c; print(c.__file__)"
```

`CMakeLists.txt` builds the static library `quicksynth_core` at `-O3` and the
nanobind module `_core` at nanobind's `-Os`. Both compile with
`-ffp-contract=off -Wall -Wextra`. Any other build of the core needs
`-ffp-contract=off` too ([Bit for bit](architecture.md#bit-for-bit)). The
option `QUICKSYNTH_TESTS`, off for the wheel, adds the C++ tests of
`tests/cpp/`.

## Wheels

`make wheel` runs `uv build`, which writes the sdist and a wheel into `dist/`.
The sdist holds the sources, the tests and `tools/`.

cibuildwheel builds the release wheels from `[tool.cibuildwheel]` in
`pyproject.toml`. It covers CPython 3.12, 3.13 and 3.14 on manylinux x86_64
and aarch64, macOS arm64 and x86_64 (deployment target 11.0), and Windows x64.
Each runner builds only its own architecture. auditwheel, delocate or
delvewheel checks each wheel and sets its platform tag, and delvewheel copies
the C++ runtime DLL into the Windows wheel. Each wheel's test runs `quicksynth`
once, then the Python tests that need neither the source tree's build nor the
paper's records.

nanobind and its tsl::robin_map are compiled into the module. Their license
texts go into the wheel's metadata, under `licenses/bundled/`.

`.github/workflows/wheels.yml` builds the sdist and the wheels, one runner per
platform, and keeps them as artifacts of the run. The weekly run builds
linux-x86_64 only. The publish job uploads a tagged release to PyPI and stays
switched off until the first release ([releasing.md](releasing.md)).

`cibuildwheel --platform linux` runs the Linux builds locally in docker.
Without docker, build for the host's glibc in a build directory of its own:

```bash
uv build --wheel -C build-dir=/tmp/qs-build
uvx --with patchelf auditwheel repair -w wheelhouse dist/*.whl
```

## Documentation site

`make docs` builds the site into `build/docs` with Sphinx, MyST and the PyData
theme, from the `docs` dependency group. Warnings fail the build.
`.github/workflows/docs.yml` builds the site on every pull request and
publishes it to GitHub Pages from `main`.

The source directory is the repository root, so the pages link to each other
by their repository paths. `docs/conf.py` lists the pages and keeps
`RESULTS.md`, which the README links to, out of the navigation. `docs/api.md`
reads the docstrings with the C++ core mocked.

The interactive pages `docs/algorithm.md` and `docs/benchmark.md` keep their
scripts and data in `docs/_static/showcase/`. `results.json` there holds the
aggregates plotted in the figures of `RESULTS.md`. The paper's figure script
writes the file and the figures, so replace them together.
`python tools/showcase_data.py` rewrites the algorithm page's data. It times
the scan, so run it on an idle machine.

## Tests

The library is tested in C++ and the package in Python. No test reads the
paper's records. [CONTRIBUTING.md](../CONTRIBUTING.md) lists the make targets
and [the CI jobs](../CONTRIBUTING.md#continuous-integration) that run them.

### C++ tests

`make test-cpp` configures `build/tests` with `-DQUICKSYNTH_TESTS=ON`, builds
the tests and runs them with ctest. The first build fetches doctest and
nlohmann_json, pinned in `tests/cpp/CMakeLists.txt`.

```bash
make test-cpp CTEST_ARGS="-L golden"   # only the replays of the recorded fixtures
make test-cpp CTEST_JOBS=4             # tests run at once (default 8)
```

Each concept has one doctest file in the directory of its step:
`tests/cpp/<step>/test_<concept>.cpp` tests
`include/quicksynth/<step>/<concept>.hpp` ([Layout](architecture.md#layout)).
ctest runs each test case on its own and labels the replays of recorded
fixtures `golden`. `tests/cpp/fixtures.hpp` reads the JSON fixtures, and
`tests/cpp/oracles.hpp` is the long-double counterpart of
`tests/python/oracles.py`.

### Python tests

The tests in `tests/python/` use the package's public API. pytest puts
`tests/python/`, `tests/fixtures/` and `tools/` on the import path.
`conftest.py` provides the repository root, a fixture reader and the
calibration targets.

| File | Tests |
| --- | --- |
| `test_api.py` | QuickSynth, its one-call forms, `prepare`, the policy, result types |
| `test_exact.py` | the exact layer, its errors, elements past 64 bits |
| `test_golden.py` | `synthesise` against `golden_synthesise.json` |
| `test_deploy.py` | edge cases: near-Clifford targets, time limit, errors, interrupts |
| `test_stopping_rule.py` | the model stopping rule through `synthesise` |
| `test_cli.py` | the `quicksynth` command |
| `test_counting.py` | the paper's counting in `counting.py` |
| `test_paper.py` | the paper's lower-bound statements, re-derived independently |
| `test_oracles.py` | the reference helpers of `oracles.py` |
| `test_fingerprint.py` | answers against `fingerprint.json` ([The fingerprint](#the-fingerprint)) |
| `test_layers.py` | the source layout, the include rule, CMake's source lists |
| `test_license.py` | the MIT banner in every tracked source file |

### Slow tests and calibration targets

Long tests carry the pytest marker and the ctest label `slow`. Run only them
with `pytest -m slow` or `ctest -L slow`. `make test-fast` leaves them out of
both suites, as `pytest -m "not slow"` and `ctest -LE slow` do.

Tests that draw Haar targets read `tests/fixtures/calibration_targets.json`,
through the `calibration_targets` fixture of `conftest.py` in Python and
`calibration_targets()` of `tests/cpp/fixtures.hpp` in C++.

## Golden fixtures

A golden fixture is a JSON file in `tests/fixtures/` that stores results of an
earlier version of the library. The tests rerun the same inputs and require the
same results, so a change that should keep every answer can show that it does.

The rules:

- A port, a refactor or a speed-up must reproduce every fixture unchanged.
- Only a change meant to alter an answer, such as a new stopping rule or a fix
  that changes a certificate, may re-record. It re-records only the runs or
  cases it changes, and its commit states which changed and why.
- No script in the tree writes the recorded or the frozen fixtures, so a
  change that re-records them brings its own recorder.
- The frozen fixtures are the reference the core was ported against and are
  not re-recorded, with the few exceptions in
  [their history](../tests/fixtures/README.md#re-recordings).
- In `golden_rotation.json` only QuickSynth's own costs may be re-recorded.

| File | Contents | Kind |
| --- | --- | --- |
| `calibration_targets.json` | 100 Haar-random targets (Euler angles and matrix) that the tests draw from | input |
| `golden_scan.json`, `golden_split.json`, `golden_combined.json` | search runs | [recorded](#recorded-fixtures) |
| `golden_words.json`, `golden_balls.json`, `golden_prefixes.json` | exact layer, balls, prefix sets | [frozen](#frozen-fixtures) |
| `golden_synthesise.json` | end-to-end `synthesise` calls | [frozen](#frozen-fixtures) |
| `golden_window.json` | model stopping rule | [reference](#the-stopping-rule-fixture) |
| `golden_rotation.json` | stage 0 on z-rotations | [reference](#the-rotation-fixture) |
| `golden_relative_norm.json` | norm equations | [reference](#the-norm-equation-fixture) |
| `fingerprint.json` | 444 fixed calls | [fingerprint](#the-fingerprint) |

[Fixture records](../tests/fixtures/README.md) tells how each fixture was
recorded and lists every change that re-recorded one.

### Recorded fixtures

`golden_scan.json` holds runs of the scale scan, `golden_split.json` runs of
the prefix split on a radius-6 ball, and `golden_combined.json` runs of the
scan followed by the split. They cover the first 12 calibration targets at ε
down to 1e-7, on one thread. The tests labelled `golden` in
`tests/cpp/search_scales/` replay every run.

- The scan runs come in four variants: a fixed throughput, `certify` (replayed
  with `rate=None` and no scale cap), `scales=2` and `rate=None`.
- Every run must return the recorded element with its cost, n_t and n_q. The
  split must also return its distance and its log's scales and leaves to the
  last bit.
- The scan's recorded `ties` are not admitted.

### Frozen fixtures

`golden_words.json`, `golden_balls.json`, `golden_prefixes.json` and
`golden_synthesise.json` hold outputs of the Python implementation that the C++
core replaced. They cover elements and words, every coset of the small balls,
the prefix sets and `synthesise` end to end.

- The C++ tests and `test_golden.py` require every value to the last bit,
  except the matrices and distances, which the old code computed with mpmath.
- The C++ tests hold those to 4e-16. `test_golden.py` holds them to 1e-15,
  because its targets also carry the host's LAPACK bits.
- `test_golden.py` replays every run of `golden_synthesise.json` with the fixed
  stopping rule, whole scales and no stage 0, which the runs predate. Off
  Linux x86-64 it also admits the transposed word for a target equal to its
  transpose.

### The stopping-rule fixture

`golden_window.json` holds the model stopping rule's values from a reference
implementation, from the census to the dynamic programme on a grid of states.
Nothing in it is timed. `tests/cpp/search_scales/test_stopping_rule.cpp`
compares the core with it. `make_window.py` rewrites it from the reference
scripts, which are not in this repository:

```bash
uv run --no-sync python tests/fixtures/make_window.py REFERENCE_DIR
```

### The rotation fixture

`golden_rotation.json` holds stage 0's references on z-rotations R_z(θ) at the
default price, at ε = 1e-3, 1e-4, 1e-6, 5.68e-8, 1e-8 and 1e-10. Its `meta`
names the angle sets, the sources and the commit. The stage 0 tests in
`tests/cpp/test_synthesise.cpp` read it, verify each answer's distance and
cost, and run the larger sections only in the `slow` tests.

| Section | Rows | Stage 0 must |
| --- | --- | --- |
| `clifford_t` | 120 angles per ε, with pygridsynth's least T-count | answer at that T-count, certified |
| `sqrt_t` | 50 angles per ε, with an external reference cost and QuickSynth's default and `fast` costs | answer at those costs or with a cheaper pure-T word; default: within 2 of the reference, below it in each ε's mean; `fast`: within half a gate in the mean |
| `zz` | the 392 ZZ angles of the benchmark workloads, with QuickSynth's costs | answer at those costs, certified on Clifford+T |
| `diagonal` | the 35 runs of `golden_synthesise.json` that stage 0 answers | answer at that cost, at most the recorded joint walk's |
| `near_exact` | 24 angles near kπ/4 at 1e-8 and 1e-10, with pygridsynth's T-count and QuickSynth's costs | on Clifford+T, answer at that T-count, certified unless a walk was cut |

### The norm-equation fixture

`golden_relative_norm.json` holds 2696 relative norm equations t t* = ξ over
Z[ω] and Z[ζ16], each with its number of solutions up to roots of unity, from
PARI/GP's `idealfactor`. Its `description` names the families.
`tests/cpp/solve_norm_equation/test_relative_norm.cpp` checks the solver's
outcome and number of solutions against it.
`tests/fixtures/check_relative_norm.py` recomputes every count with PARI/GP
through cypari2, which the package does not depend on.

## The fingerprint

`tests/fixtures/fingerprint.json` holds QuickSynth's answers to a fixed set of
calls, taken from the golden runs, the benchmark campaign, the stopping-rule
evaluation and a few chosen targets
([case sources](../tests/fixtures/README.md#fingerprint-cases)).
`tools/fingerprint.py` checks and records it:

```bash
uv run python tools/fingerprint.py check -n 6               # every case, on six processes
uv run python tools/fingerprint.py check --fast             # the fast subset
uv run python tools/fingerprint.py check --case c1t-1e-08-05
uv run python tools/fingerprint.py record -n 6              # needs the benchmark records
```

A case is one call. It stores ε, the gate set, every setting that differs
from the default, and the target as the nearest unitary the core receives, in
hexadecimal floats. The check therefore does not depend on the host's LAPACK.
The 99 cases whose target stage 0 takes run with `rotation_path: false`, since
their answers predate it. The others run with
stage 0 on, which pins that it leaves their targets alone.

`check` compares the word, the cost, n_t, n_q, the distance to the last bit,
`certified`, `stages`, the search summary and every search-log entry except its
time. It compares a call that raises by its exception's type and text, prints
each field that differs and exits with 1. A full check takes about ten minutes
of one CPU.

`test_fingerprint.py` runs each case as a test, on Linux x86-64 only, where the
fixture was recorded. The 206 `fast` cases run in every `pytest` run, and so in
CI, in a few seconds. The others carry the `slow` marker. A log entry's
`points` may move by 2 or by 5%, whichever is larger, because on the GitHub
runners they differ from the recording host by up to 3%.
