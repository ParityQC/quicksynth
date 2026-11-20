<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# Development

This guide covers building the package and the wheels, running the tests, and
the rules of the golden fixtures. Run every command from the repository root.

## Prerequisites

| For | Needed |
| --- | --- |
| the package | Python ≥ 3.12, a C++17 compiler, CMake ≥ 3.18 |
| development | [uv](https://docs.astral.sh/uv/) 0.12 or later, which reads `uv.lock` (lockfile revision 3) without rewriting it |

The core uses only the C++ standard library; nanobind, which the build
fetches, binds it to Python. The dev container (`.devcontainer/`) ships the
C++ toolchain and uv 0.12.17.

On Windows the compiler must be clang-cl, the LLVM compiler that Visual Studio
installs with its "C++ Clang tools for Windows" component: the core's
arithmetic uses `__int128`, which MSVC's own compiler does not have, and CMake
stops with a message under it. `pip install .` builds with Visual Studio's
ClangCL toolset; with another generator, set `CMAKE_CXX_COMPILER=clang-cl`.
CMake links compiler-rt's builtins library, which clang-cl names, for the
conversions between `__int128` and `double`.

## Building

```bash
make sync    # uv sync --all-extras
make cpp     # uv sync --all-extras --reinstall-package quicksynth
make wheel   # uv build: sdist and wheel into dist/
```

`make sync` installs three things: the `dev` dependency group, the extras
`plot` (matplotlib) and `cpp` (scikit-build-core,
nanobind, CMake, Ninja), and the package as an editable install:

- `python/quicksynth/` is imported in place, so an edit to the Python takes
  effect at once.
- scikit-build-core compiles `quicksynth._core` from `CMakeLists.txt` and
  installs it into the environment's site-packages. The build is Release with
  `-O3`, and the build tree is `build/<wheel tag>/`. `editable.rebuild` is off,
  so an edit under `include/` or `src/` or to `CMakeLists.txt` takes effect
  only after `make cpp`.

uv installs into `$UV_PROJECT_ENVIRONMENT` if that is set, and into `.venv/`
otherwise. Each checkout needs its own environment; two checkouts that share
one also share one `_core`. The dev container sets it to `.venv-container`, a
path relative to the project root, so every checkout and worktree in the
container gets its own. An absolute path would be shared by all of them.

```bash
make sync
uv run python -c "import quicksynth._core as c; print(c.__file__)"   # the module in use
```

The core, both the library and the module, is compiled with
`-ffp-contract=off`. The inputs of the lattice reduction are fixed sequences
of correctly rounded operations, and a fused multiply-add would change them.
Any other build of the core needs the same flag.

`make lint`, `make format` and `make typecheck` run ruff and mypy.
`uv run pre-commit install` installs the hooks of `.pre-commit-config.yaml`:
ruff check with `--fix`, ruff format, clang-format on the C++ files
(`.clang-format`), and the standard file checks. The end-of-file fixer skips
`tests/fixtures/calibration_targets.json`, which is kept byte for byte as it
was written, without a final newline. The hook pins clang-format 14.0.6, the dev
container's, because other versions lay the same `.clang-format` out
differently; `uv run pre-commit run clang-format --all-files` formats every
C++ file.

### Wheels

cibuildwheel builds the wheels from `[tool.cibuildwheel]` in `pyproject.toml`:
CPython 3.12, 3.13 and 3.14 on manylinux x86_64 and aarch64, on macOS arm64
and x86_64 (deployment target 11.0), and on Windows x64. A runner builds only
its own architecture, so each of the five platforms needs a runner of its own.
The platform's repair tool checks each wheel and sets its platform tag:
auditwheel on Linux, delocate on macOS and delvewheel on Windows, which also
copies the C++ runtime DLL into the wheel. The build installs the license texts
of nanobind and of the tsl::robin_map it includes, which every build compiles
into the module, into the wheel's metadata, under `licenses/bundled/`. The test
of each wheel runs `quicksynth` once, then the Python tests that need neither
the paper's records nor the source tree's build, against the installed wheel.

`.github/workflows/wheels.yml` runs cibuildwheel on one runner per platform,
builds the sdist, and keeps both as artifacts of the run; the weekly run builds
linux-x86_64 only. Its publish job uploads a tagged release to PyPI
([releasing.md](releasing.md)).

`cibuildwheel --platform linux` runs the Linux builds locally in docker. The
same build works without docker, for the host's glibc:

```bash
uv build --wheel -C build-dir=/tmp/qs-build
uvx --with patchelf auditwheel repair -w wheelhouse dist/*.whl
```

The build directory is a new one, which keeps the build apart from the
development build in `build/`.

`uv build` also builds the sdist, which holds the sources, the tests and
`tools/`.

### Documentation site

`make docs` builds the site into `build/docs` (`build/docs/index.html`) with
Sphinx, MyST and the PyData theme, from `uv sync --group docs`. Its source
directory is the repository root, so the Markdown files keep linking to each
other by their repository paths; `docs/conf.py` lists the pages, `docs/index.md`
is the landing page and `docs/api.md` the API reference, which reads the
docstrings with the C++ core mocked. Warnings fail the build.
`docs/how-it-works/` holds the interactive pages; their scripts and data are in
`docs/_static/showcase/`. `results.json` holds the aggregates plotted in the
figures of `RESULTS.md`, exported from the paper's figure script with them:
replace both together. `python tools/showcase_data.py` rewrites the data of
"How QuickSynth finds a circuit" (`mechanism.json`, `points.bin`) from the
library; it times the scan, so run it on an idle machine.
`.github/workflows/docs.yml` builds it on every pull request and publishes it
to GitHub Pages from `main`.

## Tests

```bash
make test-cpp         # the C++ tests: cmake -DQUICKSYNTH_TESTS=ON in build/tests, then ctest
make test-fast        # both suites without the slow tests
make test             # both suites, slow tests included
make test-sanitizers  # the C++ tests without the slow ones, under ASan and UBSan
make check            # lint, typecheck, then both suites without the slow tests: what CI runs
```

The C++ build fetches doctest 2.4.11 and nlohmann_json 3.11.3 once, by URL and
SHA256 (`tests/cpp/CMakeLists.txt`). `make test-cpp CTEST_ARGS="-L golden"`
runs only the replays of the recorded fixtures, `-LE slow` leaves out the slow
tests, and `CTEST_JOBS` sets how many run at once (8).

CI (`.github/workflows/ci.yml`) runs these checks on Linux, macOS and Windows;
[CONTRIBUTING.md](../CONTRIBUTING.md#continuous-integration) lists its jobs.

The library is tested in C++ and the package in Python. `tests/cpp/` holds one
file per concept, in the directory of the step it tests
(`tests/cpp/<step>/test_<concept>.cpp` for
`include/quicksynth/<step>/<concept>.hpp`, and `tests/cpp/test_<concept>.cpp`
for the entry's headers at the top level), written with doctest; ctest runs each
test case on its own, labelled `golden` for the replays of the recorded fixtures
and `slow` for the slow ones.
`tests/cpp/fixtures.hpp` reads the JSON fixtures, integers past 64 bits
included, and `tests/cpp/oracles.hpp` is the C++ counterpart of
`tests/python/oracles.py`, in long double.

The Python tests, in `tests/python/`, import the package's public API, and the
paper's counting, `tests/fixtures/counting.py`, which pytest puts on the import
path with `tests/python/`.
`tests/python/oracles.py` holds what the tests check the package with: the
diamond distance in double precision, the magic counts of a word, a key of an
element's Clifford coset, and the Haar-random and stored targets.
`tests/python/conftest.py` gives the repository root, a reader of the JSON
fixtures and the calibration targets.

| File | Pins |
| --- | --- |
| `test_api.py` | QuickSynth and its one-call forms, `prepare`, the policy reaching the stage settings, the result's field types |
| `test_exact.py` | the exact layer through the package: the README's example, its errors, the spelling of random words, elements past 64 bits, the T-count off Clifford+T, and the matrices of elements past 64 bits and of a word past a double's range against 50 digits (mpmath) |
| `test_golden.py` | `synthesise` against the QuickSynth runs of `golden_synthesise.json` |
| `test_deploy.py` | targets near a Clifford, Clifford+T below 1e-4 and at 1e-9, the time limit, the same word under load, settings at exactly 1e-7, error order and text, verbose lines, the global phase, an interrupted call |
| `test_stopping_rule.py` | the model stopping rule through `synthesise`: a refused price, the fallback, the default, the split's log, determinism over threads |
| `test_cli.py` | the `quicksynth` command |
| `test_counting.py` | the paper's counting, `counting.py`: the recurrence of the levels, the Poisson mean against the union bound, `ball_measure` against `golden_window.json` |
| `test_paper.py` | the statements the paper's lower bound rests on, re-derived independently of the library |
| `test_oracles.py` | the oracles of `tests/python/oracles.py` |
| `test_fingerprint.py` | QuickSynth's answers against `fingerprint.json`, field for field ([The fingerprint](#the-fingerprint)) |
| `test_layers.py` | the pipeline layout of `include/quicksynth/` and `src/`, its include rule and private headers, the C++ tests in the step of their subject, CMake building every source and test |

The tests marked `slow` are in `test_deploy.py`, `test_stopping_rule.py`,
`test_golden.py` and `test_fingerprint.py` (`-m slow`
selects them), and the C++ tests labelled `slow` are in
`tests/cpp/search_scales/test_scan.cpp`,
`tests/cpp/synthesise_rotation/test_rotation_body.cpp` and
`tests/cpp/test_synthesise.cpp` (`-L slow`).
The tests that draw Haar targets read
`tests/fixtures/calibration_targets.json`, through the `calibration_targets`
fixture of `tests/python/conftest.py` in Python and `calibration_targets()` of
`tests/cpp/fixtures.hpp` in C++.

No test reads the paper's records.

### Golden fixtures

A golden fixture is a JSON file in `tests/fixtures/` that stores results of an
earlier version of the library. The tests run the same inputs again and require
the same results, so a change that should not alter any answer can show that it
does not.

| File | Contents | Kind |
| --- | --- | --- |
| `calibration_targets.json` | 100 Haar-random targets (Euler angles and matrix) that the tests draw from | input |
| `golden_scan.json`, `golden_split.json`, `golden_combined.json` | runs of the scale scan, the prefix split and both | recorded |
| `golden_words.json`, `golden_balls.json`, `golden_prefixes.json` | the exact layer (elements, Bloch images, prices, words, spellings, distances), the cost balls and the prefix sets | frozen |
| `golden_synthesise.json` | end-to-end `synthesise` calls | frozen |
| `golden_window.json` | reference values of the model stopping rule | reference |
| `golden_rotation.json` | stage 0's references on z-rotations | reference |
| `golden_relative_norm.json` | relative norm equations t t* = ξ over Z[ω] and Z[ζ16] with their number of solutions | reference |
| `fingerprint.json` | the answers to 444 fixed calls ([The fingerprint](#the-fingerprint)) | fingerprint |

`tests/fixtures/` also holds `counting.py`, the paper's counting that
`test_counting.py` checks, and `make_window.py`, which writes
`golden_window.json`. No script in the tree re-records the recorded or the
frozen fixtures.

**Recorded fixtures.** `golden_scan.json` holds runs of the scale scan at
several prices and on Clifford+T in four variants: under a fixed throughput,
`certify` (replayed with `rate=None` and no scale cap), `scales=2` and
`rate=None`. `golden_split.json` holds runs of the prefix split on a radius-6
ball, and `golden_combined.json` runs of the scan followed by the split. The
tests labelled `golden` in `tests/cpp/` (`search_scales/test_scan.cpp`,
`test_split.cpp`, `test_scan_then_split.cpp`) run every recorded run again and
compare the results.

```bash
make test-cpp CTEST_ARGS="-L golden"   # the replays
```

- The targets are the first 12 of `tests/fixtures/calibration_targets.json`. Each file
  stores them in its `meta`, with the git head, the numpy version and the CPU
  seconds of the recording.
- Everything runs on one thread.
- The fixed-throughput variants pass 1.4e7 (√T) and 2.5e7 (T) explicitly. These
  are the values the library's reference throughput had before the
  calibration, so the fixtures do not follow the library constant.
- The targets of the scan and the QuickSynth runs are taken at ε from 1e-2 to
  1e-6, and the first 4 also at 1e-7. The split runs at 1e-3 and 1e-5 at
  (1, 3) (prefix scales 1 and 2) and at 1e-4 and 1e-6 on Clifford+T.
- The replay is strict: every run must return the recorded element, with its
  cost, n_t and n_q, and the split also its distance, recomputed from the
  element, and its composite scales and leaves, read from its log, to the last
  bit. The scan's recorded `ties` are not admitted.

No script in the tree writes them now. A port, a refactor or a speed-up must
reproduce them unchanged; only a change meant to alter what a recorded run
returns, such as a new stopping rule or a fix that changes a certificate, may
re-record the runs it changes, with a recorder of its own, stating in the commit
which runs changed and why. QuickSynth's own LLL (`lll.cpp`), which replaced
FLINT's, reduces some bases differently: 126 runs were re-recorded for their
point counts and 9 for another element of the same cost.

**Frozen fixtures.** These are `golden_words.json`, `golden_balls.json`,
`golden_prefixes.json` and `golden_synthesise.json`. They hold the outputs
of an earlier Python implementation of the exact layer, which the C++ core
replaced: elements, Bloch images, prices, words, spellings and distances; every
coset of the small balls; prefix sets and prefix radii; and QuickSynth's
`synthesise` end to end. The C++ tests
(`tests/cpp/price_candidates`, `tests/cpp/exact_synthesis`, and the Bloch,
ball and distance tests in `tests/cpp/data_structures`) and
`test_golden.py` require the same values, to the last bit of every float but
the matrices and distances. The core evaluates those in double precision rather
than as CPython and mpmath (at 40 digits) did, so the C++ tests hold them to
4e-16, and `test_golden.py`, whose targets also carry the host's LAPACK bits,
to 1e-15. `test_golden.py` requires the recorded word itself on Linux x86-64,
where the fixture was recorded; elsewhere the host's LAPACK gives the target's
nearest unitary other last bits, and a target equal to its transpose may then
answer with the transposed word (run 107 on macOS x86-64), which is admitted
there. They are the reference the core was ported against, and they are not
re-recorded, except runs 13, 53 and 57 of `golden_synthesise`: QuickSynth's
own LLL answers them with another word of the same cost, and run 13's distance
is the new word's, computed at 60 digits.
They were recorded with that Python implementation, which is not in this
repository:

- One recording wrote the exact layer, the balls and the prefix sets into one
  file, since split into the first three; the `cpu_seconds` of each is that of
  the whole recording. The ball at the price (1, 1) was recorded to radius 5
  and is kept to radius 4.
- The `synthesise` runs were taken with the Python word and check.
- The `past_64_bits` section of `golden_words`, added later the same way, holds
  three words whose coefficients pass 64 bits (74, 71 and 125 bits), with their
  Bloch images, words and prices at five prices and their T-counts.
- Runs 30-122 of `golden_synthesise` were added later the same way: among
  them auto below 1e-3, a scan at (2, 5), four workers, a second-stage answer
  near a Clifford, and Cliffords, magic gates and short words on both gate sets
  down to 1e-8.
- The file keeps the paper's 73 QuickSynth runs in their order, byte for byte
  except runs 13, 53 and 57; `meta.paper_v1_runs` gives the index of each.

`test_golden.py` replays a recorded `method="scan"`, and `"auto"` below 1e-3,
which picked the scan, through QuickSynth under the fixed stopping rule, on
whole scales (`fine_levels=False`) and without stage 0, which the runs predate.

**The stopping-rule fixture.** `tests/fixtures/make_window.py` writes
`golden_window.json` from a reference implementation of the model stopping
rule, a set of scripts that is not in this repository: the census, its minima
at nine prices, `ball_measure`, and for five settings the one-step values, the
dynamic programme, the savings and the split's saving on a grid of states.
Nothing in it is timed.

**The rotation fixture.** `golden_rotation.json` holds stage 0's references on
z-rotations R_z(θ) at the default price, ε = 1e-3, 1e-4, 1e-6, 5.68e-8, 1e-8
and 1e-10; `meta` names the angle sets, the sources and the commit.

- `clifford_t`: 120 angles (100 uniform, 20 small) at each ε with the least
  T-count of the Ross-Selinger search with its level order corrected
  (pygridsynth); stage 0 must answer it, certified.
- `sqrt_t`: 50 angles at each ε with a reference cost, the native Clifford+√T
  cost of an external synthesiser run as a black box, and QuickSynth's default
  and `fast` costs. The default must stay within 2 of the reference on each
  angle and below it in the mean of each ε, `fast` within half a gate of it in
  the mean, and both at their recorded cost.
- `zz`: the 392 ZZ angles of the trotter, qaoa and adiabatic workloads with
  QuickSynth's T-count and default and `fast` costs.
- `diagonal`: the 35 runs of `golden_synthesise.json` that stage 0 answers, at
  its cost; on Clifford+T at most the recorded joint walk's and equal to it
  where that was certified.
- `near_exact`: 24 angles θ = kπ/4 + cε (k = 0, 1, 2; c = ±3, ±10) at ε = 1e-8
  and 1e-10, whose rotation bodies hold slabs of millions of points, with the
  pygridsynth T-count as for `clifford_t` and QuickSynth's Clifford+√T default
  and `fast` costs. On Clifford+T stage 0 must answer at that T-count, certified
  unless a walk was cut.

The C++ tests of stage 0 in `tests/cpp/test_synthesise.cpp` read it: every
Clifford+T row, 20 `sqrt_t` angles at three ε, the diagonal runs and the
Clifford+T `near_exact` rows in the default run, all `sqrt_t` rows, the `zz`
angles and the Clifford+√T `near_exact` rows under `slow`. Each answer is
verified: its distance against the long-double oracle and its cost priced
again from its word. Only a change meant to alter an answer may re-record the
QuickSynth costs, stating which changed and why.

**The norm-equation fixture.** `golden_relative_norm.json` holds 2696
relative norm equations t t* = ξ over Z[ω] and Z[ζ16] with their number of
solutions up to roots of unity, from PARI/GP's `idealfactor`; its
`description` names the families. `tests/fixtures/check_relative_norm.py`
recomputes every count with PARI/GP through cypari2, which is not a dependency
of the package (last run: PARI/GP 2.17.2, cypari2 2.2.4, no count differs).
`tests/cpp/solve_norm_equation/test_relative_norm.cpp` checks the solver's
outcome and number of solutions against it.

### The fingerprint

`tests/fixtures/fingerprint.json` holds what QuickSynth answers on a fixed set
of calls, so that a change meant to keep every answer can show that it does.
`tools/fingerprint.py` checks and records it:

```bash
uv run python tools/fingerprint.py check -n 6               # every case, on six processes
uv run python tools/fingerprint.py check --fast             # the fast subset
uv run python tools/fingerprint.py check --case c1t-1e-08-05
uv run python tools/fingerprint.py record -n 6              # needs the benchmark records
```

- A case is one call: the target, stored as the nearest unitary the core
  receives in hexadecimal floats; ε; and the gate set, the price, the stopping
  rule and any other setting that differs from the default. The check passes the
  stored matrix to the core's `QuickSynth.synthesise`, so it does not depend on
  the host's LAPACK. The 99 cases whose target stage 0 takes carry
  `rotation_path: false` in their settings, since their answers predate it;
  every other case runs with stage 0 on, so the check also pins that stage 0
  leaves every other target alone.
- Per case it compares the word, the cost, n_t, n_q, the distance to the last
  bit, `certified`, `stages`, the answering stage's search summary and every
  entry of the search log except its time; a call that raises is compared by
  the type and text of its exception. `check` prints each field that differs
  and exits with 1. The libraries' versions are printed, not compared.
- The 444 cases: every QuickSynth run of `golden_synthesise.json` under the
  fixed stopping rule it was recorded with and, except its three slow runs,
  under the model stopping rule; from the benchmark campaign, under the fixed
  stopping rule, 149 stage-1 calls (up to seven per setting) and one stage-2
  call; from the evaluation of the model stopping rule, 18 calls that answered
  at stage 2 and 71 stage-1 calls; the targets near a Clifford of
  `test_deploy.py`, one of which answers at stage 3; rotations and other targets
  with a Clifford symmetry on both gate sets; exact words moved by ε/3 at
  1e-10 to 1e-13; and the first campaign call of each setting again with
  `fine_levels=False` (`whole-*`), which keeps the whole-scale search the
  default ran before the finer levels. Every case taken from a benchmark record
  returned the recorded word, cost and stages when the fixture was first
  recorded.
- `meta` names the commit of the last full recording (`record`), which a
  re-record of some cases in place keeps, and the number of cases. It was
  recorded in the dev container, and recorded again after the switch to
  QuickSynth's own LLL, which changed the per-walk point counts of 120 cases
  and the word, at equal cost, of 16. It was recorded a third time for the
  refit of the stopping rule's time model and reference throughputs, which
  changed the search of 29 model-rule cases and lowered the cost of two
  (`golden-016-model` 28 to 27, `window-13-1e-04-78` 28 to 26). When the finer
  levels became the default on Clifford+√T, the 192 cases whose search changed
  were recorded again in place, among them 32 with another word (cost 49 lower
  in all, two dearer), and the 17 `whole-*` cases were added with the answers
  recorded before. The full check takes about ten minutes of one CPU.
- `test_fingerprint.py` runs each case as a test: the 206 marked `fast`, those
  that walk the fewest leaves at stage 1 and the `whole-*` cases of those, in
  every `pytest` run, and so in CI (a few seconds); the others under the
  `slow` marker. It runs on Linux x86-64 only, where the fixture was recorded.
- As for the recorded fixtures, only a change meant to alter an answer may
  re-record the fingerprint, stating in its commit which cases changed and why.
