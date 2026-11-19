# QuickSynth

**[Preliminary results](RESULTS.md):** cost and compile time against
cyclosynth and pygridsynth, ahead of the paper.

Near-optimal synthesis of single-qubit unitaries over **Clifford+T** and
**Clifford+√T** at any magic-state price. QuickSynth is implemented in C++
with an easy to use Python interface. If you use the algorithm for any
scientific publication, please cite us.

Given a 2x2 unitary and a tolerance ε, QuickSynth returns a gate word over
`H, S, X, Y, Z, T, Q` and the inverses `s, t, q` (`Q` = √T) within ε in
diamond distance, up to a global phase, at the least magic-state cost
`t_cost·n_T + q_cost·n_√T` it can find, fast at any precision.
It works on Clifford+√T and Clifford+T, at any price.

Every answer is checked before it is returned: the word is assembled exactly
over `Z[ζ16]`, and its distance to the target, evaluated from that element in
double precision without cancellation (typically to about 1e-16, at most
4.4e-16 measured), must lie within ε − 5e-16. The same input always gives the
same word, and between two magic gates a word carries at most two Clifford
letters.

## Installation

### From PyPI

QuickSynth is not on PyPI yet; until then, install it from source (below).
Once it is released, a prebuilt wheel installs it without a compiler:

```bash
pip install quicksynth
```

Wheels are built for CPython 3.12 to 3.14 on Linux (x86_64, aarch64), macOS
(arm64, x86_64) and Windows (x64). On any other platform, `pip install
quicksynth` builds the package from the source distribution, which needs the
tools listed below.

### From source

A build from source needs Python 3.12 or later and a C++17 compiler (clang-cl
on Windows). The build backend, scikit-build-core with nanobind, fetches CMake
and Ninja; the core needs no other library. In a clone of this repository:

```bash
pip install .            # builds the native core quicksynth._core
```

For development, use [uv](https://docs.astral.sh/uv/) 0.12 or later:

```bash
uv sync --all-extras     # editable install with the dev tools
make cpp                 # rebuild the native core after editing include/ or src/
make check               # lint, mypy and both test suites without the slow tests
make demo                # examples/quickstart.py
```

The repository also has a dev container (Python 3.12, uv and a C++
toolchain) for "Reopen in Container"
([.devcontainer/README.md](.devcontainer/README.md)).
[docs/development.md](docs/development.md) covers the build, the wheels and the
tests, and [CONTRIBUTING.md](CONTRIBUTING.md) the checks a change must pass.

## Command line

```bash
quicksynth pi/8 1e-3                     # R_z(pi/8) to 1e-3, by QuickSynth
quicksynth 0.7 1e-4 --axis x --json      # R_x(0.7), machine-readable
quicksynth --unitary '[["0.6+0.8j",0],[0,1]]' 1e-3 --qasm
quicksynth pi/8 1e-6 --gateset clifford_t  # cost = T-count
```

`THETA` may be an expression in `pi`, `e` and `sqrt`. The command takes
`--price T SQRT_T` (default `1 3`), `--gateset {sqrt_t,clifford_t}`,
`--workers`, `--unitary-file` (a `.npy` file), `--time-limit`,
`--stopping-rule`, `--extra-scales` and `--point-budget`; `--json` and `--qasm`
select the output. `--help` for listing all the commands.

## Python API

```python
import numpy as np
import quicksynth as qs

r = qs.synthesise_rotation(np.pi / 8, 1e-3)    # R_z(theta) = exp(-i theta Z / 2)
U = qs.rotation(0.7, "x")                      # or any 2x2 unitary
r = qs.synthesise(U, 1e-4, q_cost=2)           # at the price T = 1, sqrt(T) = 2
r = qs.QuickSynth(gateset="clifford_t").synthesise(U, 1e-3)   # T-count
r.word, r.cost, r.n_t, r.n_q, r.distance, r.certified
r.gates        # circuit order: the first gate acts first
r.unitary      # the word's matrix
r.to_qasm()    # OpenQASM 3, sqrt(T) as p(pi/8)
r.as_dict()
```

`synthesise(target, epsilon, **overrides)` accepts any 2x2 unitary; a matrix
within 1e-6 of unitary (float32 included) is replaced by the nearest unitary.
`synthesise_rotation(theta, epsilon, axis="z", ...)` does the same for a
rotation, and `QuickSynth(config=None, **overrides)` keeps one configuration
for many calls. The most used options:

| Option | Default | Meaning |
| --- | --- | --- |
| `gateset` | `"sqrt_t"` | `"clifford_t"` restricts the words to Clifford+T |
| `t_cost`, `q_cost` | `1`, `3` | magic-state price of `T` and `√T` |
| `fast` | `False` | search less past the first answer: faster, at a slightly higher cost on Clifford+√T |
| `workers` | `1` | threads for the scan's walk |
| `time_limit` | `None` | seconds after which no further search stage starts |

`r.certified` says whether the cost is proved minimal or is an upper bound.
Every option, and how the search uses them, is in
[docs/python-api.md](docs/python-api.md).

The exact layer is available on its own, at the package's top level. A word's
element `M / √2^k` over `Z[ζ16]` is exact, in integers of any size, and its
price and a cheapest word for it are read off its normal form. That word is the
element up to a global phase `ζ^j`:

```python
from quicksynth import bloch_image, matrix_word, normal_form, price, word_unitary

element = word_unitary("HSQHTHQ")   # k = 3 and the integer coefficients of M
price(element, 1, 3)                # (7, 1, 2): the cost and n_T, n_√T
word = matrix_word(element, 1, 3)   # "HQSHSXTXHQ", checked against the element
normal_form(word)                   # "HQSHtSHQ", at most two Cliffords per gate
assert bloch_image(word_unitary(word)) == bloch_image(element)  # up to a global phase
```

## Repository layout

```
include/quicksynth/  the C++ headers, among them the bindings (the nanobind
                     module quicksynth._core)
src/                 the C++ sources, in the same directories
python/quicksynth/   the Python package: QuickSynth, the exact synthesis algorithm,
                     SynthesisConfig and the command line interface
tests/               the C++ tests (tests/cpp/), the Python tests (tests/python/)
                     and the golden fixtures both read (tests/fixtures/)
docs/                the documentation
examples/            quickstart.py
```

## Documentation

* [docs/python-api.md](docs/python-api.md): every option of `SynthesisConfig`
  and how the search works.
* [RESULTS.md](RESULTS.md): preliminary benchmark figures, ahead of the paper.
* [docs/architecture.md](docs/architecture.md): how the library fits together:
  one call mapped to its files, the layout, the boundary between Python and
  C++, where each constant lives, what keeps results bit for bit, threads,
  interrupts and errors.
* [docs/development.md](docs/development.md): building, the wheels, the tests
  and the rules of their fixtures.
* [CONTRIBUTING.md](CONTRIBUTING.md): the development setup, the checks a
  change must pass and the continuous integration.

## Citation

If you use this code, please cite the paper; `CITATION.cff` holds the
reference, and GitHub shows it under "Cite this repository". The paper, with a
detailed description of the algorithm and of the results, is coming soon.

## License

MIT; see `LICENSE.txt`.

The extension module compiles in nanobind (BSD 3-Clause) and the
tsl::robin_map it includes (MIT). Each wheel carries their license texts in its
metadata directory under `licenses/bundled/`.
