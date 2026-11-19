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

`QuickSynth(config=None, **overrides)` is the synthesiser, and
`.synthesise(target, epsilon)` accepts any 2x2 unitary; a matrix within 1e-6 of
unitary (float32 included) is replaced by the nearest unitary.
`synthesise(target, epsilon, config=None, **overrides)` is QuickSynth's one
call, `synthesise_rotation(theta, epsilon, axis="z", ...)` does the same for a
rotation, and `rotation(theta, axis)` builds its matrix. Keyword overrides are
fields of `SynthesisConfig`:

| Field | Default | Meaning |
| --- | --- | --- |
| `gateset` | `"sqrt_t"` | `"clifford_t"` restricts the words to Clifford+T |
| `t_cost`, `q_cost` | `1`, `3` | magic-state price of `T` and `√T` |
| `extra_scales` | `None` (2) | under the fixed stopping rule on whole scales, the most scales the scan walks past its first answer; no effect with `fine_levels` |
| `point_budget` | `None` (3e6) | the most estimated lattice points the scan walks per determinant class; larger scales are skipped |
| `workers` | `1` | threads for the scan's walk |
| `time_limit` | `None` | seconds after which QuickSynth starts no further escalation stage; stage 0 counts against it |
| `stopping_rule` | `"fixed"` | how far QuickSynth's first stage searches past its first answer: `"fixed"`, a fixed number of scales, or `"model"`, decided by a counting model of the gates within ε |
| `fast` | `False` | the rate rule at 3 instead of 0.3 magic states a second, in the scan and the split: on whole scales no walk past the first answer predicted to take more than about 1 s; with `fine_levels` at 6 instead of 0.2 |
| `fine_levels` | `None` (on) | on Clifford+√T, search four finer levels per denominator scale and, past the first answer, walk the scan or split level whose expected saving most exceeds 0.2 magic states a second (6 with `fast`) times its predicted time; `False` walks whole scales only, at 0.3 (3 with `fast`); no effect on Clifford+T; stages under the fixed stopping rule only |
| `rotation_path` | `True` | a target within `rotation_tolerance` · ε of a z-rotation up to Cliffords is searched first on the z-rotation path (stage 0); no answer falls through to the stages |
| `rotation_tolerance` | `0.01` | largest off-axis part, in units of ε, that stage 0 takes; it searches at ε less that part |
| `rotation_solvable` | `None` | stage 0's count rule on Clifford+√T: whole levels until this many diagonal entries had a solvable norm equation (256, or 40 with `fast`) |
| `verbose` | `False` | progress on stdout |

Stage 0 takes a target that is a z-rotation up to Cliffords, such as every
target of `synthesise_rotation`. It walks the diagonal entry u of the candidate
matrices level by level and solves the norm equation for the other entry
exactly, as in Ross and Selinger's z-rotation synthesis. On one thread at
ε = 1e-3 to 1e-10, on Clifford+T it returns the least T-count, certified for
an exact z-rotation, in about 1 ms; on Clifford+√T it stops by its count rule
after about 35 ms, or 6 ms with `fast` at a mean cost 0.2 to 0.8 higher.

QuickSynth's stopping rule decides how far its first stage searches past the
first answer. Under the fixed stopping rule, the default, on Clifford+√T the
search walks the finer levels: past the first answer it walks the scan or split
level whose expected saving most exceeds 0.2 magic states a second times its
predicted time, and stops when none does. On whole scales (Clifford+T, or
`fine_levels=False`) the scan walks at most `extra_scales` scales past its first
hit while a fixed fraction of the best cost saved is worth the time and, below
ε = 1e-5, the split walks a fixed number of composite scales. Under the model
stopping rule, which walks whole scales, a census counts the
gates of each denominator scale and cost, and the gates within ε of a target are
taken as Poisson with mean that count times the measure of the ε-ball. The scan
walks the next scale while its expected saving in magic states exceeds 0.3 times
its predicted CPU seconds, from a time model fitted on the calibration set, and
stops for good once no deeper gate can be cheaper than the answer in hand. The
split then walks each further composite scale with the prefix scale of the best
expected value. A target with a Clifford symmetry, which the model does not
describe, a price in the hundreds or more, whose model would not fit in memory,
and every later escalation stage run the fixed stopping rule. The model rule
returns cheaper circuits from ε = 1e-4 to 1e-7 at the catalysis price, at two to
three times the CPU time, and dearer ones at 1e-8.

On whole scales the split's last composite scale is the dearest walk, and its
size depends on where ε falls between two scales, so the time rises and falls
with ε: with `fine_levels=False`, at the catalysis price, 0.4 s a rotation at
ε = 1e-7 and 3 s at 2e-7. On Clifford+√T the default walks the finer levels
instead, so it stops between whole scales: on 40 random targets at 31
precisions it returns circuits 0.3 magic states cheaper on average than whole
scales, at 0.7 times their CPU time, and its time and cost vary far less from
one ε to the next. With `fast=True` a walk past the first answer must save more
than 6 magic states a second of its predicted time (3 on whole scales, the
split's from the model rule's time model), so the dearest levels are left out:
at 2e-7 about 0.4 s a rotation, at a cost about 1 higher than whole scales
without `fast`, and about 0.3 s at a cost about 2 higher with
`fine_levels=False`.

`SynthesisConfig` is frozen; `dataclasses.replace(cfg, **changes)` gives a
modified copy. `prepare(epsilon)` builds, ahead of a timed call, the cost
balls the split takes its prefixes from at `epsilon`. `certified` on the
result says whether the cost is proved minimal (on Clifford+T, where the scan
reaches its scale bound or stage 0 its level bound) or is an upper bound. The
cost balls depend only on the price; they are built in memory, within a
fraction of a second, the first time a process needs them, and nothing is
written to disk.

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
