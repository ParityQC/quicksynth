<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# QuickSynth

Near-optimal synthesis of single-qubit unitaries over **Clifford+T** and
**Clifford+√T** at any magic-state price, in a C++ core with a Python interface.

Given a 2x2 unitary and a tolerance ε, QuickSynth returns a word over
`H, S, X, Y, Z, T, Q` and its inverses `s, t, q` (`Q` = √T) within ε in
diamond distance, up to a global phase, at the least magic-state cost
`t_cost·n_T + q_cost·n_√T` it can find. Every answer is checked before it is
returned, and the same input always gives the same word.

```python
import numpy as np
import quicksynth as qs

r = qs.synthesise_rotation(np.pi / 8, 1e-3)   # R_z(pi/8) within 1e-3
r.word, r.cost, r.n_t, r.n_q
```

::::{grid} 1 1 2 3
:gutter: 3

:::{grid-item-card} Getting started
:link: ../README
:link-type: doc

Installation, the command line and the Python API.
:::

:::{grid-item-card} Algorithm
:link: algorithm
:link-type: doc

How the search finds a circuit, in interactive figures.
:::

:::{grid-item-card} Benchmark
:link: benchmark
:link-type: doc

Cost and compile time against cyclosynth and pygridsynth.
:::

:::{grid-item-card} API reference
:link: api
:link-type: doc

Every function, class and option of the Python package.
:::

:::{grid-item-card} Developer guide
:link: developers
:link-type: doc

How the library is built, tested and released.
:::
::::

```{toctree}
:hidden:

Getting started <../README>
Algorithm <algorithm>
benchmark
api
developers
```
