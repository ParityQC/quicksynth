<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# Architecture

QuickSynth is a C++17 library with a thin Python package on top.

- `include/quicksynth/` and `src/` hold the library. It contains no Python and
  builds as the static library `quicksynth_core`. Each header opens with its
  concept, the paper's label for it and its contract. This page shows how the
  parts fit together, not how each algorithm works.
- The bindings in `include/quicksynth/bindings/` and `src/bindings/` build the
  nanobind module `quicksynth._core`. It exposes the exact layer and the entry
  point, QuickSynth, with its policy.
- `python/quicksynth/` checks the arguments, replaces the target by its nearest
  unitary in numpy, calls `_core` and wraps the results.

The tests compare the library with golden fixtures, the stored results of
earlier versions ([development.md](development.md#golden-fixtures)).

## One call, step by step

A call searches over denominator scales. A target that is a z-rotation up to
Cliffords goes first to [stage 0](#stage-0), which walks only the matrix's
diagonal entry and solves for the other one. Otherwise, or if stage 0 has no
answer, the stages run. In each stage the scan picks the scales, the
enumerator searches one scale, and the split reaches deeper scales through
short prefixes. Only the winner becomes a word, and it is checked once, at the
end.

The right-hand column names each part's source under `src/` ([Layout](#layout)).

```
synthesise(U, eps)                                       synthesise.cpp
 |- target -> nearest unitary                            python/quicksynth/api.py
 |- search at eps - 1e-15                                synthesise.cpp
 |- STAGE 0, if U = C1 D C2 within the tolerance         synthesise.cpp
 |    DIAGONAL FORM: Cliffords C1, C2, diagonal D        synthesise_rotation/diagonal_form.cpp
 |    for each level (K, sub) and class:                 synthesise_rotation/rotation_search.cpp
 |      ROTATION BODY: the entries u near D              synthesise_rotation/rotation_body.cpp
 |      NORM EQUATION t t* = 2^K - u u*, exactly         solve_norm_equation/relative_norm.cpp,
 |                                                       prime_ideals.cpp, integer_factor.cpp,
 |                                                       prime_field.cpp
 |      each completion within eps: PRICE                price_candidates/price.cpp
 |    STOP: level bound (Clifford+T, certified),
 |          count rule (Clifford+sqrt(T))
 |    winner -> C1 + WORD + C2, accept,                  exact_synthesis/descent.cpp
 |              else fall through to the stages
 '- for each stage, until one answers                    policy.cpp
     scan_then_split                                     search_scales/scan_then_split.cpp
      |- SCAN over the denominator scales K              search_scales/scan.cpp
      |    for each K and each determinant class:
      |      ENUMERATOR(K, class)                        enumerate_lattice/enumerator.cpp
      |        SEARCH REGION: the scale-K points         enumerate_lattice/search_region.cpp
      |                       within eps of the target
      |        LLL-reduce its basis                      enumerate_lattice/lll.cpp
      |        Fincke-Pohst walk inside it               enumerate_lattice/walk.cpp
      |        each survivor: exact ring check, PRICE    price_candidates/price.cpp,
      |                                                  data_structures/bloch.cpp,
      |                                                  exact_synthesis/descent.cpp
      |                                                  (cost only)
      |    STOPPING RULE: walk the next scale or stop    search_scales/stopping_rule.cpp
      |    SCALE BOUND: no deeper scale can be cheaper   search_scales/scale_bound.cpp
      |                 -> stop; the answer is certified
      |- SPLIT, unless the scan certified                search_scales/split.cpp
      |    prefixes P of scale k, cheapest first         data_structures/prefix_sets.cpp
      |    x remainders R: ENUMERATOR(K - k) around P^-1 V
      |- (finer levels on Clifford+sqrt(T): SCAN and      search_scales/level_search.cpp
      |   SPLIT interleaved level by level; on by default)
      '- the cheaper of the two answers
     winner -> WORD (once)                               exact_synthesis/descent.cpp,
                                                         exact_synthesis/normal_form.cpp
            -> normal form, distance <= eps - 5e-16:     exact_synthesis/normal_form.cpp,
               accept, else the next stage               data_structures/distance.cpp
synthesise, on Clifford+sqrt(T): the Clifford+T          synthesise.cpp
     search (stage 0, the stages) below the answer's
     cost; a cheaper pure-T word replaces it
```

- The search region of ε and a scale K holds the scale-K points within ε of the
  target. The enumerator encloses it in an ellipsoid, so every candidate is a
  lattice point inside the ellipsoid.
- The enumerator finds every element of one scale and determinant class within
  ε. It LLL-reduces the region's basis, walks the points inside the ellipsoid,
  and checks and prices each survivor exactly. Later scales and the split's
  remainders reuse its cached set-up.
- The scan starts at the first scale that can hold an element within ε and
  walks on while the [stopping rule](#stopping-rule) allows. No element of scale
  K or deeper costs less than t (2K - 3) on Clifford+T, or min(t/2, q/3) (4K - 7)
  on Clifford+√T (`scale_bound.hpp`). Once that bound reaches the best cost, the
  scan stops, and its answer is optimal.
- The split factors an element of scale K as P R, with P of scale k, R of scale
  K - k and d(R, P^-1 V) = d(U, V). For the cheapest prefixes P it runs the
  enumerator at scale K - k around P^-1 V. It is a heuristic, since an element
  need not factor through a prefix of the chosen scale. Unless the scan
  certified, the model rule always splits, and the fixed rule only below
  `split_below` (1e-5) or when the scan found nothing.
- The walk prices every candidate without building its word. On Clifford+T the
  price is the T-count, half the Bloch exponent. On Clifford+√T it is the weight
  of the normal form, which the descent reads off the Bloch image.
- For the winner the descent runs again and writes the letters. `accept`
  rewrites the word in normal form, multiplies it out exactly and checks its
  distance to the target ([Bit for bit](#bit-for-bit)).

On Clifford+√T a Clifford+T search follows the answer. A pure-T word lies at
about half its T-count in scale, far deeper than a √T word of the same cost, so
the √T search stops before it. The Clifford+T search stops where t (2K - 3)
reaches the answer's cost, usually on its first scales. A cheaper pure-T word
replaces the answer with `source` `"clifford_t"`. This happens at a high √T
price such as (1, 6).

## Layout

The headers are in `include/quicksynth/` and the sources in `src/`: the
bindings, the entry, one directory per step of a call, and two foundations. A
name without an extension stands for a header and its source.

```
include/quicksynth/, src/
  bindings/                             the module quicksynth._core
    module.cpp  python_conversions  exact_bindings.cpp
    synthesis_bindings.cpp
  synthesise  policy                    1. the entry: stage 0, the stages, the
                                           policy, the cached balls
  synthesise_rotation/                  2. stage 0, the z-rotation path: targets
    diagonal_form  rotation_body           Clifford-equivalent to a z-rotation,
    rotation_search                        the diagonal entries u near one, the
                                           level search
  solve_norm_equation/                  3. the relative norm equation t t* = xi,
    relative_norm  prime_ideals            from the prime ideals above the
    integer_factor  prime_field            factors of N(xi) and the factors of
                                           the cyclotomic polynomial mod p
  search_scales/                        4. which scales to search
    scan_then_split  level_search  scan  stopping_rule  scale_bound  split
  enumerate_lattice/                    5. every element of one scale within
    enumerator  search_region  lll         epsilon
    setup  walk
  price_candidates/                     6. the cost of a candidate, without its
    price                                  word
  exact_synthesis/                      7. the descent and the word
    descent  normal_form
  data_structures/                      foundation: the data structures the
    cyclotomic_matrix  int_mp              steps share: the exact matrices over
    gate_set  clifford  su2  cost.hpp      the ring, the integers of any size,
    bloch  cost_ball  prefix_sets          the gate sets, the Cliffords, SU(2),
    distance  double_double                Price and Cost, the Bloch image, the
    number_field                           cost ball, the split's prefix sets,
                                           the final distance, double-double
                                           numbers and the cyclotomic rings
  support/                              foundation: the host's hooks, the
    host  parallel  diagnostics            search log and summary and the
    clock.hpp  lru.hpp                     verbose walk line, the clock, the
                                           thread pool, the LRU cache
```

A file includes a header as `"quicksynth/<directory>/<file>.hpp"`, or as
`"quicksynth/<file>.hpp"` at the top level. What it may include depends on its
directory:

| A file in | May include |
| --- | --- |
| `bindings/` | every directory |
| the top level (the entry) | every directory but `bindings/` |
| a step, 2 to 7 | its own directory, the later steps and the foundations |
| `data_structures/` | itself and `support/` |
| `support/` | itself |

Only the bindings include nanobind. The Pricer of `price_candidates/` calls the
descent's cost walk, and the entry calls the descent's word walk for the
winner. Since `include/` is the library's only include path, the compiler
cannot enforce these rules, so `tests/python/test_layers.py` reads the
`#include` lines and does. The C++ tests mirror the layout
([C++ tests](development.md#c-tests)).

## The boundary between Python and C++

The bindings expose the exact layer and the entry point. The exact layer
(`exact_bindings.cpp`) is `word_matrix`, `bloch_image`, `price` (with
`t_count`), `matrix_word` and `normal_form`, and the package wraps `word_matrix`
as `word_unitary`. The entry point (`synthesis_bindings.cpp`) is `QuickSynth`,
with `synthesise` and `prepare`, and its result `Synthesised`. The module also
exports the policy constants that `api.POLICY` reads, and
`reference_throughput`.

Data crosses as builtins and numpy arrays (`python_conversions.hpp`). A target
is a C-contiguous `(2, 2)` complex128 array, which the core takes to be
unitary. An exact element is `k` and the 32 hexadecimal coefficients of its
entries, row-major, each over 1, ζ, ..., ζ^7. Gate sets and stopping rules cross
as names, and absent settings as `None`. Every search and build releases the
GIL for its whole duration.

Python keeps the argument checks, `SynthesisConfig`, `SynthesisResult` with its
QASM output and the command line, and repeats no default of the core. It also
keeps `nearest_unitary`, which replaces the target by its polar factor from
numpy's SVD. Every later step starts from its last bits, which depend on the
LAPACK build, so it is not ported. `CyclotomicMatrix.to_numpy` sums an element's
matrix with the core's operations in the core's order, so its bits are the
core's.

## Constants

| Constant | Value | Role | Home |
| --- | --- | --- | --- |
| `kEpsilonMargin` | 1e-15 | taken off ε before every search | `policy.hpp` |
| `kDistanceError` | 5e-16 | the final check's rounding allowance | `data_structures/distance.hpp` |
| `kTargetUnitarity` | 1e-12 | the final distance's unitarity tolerance | `data_structures/distance.hpp` |
| `kPrefixRadius` | 6 | prefix-ball radius for ε >= 1e-7 | `policy.hpp` |
| `kSplitBelow` | 1e-5 | the fixed rule splits below it | `search_scales/scan_then_split.hpp` |
| `kMaxPrefixes` | 600 | prefixes of stage 1's split | `search_scales/scan_then_split.hpp` |
| `kDeepEpsilon` | 1e-7 | below it, deeper prefixes and larger budgets | `search_scales/scan_then_split.hpp` |
| `kDeeperEpsilon` | 1e-8 | below it, Clifford+T splits on every coset | `search_scales/scan_then_split.hpp` |
| `prefix_scale_for`, `split_budget_for` | 2, 3e7; below 1e-7: 3 (Clifford+T 4, 6 below 1e-8), 3e8 | stage 1's prefix scale and split budget | `search_scales/scan_then_split.cpp` |
| `kDeeperLastBudget` | 1.5e10 | last split budget, Clifford+T below 1e-8 | `policy.hpp` |
| `kDefaultRate`, `kFastRate` | 0.3, 3 per second | the rate rule, default and `fast` | `search_scales/stopping_rule.hpp` |
| `kFineDefaultRate`, `kFineFastRate` | 0.3, 24 per second | LevelSearch's rate, default and `fast` | `search_scales/stopping_rule.hpp` |
| `kExtraScales`, `kPointBudget` | 2, 3e6 | scales past the first answer, points per scan walk | `search_scales/scan.hpp` |
| `kSplitPointBudget` | 3e7 | the split's budget for ε >= 1e-7 | `search_scales/split.hpp` |
| `rotation_tolerance` | 1e-2 | stage 0's off-axis allowance, in units of ε | `synthesise.hpp` (`Settings`) |
| `kRotationSolvable`, `kRotationSolvableFast` | 256, 40 | stage 0's count rule, default and `fast` | `synthesise_rotation/rotation_search.hpp` |
| `MIN_EPSILON`, `UNITARITY_TOLERANCE` | 1e-13, 1e-6 | Python's argument checks | `python/quicksynth/api.py` |

Rates count magic states saved per second. The calibration and limit constants
live, with the measurements behind them, in `search_scales/stopping_rule.hpp`,
`level_search.hpp` and `split.hpp`, `enumerate_lattice/enumerator.hpp`,
`synthesise_rotation/rotation_search.hpp` and `solve_norm_equation/*.hpp`.

Three tests pin the defaults. "the defaults are the values the fixtures were
recorded at" (`tests/cpp/search_scales/test_scan_then_split.cpp`) pins the
scan's, the split's and the enumerator's. "the policy reaches the stage
settings" (`tests/cpp/test_policy.cpp`) and `test_the_policy_is_read_from_the_core`
(`tests/python/test_api.py`) pin the policy in the stages and in Python.

## The entry point

`synthesise.hpp` holds QuickSynth. On Clifford+T it takes the price with
`q = 3` (`Price::clifford_t`), which the T-only ball does not use.

### Stage 0

Stage 0 (`rotation_path`, on by default) runs first, on one thread. It takes a
target U = C1 D C2 up to phase, with Cliffords C1, C2 and D diagonal, when U's
off-axis part is at most `rotation_tolerance` times ε. A target farther off the
axis goes to the stages, which are far slower on near-diagonal targets at small
ε: tens of seconds at 1e-8, and often no answer at 1e-12.

A candidate near D is [[u, -t* w], [t, u* w]] / √2^K, and its distance depends
on u alone [Ross, Selinger, arXiv:1403.2975]. For each level the search walks
the u near D and solves the norm equation t t* = 2^K - u u* exactly.

On Clifford+T the search stops once the next level's T-count bound reaches the
best cost. On Clifford+√T that bound is too weak, so the search stops once
`rotation_solvable` u (256, or 40 with `fast`) had a solvable equation. A walk
cut at `kRotationPointLimit` or `kRotationWalkSolvable` loses the certificate.
`certified` holds only on Clifford+T, when the level bound stopped the search
with nothing cut and the target is diagonal up to rounding.

With no candidate, the call falls through to stage 1. A Clifford+√T candidate
found before a level or scale cap is an upper bound: it is the answer if the
stages find nothing cheaper or `time_limit` stops them. A candidate found when
`time_limit` stops stage 0 is the answer at once.

### The stages

When a stage returns nothing, or its answer fails the check, the next stage
runs with larger budgets and deeper prefixes. Most calls answer in the first.
`stages` (`policy.hpp`) lists them.

| Stage | Prefix scale | Prefixes | Split budget | Scan point budget |
| --- | --- | --- | --- | --- |
| 1 | `prefix_scale_for(ε)` | 600 (all on Clifford+T below 1e-8) | `split_budget_for(ε)` | 3e6, or `point_budget` |
| 2 ... | each deeper scale up to 4 (on Clifford+T 5 below 1e-7, 7 below 1e-8) | all | 3e8 | 3e6 |
| last | `prefix_scale_for(ε)` (8 on Clifford+T below 1e-8) | all | 3e9 (1.5e10 on Clifford+T below 1e-8) | 1e9 |

A call has 4 stages for ε >= 1e-7 and 3 below. Each stage splits with the
prefix sets of its scale and every shallower one, except on Clifford+T below
1e-8, where it takes only its own scale's cosets. Each stage:

1. polls the interrupt check and the time limit,
2. takes the prefix ball of `split_radius` levels, from the cache or built,
3. runs `ScanThenSplit` on the target,
4. under `verbose`, prints `stage i/n {changes}: cost C`,
5. writes the answer as a word with `Descent::word` and passes it to `accept`.

Before each stage the call checks `time_limit`. Once it has passed, the call
returns stage 0's word if there is one. Otherwise it raises `TimeoutError`, from
stage 2 on. A running stage always finishes.
`prepare(ε)` builds every ball and prefix set a call at ε uses.

### Stopping rule

The stopping rule decides how far past its first answer a search walks. Stage 1
runs the configured rule (`"fixed"` by default), and later stages the fixed
rule. Stage 1 also falls back to the fixed rule where the model rule has no
model: for a price with too many cost levels and for a target with a Clifford
symmetry (`stopping_fallback`).

On whole scales the fixed rule walks a further scale while `kScaleSaving` times
the best cost exceeds `rate` times the scale's predicted time at
`reference_throughput`. The scan walks at most `extra_scales` scales past the
first answer. Below `split_below` the split walks `kFixedSplitScales` scales.

On Clifford+√T the fixed rule searches finer levels by default, four per scale
(`fine_levels`). LevelSearch (`level_search.hpp`) interleaves the scan and the
split level by level. After the first answer, each round values the scan's next
levels and, for each of the split's prefixes, the composite levels past its own
last, at the expected saving minus the rate times the predicted time. A prefix
saves the gates it reaches first: those whose normal form starts with it and
with no longer prefix of the set. The round walks the scan's best level or the
split's batch of the `kPrefixBatch` most valuable prefixes, each to its own
level, whichever is worth more, and stops when neither has positive value.

The model rule walks whole scales. Before each it asks whether the expected
saving, from the census of normal forms, is worth `rate` times the scale's time
under `time_model`. Where this one-step test is not optimal, a dynamic
programme over the deeper scales decides. Its split starts at the first scale
the scan did not walk in full.

## Bit for bit

The same target and settings give the same word on one host and toolchain,
whatever the load and the thread count:

- The rate rule prices a scale's time at a fixed reference, never at the rate
  it measures (`reference_throughput`, `time_model`).
- The region, the basis, its centre and the LLL's Gram-Schmidt are fixed
  sequences of correctly rounded operations. Fused multiply-adds would change
  them, so both CMake targets compile with `-ffp-contract=off`.
- The walk order is fixed, and the scan keeps the first of the cheapest hits. A
  parallel walk returns the serial walk's survivors in the same order. Stage 0's
  count rule counts solutions, not time.
- Elements, prices and words are exact. The final distance is computed in
  double precision as the searches compute it. Its error, measured at most
  4.4e-16, stays below `kDistanceError` = 5e-16, which the 1e-15 margin exceeds.
- A coefficient past 64 bits enters the matrix truncated to 53 bits
  (`int_mp::truncated_frexp`).

Results still depend on the host through the target's bits (numpy's SVD), and
on the wall clock through `time_limit` alone.

## Threads, interrupts and errors

### Threads

`workers` (default 1) is a call's thread count. The threads share the scan's
walk, the set-ups of a batch and the split's batch of walks, through one
process-wide pool (`support/parallel.hpp`). A loop started while another holds
the pool runs on its caller alone. The results do not depend on `workers`, but
the CPU time does, so take CPU-time figures at one worker.

Calls from several Python threads are safe. The searches are built per call,
and cached balls are only read. One mutex guards the prefix-ball caches and is
held while a ball is built, so concurrent calls do not build the same ball
twice.

### Interrupts and progress

`support/host.hpp` declares two hooks, which the module sets on import. The
interrupt check calls `PyErr_CheckSignals()`, so a pending `KeyboardInterrupt`
leaves the call, on the main thread only. The printer sends each `verbose` line
through `builtins.print`. The check runs between stages, stage-0 walks and the
levels of a ball build or prefix collection, never inside a walk. An
interrupted build is not cached.

### Errors

| Condition | Raised in | Python sees |
| --- | --- | --- |
| A price below 1, an unknown `gateset` or `stopping_rule`, `time_limit <= 0`, `workers < 1`, `rotation_tolerance` outside [0, 1), `rotation_solvable < 1` | `config.py` | `ValueError` |
| An unknown override keyword | `dataclasses.replace` | `TypeError` |
| A target that is not numeric, not 2x2, not finite, or not unitary within 1e-6 | `api.nearest_unitary` | `ValueError` |
| `epsilon` that is not a number, not in (0, 1), or below 1e-13 | `api` | `ValueError` |
| `t_cost > 2 q_cost` on Clifford+√T; a price of 2^31 or more | `synthesise.cpp`, `cost_ball.cpp`: `std::invalid_argument` | `ValueError` |
| The time limit passed before stage 2 or a later stage, and stage 0 left no word | `SynthesisTimeout` | `TimeoutError` |
| No stage answers | `std::runtime_error` | `RuntimeError` |
| An LLL transform, its inverse or a survivor too wide for the walk's integers | `std::overflow_error` | `OverflowError` |
| A failed internal consistency check | `std::logic_error` | `RuntimeError` |

The checks raise in this order: the configuration, the target, then ε.
`test_errors_keep_their_order_and_text` (`tests/python/test_deploy.py`) pins
the order and the messages.
