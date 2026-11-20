<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# Architecture

QuickSynth is a C++17 library with a thin Python package on top.

- `include/quicksynth/` holds the library's headers and `src/` its sources.
  The library contains no Python and builds as the static library `quicksynth_core`.
  Each header opens with its concept, the paper's label for it and its contract;
  this document says how the parts fit, not how each algorithm works.
- The bindings, `include/quicksynth/bindings/` and `src/bindings/`, are the
  nanobind module `quicksynth._core`. It exposes the exact layer and the entry
  point, QuickSynth, with its policy.
- `python/quicksynth/` checks the arguments, replaces the target by its nearest
  unitary in numpy, calls `_core` and wraps the results.

The tests compare the library with golden fixtures in `tests/fixtures/`. These
JSON files store results of earlier versions: runs of the scan, the split and
`synthesise`, the exact layer, the cost balls and prefix sets, the stopping
rule's reference values, and the answers to a fixed set of calls
(`fingerprint.json`). A change that should not alter any answer must reproduce
them bit for bit. The same directory holds the Haar-random targets the tests
draw from. [development.md](development.md) covers building, the tests and the
rules of the fixtures.

## One call, step by step

A QuickSynth call is a search over denominator scales. A target that is a
z-rotation up to Cliffords first goes to stage 0, the z-rotation path, which
walks only the matrix's diagonal entry and solves for the other one. Otherwise,
or when stage 0 has no answer, the stages run: the scan decides which scales to
search; the enumerator searches one scale; the split reaches deeper scales
through short prefixes. The answer is written as a word and checked once, at
the end.

The right-hand column names the source of each part, a C++ file under `src/`.
Its directory is the step it belongs to ([Layout](#layout)): the entry at the
top level, then `synthesise_rotation/`, `solve_norm_equation/`,
`search_scales/`, `enumerate_lattice/`, `price_candidates/` and
`exact_synthesis/`; the data structures they share are in `data_structures/`.

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
synthesise, on Clifford+sqrt(T): the Clifford+T       synthesise.cpp
     search (stage 0, the stages) below the answer's
     cost; a cheaper pure-T word replaces it
```

- **Stage 0.** U = C1 D C2 up to phase, with Cliffords C1, C2 and D diagonal,
  when the Bloch image of U maps a coordinate axis to one (`diagonal_form`);
  the off-axis part of U, at most `rotation_tolerance` (1e-2) times ε, comes
  off the search's ε. A candidate near D is [[u, -t* w], [t, u* w]] / √2^K, and
  its distance depends on u alone [Ross, Selinger, arXiv:1403.2975]. The
  search walks the u of each level (scale K, sublattice δ^sub, δ = 1 - ζ) in
  the rotation body, solves t t* = 2^K - u u* exactly over Z[ω] or Z[ζ16] from
  the factors of its norm, and prices each completion through the enumerator's
  exact check. On Clifford+T it stops once the next level's T-count bound
  (T >= l - 2) reaches the best cost; an answer for an exactly diagonal target
  is then the certified optimum. On Clifford+√T that bound is about half the
  costs found, so a count rule stops it: whole levels until `rotation_solvable`
  u (256; 40 with `fast`) had a solvable equation. There the phases ζ^i of each
  t are priced only if they can beat the best (`rotation_search.hpp`). A walk
  past `kRotationPointLimit` points, or on Clifford+√T past
  `kRotationWalkSolvable` times the count rule's solvable u, is cut: the
  certificate is lost and the search goes on. Directions within a few ε of a
  multiple of π/4 hold such slabs of points (`rotation_body.hpp`). Caps on
  levels, scale and time end the search; with no candidate the call falls
  through to stage 1. A Clifford+√T candidate found before a cap is kept as an
  upper bound: the stages run, and it is the answer if they find nothing
  cheaper or `time_limit` stops them. A target off the axis by more than
  `rotation_tolerance` goes to the stages, which are far slower on
  near-diagonal targets at small ε (tens of seconds at 1e-8, often no answer
  at 1e-12) than stage 0 just inside the tolerance.
- **Search region.** For a precision epsilon and a denominator scale K, the
  search region holds the points of scale K within epsilon of the target. At
  the identity place these points lie in the cap of the unit 4-ball around the
  target; at every other conjugate place they lie on the unit sphere. The
  enumerator encloses the region in its Loewner ellipsoid, so every candidate
  is a lattice point inside the ellipsoid.
- **Enumerator.** The enumerator finds every element of one scale and one
  determinant class within epsilon. It LLL-reduces the region's basis in the
  ellipsoid's metric, walks all lattice points inside the ellipsoid, and checks
  and prices each survivor exactly in the ring. It caches the set-up of a target
  direction, so the later scales of a target and the remainders of a split
  reuse it.
- **Scan.** The scan starts at the first scale that can hold an element within
  epsilon and walks further scales while the stopping rule judges them worth
  their time. At the price (t, q), no element of scale K or deeper costs less
  than t (2K - 3) on Clifford+T or min(t/2, q/3) (4K - 7) on Clifford+sqrt(T)
  (`scale_bound.hpp`). Once the next scale's bound reaches the best cost found,
  the scan stops, and its answer is optimal.
- **Stopping rule.** The stopping rule decides how far past its first answer a
  search walks. By default, and in every later stage, the fixed stopping rule
  applies. On whole scales (Clifford+T, or `fine_levels=False`) it walks a
  further scale while its rate rule holds (`kScaleSaving` times the best cost
  exceeds `rate` times the scale's predicted time), at most `extra_scales`
  scales past the first answer; on finer levels LevelSearch decides instead. If
  the configuration asks for it, the first stage uses the model stopping rule
  instead, which walks whole scales and weighs the expected saving of a scale,
  from a census of the normal forms, against its predicted time.
- **Split.** An element of scale K factors as P R, with P of scale k and R of
  scale K - k, and d(R, P^-1 V) = d(U, V). The split takes the cheapest
  prefixes P and runs the enumerator at scale K - k around the rotated target
  P^-1 V, in a region 2^(deg k) times smaller. The split is a heuristic: an
  element need not factor through a prefix of the chosen scale. Unless the scan
  certified its answer, the split runs at every precision under the model
  stopping rule, and under the fixed stopping rule only below `split_below`
  (1e-5) or when the scan found nothing.
- **Finer levels.** On Clifford+√T under the fixed stopping rule, the scan and
  the split search four levels per denominator scale, the sublattices
  δ^sub with δ = 1 - ζ16, and LevelSearch (`level_search.hpp`) interleaves them
  level by level. This is the default; `fine_levels=False` keeps whole scales.
  Before the first answer the scan walks its levels in order and, without a
  scan answer, the split walks the levels past the scan's last; a level whose
  LLL set-ups would take longer than its walks is skipped. After the first
  answer each round prices every candidate level, scan or split, at its
  expected saving minus the rate times its predicted seconds (walks and new LLL
  set-ups), walks the candidate of largest positive value, and stops when none
  is positive. The rate is `kFineDefaultRate`, or `kFineFastRate` with `fast`;
  `extra_scales` does not apply. Candidates end at the point budgets and where
  a scale bound reaches the best cost, and split candidates after a scan
  answer come only below `split_below`.
- **Pricing.** The walk prices every candidate without building its word. On
  Clifford+T the price is the T-count, half the Bloch exponent. On
  Clifford+sqrt(T) the price is the weight of the normal form: the descent reads
  the normal form off the Bloch image, syllable by syllable, and counts the T
  and Q syllables.
- **Answer.** Only the winner becomes a word. The descent runs again, now
  writing the letters, and then appends the closing Clifford and merges
  letters. `accept` rewrites the word in normal form (`normal_form`),
  multiplies it out exactly and evaluates its diamond distance to the target
  (`precise_distance`) in double precision, from the exact element and the
  chord between the two SU(2) vectors, as the searches measure it. The
  measured error of this distance is at most 4.4e-16, so `accept` takes the
  answer when the distance is at most epsilon - 5e-16 (`kDistanceError`). The searches run at epsilon - 1e-15, a
  margin that exceeds this allowance. The normal form keeps the search's cost
  (asserted in `tests/cpp/test_synthesise.cpp`); on Clifford+T it has no `Q`.
- **Stages.** When a stage returns nothing, or its answer fails the check, the
  next stage runs with larger budgets and deeper prefixes (`stages()` in
  policy.cpp). Most calls answer in the first stage.
- **Clifford+T comparison.** A pure-T word of T-count n has scale about n/2,
  far deeper than a sqrt(T) word of its cost, so the sqrt(T) search stops
  before it. On Clifford+sqrt(T) the Clifford+T search therefore runs after the
  answer, stage 0 and then the stages, with the answer's cost as an incumbent:
  the scan and the split stop where t (2K - 3) reaches that cost, a scan
  stopped there leaves out the split and ends the search, and a cheaper pure-T
  word that passes the check replaces the answer (`source` `"clifford_t"`,
  `stages` of the Clifford+T search). At most prices the bound reaches the cost
  on the scan's first scales, and no prefix ball is built; at a high sqrt(T)
  price, such as (1, 6), the pure-T word wins.

## Layout

The headers are in `include/quicksynth/` and the sources in `src/`. A name
without an extension below stands for a header and its source.

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

A file includes its own directory, the steps after its own and the foundations,
as `"quicksynth/<directory>/<file>.hpp"` (`"quicksynth/<file>.hpp"` at the top
level). The bindings may include every directory, the entry every directory but
the bindings, a step never one before it, `data_structures/` only itself and
`support/`, and `support/` only itself. The Pricer of `price_candidates/` calls
the descent's cost walk, and the entry calls the descent's word walk for the
winner. Only the bindings include nanobind.

`include/` is the library's one include path, so the compiler does not enforce
these rules; `tests/python/test_layers.py` reads the `#include` lines and does.
The C++ tests mirror the steps:
`tests/cpp/<directory>/test_<concept>.cpp` tests
`include/quicksynth/<directory>/<concept>.hpp`, and the tests of the entry,
whose headers sit at the top level, are in `tests/cpp/` itself, beside the
helpers the tests share.

## The boundary between Python and C++

The bindings expose two things:

- **The exact layer** (`exact_bindings.cpp`): `word_matrix`, `bloch_image`,
  `price` (with `t_count` for the T-count), `matrix_word` and `normal_form`,
  which the package exposes under the same name.
- **The entry point** (`synthesis_bindings.cpp`): `QuickSynth` (`synthesise`,
  `prepare`), built from all its settings at once, with an absent value as
  `None`; its result `Synthesised`; and the policy's constants
  (`EPSILON_MARGIN`, `SPLIT_BELOW`, `RATE`, `PREFIX_RADIUS`, `PREFIX_COUNT`)
  and `reference_throughput`.

Data crosses in these forms (`python_conversions.hpp`):

- a target is a C-contiguous `(2, 2)` complex128 array, which the core takes
  to be unitary: the final distance refuses a target farther than
  `kTargetUnitarity` (1e-12, Frobenius norm of `U^† U - I`) from unitary with
  `ValueError`;
- an exact element is `k` and the 32 coefficients of the entries in
  hexadecimal (row-major, each over 1, ζ, ..., ζ^7), the form `bloch_image`,
  `price` and `matrix_word` take; `word_matrix` returns Python ints of any
  size;
- a gate set is its name, `"sqrt_t"` or `"clifford_t"`, and a stopping rule
  its name, `"model"` or `"fixed"`;
- results come back as builtins: `Synthesised`'s fields, its log as a list of
  dicts, its search summary as a dict.

Every search and build releases the GIL for its whole duration.

What stays in Python, in `python/quicksynth/`:

- **`nearest_unitary`.** The target must be a finite 2x2 matrix within 1e-6 of
  unitary (Frobenius norm of `U^† U - I`), and is replaced by its polar factor
  `W Vh` from numpy's SVD. Its last bits depend on the LAPACK build, and every
  search, distance and final check starts from them, so it is not ported.
- **`CyclotomicMatrix.to_numpy`.** An element's complex matrix, summed with
  the core's operations in the core's order (`matrix_unitary` in
  `data_structures/su2.cpp`, and `matrix_of` in `distance.cpp` past 64 bits),
  so its bits are the core's.
- The argument checks, `SynthesisConfig`, `SynthesisResult` (with QASM output),
  the command line, and the marshalling of `CyclotomicMatrix`. No default the
  core defines is repeated in Python: `api.POLICY` reads the core's constants.

## Constants

| Constant | Value | Home |
| --- | --- | --- |
| `kEpsilonMargin` | 1e-15 | `policy.hpp` |
| `kDistanceError` (the final check's allowance; `kEpsilonMargin` must exceed it) | 5e-16 | `data_structures/distance.hpp` |
| `kTargetUnitarity` (how far from unitary a target of the final distance may lie) | 1e-12 | `data_structures/distance.hpp` |
| `kPrefixRadius` | 6 | `policy.hpp` |
| `kSplitBelow` | 1e-5 | `search_scales/scan_then_split.hpp` |
| `kMaxPrefixes` (also the prefixes the split's ball must hold below `kDeepEpsilon`) | 600 | `search_scales/scan_then_split.hpp` |
| `kDeepEpsilon` (below it the split goes deeper: `prefix_scale_for`, `split_budget_for`, `split_radius`) | 1e-7 | `search_scales/scan_then_split.hpp` |
| `kDeeperEpsilon` (below it Clifford+T splits at prefix scale 6 with every coset of the scale: `split_radius` 2·scale·`t_cost`, no `max_prefixes` cap, no shallower fallback sets) | 1e-8 | `search_scales/scan_then_split.hpp` |
| `prefix_scale_for`, `split_budget_for` | 2 and 3e7; below ε = 1e-7, 3 (4 on Clifford+T, 6 below 1e-8) and 3e8 | `search_scales/scan_then_split.cpp` |
| `kDeeperLastBudget` (the last stage's split budget on Clifford+T below `kDeeperEpsilon`) | 1.5e10 | `policy.hpp` |
| `kDefaultRate`, `kFixedSplitScales`, `kScaleSaving` | 0.3 magic states a second, 3, 0.06 | `search_scales/stopping_rule.hpp` |
| `kFastRate` (`fast`: the scan's rate and the fixed split's, `split_rate`) | 3 magic states a second | `search_scales/stopping_rule.hpp` |
| `kFineDefaultRate`, `kFineFastRate` (LevelSearch's rate on Clifford+√T, where finer levels are on unless `fine_levels=False`: by default and under `fast`) | 0.2, 6 magic states a second | `search_scales/stopping_rule.hpp` |
| `walk_seconds` knots, `kFreshSetupSeconds`, `kDerivedSetupSeconds` (LevelSearch's walk and LLL set-up times) | 21 µs at 0.1 points to 123 ms at 1e6, slope 0.75 above; 200 µs, 60 µs | `search_scales/stopping_rule.cpp`, `.hpp` |
| `kLevelSplitCoverage`, `kCensusSplitRateDivisor`; `kDecayFirstSaving`, `kDecayPerLevel`, `kDecaySplitCoverage`, `kDecayScanRateDivisor`, `kDecaySplitRateDivisor` (LevelSearch's census and decay gains and rate shares) | 0.25, 4; 0.026, 0.8, 0.6, 2, 6 | `search_scales/stopping_rule.hpp` |
| `kSkipHeadroom` (LevelSearch skips a set-up-bound level only if the next is this far under the split budget) | 4 | `search_scales/level_search.hpp` |
| `kReferenceThroughputSqrtT`, `kReferenceThroughputCliffordT` (`reference_throughput`, the fixed stopping rule's) | 8.3e6 estimated points, 2.1e7 leaves a second | `search_scales/stopping_rule.hpp` |
| `kTimeModelSqrtT`, `kTimeModelCliffordT` (the model stopping rule's `TimeModel`): scan θ and o, split a and θ_split | Clifford+√T 9.61e6, 1.88e-3 s, 1.85e-4 s, 2.71e6; Clifford+T 1.77e7, 1.78e-5 s, 5.28e-5 s, 1.94e7 | `search_scales/stopping_rule.hpp` |
| `kMaxModelCosts` (above it a price gets the fixed stopping rule) | 16384 cost levels | `search_scales/stopping_rule.hpp` |
| `kExtraScales`, `kPointBudget` | 2, 3e6 | `search_scales/scan.hpp` |
| `kSplitPointBudget`, `kSplitTransformCache`, `kSplitScaleGuard` | 3e7, 8192, 64 | `search_scales/split.hpp` |
| `kTransformCache` | 4096 | `enumerate_lattice/enumerator.hpp` |
| `rotation_tolerance` (stage 0: the largest off-axis part, in units of ε) | 1e-2 | `synthesise.hpp` (`Settings`) |
| `kRotationSolvable`, `kRotationSolvableFast` (stage 0's count rule on Clifford+√T, default and `fast`) | 256, 40 | `synthesise_rotation/rotation_search.hpp`, `policy.hpp` |
| Lazy phase sweep (a phase ζ^i of a choice is priced only if phase 0's cost less this is below the best) | `2 q_cost` for odd i, `2 min(t_cost, 2 q_cost)` for i = 2; exact, no calibrated margin | `src/synthesise_rotation/rotation_search.cpp` |
| `kRotationLevelsPastHit`, `kRotationLevelsPastUseful` (stage 0 gives up this many levels past the later of the first solvable and the first useful level, and past the first useful one) | 16, 32 | `synthesise_rotation/rotation_search.hpp` |
| `kRotationPointLimit`, `kRotationWalkSolvable`, `kRotationMaxScale` (points of one walk before it is cut; on Clifford+√T, solvable u of one walk in multiples of the count rule's; the largest scale, coefficients below the 2^60 of `Enumerator::matrices`) | 2^20, 16, 118 | `synthesise_rotation/rotation_search.hpp` |
| `kChoiceCap` (norm-equation solutions per u; a cut list clears the certificate) | 64 | `solve_norm_equation/relative_norm.hpp` |
| `kRhoAttempts`, `kPrimeCacheCapacity` (Pollard-Brent attempts per cofactor; cached prime generators) | 32, 256 | `solve_norm_equation/integer_factor.hpp`, `prime_ideals.hpp` |
| `MIN_EPSILON`, `UNITARITY_TOLERANCE` | 1e-13, 1e-6 | `python/quicksynth/api.py` (the checks only) |

The C++ test "the defaults are the values the fixtures were recorded at"
(`tests/cpp/search_scales/test_scan_then_split.cpp`) pins the scan's, the
split's and the enumerator's defaults, the C++ test "the policy reaches the
stage settings" (`tests/cpp/test_policy.cpp`) the policy's in the stages, and
`test_the_policy_is_read_from_the_core` (`tests/python/test_api.py`) the
policy's as Python reads them.

## The entry point

`synthesise.hpp` holds QuickSynth. On Clifford+T it takes the price with
`q = 3` (`Price::clifford_t`), which the T-only ball does not use.

**Stage 0** (`rotation_path`, on by default) runs first, on one thread, and
builds no prefix ball. It answers with the word C1 + w + C2, w the descent's
word of the winner, through `accept`; the word is priced again and must cost
what the search found. `certified` holds on Clifford+T when the level bound
stopped the search with every candidate decided and no solution list cut, and
the target is diagonal up to rounding, and no walk was cut. Stage 0 checks
`time_limit` before each walk.

`stages(epsilon, gateset, first)` (`policy.hpp`) lists the escalation stages;
`QuickSynth::first_stage` sets stage 1 from the configuration:

| Stage | Prefix scale | Prefixes | Split budget | Scan point budget |
| --- | --- | --- | --- | --- |
| 1 | `prefix_scale_for(ε)` | 600 (all on Clifford+T below 1e-8) | `split_budget_for(ε)` | 3e6, or `point_budget` |
| 2 ... | each deeper scale up to 4 (on Clifford+T 5 below 1e-7, 7 below 1e-8) | all | 3e8 | 3e6 |
| last | `prefix_scale_for(ε)` (8 on Clifford+T below 1e-8) | all | 3e9 (1.5e10 on Clifford+T below 1e-8) | 1e9 |

A call has 4 stages for ε >= 1e-7 and 3 below. Each stage splits with the
prefix sets of its scale and every shallower one, except on Clifford+T below
1e-8, where it takes its scale's cosets only: the first answer there lies
about 5 scales deeper per decade of ε, and scale 4 at 3e8 stops short of it.
Every stage takes stage 1's rates (by default `kDefaultRate` for the scan and
the split, `kFineDefaultRate` for LevelSearch), `extra_scales` (2, or the
configuration's value; whole scales only) and `split_below` =
`kSplitBelow - kEpsilonMargin`. These settings follow the requested ε, not
the search target. Each stage:

1. polls the interrupt check and, from stage 2 on, raises `TimeoutError` once
   `time_limit` has passed since the call began (a running stage finishes); if
   stage 0 left a word, that word is returned instead, from stage 1 on;
2. takes the prefix ball of `split_radius` levels (`prefix_ball`), and builds
   it if the in-memory cache does not hold it;
3. runs `ScanThenSplit` on the target;
4. under `verbose`, prints `stage i/n {changes}: cost C`, where `changes` are
   the settings that differ from stage 1;
5. writes the answer as a word with `Descent::word` and passes it to `accept`;
   on Clifford+T, `certified` is the scan's certificate.

`prepare(ε)` builds every ball and prefix set a call at ε uses, so that a timed
call pays for none.

**Stopping rule.** QuickSynth decides the stopping rule once per call. Stage 1
runs the configuration's rule (`"fixed"` by default) and every later stage runs
the fixed rule. `stopping_model` builds the model of a (gate set, price, ε)
once and keeps the 32 most recently used. Two cases fall back to the fixed
rule: a price whose cost cap would pass `kMaxModelCosts` (16384), for which
`stopping_model` returns no model, and a target with a Clifford symmetry
(`stopping_fallback(U)`).

- Under the model rule, the scan asks `StoppingModel::walk` before each
  further scale whether the scale's expected saving, from the census N(K, c)
  and the measure of the ε-ball, is worth `rate` times its time under
  `time_model(gateset)`. Where the one-step test is not optimal, a dynamic
  programme over the deeper scales decides. The scan stops for good at
  `integer_scale_bound`, the least cost of any deeper gate. The split then
  starts at the first scale the scan did not walk in full and walks each
  composite scale with the prefix scale of the best expected value. It uses the
  stage's prefix sets and those of `stopping_extra_scales` (scales up to 3, from
  `QuickSynth::stopping_ball`).
- Under the fixed rule on whole scales, the scan walks at most `extra_scales`
  scales past the first hit, under the rate rule at `reference_throughput`, and
  below `split_below` the split walks `kFixedSplitScales` composite scales. On
  finer levels LevelSearch walks instead (above).

## Bit for bit

The same target and settings give the same word on one host and toolchain,
whatever the load and the thread count.

- **The rate rule's clock.** The rule prices a scale's time at a fixed
  reference, never at the rate it measures: `reference_throughput` under
  the fixed stopping rule and `time_model` under the
  model stopping rule, so where the scan stops depends on the target alone.
- **The lattice inputs.** The search direction, the region, the basis and its
  centre are fixed sequences of correctly rounded operations, and so is the
  LLL's floating-point Gram-Schmidt, so every host reduces the same integers
  to the same basis. Both targets compile with `-ffp-contract=off`, since an
  inline function compiled into both keeps one copy at link time. The SU(2)
  functions stay out of line, and the look-alike float routines `ldexp` and
  `pow` are never merged. The bits of the target (numpy's SVD) are those of
  the host.
- **The order.** The walk order is fixed and the scan keeps the first of the
  cheapest hits; a parallel walk returns the serial walk's survivors in the
  same order; each split walk has its own leaf cap. Stage 0 keeps the first of
  its cheapest candidates, and its count rule counts solutions, not time.
- **Exact arithmetic.** Elements, prices, words and spellings are exact, and
  the final distance is the diamond distance in double precision, computed
  as the searches compute it, from the exact element. Its error is absolute, typically about 1e-16 and
  measured at most 4.4e-16 at every distance; the final check allows
  `kDistanceError` = 5e-16 for it, and the searches' margin of 1e-15 exceeds
  that (`data_structures/distance.hpp`). A coefficient past 64 bits enters the
  matrix truncated to 53 bits (`int_mp::truncated_frexp`) and through `ldexp`,
  and past 1000 bits every coefficient is scaled by a power of two, exactly,
  so a word of any length gives a finite matrix.

Results still depend on the wall clock through `time_limit` alone.

## Threads, interrupts and errors

**Threads.** `workers` (default 1) is a call's thread count. It splits the
scan's walk over its top coordinates (`walk`), the set-ups of a batch
(`SetupCache`) and the split's batch of walks (`walk_batch`). Ball builds and
prefix collection run on the calling thread.
All of it runs through `parallel_for` on one process-wide pool
(`support/parallel.hpp`): a loop started while another holds the pool, nested
or from a concurrent call, runs on its caller alone, and the first exception of
a loop is rethrown once every worker has returned. Each part returns what the
serial code returns. The CPU time does depend on `workers`, so CPU-time figures
are taken at one worker.

Calls from several Python threads are safe. The searches are built per call;
cached balls are only read. The prefix-ball caches of `synthesise.cpp` sit
behind one mutex, held while a ball is built, so a ball is built at most once
per key while it stays cached. `CostBall::prefix_sets` has a mutex of its own and
collects without it, since the collection polls the interrupt check, which
takes the GIL.

**Interrupts and progress.** `support/host.hpp` declares two hooks, which the
module sets on import:

- The interrupt check takes the GIL and calls `PyErr_CheckSignals()`; a pending
  `KeyboardInterrupt` throws and leaves the call. It is polled before each
  QuickSynth stage, each walk of stage 0, each level of a ball being built,
  each radius `prefix_radius` tries and each level of a prefix collection, not
  inside a walk, a search or one level's build. A stopped step leaves nothing
  behind: an interrupted build is not cached, interrupted prefix sets are not
  kept. Signal handlers run on the main thread only, so a call from another
  thread is not interrupted.
- The printer sends each `verbose` line through `builtins.print`: one per walk
  of the scan and one per QuickSynth stage. A pending signal handler can run in
  `print`, so under `verbose` an interrupt can also leave from the printer.

**Errors.**

| Condition | Raised in | Python sees |
| --- | --- | --- |
| A price below 1, an unknown `gateset` or `stopping_rule`, `time_limit <= 0`, `workers < 1` | `config.py` | `ValueError` |
| An unknown override keyword | `dataclasses.replace` | `TypeError` |
| A target that is not numeric, not 2x2, not finite, or not unitary within 1e-6 | `api.nearest_unitary` | `ValueError` |
| `epsilon` that is not a number, not in (0, 1), or below 1e-13 | `api` | `ValueError` |
| `t_cost > 2 q_cost` on Clifford+√T; a price of 2^31 or more | `synthesise.cpp`, `cost_ball.cpp`: `std::invalid_argument` | `ValueError` |
| The time limit passed before stage 2 or a later stage, and stage 0 left no word | `SynthesisTimeout` | `TimeoutError` |
| No stage answers | `std::runtime_error` | `RuntimeError` |
| An LLL transform, its inverse or a survivor too wide for the walk's integers | `std::overflow_error` | `OverflowError` |
| A failed internal consistency check | `std::logic_error` | `RuntimeError` |

The checks raise in this order: the configuration, the target, then ε;
`test_errors_keep_their_order_and_text` (`tests/python/test_deploy.py`) pins
the order and the messages.

## Build

`CMakeLists.txt` defines `quicksynth_core`, the static library (one source list
per step, compiled with `-O3`), and `_core`, the nanobind module (`NB_STATIC`,
compiled with `-Os` by nanobind), both with `-ffp-contract=off -Wall -Wextra`.
scikit-build-core runs it on `pip install .`, `uv sync` and `uv build`. The
editable install imports `python/quicksynth/` in place and installs the
compiled module into site-packages; it does not rebuild on import, so after an
edit under `include/` or `src/`, `make cpp` (`uv sync --all-extras
--reinstall-package quicksynth`) rebuilds it in `build/<wheel tag>/`. The option
`QUICKSYNTH_TESTS` (off for the wheel) adds the C++ tests of `tests/cpp/`.
