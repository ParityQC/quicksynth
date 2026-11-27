<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# Fixture records

How the golden fixtures in this directory were recorded, and every change that
re-recorded one. [Golden fixtures](../../docs/development.md#golden-fixtures)
explains what each holds and how the tests replay it.

## Recorded fixtures

`golden_scan.json`, `golden_split.json` and `golden_combined.json`:

- The targets are the first 12 of `calibration_targets.json`. Each file stores
  them in its `meta`, with the git head, the numpy version and the CPU seconds
  of the recording.
- The scan and combined runs take ε from 1e-2 to 1e-6, and the first 4 targets
  also 1e-7. The split runs at 1e-3 and 1e-5 at the price (1, 3), with prefix
  scales 1 and 2, and at 1e-4 and 1e-6 on Clifford+T.
- The fixed-throughput variants pass 1.4e7 (√T) and 2.5e7 (T) explicitly.
  These were the library's reference throughputs before the calibration, so
  the fixtures do not follow the library constant.

## Frozen fixtures

`golden_words.json`, `golden_balls.json`, `golden_prefixes.json` and
`golden_synthesise.json` came from the Python implementation of the exact
layer, which computed matrices and distances with CPython and mpmath at 40
digits. That implementation is not in this repository.

- One recording wrote the exact layer, the balls and the prefix sets into one
  file, since split into the first three. The `cpu_seconds` of each is that of
  the whole recording.
- The ball at the price (1, 1) was recorded to radius 5 and is kept to
  radius 4.
- The `past_64_bits` section of `golden_words.json` was added later the same
  way. It holds three words whose coefficients pass 64 bits (74, 71 and 125
  bits), with their Bloch images, words and prices at five prices, and their
  T-counts.
- The `synthesise` runs used the Python word and check. The paper's runs
  30-122 were added later the same way. They include auto below 1e-3, a scan
  at (2, 5), four workers, a second-stage answer near a Clifford, and
  Cliffords, magic gates and short words on both gate sets down to 1e-8.
- `golden_synthesise.json` keeps the paper's 73 QuickSynth runs in their
  order, byte for byte except the runs at file indices 9, 13, 53, 55 and 57
  ([Re-recordings](#re-recordings)). `meta.paper_v1_runs` gives each run's
  number in the paper.

On macOS x86-64 the paper's run 107 (file index 63) answers with the
transposed word, which `test_golden.py` admits off Linux x86-64.

## Norm-equation counts

The last run of `check_relative_norm.py` used PARI/GP 2.17.2 and cypari2 2.2.4,
and no count of `golden_relative_norm.json` differed.

## Fingerprint cases

Where the cases of `fingerprint.json` come from:

| Cases | Source |
| --- | --- |
| 143 `golden-*` | every QuickSynth run of `golden_synthesise.json` under its recorded fixed stopping rule and, except its three slow runs, under the model stopping rule |
| 150 `c1-*`, `c1t-*` | the benchmark campaign, under the fixed stopping rule: 149 stage-1 calls (up to seven per setting) and one stage-2 call |
| 89 `window-*` | the evaluation of the model stopping rule: 18 calls that answered at stage 2 and 71 stage-1 calls |
| 3 `deploy-*` | the near-Clifford targets of `test_deploy.py`, one answered at stage 3 |
| 32 `symmetric-*` | rotations and other targets with a Clifford symmetry, on both gate sets |
| 10 `deep-*` | exact words moved by ε/3, at 1e-10 to 1e-13 |
| 17 `whole-*` | the first campaign call of each setting with `fine_levels=False`, the whole-scale search that preceded the finer levels |

Every case taken from a benchmark record returned the recorded word, cost and
stages when the fixture was first recorded, in the dev container. `meta` names
the number of cases and the commit of the last full recording (`record`). A
re-record of some cases in place keeps that commit.

## Re-recordings

Every change that re-recorded a fixture, oldest first. Point counts are the
per-walk point counts of the search log. Another element or word has the same
cost unless the row gives costs.

| Change | Recorded fixtures | `golden_synthesise.json` | `fingerprint.json` |
| --- | --- | --- | --- |
| QuickSynth's own LLL (`lll.cpp`) replaced FLINT's and reduces some bases differently | 126 runs for their point counts, 9 for another element | file indices 13, 53 and 57: another word; run 13's distance is the new word's, at 60 digits | full re-record: point counts of 120 cases, another word in 16 |
| Refit of the stopping rule's time model and reference throughputs | none | none | full re-record: the search of 29 model-rule cases changed, and two cost less (`golden-016-model` 28 to 27, `window-13-1e-04-78` 28 to 26) |
| Finer levels became the default on Clifford+√T | none | none | the 192 cases whose search changed, in place: 32 with another word (cost 49 lower in all, two dearer); the 17 `whole-*` cases added with the answers recorded before |
| Fix of the full test run (#4): one case's point counts had moved by up to 9%, past the test's 5% tolerance | none | none | `symmetric-rz0.1-sqrt_t-1e-07`'s point counts, in place; its answer is unchanged |
| Tighter ellipsoid around the cap (#9, `search_region.cpp`): walks visit a median of 8.7 times fewer leaves | 776 runs for their point counts (704 of the scan, all 72 of the split), 107 for another element (78 of the scan, 1 of the split, 28 combined) | file indices 9, 13, 53 and 55: another word, the distances the new words', at 60 digits | full re-record: point counts of 291 cases, another word in 37; no cost, scale, estimate or certificate changed |
| Double-double set-up of the cap below ε ≈ 3.5e-7 (#14, `search_region.cpp`): walks at 1e-7 visit 9.1 times fewer leaves | the 48 scan runs at 1e-7 for their point counts; no element changed | none | the 91 cases at 1e-7 to 1e-13 whose point counts changed, in place: 3 with another word |
| Per-prefix pricing of the split in LevelSearch (#18, `stopping_rule.cpp`, `level_search.cpp`): after the first answer the split walks batches of prefixes, each to its own level | none | none | the 4 fast cases whose search changed, in place: `golden-119-fixed`, `golden-119-model`, `deep-HTSHTHT-sqrt_t-1e-13` and `deep-HQSHTQHSQ-sqrt_t-1e-12`, one split walk a scale shallower on a batch of prefixes; no answer changed. Recorded with `tools/fingerprint.py update` on macOS arm64, where these differences equal those on Linux x86-64; the slow cases are still to be re-recorded on Linux |
