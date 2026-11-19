# Python API reference

The README shows the calls; this page gives every field of `SynthesisConfig`
and how the search uses them.

## `SynthesisConfig`

Every keyword override of `QuickSynth`, `synthesise` and `synthesise_rotation`
is a field of `SynthesisConfig`:

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

`SynthesisConfig` is frozen; `dataclasses.replace(cfg, **changes)` gives a
modified copy. `prepare(epsilon)` builds, ahead of a timed call, the cost
balls the split takes its prefixes from at `epsilon`. `certified` on the
result says whether the cost is proved minimal (on Clifford+T, where the scan
reaches its scale bound or stage 0 its level bound) or is an upper bound. The
cost balls depend only on the price; they are built in memory, within a
fraction of a second, the first time a process needs them, and nothing is
written to disk.

## How the search works

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
