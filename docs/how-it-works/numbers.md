<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# QuickSynth in numbers

```{admonition} Preliminary results
:class: note

These numbers come from the benchmark campaign behind the paper, which is
coming soon, and may still change. They are the values plotted in the
[preliminary results](../../RESULTS.md). Every method runs on one thread per
rotation.
```

A rotation's **cost** is the number of magic states its circuit consumes, at a
price `(t, q)` per `T` and per `√T` gate. The prices compared are the gate count
`(1,1)`, power `(2,3)`, unfolding `(2,5)` and catalysis `(1,3)`; on Clifford+T the
cost is the T-count.

```{raw} html
<div class="qs-tiles" id="qs-tiles" aria-live="polite"></div>
```

## Cost against compile time

Each workload is compiled one rotation after another: a million Haar-random
rotations at ε = 10<sup>-7</sup>, or the distinct rotations of a Trotter, QAOA
or adiabatic-evolution circuit. Lower is cheaper, left is faster. On Clifford+√T
the small diamonds are cyclosynth at its other measured deadlines.

```{raw} html
<div class="qs-figure" id="qs-compile">
  <div class="qs-controls" role="group" aria-label="Gate set">
    <button type="button" data-gate-set="clifford_sqrt_t" aria-pressed="true">Clifford+√T, catalysis (1,3)</button>
    <button type="button" data-gate-set="clifford_t" aria-pressed="false">Clifford+T, T-count</button>
  </div>
  <div class="qs-controls" role="group" aria-label="Workload">
    <button type="button" data-column="haar" aria-pressed="true">10<sup>6</sup> random rotations</button>
    <button type="button" data-column="trotter" aria-pressed="false">Trotter</button>
    <button type="button" data-column="qaoa" aria-pressed="false">QAOA</button>
    <button type="button" data-column="adiabatic" aria-pressed="false">Adiabatic</button>
  </div>
  <div class="qs-legend"></div>
  <div class="qs-chart"></div>
  <p class="qs-caption"></p>
  <details class="qs-table"><summary>Show the numbers</summary><div></div></details>
</div>
```

## Cost against precision

The mean cost above the lower bound (top) and the CPU seconds per rotation
(bottom) over 40 Haar-random targets, from ε = 10<sup>-3</sup> to 10<sup>-8</sup>
on Clifford+√T, the finest precision cyclosynth's Clifford+√T search accepts,
and to 10<sup>-9</sup> on Clifford+T. The squares are the certified optimum of
the same targets, where an exhaustive search certifies all of them.

```{raw} html
<div class="qs-figure" id="qs-sweep">
  <div class="qs-controls" role="group" aria-label="Price">
    <button type="button" data-sweep="0" aria-pressed="false">gate count (1,1)</button>
    <button type="button" data-sweep="1" aria-pressed="false">power (2,3)</button>
    <button type="button" data-sweep="2" aria-pressed="false">unfolding (2,5)</button>
    <button type="button" data-sweep="3" aria-pressed="true">catalysis (1,3)</button>
    <button type="button" data-sweep="4" aria-pressed="false">Clifford+T, T-count</button>
  </div>
  <div class="qs-legend"></div>
  <div class="qs-chart"></div>
  <p class="qs-caption"></p>
  <details class="qs-table"><summary>Show the numbers</summary><div></div></details>
</div>
```

## Why √T saves cost

Every bit of precision costs a fixed amount on average: the lower bound on the
mean least cost grows linearly in log<sub>2</sub>(1/ε). Its slope follows from
how fast the number of distinct circuits grows with their cost. With √T there
are more circuits of each cost, so each bit is cheaper, most of all at the
gate-count price.

```{raw} html
<div class="qs-figure" id="qs-growth">
  <div class="qs-legend"></div>
  <div class="qs-chart qs-pair"></div>
  <p class="qs-caption"></p>
  <details class="qs-table"><summary>Show the numbers</summary><div></div></details>
</div>
<noscript><p>The interactive figures need JavaScript; the
<a href="../../RESULTS.html">preliminary results</a> show the same data.</p></noscript>
<script src="https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js"
        integrity="sha384-CjloA8y00+1SDAUkjs099PVfnY2KmDC2BZnws9kh8D/lX1s46w6EPhpXdqMfjK6i"
        crossorigin="anonymous"></script>
<script src="../../_static/showcase/numbers.js"></script>
```
