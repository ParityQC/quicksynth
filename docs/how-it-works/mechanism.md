<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# How QuickSynth finds a circuit

A quantum computer applies a rotation as a circuit of gates. The Clifford gates
`H` and `S` are cheap; the magic gates `T` and `√T` are expensive, each consuming
a magic state. QuickSynth looks for the cheapest circuit within a distance ε of
the rotation. These figures show the idea behind the search, with data from the
library itself.

## Every circuit is a point

Each dot is a state that some circuit of at most the chosen cost prepares from
`|0⟩`; darker dots are cheaper. The circle around the target holds the states
within ε of it, and the ringed dot is the cheapest circuit inside: the answer.
Raise the cost and the sphere fills up, and with `√T` it fills much faster, so
the cheapest dot inside the circle is cheaper. Drag a sphere to turn it, and
click it to move the target.

```{raw} html
<div class="qs-figure" id="qs-sphere">
  <div class="qs-sliders">
    <label>circuits up to cost <output data-out="cost"></output>
      <input type="range" data-input="cost" min="0" max="9" step="1" value="9"></label>
    <label>distance ε <output data-out="eps"></output>
      <input type="range" data-input="eps" min="0" max="6" step="1" value="3"></label>
    <div class="qs-controls"><button type="button" data-action="target">New target</button></div>
  </div>
  <div class="qs-ramp" aria-hidden="true"></div>
  <div class="qs-chart qs-spheres"></div>
  <p class="qs-caption"></p>
  <details class="qs-table"><summary>Show the numbers</summary><div></div></details>
</div>
```

QuickSynth matches the whole rotation, not only the state it prepares: a point in
the three-dimensional space of rotations instead of on the two-dimensional sphere.
The picture is the same one dimension up.

## Scale by scale

The matrix of every such circuit has entries of the form *x* / √2<sup>K</sup>, with
*x* in a ring of algebraic integers; *K* is the circuit's **scale**. At one scale,
the circuits near the target are the points of a lattice in a small region, and
QuickSynth walks all of them exactly. Cheap circuits have small scales, so the
search goes scale by scale. Two numbers decide when it stops: the cheapest cost
found so far, and a lower bound on the cost of every circuit at the next scale or
deeper. Step through a search for one random rotation.

```{raw} html
<div class="qs-figure" id="qs-scan">
  <div class="qs-controls" role="group" aria-label="Gate set">
    <button type="button" data-gate-set="clifford_t" aria-pressed="true">Clifford+T, T-count</button>
    <button type="button" data-gate-set="clifford_sqrt_t" aria-pressed="false">Clifford+√T, catalysis (1,3)</button>
  </div>
  <div class="qs-controls" role="group" aria-label="Precision">
    <button type="button" data-epsilon="0.001" aria-pressed="true">ε = 10<sup>-3</sup></button>
    <button type="button" data-epsilon="1e-05" aria-pressed="false">ε = 10<sup>-5</sup></button>
  </div>
  <div class="qs-controls qs-steps" role="group" aria-label="Steps">
    <button type="button" data-step="back">Back</button>
    <button type="button" data-step="next">Next scale</button>
    <button type="button" data-step="all">Show all</button>
  </div>
  <p class="qs-narration" aria-live="polite"></p>
  <div class="qs-legend"></div>
  <div class="qs-chart"></div>
  <p class="qs-caption"></p>
  <details class="qs-table"><summary>Show the numbers</summary><div></div></details>
</div>
```

On Clifford+T the bound is tight, and the search usually ends with a
certificate: no cheaper circuit exists. On Clifford+√T the bound is looser, and
each scale holds about 256 times as many points as the last, so the search stops
when the next scale would cost more time than it is likely to save.

## Finding the points

At one scale, the circuits within ε of the target are points of a lattice: the
integer coefficients of the circuit's first column, 8 of them on Clifford+T and 16
on Clifford+√T. They lie in a thin search region, which QuickSynth encloses in an
ellipsoid. Its walk (the Fincke–Pohst method) fixes one coordinate at a time: each
choice of the last coordinate leaves an interval for the next, so it only ever
visits points inside the ellipsoid. Every point it finds is then checked exactly:
is it a unitary within ε of the target, and what does its circuit cost?

How much work the walk does depends on the basis it describes the lattice in. A
skewed basis makes it try many rows that hold no point; the basis reduced by the
LLL algorithm, in the ellipsoid's own metric, lines the rows up with the region.
The figure shows the same idea in two dimensions.

```{raw} html
<div class="qs-figure" id="qs-walk">
  <div class="qs-controls" role="group" aria-label="Basis">
    <button type="button" data-basis="skewed" aria-pressed="true">Skewed basis</button>
    <button type="button" data-basis="reduced" aria-pressed="false">LLL-reduced basis</button>
  </div>
  <div class="qs-sliders">
    <label>scale <output data-out="scale"></output>
      <input type="range" data-input="scale" min="0" max="4" step="1" value="2"></label>
    <div class="qs-controls">
      <button type="button" data-action="walk">Walk</button>
      <button type="button" data-action="finish">Show the result</button>
    </div>
  </div>
  <p class="qs-narration" aria-live="polite"></p>
  <div class="qs-legend"></div>
  <div class="qs-chart"></div>
  <p class="qs-caption"></p>
</div>
```

## Reaching deeper: the split

Every circuit of scale *K* is a product *P R* of a prefix *P* of a small scale *k*
and a remainder *R* of scale *K* − *k*. The split takes the cheapest prefixes and,
for each, searches only for the remainder, near *P*<sup>−1</sup> times the target:
a far smaller search. It can miss a circuit that does not start with one of the
chosen prefixes, so it adds to the scan rather than replacing it. Below about
ε = 10<sup>-5</sup> it reaches scales the scan cannot afford.

```{raw} html
<div class="qs-figure" id="qs-split">
  <div class="qs-controls" role="group" aria-label="Precision">
    <button type="button" data-epsilon="1e-07" aria-pressed="true">ε = 10<sup>-7</sup></button>
    <button type="button" data-epsilon="1e-08" aria-pressed="false">ε = 10<sup>-8</sup></button>
  </div>
  <div class="qs-legend"></div>
  <div class="qs-chart"></div>
  <p class="qs-caption"></p>
  <details class="qs-table"><summary>Show the numbers</summary><div></div></details>
</div>
```

## The price decides the winner

The cost of a circuit with *n*<sub>T</sub> `T` gates and *n*<sub>√T</sub> `√T`
gates is *t* *n*<sub>T</sub> + *q* *n*<sub>√T</sub>, and the prices *t* and *q*
depend on how the magic states are made. Each dot is the answer QuickSynth gives
for the rotation of the scale-by-scale figure, at ε = 10<sup>-4</sup>, at one
price. The line holds
every circuit of the winner's cost at the chosen price; nothing found lies below
it. Make `√T` dearer and the line tilts, and the winner moves towards fewer `√T`
gates.

```{raw} html
<div class="qs-figure" id="qs-price">
  <div class="qs-sliders">
    <label>price of √T against T, r = q/t: <output data-out="price"></output>
      <input type="range" data-input="price" min="0" max="6" step="1" value="5"></label>
  </div>
  <div class="qs-legend"></div>
  <div class="qs-chart"></div>
  <p class="qs-word"></p>
  <p class="qs-caption"></p>
  <details class="qs-table"><summary>Show the numbers</summary><div></div></details>
</div>
<noscript><p>The interactive figures need JavaScript.</p></noscript>
<script src="https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js"
        integrity="sha384-CjloA8y00+1SDAUkjs099PVfnY2KmDC2BZnws9kh8D/lX1s46w6EPhpXdqMfjK6i"
        crossorigin="anonymous"></script>
<script src="../../_static/showcase/common.js"></script>
<script src="../../_static/showcase/mechanism.js"></script>
```
