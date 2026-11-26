/* Copyright (c) 2026 Parity Quantum Computing GmbH
   SPDX-License-Identifier: MIT */

/* Benchmark: headline tiles and the three interactive figures of
   docs/benchmark.md, drawn with D3 from results.json, the
   aggregates plotted in RESULTS.md's figures. */

(function () {
  "use strict";

  const DATA_URL = new URL("results.json", document.currentScript.src);

  const METHODS = {
    quicksynth: { label: "QuickSynth", color: "--qs-quicksynth", symbol: "circle", fill: true },
    quicksynth_fast: {
      label: "QuickSynth, fast",
      color: "--qs-quicksynth",
      symbol: "circle",
      fill: false,
      dash: "5 3",
    },
    cyclosynth: { label: "cyclosynth", color: "--qs-cyclosynth", symbol: "diamond", fill: true },
    pygridsynth: { label: "pygridsynth", color: "--qs-pygridsynth", symbol: "triangle", fill: true },
    optimum: { label: "certified optimum", color: "--pst-color-text-base", symbol: "square", fill: true },
  };
  const SYMBOLS = {
    circle: d3.symbolCircle,
    diamond: d3.symbolDiamond,
    triangle: d3.symbolTriangle,
    square: d3.symbolSquare,
  };

  const { sup, eps, signed, duration, seconds, color, el, tooltip, table, axisStyle, pressed, live } =
    window.qsShowcase;

  // --- shared pieces -----------------------------------------------------------------

  function symbolPath(kind, size) {
    return d3.symbol(SYMBOLS[kind], size)();
  }

  /** A marker for a method: filled, or hollow with a coloured ring. */
  function marker(g, method, size) {
    const m = METHODS[method];
    const c = color(m.color);
    return g
      .append("path")
      .attr("d", symbolPath(m.symbol, size))
      .attr("fill", m.fill ? c : color("--pst-color-background"))
      .attr("stroke", m.fill ? color("--pst-color-background") : c)
      .attr("stroke-width", m.fill ? 1 : 1.8)
      .classed("qs-lead", method.startsWith("quicksynth"));
  }

  /** Legend entries: [{method, line}] with a marker, a line or both. */
  function legend(node, entries) {
    node.replaceChildren();
    for (const e of entries) {
      const span = el("span");
      const svg = d3.create("svg").attr("width", 24).attr("height", 12).attr("aria-hidden", "true");
      const c = color(e.color || METHODS[e.method].color);
      if (e.line) {
        svg
          .append("line")
          .attr("x1", 1)
          .attr("x2", 23)
          .attr("y1", 6)
          .attr("y2", 6)
          .attr("stroke", c)
          .attr("stroke-width", 2)
          .attr("stroke-dasharray", e.dash || (e.method && METHODS[e.method].dash) || null);
      }
      if (e.method && e.marker !== false) {
        marker(svg.append("g").attr("transform", "translate(12,6)"), e.method, 50);
      }
      span.append(svg.node(), document.createTextNode(e.label || METHODS[e.method].label));
      node.append(span);
    }
  }

  function hitTargets(g, points, tip, describe, focusable = true) {
    const sel = g
      .selectAll("g.point")
      .data(points)
      .join("g")
      .attr("class", "point")
      .attr("transform", (p) => `translate(${p.x},${p.y})`);
    sel
      .append("circle")
      .attr("class", "hit")
      .attr("r", 13)
      .attr("tabindex", focusable ? 0 : null)
      .attr("role", focusable ? "img" : null)
      .attr("aria-label", (p) => (focusable ? describe(p).aria : null));
    sel
      .append("circle")
      .attr("class", "ring")
      .attr("r", 8)
      .attr("fill", "none")
      .attr("stroke", "none")
      .attr("pointer-events", "none");
    const on = (event, p) => {
      const node = event.currentTarget.parentNode;
      d3.select(node).select(".ring").classed("active", true);
      const d = describe(p);
      tip.show(p.x + p.dx, p.y + p.dy, d.head, d.rows);
    };
    const off = (event) => {
      d3.select(event.currentTarget.parentNode).select(".ring").classed("active", false);
      tip.hide();
    };
    sel
      .select(".hit")
      .on("pointerenter", on)
      .on("focus", on)
      .on("pointerleave", off)
      .on("blur", off);
  }

  // --- tiles -------------------------------------------------------------------------

  function tiles(data) {
    const node = document.getElementById("qs-tiles");
    if (!node) return;
    const haar = data.compile.find((p) => p.gate_set === "clifford_sqrt_t" && p.column === "haar");
    const per = (name) => (haar.points[name].hours * 3600) / haar.rotations;
    const cost = (name) => haar.points[name].cost;
    const over = (a, b) => 100 * (cost(a) / cost(b) - 1);
    // the precisions where the fast setting is faster than cyclosynth, at the
    // price where it is so at the fewest
    let faster = Infinity;
    let total = 0;
    for (const s of data.sweep.filter((s) => s.series.quicksynth_fast)) {
      const cyclo = new Map(s.series.cyclosynth.map((p) => [p.epsilon, p.cpu_seconds]));
      const both = s.series.quicksynth_fast.filter((p) => cyclo.has(p.epsilon));
      faster = Math.min(faster, both.filter((p) => p.cpu_seconds < cyclo.get(p.epsilon)).length);
      total = Math.max(total, both.length);
    }
    const savings = data.growth.series
      .filter((s) => s.gate_set === "clifford_sqrt_t")
      .map((s) => s.saving_percent);
    const items = [
      [
        `${per("quicksynth_fast").toFixed(2)} s`,
        `per rotation, wall-clock, with the fast setting on Clifford+√T at ε = 10⁻⁷; cyclosynth takes ${per(
          "cyclosynth"
        ).toFixed(2)} s, ${(per("cyclosynth") / per("quicksynth_fast")).toFixed(0)} times as long`,
      ],
      [
        `+${over("cyclosynth", "quicksynth").toFixed(1)}%`,
        `magic-state cost of cyclosynth's circuits over QuickSynth's (+${over(
          "cyclosynth",
          "quicksynth_fast"
        ).toFixed(1)}% over the fast setting's)`,
      ],
      [
        `${faster} of ${total}`,
        "precisions from 10⁻³ to 10⁻⁸ where the fast setting is faster than cyclosynth, at each of the four prices",
      ],
      [
        `${Math.round(Math.min(...savings))}–${Math.round(Math.max(...savings))}%`,
        "less cost per bit of precision with √T than with Clifford+T alone, by price",
      ],
    ];
    node.replaceChildren(
      ...items.map(([v, l]) => {
        const t = el("div", { class: "qs-tile" });
        t.append(el("div", { class: "qs-tile-value" }, v), el("div", { class: "qs-tile-label" }, l));
        return t;
      })
    );
  }

  // --- cost against compile time -----------------------------------------------------

  const COMPILE_TITLE = { haar: "10⁶ random rotations", trotter: "Trotter", qaoa: "QAOA", adiabatic: "Adiabatic" };
  // in hours; about evenly spaced on a log axis
  const DURATION_TICKS = [
    [1 / 3600, "1 s"],
    [10 / 3600, "10 s"],
    [1 / 60, "1 min"],
    [10 / 60, "10 min"],
    [1, "1 h"],
    [6, "6 h"],
    [24, "1 day"],
    [168, "1 week"],
    [720, "1 month"],
    [8766, "1 year"],
  ];

  function compileFigure(data) {
    const root = document.getElementById("qs-compile");
    if (!root) return null;
    const state = { gate: "clifford_sqrt_t", column: "haar" };
    const chart = root.querySelector(".qs-chart");
    const tip = tooltip(chart);

    function render() {
      pressed(root, "data-gate-set", state.gate);
      pressed(root, "data-column", state.column);
      const panel = data.compile.find((p) => p.gate_set === state.gate && p.column === state.column);
      const names = ["pygridsynth", "cyclosynth", "quicksynth", "quicksynth_fast"].filter((n) => panel.points[n]);
      const entries = names.slice().reverse().map((n) => ({ method: n }));
      if (panel.front.length) {
        entries.push({ method: "cyclosynth", label: "cyclosynth, other deadlines" });
      }
      if (panel.bound !== null) entries.push({ color: "--qs-bound", line: true, dash: "4 3", label: "lower bound" });
      legend(root.querySelector(".qs-legend"), entries);

      const width = chart.clientWidth;
      const height = Math.max(260, Math.min(400, width * 0.55));
      const m = { top: 14, right: 20, bottom: 46, left: 54 };
      const xs = [...Object.values(panel.points).map((p) => p.hours), ...panel.front.map((p) => p.hours)];
      const ys = [
        ...Object.values(panel.points).flatMap((p) => [p.cost - p.se, p.cost + p.se]),
        ...panel.front.map((p) => p.cost),
      ];
      if (panel.bound !== null) ys.push(panel.bound);
      const x = d3
        .scaleLog()
        .domain([d3.min(xs) / 3, d3.max(xs) * 3])
        .range([m.left, width - m.right]);
      const [ylo, yhi] = d3.extent(ys);
      const pad = (yhi - ylo || 1) * 0.15;
      const y = d3
        .scaleLinear()
        .domain([ylo - pad, yhi + pad])
        .nice()
        .range([height - m.bottom, m.top]);

      const svg = d3
        .create("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("role", "img")
        .attr(
          "aria-label",
          `Mean cost against compile time, ${COMPILE_TITLE[state.column]}, ${
            state.gate === "clifford_t" ? "Clifford+T" : "Clifford+√T"
          }`
        );
      const [d0, d1] = x.domain();
      const ticks = DURATION_TICKS.filter(([h]) => h >= d0 && h <= d1);
      svg
        .append("g")
        .attr("class", "grid")
        .selectAll("line")
        .data(y.ticks(5))
        .join("line")
        .attr("x1", m.left)
        .attr("x2", width - m.right)
        .attr("y1", (v) => y(v))
        .attr("y2", (v) => y(v));
      svg
        .append("g")
        .attr("transform", `translate(0,${height - m.bottom})`)
        .call(
          d3
            .axisBottom(x)
            .tickValues(ticks.map((t) => t[0]))
            .tickFormat((h) => ticks.find((t) => t[0] === h)[1])
        )
        .call(axisStyle);
      svg.append("g").attr("transform", `translate(${m.left},0)`).call(d3.axisLeft(y).ticks(5)).call(axisStyle);
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("x", (m.left + width - m.right) / 2)
        .attr("y", height - 8)
        .attr("text-anchor", "middle")
        .text(
          width < 560
            ? "time to compile (log scale; left is faster)"
            : "time to compile, one rotation after another (log scale; left is faster)"
        );
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("transform", `translate(14,${(m.top + height - m.bottom) / 2}) rotate(-90)`)
        .attr("text-anchor", "middle")
        .text(state.gate === "clifford_t" ? "mean T-count per rotation" : "mean cost per rotation");

      if (panel.bound !== null) {
        svg
          .append("line")
          .attr("x1", m.left)
          .attr("x2", width - m.right)
          .attr("y1", y(panel.bound))
          .attr("y2", y(panel.bound))
          .attr("stroke", color("--qs-bound"))
          .attr("stroke-dasharray", "4 3");
        svg
          .append("text")
          .attr("class", "note")
          .attr("x", m.left + 6)
          .attr("y", y(panel.bound) - 5)
          .text("lower bound");
      }

      // arrows from QuickSynth to cyclosynth's time at QuickSynth's cost
      const defs = svg.append("defs");
      defs
        .append("marker")
        .attr("id", "qs-arrow")
        .attr("viewBox", "0 0 10 10")
        .attr("refX", 9)
        .attr("refY", 5)
        .attr("markerWidth", 7)
        .attr("markerHeight", 7)
        .attr("orient", "auto")
        .append("path")
        .attr("d", "M0,0 L10,5 L0,10 z")
        .attr("fill", color("--pst-color-text-muted"));
      for (const n of ["quicksynth", "quicksynth_fast"]) {
        const p = panel.points[n];
        const to = panel.matched_hours[n];
        if (!p || !to || to <= p.hours * 1.05) continue;
        const x0 = x(p.hours) + 9;
        const x1 = x(to) - 6;
        svg
          .append("line")
          .attr("x1", x0)
          .attr("x2", x1)
          .attr("y1", y(p.cost))
          .attr("y2", y(p.cost))
          .attr("stroke", color("--pst-color-text-muted"))
          .attr("stroke-width", 1)
          .attr("marker-end", "url(#qs-arrow)");
        const label = `×${(to / p.hours).toFixed(to / p.hours < 10 ? 1 : 0)}`;
        const w = 16 + 7.5 * label.length;
        // on the arrow when it is long enough, above it otherwise
        const py = x1 - x0 > w + 16 ? y(p.cost) : y(p.cost) - 15;
        const pill = svg.append("g").attr("class", "qs-pill").attr("transform", `translate(${(x0 + x1) / 2},${py})`);
        pill.append("rect").attr("x", -w / 2).attr("y", -10).attr("width", w).attr("height", 20).attr("rx", 10);
        pill.append("text").attr("text-anchor", "middle").attr("dy", "0.35em").text(label);
      }

      // cyclosynth's other deadlines: filled on one thread, hollow on several
      const fg = svg.append("g");
      const front = panel.front.map((p) => ({ ...p, x: x(p.hours), y: y(p.cost), dx: 0, dy: 0 }));
      for (const p of front) {
        const multi = /(\d+) threads/.exec(p.setting);
        const g = fg.append("g").attr("transform", `translate(${p.x},${p.y})`);
        const c = color("--qs-cyclosynth");
        g.append("path")
          .attr("d", symbolPath("diamond", 28))
          .attr("fill", multi && +multi[1] > 1 ? color("--pst-color-background") : c)
          .attr("stroke", c)
          .attr("stroke-width", 1)
          .attr("opacity", 0.8);
      }

      // the methods, with error bars of one standard error
      const pts = names.map((n) => ({
        name: n,
        ...panel.points[n],
        x: x(panel.points[n].hours),
        y: y(panel.points[n].cost),
        dx: 0,
        dy: 0,
      }));
      for (const p of pts) {
        if (p.se > 0 && p.name !== "quicksynth_fast") {
          svg
            .append("line")
            .attr("x1", p.x)
            .attr("x2", p.x)
            .attr("y1", y(p.cost - p.se))
            .attr("y2", y(p.cost + p.se))
            .attr("stroke", color(METHODS[p.name].color))
            .attr("stroke-width", 1.2);
        }
        marker(svg.append("g").attr("transform", `translate(${p.x},${p.y})`), p.name, 110);
        if (p.name !== "cyclosynth" && Math.abs(p.percent) >= 0.5) {
          svg
            .append("text")
            .attr("class", "label")
            .attr("x", p.x + 9)
            .attr("y", p.y + 17)
            .text(signed(p.percent));
        }
      }

      const per = (h) => (state.column === "haar" ? ` (${seconds((h * 3600) / panel.rotations)} per rotation)` : "");
      const hits = svg.append("g");
      hitTargets(hits, [...front.map((p) => ({ ...p, front: true })), ...pts], tip, (p) => {
        const head = p.front ? `cyclosynth, ${p.setting}` : METHODS[p.name].label;
        const rows = [
          {
            value: p.se ? `${p.cost.toFixed(2)} ± ${p.se.toFixed(2)}` : p.cost.toFixed(2),
            label: state.gate === "clifford_t" ? "mean T-count" : "mean cost",
          },
          { value: duration(p.hours), label: `to compile${per(p.hours)}` },
        ];
        if (!p.front && p.name !== "cyclosynth") {
          rows.push({
            value: signed(p.percent),
            label: p.name === "pygridsynth" ? "cost against QuickSynth" : "cost against cyclosynth",
          });
        }
        return { head, rows, aria: `${head}: ${rows.map((r) => `${r.value} ${r.label}`).join(", ")}` };
      });
      chart.replaceChildren(svg.node(), chart.querySelector(".qs-tooltip") || el("div"));
      tip.hide();

      // caption and table
      const parts = ["quicksynth", "quicksynth_fast"]
        .filter((n) => panel.points[n])
        .map((n) => {
          const p = panel.points[n];
          const r = panel.matched_hours[n] / p.hours;
          const who = n === "quicksynth_fast" ? "the fast setting" : "QuickSynth";
          const cost = Math.abs(p.percent) < 0.05 ? "the same cost" : `${signed(p.percent)} cost`;
          return `${who} compiles it ${r.toFixed(1)} times as fast as cyclosynth at its default deadline, at ${cost}`;
        });
      const caption = `${COMPILE_TITLE[state.column]} at ε = ${eps(panel.epsilon)}: ${parts.join("; ")}.`;
      root.querySelector(".qs-caption").textContent =
        caption + (panel.points.pygridsynth ? ` pygridsynth's percentage compares it with QuickSynth.` : "");
      table(
        root.querySelector(".qs-table div"),
        ["method", "mean cost", "standard error", "time to compile", "cost against reference"],
        [
          ...pts.map((p) => [
            METHODS[p.name].label,
            p.cost.toFixed(2),
            p.se.toFixed(2),
            duration(p.hours),
            p.name === "cyclosynth" ? "reference" : signed(p.percent),
          ]),
          ...front.map((p) => [`cyclosynth, ${p.setting}`, p.cost.toFixed(2), "", duration(p.hours), ""]),
        ]
      );
    }

    root.addEventListener("click", (e) => {
      const b = e.target.closest("button");
      if (!b) return;
      if (b.dataset.gateSet) state.gate = b.dataset.gateSet;
      if (b.dataset.column) state.column = b.dataset.column;
      render();
    });
    return { root, render };
  }

  // --- cost against precision --------------------------------------------------------

  function sweepFigure(data) {
    const root = document.getElementById("qs-sweep");
    if (!root) return null;
    const state = { index: 3 };
    const chart = root.querySelector(".qs-chart");
    const tip = tooltip(chart);

    function render() {
      pressed(root, "data-sweep", state.index);
      const s = data.sweep[state.index];
      const ct = s.gate_set === "clifford_t";
      const names = ["quicksynth", "quicksynth_fast", "cyclosynth", "pygridsynth"].filter((n) => s.series[n]);
      const entries = names.map((n) => ({ method: n, line: true, marker: true }));
      if (s.optimum.length) entries.push({ method: "optimum" });
      entries.push({ color: "--qs-bound", line: true, dash: "4 3", label: "lower bound" });
      legend(root.querySelector(".qs-legend"), entries);

      const width = chart.clientWidth;
      const m = { top: 12, right: 20, bottom: 42, left: 58 };
      const h1 = Math.max(170, Math.min(240, width * 0.3));
      const h2 = h1;
      const gap = 18;
      const height = m.top + h1 + gap + h2 + m.bottom;
      const all = [...new Set(names.flatMap((n) => s.series[n].map((p) => p.epsilon)))].sort((a, b) => b - a);
      const finest = d3.min(all);
      const x = d3
        .scaleLog()
        .domain([10 ** -2.75, finest * 10 ** -0.25])
        .range([m.left, width - m.right]);
      // cost above the bound: pygridsynth is off this scale (noted instead)
      const costNames = names.filter((n) => n !== "pygridsynth");
      const above = [
        0,
        ...costNames.flatMap((n) => s.series[n].flatMap((p) => [p.above_bound - (p.se || 0), p.above_bound + (p.se || 0)])),
        ...s.optimum.map((p) => p.above_bound),
      ];
      const [alo, ahi] = d3.extent(above);
      const yc = d3
        .scaleLinear()
        .domain([alo - 0.05 * (ahi - alo), ahi + 0.12 * (ahi - alo)])
        .nice()
        .range([m.top + h1, m.top]);
      const cpus = names.flatMap((n) => s.series[n].map((p) => p.cpu_seconds).filter((v) => v > 0));
      const yt = d3
        .scaleLog()
        .domain([d3.min(cpus) / 1.6, d3.max(cpus) * 2])
        .range([m.top + h1 + gap + h2, m.top + h1 + gap]);

      const svg = d3
        .create("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("role", "img")
        .attr("aria-label", `Mean cost above the lower bound and CPU time against precision, ${s.price_name}`);
      const decades = d3.range(3, 10).map((k) => 10 ** -k).filter((e) => e >= finest * 0.99);

      for (const [yScale, ticks] of [
        [yc, yc.ticks(4)],
        [yt, yt.ticks(4).filter((v) => Math.abs(Math.log10(v) % 1) < 1e-9)],
      ]) {
        svg
          .append("g")
          .attr("class", "grid")
          .selectAll("line")
          .data(ticks)
          .join("line")
          .attr("x1", m.left)
          .attr("x2", width - m.right)
          .attr("y1", (v) => yScale(v))
          .attr("y2", (v) => yScale(v));
      }
      svg
        .append("g")
        .attr("transform", `translate(0,${height - m.bottom})`)
        .call(d3.axisBottom(x).tickValues(decades).tickFormat(eps))
        .call(axisStyle);
      svg
        .append("g")
        .attr("transform", `translate(0,${m.top + h1})`)
        .call(d3.axisBottom(x).tickValues(decades).tickFormat(() => ""))
        .call(axisStyle);
      svg.append("g").attr("transform", `translate(${m.left},0)`).call(d3.axisLeft(yc).ticks(4)).call(axisStyle);
      const cpuTicks = yt.ticks(4).filter((v) => Math.abs(Math.log10(v) % 1) < 1e-9);
      svg
        .append("g")
        .attr("transform", `translate(${m.left},0)`)
        .call(
          d3
            .axisLeft(yt)
            .tickValues(cpuTicks)
            .tickFormat((v) => (v >= 1 ? `${v} s` : `${Math.round(v * 1000)} ms`))
        )
        .call(axisStyle);
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("x", (m.left + width - m.right) / 2)
        .attr("y", height - 6)
        .attr("text-anchor", "middle")
        .text("precision ε (log scale; finer to the right)");
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("transform", `translate(14,${m.top + h1 / 2}) rotate(-90)`)
        .attr("text-anchor", "middle")
        .text(ct ? "T-count above bound" : "cost above bound");
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("transform", `translate(14,${m.top + h1 + gap + h2 / 2}) rotate(-90)`)
        .attr("text-anchor", "middle")
        .text("CPU per rotation");

      // the bound is the zero of the top panel
      svg
        .append("line")
        .attr("x1", m.left)
        .attr("x2", width - m.right)
        .attr("y1", yc(0))
        .attr("y2", yc(0))
        .attr("stroke", color("--qs-bound"))
        .attr("stroke-dasharray", "4 3");

      const isDecade = (e) => Math.abs(Math.log10(e) - Math.round(Math.log10(e))) < 1e-6;
      for (const n of names) {
        const c = color(METHODS[n].color);
        const pts = s.series[n];
        const dash = METHODS[n].dash || null;
        const lead = n.startsWith("quicksynth");
        if (n !== "pygridsynth") {
          if (pts[0].se !== undefined && !ct && n !== "quicksynth_fast") {
            svg
              .append("path")
              .attr("fill", c)
              .attr("fill-opacity", 0.14)
              .attr(
                "d",
                d3
                  .area()
                  .x((p) => x(p.epsilon))
                  .y0((p) => yc(p.above_bound - p.se))
                  .y1((p) => yc(p.above_bound + p.se))(pts)
              );
          }
          svg
            .append("path")
            .attr("fill", "none")
            .attr("stroke", c)
            .attr("stroke-width", lead ? 2.6 : 2)
            .classed("qs-lead", lead)
            .attr("stroke-dasharray", ct && n === "cyclosynth" ? "6 4" : dash)
            .attr(
              "d",
              d3
                .line()
                .x((p) => x(p.epsilon))
                .y((p) => yc(p.above_bound))(pts)
            );
          for (const p of pts.filter((p) => isDecade(p.epsilon))) {
            marker(svg.append("g").attr("transform", `translate(${x(p.epsilon)},${yc(p.above_bound)})`), n, 34);
          }
        }
        const cpu = pts.filter((p) => p.cpu_seconds > 0);
        svg
          .append("path")
          .attr("fill", "none")
          .attr("stroke", c)
          .attr("stroke-width", lead ? 2.6 : 2)
          .classed("qs-lead", lead)
          .attr("stroke-dasharray", dash)
          .attr(
            "d",
            d3
              .line()
              .x((p) => x(p.epsilon))
              .y((p) => yt(p.cpu_seconds))(cpu)
          );
        for (const p of cpu.filter((p) => isDecade(p.epsilon))) {
          marker(svg.append("g").attr("transform", `translate(${x(p.epsilon)},${yt(p.cpu_seconds)})`), n, 34);
        }
      }
      for (const p of s.optimum) {
        marker(svg.append("g").attr("transform", `translate(${x(p.epsilon)},${yc(p.above_bound)})`), "optimum", 46);
      }

      // notes: cyclosynth's CPU against QuickSynth's, pygridsynth's cost
      const notes = [];
      const ours = s.series.quicksynth_fast || s.series.quicksynth;
      if (s.series.cyclosynth && ours) {
        const c = new Map(s.series.cyclosynth.map((p) => [p.epsilon, p.cpu_seconds]));
        const r = ours.filter((p) => c.has(p.epsilon)).map((p) => c.get(p.epsilon) / p.cpu_seconds);
        const who = s.series.quicksynth_fast ? "the fast setting's" : "QuickSynth's";
        notes.push(`cyclosynth takes ${d3.min(r).toFixed(1)} to ${d3.max(r).toFixed(1)} times ${who} CPU time.`);
      }
      if (s.series.pygridsynth && s.series.quicksynth) {
        const b = new Map(s.bound.map((p) => [p.epsilon, p.cost]));
        const q = new Map(s.series.quicksynth.map((p) => [p.epsilon, p.above_bound + b.get(p.epsilon)]));
        const r = s.series.pygridsynth
          .filter((p) => q.has(p.epsilon) && p.above_bound !== null)
          .map((p) => (p.above_bound + b.get(p.epsilon)) / q.get(p.epsilon));
        notes.push(
          `pygridsynth's T-count, off the scale above, is ${d3.min(r).toFixed(1)} to ${d3.max(r).toFixed(1)} times QuickSynth's.`
        );
      }
      if (ct) notes.push("On Clifford+T, cyclosynth (dashed) returns the same T-counts as QuickSynth.");
      root.querySelector(".qs-caption").textContent = `${s.targets} targets. ${notes.join(" ")}`;

      // crosshair: snaps to the nearest precision, lists every series there
      const cross = svg.append("g").attr("pointer-events", "none").style("display", "none");
      cross
        .append("line")
        .attr("y1", m.top)
        .attr("y2", height - m.bottom)
        .attr("stroke", color("--pst-color-text-muted"))
        .attr("stroke-width", 1);
      const overlay = svg
        .append("rect")
        .attr("x", m.left)
        .attr("y", m.top)
        .attr("width", width - m.left - m.right)
        .attr("height", height - m.top - m.bottom)
        .attr("fill", "transparent")
        .attr("tabindex", 0)
        .attr("aria-label", "Precision readout: use the arrow keys to move");
      let index = 0;
      const show = (i) => {
        index = Math.max(0, Math.min(all.length - 1, i));
        const e = all[index];
        cross.style("display", null).attr("transform", `translate(${x(e)},0)`);
        const rows = [];
        for (const n of names) {
          const p = s.series[n].find((q) => q.epsilon === e);
          if (!p) continue;
          rows.push({
            value: p.above_bound === null ? "–" : p.above_bound.toFixed(2),
            label: `${METHODS[n].label}, ${seconds(p.cpu_seconds)} CPU`,
            color: color(METHODS[n].color),
            dash: !!METHODS[n].dash,
          });
        }
        const o = s.optimum.find((q) => q.epsilon === e);
        if (o) rows.push({ value: o.above_bound.toFixed(2), label: "certified optimum" });
        tip.show(x(e), m.top + h1 / 2, `ε = ${eps(e)}, ${ct ? "T-count" : "cost"} above bound`, rows);
      };
      const nearest = (px) => {
        const le = Math.log10(x.invert(px));
        return d3.minIndex(all, (e) => Math.abs(Math.log10(e) - le));
      };
      overlay
        .on("pointermove", (event) => show(nearest(d3.pointer(event)[0])))
        .on("pointerleave blur", () => {
          cross.style("display", "none");
          tip.hide();
        })
        .on("focus", () => show(index))
        .on("keydown", (event) => {
          if (event.key === "ArrowRight") show(index + 1);
          else if (event.key === "ArrowLeft") show(index - 1);
          else return;
          event.preventDefault();
        });

      chart.replaceChildren(svg.node(), chart.querySelector(".qs-tooltip") || el("div"));
      tip.hide();

      table(
        root.querySelector(".qs-table div"),
        ["ε", "lower bound", ...names.flatMap((n) => [`${METHODS[n].label}: above bound`, `${METHODS[n].label}: CPU`]), "optimum: above bound"],
        all.map((e) => {
          const row = [eps(e), (s.bound.find((p) => p.epsilon === e) || {}).cost?.toFixed(2) ?? ""];
          for (const n of names) {
            const p = s.series[n].find((q) => q.epsilon === e);
            row.push(p && p.above_bound !== null ? p.above_bound.toFixed(2) : "", p ? seconds(p.cpu_seconds) : "");
          }
          const o = s.optimum.find((q) => q.epsilon === e);
          row.push(o ? o.above_bound.toFixed(2) : "");
          return row;
        })
      );
    }

    root.addEventListener("click", (e) => {
      const b = e.target.closest("button");
      if (!b || b.dataset.sweep === undefined) return;
      state.index = +b.dataset.sweep;
      render();
    });
    return { root, render };
  }

  // --- growth of the cost ball and the slope of the bound ------------------------------

  function growthFigure(data) {
    const root = document.getElementById("qs-growth");
    if (!root) return null;
    const chart = root.querySelector(".qs-chart");
    const tip = tooltip(chart);
    const g = data.growth;
    const sqrt = g.series.filter((s) => s.gate_set === "clifford_sqrt_t");
    const colorOf = (s) => (s.gate_set === "clifford_t" ? "--qs-clifford-t" : `--qs-price-${sqrt.indexOf(s)}`);
    const name = (s) => (s.gate_set === "clifford_t" ? "Clifford+T, T-count" : `${s.price_name} (${s.price.join(",")})`);

    function render() {
      legend(
        root.querySelector(".qs-legend"),
        g.series.map((s) => ({ color: colorOf(s), line: true, label: name(s) }))
      );
      const full = chart.clientWidth;
      const pair = full >= 640;
      chart.style.gridTemplateColumns = pair ? "1fr 1fr" : "1fr";
      const width = pair ? (full - 16) / 2 : full;
      const height = Math.max(240, Math.min(320, width * 0.75));
      const m = { top: 12, right: 44, bottom: 44, left: 56 };

      // (a) the number of Clifford cosets of each cost, in units of the T price
      const a = d3.create("svg").attr("width", width).attr("height", height).attr("role", "img");
      a.attr("aria-label", "Number of distinct circuits against their cost, by price");
      const x = d3.scaleLinear().domain([0, 20]).range([m.left, width - m.right]);
      const y = d3
        .scaleLog()
        .domain([1, 1e10])
        .range([height - m.bottom, m.top]);
      a.append("g")
        .attr("class", "grid")
        .selectAll("line")
        .data([1, 1e2, 1e4, 1e6, 1e8, 1e10])
        .join("line")
        .attr("x1", m.left)
        .attr("x2", width - m.right)
        .attr("y1", (v) => y(v))
        .attr("y2", (v) => y(v));
      a.append("g")
        .attr("transform", `translate(0,${height - m.bottom})`)
        .call(d3.axisBottom(x).ticks(5))
        .call(axisStyle);
      a.append("g")
        .attr("transform", `translate(${m.left},0)`)
        .call(
          d3
            .axisLeft(y)
            .tickValues([1, 1e2, 1e4, 1e6, 1e8, 1e10])
            .tickFormat((v) => (v === 1 ? "1" : `10${sup(Math.round(Math.log10(v)))}`))
        )
        .call(axisStyle);
      a.append("text")
        .attr("class", "axis-title")
        .attr("x", (m.left + width - m.right) / 2)
        .attr("y", height - 6)
        .attr("text-anchor", "middle")
        .text("cost, in units of the T price t");
      a.append("text")
        .attr("class", "axis-title")
        .attr("transform", `translate(14,${(m.top + height - m.bottom) / 2}) rotate(-90)`)
        .attr("text-anchor", "middle")
        .text("distinct circuits of that cost");
      const pointsA = [];
      for (const s of g.series) {
        const c = color(colorOf(s));
        const xEnd = Math.min(20, Math.log(1e10 / 2.2 / s.eta) / Math.log(s.growth_per_t));
        const xsLine = d3.range(0, xEnd + 1e-9, xEnd / 80);
        a.append("path")
          .attr("fill", "none")
          .attr("stroke", c)
          .attr("stroke-width", 1.5)
          .attr(
            "d",
            d3
              .line()
              .x((v) => x(v))
              .y((v) => y(Math.max(1, s.eta * s.growth_per_t ** v)))(xsLine)
          );
        a.append("text")
          .attr("class", "label")
          .attr("x", x(xEnd) + 4)
          .attr("y", y(s.eta * s.growth_per_t ** xEnd))
          .attr("dominant-baseline", "middle")
          .text(`×${s.growth_per_t.toPrecision(3)}`);
        for (const [cx, count] of s.counts) {
          a.append("circle")
            .attr("cx", x(cx))
            .attr("cy", y(count))
            .attr("r", 3)
            .attr("fill", c)
            .attr("stroke", color("--pst-color-background"))
            .attr("stroke-width", 0.6);
          pointsA.push({ s, cost: cx, count, x: x(cx), y: y(count), dx: 0, dy: 0 });
        }
      }
      hitTargets(a.append("g"), pointsA, tip, (p) => {
        const rows = [
          { value: p.count.toLocaleString("en"), label: "distinct circuits" },
          { value: `${p.cost}`, label: "cost, in units of t" },
        ];
        return { head: name(p.s), rows, aria: "" };
      }, false);

      // (b) the slope of the bound against the price ratio r = q/t
      const b = d3.create("svg").attr("width", width).attr("height", height).attr("role", "img");
      b.attr("aria-label", "Slope of the lower bound per bit of precision against the price ratio");
      const xr = d3
        .scaleLog()
        .domain([0.5, 12])
        .range([m.left, width - m.right]);
      const ys = d3
        .scaleLinear()
        .domain([0, 3.35])
        .range([height - m.bottom, m.top]);
      b.append("g")
        .attr("class", "grid")
        .selectAll("line")
        .data(ys.ticks(4))
        .join("line")
        .attr("x1", m.left)
        .attr("x2", width - m.right)
        .attr("y1", (v) => ys(v))
        .attr("y2", (v) => ys(v));
      b.append("g")
        .attr("transform", `translate(0,${height - m.bottom})`)
        .call(d3.axisBottom(xr).tickValues([0.5, 1, 2, 3, 5, 10]).tickFormat(d3.format("~g")))
        .call(axisStyle);
      b.append("g").attr("transform", `translate(${m.left},0)`).call(d3.axisLeft(ys).ticks(4)).call(axisStyle);
      b.append("text")
        .attr("class", "axis-title")
        .attr("x", (m.left + width - m.right) / 2)
        .attr("y", height - 6)
        .attr("text-anchor", "middle")
        .text("price ratio r = q/t (log scale)");
      b.append("text")
        .attr("class", "axis-title")
        .attr("transform", `translate(14,${(m.top + height - m.bottom) / 2}) rotate(-90)`)
        .attr("text-anchor", "middle")
        .text("cost per bit of precision (units of t)");
      b.append("line")
        .attr("x1", m.left)
        .attr("x2", width - m.right)
        .attr("y1", ys(g.clifford_t_slope))
        .attr("y2", ys(g.clifford_t_slope))
        .attr("stroke", color("--qs-clifford-t"))
        .attr("stroke-width", 1.5);
      b.append("text")
        .attr("class", "note")
        .attr("x", m.left + 6)
        .attr("y", ys(g.clifford_t_slope) - 6)
        .text(`Clifford+T: ${g.clifford_t_slope} per bit`);
      b.append("path")
        .attr("fill", "none")
        .attr("stroke", color("--pst-color-text-muted"))
        .attr("stroke-width", 1.5)
        .attr(
          "d",
          d3
            .line()
            .x((p) => xr(p[0]))
            .y((p) => ys(p[1]))(g.slope_curve)
        );
      const pointsB = [];
      for (const s of sqrt) {
        const r = s.price[1] / s.price[0];
        const c = color(colorOf(s));
        b.append("line")
          .attr("x1", xr(r))
          .attr("x2", xr(r))
          .attr("y1", ys(s.slope_per_bit))
          .attr("y2", ys(g.clifford_t_slope))
          .attr("stroke", c)
          .attr("stroke-dasharray", "2 2");
        // below right of the point: the curve rises, so that corner stays free
        b.append("text")
          .attr("class", "label")
          .attr("x", xr(r) + 7)
          .attr("y", ys(s.slope_per_bit) + 14)
          .attr("dominant-baseline", "middle")
          .text(`−${Math.round(s.saving_percent)}%`);
        b.append("circle")
          .attr("cx", xr(r))
          .attr("cy", ys(s.slope_per_bit))
          .attr("r", 5)
          .attr("fill", c)
          .attr("stroke", color("--pst-color-background"))
          .attr("stroke-width", 1);
        pointsB.push({ s, r, x: xr(r), y: ys(s.slope_per_bit), dx: pair ? width + 16 : 0, dy: pair ? 0 : height + 16 });
      }
      hitTargets(b.append("g"), pointsB, tip, (p) => {
        const rows = [
          { value: p.s.slope_per_bit.toFixed(2), label: "t per bit of precision" },
          { value: `−${p.s.saving_percent.toFixed(1)}%`, label: "against Clifford+T" },
        ];
        return { head: name(p.s), rows, aria: `${name(p.s)}: ${rows.map((r) => `${r.value} ${r.label}`).join(", ")}` };
      });

      chart.replaceChildren(a.node(), b.node(), chart.querySelector(".qs-tooltip") || el("div"));
      tip.hide();
      root.querySelector(".qs-caption").textContent =
        "Left: the number of distinct circuits (Clifford cosets) of each cost; the lines are their exponential growth, labelled with the growth per unit of the T price. Right: the slope of the lower bound per bit of precision; the percentages are the saving against Clifford+T.";
      table(
        root.querySelector(".qs-table div"),
        ["price", "growth per t", "cost per bit (t)", "against Clifford+T"],
        g.series.map((s) => [
          name(s),
          `×${s.growth_per_t.toPrecision(3)}`,
          s.slope_per_bit.toFixed(2),
          s.gate_set === "clifford_t" ? "reference" : `−${s.saving_percent.toFixed(1)}%`,
        ])
      );
    }
    return { root, render };
  }

  // --- start ---------------------------------------------------------------------------

  fetch(DATA_URL)
    .then((r) => {
      if (!r.ok) throw new Error(`${r.status} ${DATA_URL}`);
      return r.json();
    })
    .then((data) => {
      tiles(data);
      const figures = [compileFigure(data), sweepFigure(data), growthFigure(data)].filter(Boolean);
      live(() => figures.forEach((f) => f.render()), figures[0].root);
    })
    .catch((err) => {
      for (const id of ["qs-compile", "qs-sweep", "qs-growth"]) {
        const n = document.getElementById(id);
        if (n) n.querySelector(".qs-caption").textContent = `The figure could not load its data (${err.message}).`;
      }
    });
})();
