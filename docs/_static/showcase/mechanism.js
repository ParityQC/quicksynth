/* Copyright (c) 2026 Parity Quantum Computing GmbH
   SPDX-License-Identifier: MIT */

/* Algorithm (docs/algorithm.md): the Bloch-sphere picture of the circuits
   up to a cost (points.bin) and the scan, scale by scale (mechanism.json).
   Both files are exported from the library itself. */

(function () {
  "use strict";

  const { sup, eps, seconds, color, el, tooltip, table, axisStyle, pressed, live } = window.qsShowcase;
  const BASE = document.currentScript.src;
  const fmt = d3.format(",");

  /** A count as 1,234 or 8.8×10⁷. */
  function count(n) {
    if (n < 1e6) return fmt(Math.round(n));
    const k = Math.floor(Math.log10(n));
    return `${(n / 10 ** k).toFixed(1)}×10${sup(k)}`;
  }

  function gates(cost, nq, price) {
    const nt = (cost - price[1] * nq) / price[0];
    const parts = [];
    if (nt) parts.push(`${nt} T`);
    if (nq) parts.push(`${nq} √T`);
    return parts.length ? parts.join(" + ") : "Cliffords only";
  }

  // --- figure 1: every circuit is a point ---------------------------------------------

  const EPS_STEPS = [0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2];

  // rotations of a 3-vector
  const rotX = (p, a) => [p[0], p[1] * Math.cos(a) - p[2] * Math.sin(a), p[1] * Math.sin(a) + p[2] * Math.cos(a)];
  const rotZ = (p, a) => [p[0] * Math.cos(a) - p[1] * Math.sin(a), p[0] * Math.sin(a) + p[1] * Math.cos(a), p[2]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const norm = (a) => {
    const n = Math.hypot(...a);
    return a.map((v) => v / n);
  };
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

  function sphereFigure(meta, buffer) {
    const root = document.getElementById("qs-sphere");
    if (!root) return null;
    const view = new DataView(buffer);
    const sets = meta.sets.map((s) => {
      const xyz = new Float32Array(3 * s.count);
      const cost = new Uint8Array(s.count);
      const nq = new Uint8Array(s.count);
      for (let i = 0; i < s.count; i++) {
        const o = (s.offset + i) * meta.record_bytes;
        xyz[3 * i] = view.getInt16(o, true) / meta.scale;
        xyz[3 * i + 1] = view.getInt16(o + 2, true) / meta.scale;
        xyz[3 * i + 2] = view.getInt16(o + 4, true) / meta.scale;
        cost[i] = view.getUint8(o + 6);
        nq[i] = view.getUint8(o + 7);
      }
      return { ...s, xyz, cost, nq };
    });
    const maxCost = d3.max(sets, (s) => s.max_cost);
    const costInput = root.querySelector('[data-input="cost"]');
    const epsInput = root.querySelector('[data-input="eps"]');
    costInput.max = String(maxCost);

    // a fixed first target, then random ones on demand
    const state = { yaw: -0.6, pitch: 0.42, target: norm([0.62, -0.35, 0.7]) };
    const rand = d3.randomNormal.source(d3.randomLcg(0.42))();
    const panes = [];

    function winner(s, budget, cosr) {
      let best = -1;
      let bestDot = -2;
      let inside = 0;
      for (let i = 0; i < s.count; i++) {
        if (s.cost[i] > budget) break; // sorted by cost
        const d = s.xyz[3 * i] * state.target[0] + s.xyz[3 * i + 1] * state.target[1] + s.xyz[3 * i + 2] * state.target[2];
        if (d < cosr) continue;
        inside++;
        if (best < 0 || (s.cost[i] === s.cost[best] && d > bestDot)) {
          best = i;
          bestDot = d;
        }
      }
      return { best, inside };
    }

    function toView(p) {
      return rotX(rotZ(p, state.yaw), state.pitch);
    }

    function fromScreen(pane, ev) {
      const rect = pane.canvas.getBoundingClientRect();
      const x = (ev.clientX - rect.left - pane.cx) / pane.r;
      const z = -(ev.clientY - rect.top - pane.cy) / pane.r;
      const rr = x * x + z * z;
      if (rr > 1) return null;
      const v = [x, -Math.sqrt(1 - rr), z];
      return rotZ(rotX(v, -state.pitch), -state.yaw);
    }

    function drawPane(pane, s, budget, theta) {
      const { ctx, size, cx, cy, r } = pane;
      const dpr = window.devicePixelRatio || 1;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, size, size);
      const lo = color("--qs-seq-0");
      const hi = color("--qs-seq-1");
      const ramp = d3.interpolateLab(lo, hi);
      const ink = color("--pst-color-text-base");
      const muted = color("--pst-color-text-muted");
      const border = color("--pst-color-border");

      ctx.beginPath();
      ctx.arc(cx, cy, r, 0, 2 * Math.PI);
      ctx.fillStyle = color("--pst-color-on-background");
      ctx.fill();
      ctx.strokeStyle = border;
      ctx.lineWidth = 1;
      ctx.stroke();

      // a circle on the sphere around axis a at angle t, front half only
      const ring = (a, t, style, width, dash) => {
        const u = norm(Math.abs(a[2]) < 0.9 ? cross(a, [0, 0, 1]) : cross(a, [1, 0, 0]));
        const w = cross(a, u);
        ctx.strokeStyle = style;
        ctx.lineWidth = width;
        ctx.setLineDash(dash || []);
        ctx.beginPath();
        let pen = false;
        for (let k = 0; k <= 120; k++) {
          const f = (2 * Math.PI * k) / 120;
          const p = [0, 1, 2].map(
            (j) => Math.cos(t) * a[j] + Math.sin(t) * (Math.cos(f) * u[j] + Math.sin(f) * w[j])
          );
          const v = toView(p);
          if (v[1] > 0) {
            pen = false;
            continue;
          }
          const X = cx + r * v[0];
          const Y = cy - r * v[2];
          if (pen) ctx.lineTo(X, Y);
          else ctx.moveTo(X, Y);
          pen = true;
        }
        ctx.stroke();
        ctx.setLineDash([]);
      };
      ring([0, 0, 1], Math.PI / 2, border, 1);
      ring([1, 0, 0], Math.PI / 2, border, 1);

      // points, the expensive ones first so the cheap ones stay on top
      const px = Math.max(1.4, Math.min(2.4, r / 110));
      let i = s.count - 1;
      while (i >= 0 && s.cost[i] > budget) i--;
      let current = -1;
      for (; i >= 0; i--) {
        const c = s.cost[i];
        if (c !== current) {
          current = c;
          ctx.fillStyle = ramp(maxCost ? c / maxCost : 0);
        }
        const v = toView([s.xyz[3 * i], s.xyz[3 * i + 1], s.xyz[3 * i + 2]]);
        if (v[1] > 0) continue;
        ctx.fillRect(cx + r * v[0] - px / 2, cy - r * v[2] - px / 2, px, px);
      }

      // axis labels
      ctx.font = `12px ${getComputedStyle(document.body).fontFamily}`;
      ctx.fillStyle = muted;
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      for (const [p, label] of [
        [[0, 0, 1], "|0⟩"],
        [[0, 0, -1], "|1⟩"],
        [[1, 0, 0], "|+⟩"],
        [[-1, 0, 0], "|−⟩"],
      ]) {
        const v = toView(p);
        if (v[1] > 0.2) continue;
        ctx.fillText(label, cx + (r + 16) * v[0], cy - (r + 16) * v[2]);
      }

      // the target, its circle and the winner
      ring(state.target, theta, ink, 1.5);
      const tv = toView(state.target);
      if (tv[1] <= 0) {
        ctx.strokeStyle = ink;
        ctx.lineWidth = 1.5;
        const X = cx + r * tv[0];
        const Y = cy - r * tv[2];
        ctx.beginPath();
        ctx.moveTo(X - 4, Y);
        ctx.lineTo(X + 4, Y);
        ctx.moveTo(X, Y - 4);
        ctx.lineTo(X, Y + 4);
        ctx.stroke();
      }
      if (pane.win.best >= 0) {
        const b = pane.win.best;
        const v = toView([s.xyz[3 * b], s.xyz[3 * b + 1], s.xyz[3 * b + 2]]);
        if (v[1] <= 0) {
          const X = cx + r * v[0];
          const Y = cy - r * v[2];
          ctx.beginPath();
          ctx.arc(X, Y, 3.5, 0, 2 * Math.PI);
          ctx.fillStyle = ramp(s.cost[b] / maxCost);
          ctx.fill();
          ctx.beginPath();
          ctx.arc(X, Y, 7, 0, 2 * Math.PI);
          ctx.strokeStyle = ink;
          ctx.lineWidth = 2;
          ctx.stroke();
        }
      }
    }

    function draw() {
      const budget = +costInput.value;
      const e = EPS_STEPS[+epsInput.value];
      const theta = 2 * Math.asin(e);
      const cosr = Math.cos(theta);
      root.querySelector('[data-out="cost"]').textContent = budget;
      root.querySelector('[data-out="eps"]').textContent = e;
      const summary = [];
      panes.forEach((pane, k) => {
        const s = sets[k];
        pane.win = winner(s, budget, cosr);
        drawPane(pane, s, budget, theta);
        const circuits = d3.sum(s.circuits.slice(0, budget + 1));
        const b = pane.win.best;
        const answer =
          b < 0
            ? "none yet: raise the cost or ε"
            : `cost ${s.cost[b]} (${gates(s.cost[b], s.nq[b], s.price)})`;
        pane.readout.replaceChildren(
          el("b", {}, s.label),
          el("br"),
          document.createTextNode(`${fmt(circuits)} circuits, ${pane.win.inside} inside the circle`),
          el("br"),
          document.createTextNode("cheapest inside: "),
          el("b", {}, answer)
        );
        pane.canvas.setAttribute(
          "aria-label",
          `${s.label}: ${fmt(circuits)} circuits up to cost ${budget}; ${pane.win.inside} within ε = ${e} of the target; cheapest inside: ${answer}.`
        );
        summary.push({ s, b });
      });
      const [t, q] = summary;
      let caption = `Up to cost ${budget}, ε = ${e}. `;
      if (t.b >= 0 && q.b >= 0) {
        const d = t.s.cost[t.b] - q.s.cost[q.b];
        caption +=
          d > 0
            ? `With √T the cheapest circuit within ε costs ${d} less.`
            : d === 0
              ? "Here both gate sets find a circuit of the same cost; move the target or lower ε."
              : "Here Clifford+T happens to be cheaper; move the target.";
      } else {
        caption += "Raise the cost or ε until a dot lies inside the circle.";
      }
      root.querySelector(".qs-caption").textContent =
        caption + " Circuits are counted up to Clifford gates at the end, which cost nothing.";
    }

    function render() {
      const holder = root.querySelector(".qs-spheres");
      const full = holder.clientWidth;
      const pair = full >= 560;
      holder.style.gridTemplateColumns = pair ? "1fr 1fr" : "1fr";
      const size = Math.min(360, pair ? (full - 16) / 2 : full);
      if (!panes.length) {
        for (let k = 0; k < sets.length; k++) {
          const box = el("div", { class: "qs-sphere-pane" });
          const canvas = el("canvas", { role: "img", tabindex: "0" });
          const readout = el("p", { class: "qs-readout" });
          box.append(canvas, readout);
          holder.append(box);
          panes.push({ canvas, readout, ctx: canvas.getContext("2d"), win: { best: -1, inside: 0 } });
          attach(panes[k]);
        }
      }
      const dpr = window.devicePixelRatio || 1;
      for (const pane of panes) {
        pane.size = size;
        pane.r = size / 2 - 26;
        pane.cx = size / 2;
        pane.cy = size / 2;
        pane.canvas.width = Math.round(size * dpr);
        pane.canvas.height = Math.round(size * dpr);
        pane.canvas.style.width = `${size}px`;
        pane.canvas.style.height = `${size}px`;
      }
      const key = root.querySelector(".qs-ramp");
      const ramp = d3.interpolateLab(color("--qs-seq-0"), color("--qs-seq-1"));
      key.replaceChildren(document.createTextNode("cheap"));
      for (let c = 0; c <= maxCost; c++) {
        const sw = el("span", { class: "qs-ramp-step", title: `cost ${c}` });
        sw.style.background = ramp(c / maxCost);
        key.append(sw);
      }
      key.append(document.createTextNode(`expensive (cost 0 to ${maxCost})`));
      draw();
      table(
        root.querySelector(".qs-table div"),
        ["cost", ...sets.map((s) => `${s.label}: circuits`), ...sets.map((s) => `${s.label}: states`)],
        d3.range(maxCost + 1).map((c) => [
          String(c),
          ...sets.map((s) => fmt(s.circuits[c] ?? 0)),
          ...sets.map((s) => fmt(s.cost.filter((v) => v === c).length)),
        ])
      );
    }

    let frame = 0;
    const redraw = () => {
      if (frame) return;
      frame = requestAnimationFrame(() => {
        frame = 0;
        draw();
      });
    };

    function attach(pane) {
      let start = null;
      pane.canvas.addEventListener("pointerdown", (ev) => {
        start = { x: ev.clientX, y: ev.clientY, yaw: state.yaw, pitch: state.pitch, moved: false };
        pane.canvas.setPointerCapture(ev.pointerId);
      });
      pane.canvas.addEventListener("pointermove", (ev) => {
        if (!start) return;
        const dx = ev.clientX - start.x;
        const dy = ev.clientY - start.y;
        if (Math.abs(dx) + Math.abs(dy) > 3) start.moved = true;
        if (!start.moved) return;
        state.yaw = start.yaw + dx / pane.r;
        state.pitch = Math.max(-1.5, Math.min(1.5, start.pitch + dy / pane.r));
        redraw();
      });
      pane.canvas.addEventListener("pointerup", (ev) => {
        if (start && !start.moved) {
          const p = fromScreen(pane, ev);
          if (p) {
            state.target = norm(p);
            redraw();
          }
        }
        start = null;
      });
      pane.canvas.addEventListener("keydown", (ev) => {
        const step = { ArrowLeft: [-0.1, 0], ArrowRight: [0.1, 0], ArrowUp: [0, -0.1], ArrowDown: [0, 0.1] }[ev.key];
        if (!step) return;
        ev.preventDefault();
        state.yaw += step[0];
        state.pitch = Math.max(-1.5, Math.min(1.5, state.pitch + step[1]));
        redraw();
      });
    }

    costInput.addEventListener("input", redraw);
    epsInput.addEventListener("input", redraw);
    root.querySelector('[data-action="target"]').addEventListener("click", () => {
      // a random direction, turned to face the viewer
      state.target = norm([rand(), rand(), rand()]);
      const v = toView(state.target);
      if (v[1] > 0) state.target = state.target.map((x) => -x);
      redraw();
    });
    return { root, render };
  }

  // --- figure 2: the scan, scale by scale -----------------------------------------------

  function scanFigure(runs) {
    const root = document.getElementById("qs-scan");
    if (!root) return null;
    const state = { gate_set: "clifford_t", epsilon: 1e-3, step: 0 };
    const tip = tooltip(root.querySelector(".qs-chart"));
    const run = () => runs.find((r) => r.gate_set === state.gate_set && Math.abs(r.epsilon - state.epsilon) < 1e-12);

    function narrate(r) {
      const sc = r.scales;
      const n = sc.length;
      if (state.step < n) {
        const s = sc[state.step];
        const prev = state.step ? sc[state.step - 1].best : null;
        let text = `Scale ${s.scale}: `;
        if (s.points === 0) text += "no lattice point lies close enough to the target. ";
        else text += `${count(s.points)} lattice points walked in ${seconds(s.seconds)}. `;
        if (s.best === null) text += s.points ? "None is a circuit within ε. " : "";
        else if (prev === null) text += `The first circuit within ε costs ${s.best}. `;
        else if (s.best < prev) text += `A cheaper circuit: cost ${s.best}. `;
        else text += `Nothing cheaper than ${s.best}. `;
        const nb = state.step + 1 < n ? sc[state.step + 1].bound : r.next.bound;
        text += `Every circuit at scale ${s.scale + 1} or deeper costs at least ${nb}.`;
        return text;
      }
      if (r.certified) {
        return `Stop. The bound for scale ${r.next.scale}, ${r.next.bound}, exceeds the best cost, ${r.cost}: no deeper scale can hold a cheaper circuit, so the answer is certified optimal. Scale ${r.next.scale} would have walked about ${count(r.next.points)} points.`;
      }
      return `Stop. A circuit at scale ${r.next.scale} might still cost as little as ${r.next.bound}, below ${r.cost}, but that scale would walk about ${count(r.next.points)} points, about ${seconds(r.next.seconds)} against ${seconds(d3.sum(sc, (s) => s.seconds))} so far. The stopping rule judges it not worth its time; the answer, cost ${r.cost}, is not certified.`;
    }

    function render() {
      const r = run();
      const sc = r.scales;
      const n = sc.length;
      state.step = Math.max(0, Math.min(n, state.step));
      pressed(root, "data-gate-set", state.gate_set);
      for (const b of root.querySelectorAll("[data-epsilon]")) {
        b.setAttribute("aria-pressed", String(Math.abs(+b.dataset.epsilon - state.epsilon) < 1e-12));
      }
      root.querySelector('[data-step="back"]').disabled = state.step === 0;
      root.querySelector('[data-step="next"]').disabled = state.step === n;
      root.querySelector('[data-step="all"]').disabled = state.step === n;
      root.querySelector('[data-step="next"]').textContent = state.step === n - 1 ? "Decide" : "Next scale";
      root.querySelector(".qs-narration").textContent = narrate(r);

      const shown = sc.slice(0, Math.min(state.step + 1, n));
      const done = state.step === n;
      const ink = color("--pst-color-text-base");
      const qs = color("--qs-quicksynth");
      const bnd = color("--qs-bound");
      const muted = color("--pst-color-text-muted");

      keyLegend(root.querySelector(".qs-legend"), [
        { label: "cheapest circuit found", color: qs, kind: "dot" },
        { label: "lower bound at this scale and deeper", color: bnd, kind: "line", dash: "4 3" },
        { label: "lattice points walked", color: muted, kind: "bar", opacity: 0.55 },
      ]);

      const chart = root.querySelector(".qs-chart");
      const width = chart.clientWidth;
      const m = { top: 12, right: 16, bottom: 40, left: 58 };
      const h1 = 190;
      const h2 = 150;
      const gap = 22;
      const height = m.top + h1 + gap + h2 + m.bottom;
      const scales = [...sc.map((s) => s.scale), r.next.scale];
      const x = d3
        .scaleBand()
        .domain(scales)
        .range([m.left, width - m.right])
        .padding(0.35);
      const xc = (k) => x(k) + x.bandwidth() / 2;

      const costs = [...sc.map((s) => s.bound), r.next.bound, ...sc.map((s) => s.best).filter((v) => v !== null)];
      const y1 = d3
        .scaleLinear()
        .domain([d3.min(costs) - 2, d3.max(costs) + 2])
        .nice()
        .range([m.top + h1, m.top]);
      const maxPts = Math.max(10, r.next.points, d3.max(sc, (s) => s.points));
      const y2 = d3
        .scaleLog()
        .domain([1, maxPts * 1.5])
        .range([m.top + h1 + gap + h2, m.top + h1 + gap])
        .clamp(true);

      const svg = d3
        .create("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("role", "img")
        .attr("aria-label", `The scan on ${r.gate_set === "clifford_t" ? "Clifford+T" : "Clifford+√T"} at ε = ${eps(r.epsilon)}, scale by scale`);

      // grids and axes
      svg
        .append("g")
        .attr("class", "grid")
        .selectAll("line")
        .data(y1.ticks(5))
        .join("line")
        .attr("x1", m.left)
        .attr("x2", width - m.right)
        .attr("y1", (v) => y1(v))
        .attr("y2", (v) => y1(v));
      svg.append("g").attr("transform", `translate(${m.left},0)`).call(d3.axisLeft(y1).ticks(5)).call(axisStyle);
      const decades = d3
        .range(0, Math.floor(Math.log10(maxPts * 1.5)) + 1, Math.log10(maxPts) > 6 ? 2 : 1)
        .map((k) => 10 ** k);
      svg
        .append("g")
        .attr("class", "grid")
        .selectAll("line")
        .data(decades)
        .join("line")
        .attr("x1", m.left)
        .attr("x2", width - m.right)
        .attr("y1", (v) => y2(v))
        .attr("y2", (v) => y2(v));
      svg
        .append("g")
        .attr("transform", `translate(${m.left},0)`)
        .call(
          d3
            .axisLeft(y2)
            .tickValues(decades)
            .tickFormat((v) => (v === 1 ? "1" : `10${sup(Math.round(Math.log10(v)))}`))
        )
        .call(axisStyle);
      svg
        .append("g")
        .attr("transform", `translate(0,${m.top + h1 + gap + h2})`)
        .call(d3.axisBottom(x).tickSizeOuter(0))
        .call(axisStyle);
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("x", (m.left + width - m.right) / 2)
        .attr("y", height - 6)
        .attr("text-anchor", "middle")
        .text("scale K");
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("transform", `translate(14,${m.top + h1 / 2}) rotate(-90)`)
        .attr("text-anchor", "middle")
        .text(r.gate_set === "clifford_t" ? "T-count" : "cost");
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("transform", `translate(14,${m.top + h1 + gap + h2 / 2}) rotate(-90)`)
        .attr("text-anchor", "middle")
        .text("points walked");

      // bars: points walked, and the next scale as an outline once the scan decides
      const bars = svg.append("g");
      const base = y2(1);
      for (const s of shown) {
        if (s.points > 0) {
          bars
            .append("rect")
            .attr("x", x(s.scale))
            .attr("width", x.bandwidth())
            .attr("y", y2(Math.max(1, s.points)))
            .attr("height", Math.max(0, base - y2(Math.max(1, s.points))))
            .attr("rx", 3)
            .attr("fill", muted)
            .attr("fill-opacity", 0.55);
        }
        bars
          .append("text")
          .attr("class", "label")
          .attr("x", xc(s.scale))
          .attr("y", y2(Math.max(1, s.points)) - 5)
          .attr("text-anchor", "middle")
          .attr("font-size", 11)
          .text(count(s.points));
      }
      if (done) {
        bars
          .append("rect")
          .attr("x", x(r.next.scale))
          .attr("width", x.bandwidth())
          .attr("y", y2(r.next.points))
          .attr("height", base - y2(r.next.points))
          .attr("rx", 3)
          .attr("fill", "none")
          .attr("stroke", muted)
          .attr("stroke-dasharray", "4 3");
        bars
          .append("text")
          .attr("class", "note")
          .attr("x", xc(r.next.scale))
          .attr("y", y2(r.next.points) - 5)
          .attr("text-anchor", "middle")
          .attr("font-size", 11)
          .text(`≈${count(r.next.points)}`);
      }

      // the bound: for scale K and deeper, revealed one scale ahead
      const boundPts = [...sc.map((s) => [s.scale, s.bound]), [r.next.scale, r.next.bound]].slice(
        0,
        Math.min(state.step + 2, n + 1)
      );
      svg
        .append("path")
        .attr("fill", "none")
        .attr("stroke", bnd)
        .attr("stroke-width", 2)
        .attr("stroke-dasharray", "4 3")
        .attr(
          "d",
          d3
            .line()
            .x((p) => xc(p[0]))
            .y((p) => y1(p[1]))(boundPts)
        );
      svg
        .append("g")
        .selectAll("circle")
        .data(boundPts)
        .join("circle")
        .attr("cx", (p) => xc(p[0]))
        .attr("cy", (p) => y1(p[1]))
        .attr("r", 3)
        .attr("fill", bnd);

      // the best cost so far
      const bestPts = shown.filter((s) => s.best !== null).map((s) => [s.scale, s.best]);
      if (bestPts.length) {
        const ext = done ? [...bestPts, [r.next.scale, r.cost]] : bestPts;
        svg
          .append("path")
          .attr("fill", "none")
          .attr("stroke", qs)
          .attr("stroke-width", 2)
          .attr(
            "d",
            d3
              .line()
              .curve(d3.curveStepAfter)
              .x((p) => xc(p[0]))
              .y((p) => y1(p[1]))(ext)
          );
        svg
          .append("g")
          .selectAll("circle")
          .data(bestPts)
          .join("circle")
          .attr("cx", (p) => xc(p[0]))
          .attr("cy", (p) => y1(p[1]))
          .attr("r", 5)
          .attr("fill", qs)
          .attr("stroke", color("--pst-color-background"))
          .attr("stroke-width", 1.5);
        const lastBest = bestPts[bestPts.length - 1];
        svg
          .append("text")
          .attr("class", "label")
          .attr("x", xc(lastBest[0]) - 9)
          .attr("y", y1(lastBest[1]) - 9)
          .attr("text-anchor", "end")
          .text(`cost ${lastBest[1]}`);
      }
      if (done) {
        const X = xc(r.next.scale);
        svg
          .append("circle")
          .attr("cx", X)
          .attr("cy", y1(r.next.bound))
          .attr("r", 8)
          .attr("fill", "none")
          .attr("stroke", ink)
          .attr("stroke-width", 1.5);
        svg
          .append("text")
          .attr("class", "label halo")
          .attr("x", X)
          .attr("y", y1(r.next.bound) + 24)
          .attr("text-anchor", "end")
          .attr("stroke", color("--pst-color-background"))
          .text(r.certified ? "optimal" : "not certified");
      }

      // hover: one column per scale
      const cols = svg.append("g");
      const rows = (k) => {
        const s = sc.find((v) => v.scale === k);
        if (!s) return { head: `scale ${k}, not walked`, rows: [{ value: `≈${count(r.next.points)}`, label: "points, predicted" }, { value: String(r.next.bound), label: "lower bound" }] };
        return {
          head: `scale ${k}`,
          rows: [
            { value: count(s.points), label: "points walked" },
            { value: seconds(s.seconds), label: "time" },
            { value: s.best === null ? "–" : String(s.best), label: "cheapest so far" },
            { value: String(s.bound), label: "lower bound, this scale and deeper" },
          ],
        };
      };
      const visible = done ? scales : shown.map((s) => s.scale);
      cols
        .selectAll("rect")
        .data(visible)
        .join("rect")
        .attr("x", (k) => x(k) - (x.step() * x.paddingInner()) / 2)
        .attr("width", x.step())
        .attr("y", m.top)
        .attr("height", h1 + gap + h2)
        .attr("fill", "transparent")
        .attr("tabindex", 0)
        .attr("aria-label", (k) => {
          const d = rows(k);
          return `${d.head}: ${d.rows.map((v) => `${v.label} ${v.value}`).join(", ")}`;
        })
        .on("pointerenter focus", (ev, k) => {
          const d = rows(k);
          tip.show(xc(k), m.top + 40, d.head, d.rows);
        })
        .on("pointerleave blur", () => tip.hide());

      chart.replaceChildren(svg.node(), chart.querySelector(".qs-tooltip") || el("div"));
      tip.hide();

      root.querySelector(".qs-caption").textContent = done
        ? `The answer: cost ${r.cost}, ${gates(r.cost, r.n_q, r.price)}, found in ${seconds(r.seconds)}: ${r.word}`
        : `A Haar-random rotation at ε = ${eps(r.epsilon)}, whole scales, one thread.`;
      table(
        root.querySelector(".qs-table div"),
        ["scale", "points walked", "time", "cheapest so far", "lower bound, this scale and deeper"],
        [
          ...sc.map((s) => [String(s.scale), fmt(s.points), seconds(s.seconds), s.best === null ? "–" : String(s.best), String(s.bound)]),
          [String(r.next.scale), `≈${count(r.next.points)} (not walked)`, `≈${seconds(r.next.seconds)}`, "", String(r.next.bound)],
        ]
      );
    }

    root.addEventListener("click", (ev) => {
      const b = ev.target.closest("button");
      if (!b) return;
      const n = run().scales.length;
      if (b.dataset.gateSet) {
        state.gate_set = b.dataset.gateSet;
        state.step = 0;
      } else if (b.dataset.epsilon) {
        state.epsilon = +b.dataset.epsilon;
        state.step = 0;
      } else if (b.dataset.step === "next") state.step++;
      else if (b.dataset.step === "back") state.step--;
      else if (b.dataset.step === "all") state.step = n;
      else return;
      render();
    });
    return { root, render };
  }

  // --- legend entries: {label, color, kind: dot | ring | line | bar | outline, dash} ------

  function keyLegend(node, entries) {
    node.replaceChildren();
    for (const e of entries) {
      const span = el("span");
      const svg = d3.create("svg").attr("width", 24).attr("height", 12).attr("aria-hidden", "true");
      if (e.kind === "dot") svg.append("circle").attr("cx", 12).attr("cy", 6).attr("r", 4).attr("fill", e.color);
      if (e.kind === "ring")
        svg.append("circle").attr("cx", 12).attr("cy", 6).attr("r", 4).attr("fill", "none").attr("stroke", e.color).attr("stroke-width", 1.5);
      if (e.kind === "line")
        svg
          .append("line")
          .attr("x1", 1)
          .attr("x2", 23)
          .attr("y1", 6)
          .attr("y2", 6)
          .attr("stroke", e.color)
          .attr("stroke-width", 2)
          .attr("stroke-dasharray", e.dash || null);
      if (e.kind === "bar" || e.kind === "outline")
        svg
          .append("rect")
          .attr("x", 7)
          .attr("y", 1)
          .attr("width", 10)
          .attr("height", 10)
          .attr("rx", 2)
          .attr("fill", e.kind === "bar" ? e.color : "none")
          .attr("fill-opacity", e.opacity ?? 1)
          .attr("stroke", e.kind === "outline" ? e.color : "none")
          .attr("stroke-dasharray", e.kind === "outline" ? "3 2" : null);
      span.append(svg.node(), document.createTextNode(e.label));
      node.append(span);
    }
  }

  // --- figure 3: the walk, in two dimensions ---------------------------------------------

  // A cartoon: the lattice Z^2 at spacing 0.9 / sqrt(2)^scale, a thin search region and
  // the ellipse around it. The walk fixes the last coordinate, then walks an interval.
  const WALK = {
    center: [0.13, -0.07],
    axes: [3.2, 0.55], // semi-axes of the ellipse
    angle: 0.42,
    region: [0.86, 0.62], // the search region, as fractions of the ellipse's semi-axes
    skewed: [
      [1, 2],
      [3, 7],
    ], // columns b1 = (1, 3), b2 = (2, 7): determinant 1, b1 across the region
  };

  function walkFigure() {
    const root = document.getElementById("qs-walk");
    if (!root) return null;
    const state = { basis: "skewed", scale: 2, shown: Infinity, timer: null };
    const scaleInput = root.querySelector('[data-input="scale"]');

    const c = Math.cos(WALK.angle);
    const s = Math.sin(WALK.angle);
    // A = R diag(1/a^2, 1/b^2) R^T
    const [a, b] = WALK.axes;
    const A = [
      [(c * c) / a ** 2 + (s * s) / b ** 2, c * s * (1 / a ** 2 - 1 / b ** 2)],
      [c * s * (1 / a ** 2 - 1 / b ** 2), (s * s) / a ** 2 + (c * c) / b ** 2],
    ];
    const form = (u, v) => u[0] * (A[0][0] * v[0] + A[0][1] * v[1]) + u[1] * (A[1][0] * v[0] + A[1][1] * v[1]);
    const inRegion = (p) => {
      const d = [p[0] - WALK.center[0], p[1] - WALK.center[1]];
      const u = c * d[0] + s * d[1];
      const v = -s * d[0] + c * d[1];
      return (u / (a * WALK.region[0])) ** 2 + (v / (b * WALK.region[1])) ** 2 <= 1;
    };

    /** Integer basis (columns) of Z^2: skewed, or Lagrange-reduced in the metric A. */
    function basisOf(kind, h) {
      let b1 = [WALK.skewed[0][0], WALK.skewed[1][0]];
      let b2 = [WALK.skewed[0][1], WALK.skewed[1][1]];
      if (kind === "reduced") {
        const n = (v) => form([v[0] * h, v[1] * h], [v[0] * h, v[1] * h]);
        for (let it = 0; it < 50; it++) {
          if (n(b2) < n(b1)) [b1, b2] = [b2, b1];
          const mu = Math.round(form([b1[0] * h, b1[1] * h], [b2[0] * h, b2[1] * h]) / n(b1));
          if (mu === 0) break;
          b2 = [b2[0] - mu * b1[0], b2[1] - mu * b1[1]];
        }
      }
      return [b1, b2];
    }

    /** The walk: rows (fixed z2) in order, each with its interval and the points on it. */
    function walk(kind, scale) {
      const h = 0.9 / Math.SQRT2 ** scale;
      const [b1, b2] = basisOf(kind, h);
      const B = [
        [b1[0] * h, b2[0] * h],
        [b1[1] * h, b2[1] * h],
      ];
      const col = (j) => [B[0][j], B[1][j]];
      const Q11 = form(col(0), col(0));
      const Q12 = form(col(0), col(1));
      const Q22 = form(col(1), col(1));
      const det = B[0][0] * B[1][1] - B[0][1] * B[1][0];
      const zc = [
        (B[1][1] * WALK.center[0] - B[0][1] * WALK.center[1]) / det,
        (-B[1][0] * WALK.center[0] + B[0][0] * WALK.center[1]) / det,
      ];
      const D = Q22 - (Q12 * Q12) / Q11;
      const span = 1 / Math.sqrt(D);
      const rows = [];
      for (let z2 = Math.ceil(zc[1] - span); z2 <= Math.floor(zc[1] + span); z2++) {
        const d2 = z2 - zc[1];
        const mid = zc[0] - (Q12 / Q11) * d2;
        const half = Math.sqrt(Math.max(0, (1 - D * d2 * d2) / Q11));
        const at = (z1) => [B[0][0] * z1 + B[0][1] * z2, B[1][0] * z1 + B[1][1] * z2];
        const points = [];
        for (let z1 = Math.ceil(mid - half); z1 <= Math.floor(mid + half); z1++) {
          const p = at(z1);
          points.push({ p, inside: inRegion(p) });
        }
        rows.push({ from: at(mid - half), to: at(mid + half), points });
      }
      return { h, b1: col(0), b2: col(1), rows };
    }

    function counts(w) {
      const pts = w.rows.flatMap((r) => r.points);
      return {
        rows: w.rows.length,
        empty: w.rows.filter((r) => !r.points.length).length,
        points: pts.length,
        inside: pts.filter((p) => p.inside).length,
      };
    }

    function render() {
      stop();
      draw();
    }

    function draw() {
      state.scale = +scaleInput.value;
      root.querySelector('[data-out="scale"]').textContent = state.scale;
      pressed(root, "data-basis", state.basis);
      const w = walk(state.basis, state.scale);
      const other = counts(walk(state.basis === "skewed" ? "reduced" : "skewed", state.scale));
      const all = counts(w);
      const ink = color("--pst-color-text-base");
      const muted = color("--pst-color-text-muted");
      const qs = color("--qs-quicksynth");
      const border = color("--pst-color-border");
      const warn = color("--qs-cyclosynth");

      keyLegend(root.querySelector(".qs-legend"), [
        { label: "search region", color: qs, kind: "bar", opacity: 0.25 },
        { label: "ellipse around it", color: ink, kind: "line", dash: "5 3" },
        { label: "row tried", color: muted, kind: "line" },
        { label: "empty row", color: warn, kind: "line", dash: "2 3" },
        { label: "circuit within ε", color: qs, kind: "dot" },
        { label: "point rejected", color: muted, kind: "ring" },
      ]);

      const chart = root.querySelector(".qs-chart");
      const width = chart.clientWidth;
      const view = [-4.2, 4.2, -2.3, 2.3];
      const height = Math.round(Math.min(420, (width * (view[3] - view[2])) / (view[1] - view[0])));
      const unit = Math.min(width / (view[1] - view[0]), height / (view[3] - view[2]));
      const X = (p) => width / 2 + unit * p[0];
      const Y = (p) => height / 2 - unit * p[1];

      const svg = d3
        .create("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("role", "img")
        .style("overflow", "hidden")
        .attr(
          "aria-label",
          `Two-dimensional walk with the ${state.basis} basis at scale ${state.scale}: ${all.rows} rows tried, ${all.empty} empty, ${all.points} points, ${all.inside} inside the search region.`
        );

      // the lattice
      const lat = [];
      const span = Math.ceil(5 / w.h);
      for (let i = -span; i <= span; i++)
        for (let j = -span; j <= span; j++) {
          const p = [i * w.h, j * w.h];
          if (Math.abs(p[0]) < 4.3 && Math.abs(p[1]) < 2.4) lat.push(p);
        }
      svg
        .append("g")
        .selectAll("circle")
        .data(lat)
        .join("circle")
        .attr("cx", X)
        .attr("cy", Y)
        .attr("r", 1.4)
        .attr("fill", border);

      // region and ellipse
      const deg = (WALK.angle * 180) / Math.PI;
      const C = WALK.center;
      svg
        .append("ellipse")
        .attr("cx", X(C))
        .attr("cy", Y(C))
        .attr("rx", unit * a * WALK.region[0])
        .attr("ry", unit * b * WALK.region[1])
        .attr("transform", `rotate(${-deg},${X(C)},${Y(C)})`)
        .attr("fill", qs)
        .attr("fill-opacity", 0.18);
      svg
        .append("ellipse")
        .attr("cx", X(C))
        .attr("cy", Y(C))
        .attr("rx", unit * a)
        .attr("ry", unit * b)
        .attr("transform", `rotate(${-deg},${X(C)},${Y(C)})`)
        .attr("fill", "none")
        .attr("stroke", ink)
        .attr("stroke-width", 1.2)
        .attr("stroke-dasharray", "5 3");

      // rows and points, up to the animation's position
      let left = state.shown;
      let lastRow = null;
      const rows = svg.append("g");
      const pts = svg.append("g");
      for (const r of w.rows) {
        if (left <= 0) break;
        left--;
        lastRow = r;
        rows
          .append("line")
          .attr("x1", X(r.from))
          .attr("y1", Y(r.from))
          .attr("x2", X(r.to))
          .attr("y2", Y(r.to))
          .attr("stroke", r.points.length ? muted : warn)
          .attr("stroke-width", r.points.length ? 1.5 : 1.5)
          .attr("stroke-dasharray", r.points.length ? null : "2 3")
          .attr("stroke-linecap", "round");
        for (const q of r.points) {
          if (left <= 0) break;
          left--;
          pts
            .append("circle")
            .attr("cx", X(q.p))
            .attr("cy", Y(q.p))
            .attr("r", 4.5)
            .attr("fill", q.inside ? qs : color("--pst-color-background"))
            .attr("stroke", q.inside ? color("--pst-color-background") : muted)
            .attr("stroke-width", 1.5);
        }
      }

      // the basis vectors, from a lattice point low in the view; the svg clips them
      const o = [Math.round(-0.9 / w.h) * w.h, Math.round(-1.7 / w.h) * w.h];
      const defs = svg.append("defs");
      defs
        .append("marker")
        .attr("id", "qs-walk-arrow")
        .attr("viewBox", "0 0 10 10")
        .attr("refX", 9)
        .attr("refY", 5)
        .attr("markerWidth", 7)
        .attr("markerHeight", 7)
        .attr("orient", "auto-start-reverse")
        .append("path")
        .attr("d", "M0,0L10,5L0,10z")
        .attr("fill", ink);
      [
        [w.b1, "b₁"],
        [w.b2, "b₂"],
      ].forEach(([v, label]) => {
        const end = [o[0] + v[0], o[1] + v[1]];
        svg
          .append("line")
          .attr("x1", X(o))
          .attr("y1", Y(o))
          .attr("x2", X(end))
          .attr("y2", Y(end))
          .attr("stroke", ink)
          .attr("stroke-width", 2)
          .attr("marker-end", "url(#qs-walk-arrow)");
        svg
          .append("text")
          .attr("class", "label halo")
          .attr("x", X(end) + 6)
          .attr("y", Y(end) - 6)
          .attr("stroke", color("--pst-color-background"))
          .text(label);
      });

      chart.replaceChildren(svg.node());

      const done = state.shown >= Infinity || !lastRow || lastRow === w.rows[w.rows.length - 1];
      const shownRows = w.rows.slice(0, w.rows.indexOf(lastRow) + 1);
      const sc = counts({ rows: shownRows });
      root.querySelector(".qs-narration").textContent =
        `Rows tried: ${sc.rows}${sc.empty ? `, ${sc.empty} of them empty` : ""}. ` +
        `Points visited: ${sc.points}; ${sc.inside} lie in the search region and are circuits within ε, ` +
        `${sc.points - sc.inside} are rejected by the exact check.` +
        (done ? "" : " …");
      root.querySelector(".qs-caption").textContent =
        `At scale ${state.scale} the ${state.basis} basis tries ${all.rows} rows to find ${all.points} points; ` +
        `the ${state.basis === "skewed" ? "reduced" : "skewed"} basis tries ${other.rows} for the same points. ` +
        "Each scale doubles the points here; in QuickSynth's 8 or 16 dimensions a scale multiplies them by 16 or 256, " +
        "and a skewed basis multiplies the empty rows at every coordinate.";
    }

    function stop() {
      if (state.timer) state.timer.stop();
      state.timer = null;
    }

    function animate() {
      stop();
      const w = walk(state.basis, +scaleInput.value);
      const total = w.rows.length + d3.sum(w.rows, (r) => r.points.length);
      state.shown = 0;
      draw();
      const per = Math.max(60, Math.min(220, 6000 / total));
      state.timer = d3.interval(() => {
        state.shown++;
        draw();
        if (state.shown >= total) {
          stop();
          state.shown = Infinity;
        }
      }, per);
    }

    root.addEventListener("click", (ev) => {
      const btn = ev.target.closest("button");
      if (!btn) return;
      if (btn.dataset.basis) {
        state.basis = btn.dataset.basis;
        state.shown = Infinity;
        render();
      } else if (btn.dataset.action === "walk") animate();
      else if (btn.dataset.action === "finish") {
        state.shown = Infinity;
        render();
      }
    });
    scaleInput.addEventListener("input", () => {
      state.shown = Infinity;
      render();
    });
    return { root, render };
  }

  // --- figure 4: the split --------------------------------------------------------------

  function splitFigure(split) {
    const root = document.getElementById("qs-split");
    if (!root) return null;
    const state = { epsilon: 1e-7 };
    const tip = tooltip(root.querySelector(".qs-chart"));

    function render() {
      const r = split.runs.find((v) => Math.abs(v.epsilon - state.epsilon) < 1e-15);
      for (const b of root.querySelectorAll("[data-epsilon]")) {
        b.setAttribute("aria-pressed", String(Math.abs(+b.dataset.epsilon - state.epsilon) < 1e-15));
      }
      const ink = color("--pst-color-text-base");
      const muted = color("--pst-color-text-muted");
      const qs = color("--qs-quicksynth");
      keyLegend(root.querySelector(".qs-legend"), [
        { label: "scan: points walked", color: muted, kind: "bar", opacity: 0.55 },
        { label: "scan alone would walk", color: muted, kind: "outline" },
        { label: "split: points walked", color: qs, kind: "bar", opacity: 0.8 },
      ]);

      const scales = [...new Set(r.walks.map((w) => w.scale))].sort((u, v) => u - v);
      const byScale = (k, kind) => r.walks.find((w) => w.scale === k && w.search === kind);
      const chart = root.querySelector(".qs-chart");
      const width = chart.clientWidth;
      const m = { top: 14, right: 16, bottom: 40, left: 58 };
      const h1 = 150;
      const h2 = 190;
      const gap = 22;
      const height = m.top + h1 + gap + h2 + m.bottom;
      const x = d3
        .scaleBand()
        .domain(scales)
        .range([m.left, width - m.right])
        .padding(0.25);
      const sub = d3
        .scaleBand()
        .domain(["scan", "split"])
        .range([0, x.bandwidth()])
        .padding(0.08);
      const best = r.walks.filter((w) => w.best !== null);
      const costs = best.map((w) => w.best);
      const y1 = d3
        .scaleLinear()
        .domain([d3.min(costs) - 4, d3.max(costs) + 4])
        .nice()
        .range([m.top + h1, m.top]);
      const maxPts = d3.max(r.walks, (w) => Math.max(w.points, w.scan_points || 0));
      const top = Math.floor(Math.log10(maxPts * 3));
      const y2 = d3
        .scaleLog()
        .domain([1, 10 ** top * 1.2])
        .range([m.top + h1 + gap + h2, m.top + h1 + gap])
        .clamp(true);
      const decades = d3.range(0, top + 1, top > 7 ? 2 : 1).map((k) => 10 ** k);

      const svg = d3
        .create("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("role", "img")
        .attr("aria-label", `Scan and split on Clifford+√T at ε = ${eps(r.epsilon)}`);
      svg
        .append("g")
        .attr("class", "grid")
        .selectAll("line")
        .data([...y1.ticks(4).map(y1), ...decades.map(y2)])
        .join("line")
        .attr("x1", m.left)
        .attr("x2", width - m.right)
        .attr("y1", (v) => v)
        .attr("y2", (v) => v);
      svg.append("g").attr("transform", `translate(${m.left},0)`).call(d3.axisLeft(y1).ticks(4)).call(axisStyle);
      svg
        .append("g")
        .attr("transform", `translate(${m.left},0)`)
        .call(
          d3
            .axisLeft(y2)
            .tickValues(decades)
            .tickFormat((v) => (v === 1 ? "1" : `10${sup(Math.round(Math.log10(v)))}`))
        )
        .call(axisStyle);
      svg
        .append("g")
        .attr("transform", `translate(0,${m.top + h1 + gap + h2})`)
        .call(d3.axisBottom(x).tickSizeOuter(0))
        .call(axisStyle);
      for (const [t, yy] of [
        ["cheapest so far", m.top + h1 / 2],
        ["points walked", m.top + h1 + gap + h2 / 2],
      ]) {
        svg
          .append("text")
          .attr("class", "axis-title")
          .attr("transform", `translate(14,${yy}) rotate(-90)`)
          .attr("text-anchor", "middle")
          .text(t);
      }
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("x", (m.left + width - m.right) / 2)
        .attr("y", height - 6)
        .attr("text-anchor", "middle")
        .text("scale K");

      // bars
      const base = y2(1);
      const bars = svg.append("g");
      for (const k of scales) {
        const sc = byScale(k, "scan");
        const sp = byScale(k, "split");
        const bar = (kind, n, style) => {
          const bx = x(k) + sub(kind);
          const top = y2(Math.max(1, n));
          if (n > 0)
            bars
              .append("rect")
              .attr("x", bx)
              .attr("width", sub.bandwidth())
              .attr("y", top)
              .attr("height", Math.max(0, base - top))
              .attr("rx", 3)
              .attr("fill", style === "outline" ? "none" : style === "split" ? qs : muted)
              .attr("fill-opacity", style === "split" ? 0.8 : 0.55)
              .attr("stroke", style === "outline" ? muted : "none")
              .attr("stroke-dasharray", style === "outline" ? "4 3" : null);
          bars
            .append("text")
            .attr("class", style === "outline" ? "note halo" : "label halo")
            .attr("stroke", color("--pst-color-background"))
            .attr("x", bx + sub.bandwidth() / 2)
            .attr("y", top - 5)
            .attr("text-anchor", "middle")
            .attr("font-size", 10.5)
            .text(style === "outline" ? `≈${count(n)}` : count(n));
        };
        if (sc) bar("scan", sc.points, "scan");
        else if (sp) bar("scan", sp.scan_points, "outline");
        if (sp) bar("split", sp.points, "split");
      }

      // cheapest so far, in the order of the walks
      let run = null;
      const path = [];
      for (const w of r.walks) {
        if (w.best !== null && (run === null || w.best < run)) run = w.best;
        if (run !== null) path.push({ w, cost: run });
      }
      const px = (w) => x(w.scale) + sub(w.search) + sub.bandwidth() / 2;
      svg
        .append("path")
        .attr("fill", "none")
        .attr("stroke", ink)
        .attr("stroke-width", 1.5)
        .attr(
          "d",
          d3
            .line()
            .curve(d3.curveStepAfter)
            .x((d) => px(d.w))
            .y((d) => y1(d.cost))(path)
        );
      svg
        .append("g")
        .selectAll("circle")
        .data(path.filter((d) => d.w.best !== null))
        .join("circle")
        .attr("cx", (d) => px(d.w))
        .attr("cy", (d) => y1(d.w.best))
        .attr("r", 5)
        .attr("fill", (d) => (d.w.search === "split" ? qs : muted))
        .attr("stroke", color("--pst-color-background"))
        .attr("stroke-width", 1.5);
      if (path.length) {
        const last = path[path.length - 1];
        svg
          .append("text")
          .attr("class", "label halo")
          .attr("x", px(last.w) - 12)
          .attr("y", y1(last.cost) + 4)
          .attr("text-anchor", "end")
          .attr("stroke", color("--pst-color-background"))
          .text(`cost ${last.cost}`);
      }

      // hover per scale
      const describe = (k) => {
        const sc = byScale(k, "scan");
        const sp = byScale(k, "split");
        const rows = [];
        if (sc) rows.push({ value: count(sc.points), label: "points walked by the scan", color: muted });
        else if (sp) rows.push({ value: `≈${count(sp.scan_points)}`, label: "points the scan alone would walk", color: muted });
        if (sp)
          rows.push({
            value: count(sp.points),
            label: `points walked by the split: prefixes of scale ${sp.prefix_scale}, remainders of scale ${k - sp.prefix_scale}`,
            color: qs,
          });
        const found = [sc, sp].filter((w) => w && w.best !== null).map((w) => w.best);
        if (found.length) rows.push({ value: String(d3.min(found)), label: "cheapest so far" });
        return { head: `scale ${k}`, rows };
      };
      svg
        .append("g")
        .selectAll("rect")
        .data(scales)
        .join("rect")
        .attr("x", (k) => x(k) - (x.step() * x.paddingInner()) / 2)
        .attr("width", x.step())
        .attr("y", m.top)
        .attr("height", h1 + gap + h2)
        .attr("fill", "transparent")
        .attr("tabindex", 0)
        .attr("aria-label", (k) => {
          const d = describe(k);
          return `${d.head}: ${d.rows.map((v) => `${v.label} ${v.value}`).join(", ")}`;
        })
        .on("pointerenter focus", (ev, k) => {
          const d = describe(k);
          tip.show(x(k) + x.bandwidth() / 2, m.top + 30, d.head, d.rows);
        })
        .on("pointerleave blur", () => tip.hide());
      chart.replaceChildren(svg.node(), chart.querySelector(".qs-tooltip") || el("div"));
      tip.hide();

      const deepest = r.walks.filter((w) => w.search === "split").pop();
      const scanBest = d3.min(r.walks.filter((w) => w.search === "scan" && w.best !== null), (w) => w.best);
      root.querySelector(".qs-caption").textContent =
        `A Haar-random rotation at ε = ${eps(r.epsilon)}. At scale ${deepest.scale} the split walks ` +
        `${count(deepest.points)} points where the scan alone would walk about ${count(deepest.scan_points)}, ` +
        `${count(Math.round(deepest.scan_points / deepest.points))} times as many. ` +
        (scanBest === undefined
          ? `The scan found no circuit within ε; the split's answer costs ${r.cost}.`
          : `The split's answer costs ${r.cost}, against ${scanBest} from the scan.`);
      table(
        root.querySelector(".qs-table div"),
        ["search", "scale", "prefix scale", "points walked", "points the scan alone would walk", "cheapest after it"],
        r.walks.map((w) => [
          w.search,
          String(w.scale),
          w.prefix_scale ?? "",
          fmt(w.points),
          w.scan_points ? `≈${count(w.scan_points)}` : "",
          w.best === null ? "–" : String(w.best),
        ])
      );
    }

    root.addEventListener("click", (ev) => {
      const b = ev.target.closest("button[data-epsilon]");
      if (!b) return;
      state.epsilon = +b.dataset.epsilon;
      render();
    });
    return { root, render };
  }

  // --- figure 5: the price decides the winner --------------------------------------------

  const PRICE_NAMES = { "1,1": "gate count", "2,3": "power", "2,5": "unfolding", "1,3": "catalysis" };

  function priceFigure(prices) {
    const root = document.getElementById("qs-price");
    if (!root) return null;
    const input = root.querySelector('[data-input="price"]');
    const sqrt = prices.rows.filter((r) => r.price[1] !== null);
    const cliffordT = prices.rows.find((r) => r.price[1] === null);
    input.max = String(sqrt.length - 1);
    const tip = tooltip(root.querySelector(".qs-chart"));
    // distinct circuits by (n_T, n_sqrtT), with the prices each wins
    const dots = [];
    for (const r of prices.rows) {
      let d = dots.find((v) => v.n_t === r.n_t && v.n_q === r.n_q);
      if (!d) dots.push((d = { n_t: r.n_t, n_q: r.n_q, word: r.word, wins: [] }));
      d.wins.push(r.price);
    }
    const priceLabel = (p) =>
      p[1] === null ? "Clifford+T" : `(${p[0]},${p[1]})${PRICE_NAMES[p.join(",")] ? ` ${PRICE_NAMES[p.join(",")]}` : ""}`;

    function wordNode(word) {
      const node = el("code", { class: "qs-gates" });
      for (const ch of word) {
        const cls = "Tt".includes(ch) ? "qs-gate-t" : "Qq".includes(ch) ? "qs-gate-q" : "qs-gate-c";
        node.append(el("span", { class: cls }, ch));
      }
      return node;
    }

    function render() {
      const row = sqrt[+input.value];
      const [t, q] = row.price;
      root.querySelector('[data-out="price"]').textContent = `${+(q / t).toFixed(2)}, price ${priceLabel(row.price)}`;
      const ink = color("--pst-color-text-base");
      const muted = color("--pst-color-text-muted");
      const qs = color("--qs-quicksynth");
      keyLegend(root.querySelector(".qs-legend"), [
        { label: "QuickSynth's answer at some price", color: qs, kind: "dot" },
        { label: "the answer at the chosen price", color: ink, kind: "ring" },
        { label: "circuits of the winner's cost", color: ink, kind: "line" },
        { label: "cheaper", color: qs, kind: "bar", opacity: 0.15 },
      ]);

      const chart = root.querySelector(".qs-chart");
      const width = chart.clientWidth;
      const m = { top: 14, right: 20, bottom: 44, left: 52 };
      const height = Math.round(Math.min(380, Math.max(260, width * 0.55)));
      const xMax = d3.max(dots, (d) => d.n_t) + 3;
      const yMax = d3.max(dots, (d) => d.n_q) + 3;
      const x = d3.scaleLinear().domain([0, xMax]).range([m.left, width - m.right]);
      const y = d3.scaleLinear().domain([0, yMax]).range([height - m.bottom, m.top]);
      const svg = d3
        .create("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("role", "img")
        .attr("aria-label", `Answers by number of T and √T gates; at price (${t},${q}) the winner has ${row.n_t} T and ${row.n_q} √T gates, cost ${row.cost}.`);
      svg
        .append("g")
        .attr("class", "grid")
        .selectAll("line")
        .data(y.ticks(5))
        .join("line")
        .attr("x1", m.left)
        .attr("x2", width - m.right)
        .attr("y1", y)
        .attr("y2", y);
      svg.append("g").attr("transform", `translate(${m.left},0)`).call(d3.axisLeft(y).ticks(5)).call(axisStyle);
      svg
        .append("g")
        .attr("transform", `translate(0,${height - m.bottom})`)
        .call(d3.axisBottom(x).ticks(7))
        .call(axisStyle);
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("x", (m.left + width - m.right) / 2)
        .attr("y", height - 8)
        .attr("text-anchor", "middle")
        .text("T gates");
      svg
        .append("text")
        .attr("class", "axis-title")
        .attr("transform", `translate(14,${(m.top + height - m.bottom) / 2}) rotate(-90)`)
        .attr("text-anchor", "middle")
        .text("√T gates");

      // t x + q y = cost, clipped to the plot; the cheaper side shaded
      const c = row.cost;
      const clip = svg.append("defs").append("clipPath").attr("id", "qs-price-clip");
      clip
        .append("rect")
        .attr("x", m.left)
        .attr("y", m.top)
        .attr("width", width - m.left - m.right)
        .attr("height", height - m.top - m.bottom);
      const g = svg.append("g").attr("clip-path", "url(#qs-price-clip)");
      const ends = [
        [0, c / q],
        [c / t, 0],
      ];
      g.append("path")
        .attr("d", `M${x(0)},${y(0)}L${x(ends[0][0])},${y(ends[0][1])}L${x(ends[1][0])},${y(ends[1][1])}Z`)
        .attr("fill", qs)
        .attr("fill-opacity", 0.12);
      g.append("line")
        .attr("x1", x(ends[0][0]))
        .attr("y1", y(ends[0][1]))
        .attr("x2", x(ends[1][0]))
        .attr("y2", y(ends[1][1]))
        .attr("stroke", ink)
        .attr("stroke-width", 1.5);

      const pts = dots.map((d) => ({ ...d, x: x(d.n_t), y: y(d.n_q), cost: t * d.n_t + q * d.n_q }));
      const win = pts.find((d) => d.n_t === row.n_t && d.n_q === row.n_q);
      svg
        .append("g")
        .selectAll("circle")
        .data(pts)
        .join("circle")
        .attr("cx", (d) => d.x)
        .attr("cy", (d) => d.y)
        .attr("r", 5.5)
        .attr("fill", qs)
        .attr("stroke", color("--pst-color-background"))
        .attr("stroke-width", 1.5);
      svg
        .append("circle")
        .attr("cx", win.x)
        .attr("cy", win.y)
        .attr("r", 10)
        .attr("fill", "none")
        .attr("stroke", ink)
        .attr("stroke-width", 2);
      svg
        .append("g")
        .selectAll("text")
        .data(pts)
        .join("text")
        .attr("class", "label halo")
        .attr("stroke", color("--pst-color-background"))
        .attr("x", (d) => (d.x > width - 120 ? d.x - 12 : d.x + 12))
        .attr("y", (d) => d.y - 10)
        .attr("text-anchor", (d) => (d.x > width - 120 ? "end" : "start"))
        .text((d) => (d.wins.some((p) => p[1] === null) && d.wins.length === 1 ? `Clifford+T, ${d.cost}` : String(d.cost)));
      svg
        .append("g")
        .selectAll("circle")
        .data(pts)
        .join("circle")
        .attr("class", "hit")
        .attr("cx", (d) => d.x)
        .attr("cy", (d) => d.y)
        .attr("r", 13)
        .attr("tabindex", 0)
        .attr("aria-label", (d) => `${d.n_t} T, ${d.n_q} √T: cost ${d.cost} at (${t},${q}); answer at ${d.wins.map(priceLabel).join(", ")}`)
        .on("pointerenter focus", (ev, d) =>
          tip.show(d.x, d.y, `${d.n_t} T + ${d.n_q} √T`, [
            { value: String(d.cost), label: `cost at (${t},${q})` },
            { value: d.wins.map(priceLabel).join(", "), label: "" },
          ])
        )
        .on("pointerleave blur", () => tip.hide());
      chart.replaceChildren(svg.node(), chart.querySelector(".qs-tooltip") || el("div"));
      tip.hide();

      const wp = root.querySelector(".qs-word");
      wp.replaceChildren(
        document.createTextNode(`The answer at (${t},${q}): ${row.n_t} T + ${row.n_q} √T = cost ${row.cost}. `),
        wordNode(row.word)
      );
      root.querySelector(".qs-caption").textContent =
        `The numbers are each answer's cost at (${t},${q}). The Clifford+T answer, ${cliffordT.n_t} T gates, is ` +
        `certified optimal over Clifford+T; at every price shown a circuit with √T gates is cheaper.`;
      table(
        root.querySelector(".qs-table div"),
        ["T", "√T", ...sqrt.map((r) => `cost at (${r.price})`), "answer at"],
        dots.map((d) => [
          String(d.n_t),
          String(d.n_q),
          ...sqrt.map((r) => String(r.price[0] * d.n_t + r.price[1] * d.n_q)),
          d.wins.map(priceLabel).join(", "),
        ])
      );
    }

    input.addEventListener("input", render);
    return { root, render };
  }

  // --- start ---------------------------------------------------------------------------

  const get = (name, kind) =>
    fetch(new URL(name, BASE)).then((r) => {
      if (!r.ok) throw new Error(`${r.status} ${name}`);
      return r[kind]();
    });
  Promise.all([get("mechanism.json", "json"), get("points.bin", "arrayBuffer")])
    .then(([meta, buffer]) => {
      const figures = [
        sphereFigure(meta.sphere, buffer),
        scanFigure(meta.scan),
        walkFigure(),
        splitFigure(meta.split),
        priceFigure(meta.prices),
      ].filter(Boolean);
      live(() => figures.forEach((f) => f.render()), figures[0].root);
    })
    .catch((err) => {
      for (const id of ["qs-sphere", "qs-scan", "qs-split", "qs-price"]) {
        const n = document.getElementById(id);
        if (n) n.querySelector(".qs-caption").textContent = `The figure could not load its data (${err.message}).`;
      }
    });
})();
