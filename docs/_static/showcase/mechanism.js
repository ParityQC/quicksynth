/* Copyright (c) 2026 Parity Quantum Computing GmbH
   SPDX-License-Identifier: MIT */

/* "How QuickSynth finds a circuit": the Bloch-sphere picture of the circuits
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

      // legend
      const leg = root.querySelector(".qs-legend");
      leg.replaceChildren();
      for (const [label, c, dash, kind] of [
        ["cheapest circuit found", qs, null, "dot"],
        ["lower bound at this scale and deeper", bnd, "4 3", "line"],
        ["lattice points walked", muted, null, "bar"],
      ]) {
        const span = el("span");
        const svg = d3.create("svg").attr("width", 24).attr("height", 12).attr("aria-hidden", "true");
        if (kind === "dot") svg.append("circle").attr("cx", 12).attr("cy", 6).attr("r", 4).attr("fill", c);
        if (kind === "line")
          svg.append("line").attr("x1", 1).attr("x2", 23).attr("y1", 6).attr("y2", 6).attr("stroke", c).attr("stroke-width", 2).attr("stroke-dasharray", dash);
        if (kind === "bar") svg.append("rect").attr("x", 7).attr("y", 1).attr("width", 10).attr("height", 10).attr("rx", 2).attr("fill", c).attr("fill-opacity", 0.55);
        span.append(svg.node(), document.createTextNode(label));
        leg.append(span);
      }

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

  // --- start ---------------------------------------------------------------------------

  const get = (name, kind) =>
    fetch(new URL(name, BASE)).then((r) => {
      if (!r.ok) throw new Error(`${r.status} ${name}`);
      return r[kind]();
    });
  Promise.all([get("mechanism.json", "json"), get("points.bin", "arrayBuffer")])
    .then(([meta, buffer]) => {
      const figures = [sphereFigure(meta.sphere, buffer), scanFigure(meta.scan)].filter(Boolean);
      live(() => figures.forEach((f) => f.render()), figures[0].root);
    })
    .catch((err) => {
      for (const id of ["qs-sphere", "qs-scan"]) {
        const n = document.getElementById(id);
        if (n) n.querySelector(".qs-caption").textContent = `The figure could not load its data (${err.message}).`;
      }
    });
})();
