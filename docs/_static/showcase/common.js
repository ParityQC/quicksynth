/* Copyright (c) 2026 Parity Quantum Computing GmbH
   SPDX-License-Identifier: MIT */

/* Helpers shared by the interactive figures of "How it works". */

window.qsShowcase = (function () {
  "use strict";

  // --- formatting --------------------------------------------------------------------

  const SUPERSCRIPT = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };

  function sup(n) {
    return String(n)
      .split("")
      .map((c) => SUPERSCRIPT[c] ?? c)
      .join("");
  }

  /** 1e-7 as "10⁻⁷", 5.68e-8 as "5.7×10⁻⁸". */
  function eps(e) {
    const k = Math.floor(Math.log10(e) + 1e-9);
    const m = e / 10 ** k;
    return Math.abs(m - 1) < 1e-6 ? `10${sup(k)}` : `${m.toFixed(1)}×10${sup(k)}`;
  }

  function signed(p, digits = 1) {
    const s = Math.abs(p).toFixed(digits);
    return p < 0 ? `−${s}%` : `+${s}%`;
  }

  /** A duration given in hours, in its natural unit. */
  function duration(hours) {
    const s = hours * 3600;
    if (s < 60) return `${s.toPrecision(2)} s`;
    if (s < 3600) return `${(s / 60).toPrecision(2)} min`;
    if (hours < 48) return `${hours.toPrecision(3)} h`;
    return `${(hours / 24).toPrecision(3)} days`;
  }

  function seconds(s) {
    return s < 1 ? `${+(s * 1000).toPrecision(2)} ms` : `${+s.toPrecision(3)} s`;
  }

  function color(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  }

  function el(tag, attrs = {}, text) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text !== undefined) e.textContent = text;
    return e;
  }

  // --- shared pieces -----------------------------------------------------------------

  /** One tooltip per chart; rows of {value, label, color, dash}. */
  function tooltip(chart) {
    let tip = chart.querySelector(":scope > .qs-tooltip");
    if (!tip) {
      tip = el("div", { class: "qs-tooltip", role: "status" });
      tip.hidden = true;
      chart.append(tip);
    }
    return {
      show(x, y, head, rows) {
        tip.replaceChildren(el("div", { class: "qs-tip-head" }, head));
        for (const r of rows) {
          const row = el("div", { class: "qs-tip-row" });
          if (r.color) {
            const key = el("span", { class: "qs-tip-key" });
            key.style.borderTopColor = r.color;
            if (r.dash) key.style.borderTopStyle = "dashed";
            row.append(key);
          }
          row.append(el("b", {}, r.value));
          if (r.label) row.append(el("span", {}, r.label));
          tip.append(row);
        }
        tip.hidden = false;
        const w = chart.clientWidth;
        const tw = tip.offsetWidth;
        tip.style.left = `${Math.min(Math.max(0, x + 14), w - tw)}px`;
        tip.style.top = `${Math.max(0, y - tip.offsetHeight - 10)}px`;
      },
      hide() {
        tip.hidden = true;
      },
    };
  }

  function table(node, head, rows) {
    const t = el("table", { class: "table" });
    const tr = el("tr");
    for (const h of head) tr.append(el("th", {}, h));
    t.append(el("thead"));
    t.firstChild.append(tr);
    const body = el("tbody");
    for (const r of rows) {
      const row = el("tr");
      for (const v of r) row.append(el("td", {}, v));
      body.append(row);
    }
    t.append(body);
    node.replaceChildren(t);
  }

  function axisStyle(g) {
    g.attr("class", "axis");
    g.selectAll(".tick line").attr("stroke-opacity", 0.8);
  }

  function pressed(group, attr, value) {
    for (const b of group.querySelectorAll(`[${attr}]`)) {
      b.setAttribute("aria-pressed", String(b.getAttribute(attr) === String(value)));
    }
  }


  /** Render now, and again on a theme switch (colours) or a width change of root. */
  function live(renderAll, root) {
    renderAll();
    new MutationObserver(renderAll).observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["data-theme"],
    });
    let width = 0;
    let timer;
    new ResizeObserver((entries) => {
      const w = Math.round(entries[0].contentRect.width);
      if (w === width) return;
      width = w;
      clearTimeout(timer);
      timer = setTimeout(renderAll, 120);
    }).observe(root);
  }

  return { sup, eps, signed, duration, seconds, color, el, tooltip, table, axisStyle, pressed, live };
})();
