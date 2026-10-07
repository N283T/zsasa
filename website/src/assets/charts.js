// Benchmark charts. Each <figure class="chart"> carries a JSON spec exported by
// website/scripts/export_benchmarks.py; this file draws it as SVG.
// Colours and fonts come from site.css, so a theme switch needs no redraw.
(() => {
  const NS = "http://www.w3.org/2000/svg";
  const DASH = { dash: "6 4", dot: "1.5 4.5" };
  const TAU = Math.PI * 2;

  const svg = (name, attrs, parent) => {
    const e = document.createElementNS(NS, name);
    for (const k in attrs) if (attrs[k] != null) e.setAttribute(k, attrs[k]);
    parent?.append(e);
    return e;
  };
  const h = (tag, cls, text) => {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text != null) e.textContent = text;
    return e;
  };
  const text = (parent, x, y, str, attrs) => {
    const t = svg("text", { x, y, ...attrs }, parent);
    t.textContent = str;
    return t;
  };
  // Format with fixed decimals; values that round to zero never print as "-0".
  const num = (v, d = 0) => {
    const r = Number(v.toFixed(d));
    return (r === 0 ? 0 : r).toLocaleString("en", { minimumFractionDigits: d, maximumFractionDigits: d });
  };

  // ----- scales -----

  function niceStep(span, count) {
    const raw = span / Math.max(1, count), p = 10 ** Math.floor(Math.log10(raw)), m = raw / p;
    return (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * p;
  }

  function linear(min, max, r0, r1, count) {
    if (min === max) { min -= 1; max += 1; }
    const step = niceStep(max - min, count);
    const d0 = Math.floor(min / step + 1e-9) * step, d1 = Math.ceil(max / step - 1e-9) * step;
    const f = (v) => r0 + ((v - d0) / (d1 - d0)) * (r1 - r0);
    f.ticks = [];
    for (let v = d0; v <= d1 + step / 2; v += step) f.ticks.push(Math.abs(v) < step / 1e6 ? 0 : v);
    f.digits = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
    return f;
  }

  function log(min, max, r0, r1) {
    const l0 = Math.log10(min), l1 = Math.log10(max), pad = (l1 - l0) * 0.07 || 0.2;
    const d0 = l0 - pad, d1 = l1 + pad;
    const f = (v) => r0 + ((Math.log10(v) - d0) / (d1 - d0)) * (r1 - r0);
    const mult = d1 - d0 > 2.4 ? [1] : d1 - d0 > 1.1 ? [1, 3] : [1, 2, 5];
    f.ticks = [];
    for (let e = Math.floor(d0); e <= Math.ceil(d1); e++)
      for (const m of mult) {
        const v = m * 10 ** e;
        if (Math.log10(v) >= d0 && Math.log10(v) <= d1) f.ticks.push(v);
      }
    f.digits = null;
    return f;
  }

  function points(values, r0, r1) {
    const inset = Math.min(48, (r1 - r0) / (values.length * 2));
    const f = (v) => {
      const i = values.indexOf(v);
      return values.length === 1 ? (r0 + r1) / 2 : r0 + inset + (i / (values.length - 1)) * (r1 - r0 - 2 * inset);
    };
    f.ticks = values;
    f.digits = 0;
    return f;
  }

  const tickLabel = (v, scale) => {
    if (scale.digits != null) return num(v, scale.digits);
    return v >= 1 ? num(v, 0) : String(v);
  };

  // ----- marks -----

  function marker(parent, s, x, y, r = 4.5) {
    const cls = `m t-${s.tool}${s.hollow ? " hollow" : ""}`;
    if (s.shape === "diamond")
      return svg("rect", { class: cls, x: x - r * 0.9, y: y - r * 0.9, width: r * 1.8, height: r * 1.8, transform: `rotate(45 ${x} ${y})` }, parent);
    if (s.shape === "square") return svg("rect", { class: cls, x: x - r * 0.85, y: y - r * 0.85, width: r * 1.7, height: r * 1.7 }, parent);
    if (s.shape === "triangle")
      return svg("path", { class: cls, d: `M${x},${y - r * 1.15}L${x + r * 1.1},${y + r * 0.8}L${x - r * 1.1},${y + r * 0.8}Z` }, parent);
    return svg("circle", { class: cls, cx: x, cy: y, r }, parent);
  }

  function key(s, withLine) {
    const k = svg("svg", { class: "chart__key", viewBox: "0 0 26 12", width: 26, height: 12, "aria-hidden": "true" });
    if (withLine) svg("line", { class: `ln t-${s.tool}`, x1: 0, x2: 26, y1: 6, y2: 6, "stroke-dasharray": DASH[s.dash] }, k);
    marker(k, s, 13, 6, 4);
    return k;
  }

  // ----- frame -----

  function frame(c, H, m, xs, ys, xDef, yDef, opts = {}) {
    const g = svg("g", { class: "axis" }, c.svg);
    for (const v of ys.ticks) {
      const y = ys(v);
      svg("line", { class: v === 0 && opts.zero ? "zero" : "grid", x1: m.l, x2: c.W - m.r, y1: y, y2: y }, g);
      text(g, m.l - 8, y + 3.5, tickLabel(v, ys), { "text-anchor": "end" });
    }
    for (const v of xs.ticks) {
      const x = xs(v);
      if (opts.xGrid) svg("line", { class: "grid", x1: x, x2: x, y1: m.t, y2: H - m.b }, g);
      else svg("line", { class: "tick", x1: x, x2: x, y1: H - m.b, y2: H - m.b + 4 }, g);
      text(g, x, H - m.b + 17, tickLabel(v, xs), { "text-anchor": "middle" });
    }
    svg("line", { class: "base", x1: m.l, x2: c.W - m.r, y1: H - m.b, y2: H - m.b }, g);
    text(g, (m.l + c.W - m.r) / 2, H - 6, xDef.label, { class: "title", "text-anchor": "middle" });
    text(g, 12, (m.t + H - m.b) / 2, yDef.label, { class: "title", "text-anchor": "middle", transform: `rotate(-90 12 ${(m.t + H - m.b) / 2})` });
  }

  const leftMargin = (scale) => 30 + Math.max(...scale.ticks.map((v) => tickLabel(v, scale).length)) * 6.7;

  function whisker(g, s, x0, y0, x1, y1) {
    const w = svg("g", { class: `wh t-${s.tool}` }, g), cap = 3;
    svg("line", { x1: x0, y1: y0, x2: x1, y2: y1 }, w);
    if (y0 === y1) {
      svg("line", { x1: x0, x2: x0, y1: y0 - cap, y2: y0 + cap }, w);
      svg("line", { x1, x2: x1, y1: y0 - cap, y2: y0 + cap }, w);
    } else {
      svg("line", { x1: x0 - cap, x2: x0 + cap, y1: y0, y2: y0 }, w);
      svg("line", { x1: x0 - cap, x2: x0 + cap, y1, y2: y1 }, w);
    }
  }

  // ----- tooltip -----

  function showTip(c, x, y, title, rows, foot) {
    const tip = c.tip;
    tip.replaceChildren();
    if (title) tip.append(h("div", "chart__tip-t", title));
    for (const r of rows) {
      const row = h("div", "chart__tip-r");
      if (r.s) row.append(key(r.s, true));
      row.append(h("b", null, r.value), h("span", null, r.label));
      tip.append(row);
    }
    if (foot) tip.append(h("div", "chart__tip-f", foot));
    tip.hidden = false;
    const w = tip.offsetWidth, ht = tip.offsetHeight;
    tip.style.left = Math.max(0, Math.min(c.W - w, x + 14 + w > c.W ? x - w - 14 : x + 14)) + "px";
    tip.style.top = Math.max(0, y - ht - 10 < 0 ? y + 14 : y - ht - 10) + "px";
  }

  const hideTip = (c) => (c.tip.hidden = true);

  const short = (def) => def.label.split(" (")[0].toLowerCase();
  const val = (v, def) => `${num(v, def.digits ?? 0)}${def.unit ? " " + def.unit : ""}`;

  // ----- bar -----

  function bar(c) {
    const { rows } = c.view, def = c.view, rowH = 30, m = { t: 6, b: 40 };
    const labelW = Math.min(c.W * 0.44, Math.max(...rows.map((r) => r.name.length)) * 7.1 + 18);
    const valueW = 22 + num(Math.max(...rows.map((r) => r.v)), def.digits).length * 7;
    const H = m.t + rows.length * rowH + m.b;
    c.svg.setAttribute("height", H);
    const xs = linear(0, Math.max(...rows.map((r) => r.v + (r.e || 0))), labelW, c.W - valueW, (c.W - labelW - valueW) / 90);
    const g = svg("g", { class: "axis" }, c.svg);
    for (const v of xs.ticks) {
      svg("line", { class: "grid", x1: xs(v), x2: xs(v), y1: m.t, y2: H - m.b }, g);
      text(g, xs(v), H - m.b + 17, tickLabel(v, xs), { "text-anchor": "middle" });
    }
    text(g, (labelW + c.W - valueW) / 2, H - 6, def.label, { class: "title", "text-anchor": "middle" });
    rows.forEach((r, i) => {
      const y = m.t + i * rowH + rowH / 2;
      const row = svg("g", { class: "row", tabindex: 0, role: "img", "aria-label": `${r.name}: ${val(r.v, def)}` }, c.svg);
      svg("rect", { class: "hit", x: 0, y: y - rowH / 2, width: c.W, height: rowH }, row);
      text(row, labelW - 10, y + 4, r.name, { class: `name${r.tool === "zsasa" ? " strong" : ""}`, "text-anchor": "end" });
      svg("rect", { class: `bar t-${r.tool}`, x: xs(0), y: y - 7, width: Math.max(1, xs(r.v) - xs(0)), height: 14 }, row);
      if (r.e) whisker(row, { tool: "ink" }, xs(r.v - r.e), y, xs(r.v + r.e), y);
      text(row, xs(r.v + (r.e || 0)) + 8, y + 4, num(r.v, def.digits), { class: "value" });
      const show = () => showTip(c, xs(r.v), y, r.name, [{ value: val(r.v, def) + (r.e ? ` ± ${num(r.e, def.digits)}` : ""), label: "" }], r.note);
      row.addEventListener("pointerenter", show);
      row.addEventListener("focus", show);
      row.addEventListener("pointerleave", () => hideTip(c));
      row.addEventListener("blur", () => hideTip(c));
    });
    return {
      head: ["Tool", def.label, "Std. dev."],
      rows: rows.map((r) => [r.name, num(r.v, def.digits), r.e ? num(r.e, def.digits) : "–"]),
    };
  }

  // ----- scatter -----

  function scatter(c) {
    const pts = c.view.points, xDef = c.spec.x, yDef = c.spec.y, H = 400;
    c.svg.setAttribute("height", H);
    const ext = (k, e) => [Math.min(...pts.map((p) => p[k] - (p[e] || 0))), Math.max(...pts.map((p) => p[k] + (p[e] || 0)))];
    const [x0, x1] = ext("x", "ex"), [y0, y1] = ext("y", "ey");
    const m = { t: 14, r: 18, b: 44, l: 0 };
    const mk = (def, a, b, r0, r1, n) => (def.scale === "log" ? log(a, b, r0, r1) : linear(a - (b - a) * 0.08, b + (b - a) * 0.08, r0, r1, n));
    let ys = mk(yDef, y0, y1, H - m.b, m.t, 6);
    m.l = leftMargin(ys);
    const xs = mk(xDef, x0, x1, m.l, c.W - m.r, (c.W - m.l) / 110);
    frame(c, H, m, xs, ys, xDef, yDef, { xGrid: true });

    const g = svg("g", null, c.svg), placed = [], boxes = [];
    const P = pts.map((p) => ({ p, x: xs(p.x), y: ys(p.y) }));
    for (const q of P) {
      if (q.p.ex) whisker(g, q.p, xs(Math.max(q.p.x - q.p.ex, 1e-9)), q.y, xs(q.p.x + q.p.ex), q.y);
      if (q.p.ey) whisker(g, q.p, q.x, ys(Math.max(q.p.y - q.p.ey, 1e-9)), q.x, ys(q.p.y + q.p.ey));
      boxes.push([q.x - 8, q.y - 8, q.x + 8, q.y + 8]);
    }
    const hit = (a, b) => a[0] < b[2] && a[2] > b[0] && a[1] < b[3] && a[3] > b[1];
    for (const q of P) {
      const w = q.p.name.length * 6.5 + 4;
      const options = [
        [q.x + 11, q.y + 4, "start", [q.x + 9, q.y - 8, q.x + 11 + w, q.y + 8]],
        [q.x - 11, q.y + 4, "end", [q.x - 11 - w, q.y - 8, q.x - 9, q.y + 8]],
        [q.x, q.y - 12, "middle", [q.x - w / 2, q.y - 24, q.x + w / 2, q.y - 9]],
        [q.x, q.y + 21, "middle", [q.x - w / 2, q.y + 9, q.x + w / 2, q.y + 24]],
      ];
      const fit = options.find(([, , , b]) => b[0] >= m.l && b[2] <= c.W - 2 && b[1] >= 0 && b[3] <= H - m.b && !placed.some((o) => hit(b, o)) && !boxes.some((o, i) => P[i] !== q && hit(b, o)));
      if (fit) {
        placed.push(fit[3]);
        text(g, fit[0], fit[1], q.p.name, { class: `lbl${q.p.tool === "zsasa" ? " strong" : ""}`, "text-anchor": fit[2] });
      }
    }
    const marks = P.map((q) => marker(g, q.p, q.x, q.y, 5.5));
    const focus = (q, i) => {
      marks.forEach((e, k) => e.classList.toggle("on", k === i));
      const d = (v, e, def) => val(v, def) + (e ? ` ± ${num(e, def.digits ?? 0)}` : "");
      showTip(c, q.x, q.y, q.p.name, [
        { value: d(q.p.y, q.p.ey, yDef), label: short(yDef) },
        { value: d(q.p.x, q.p.ex, xDef), label: short(xDef) },
      ]);
    };
    c.svg.addEventListener("pointermove", (e) => {
      const r = c.svg.getBoundingClientRect(), px = e.clientX - r.left, py = e.clientY - r.top;
      let best = -1, bd = 40 * 40;
      P.forEach((q, i) => {
        const d = (q.x - px) ** 2 + (q.y - py) ** 2;
        if (d < bd) { bd = d; best = i; }
      });
      if (best < 0) { hideTip(c); marks.forEach((e2) => e2.classList.remove("on")); } else focus(P[best], best);
    });
    c.svg.addEventListener("pointerleave", () => { hideTip(c); marks.forEach((e) => e.classList.remove("on")); });
    return {
      head: ["Tool", yDef.label, xDef.label],
      rows: pts.map((p) => [p.name, num(p.y, yDef.digits), num(p.x, xDef.digits)]),
    };
  }

  // ----- line / band -----

  function line(c, isBand) {
    const series = c.view.series, xDef = c.spec.x, yDef = c.view.y || c.spec.y, H = 360;
    c.svg.setAttribute("height", H);
    if (!series.length) return null;
    const lo = (p) => (isBand ? p[2] : p[1] - (p[2] || 0)), hi = (p) => (isBand ? p[3] : p[1] + (p[2] || 0));
    const all = series.flatMap((s) => s.pts);
    let y0 = Math.min(...all.map(lo)), y1 = Math.max(...all.map(hi));
    const labels = c.W >= 560;
    const m = { t: 14, b: 44, l: 0, r: labels ? Math.min(c.W * 0.34, Math.max(...series.map((s) => s.name.length)) * 6.5 + 22) : 18 };
    let ys;
    if (yDef.scale === "log") ys = log(Math.max(y0, 1e-9), y1, H - m.b, m.t);
    else {
      if (isBand) { y0 = Math.min(y0, 0); y1 = Math.max(y1, 0); } else y0 = 0;
      ys = linear(y0 - (isBand ? (y1 - y0) * 0.06 : 0), y1 + (y1 - y0) * 0.06, H - m.b, m.t, 6);
    }
    m.l = leftMargin(ys);
    const xv = [...new Set(all.map((p) => p[0]))].sort((a, b) => a - b);
    const xs = xDef.scale === "log" ? log(xv[0], xv.at(-1), m.l, c.W - m.r) : xDef.scale === "linear" ? linear(Math.min(0, xv[0]), xv.at(-1), m.l, c.W - m.r, 6) : points(xDef.values || xv, m.l, c.W - m.r);
    if (xDef.scale === "linear" && xDef.values) xs.ticks = xDef.values;
    frame(c, H, m, xs, ys, xDef, yDef, { zero: isBand });

    const g = svg("g", null, c.svg);
    if (isBand)
      for (const s of series)
        svg("path", { class: `area t-${s.tool}`, d: "M" + s.pts.map((p) => `${xs(p[0])},${ys(p[3])}`).join("L") + "L" + [...s.pts].reverse().map((p) => `${xs(p[0])},${ys(p[2])}`).join("L") + "Z" }, g);
    for (const s of series) {
      svg("polyline", { class: `ln t-${s.tool}`, "stroke-dasharray": DASH[s.dash], points: s.pts.map((p) => `${xs(p[0])},${ys(p[1])}`).join(" ") }, g);
      if (!isBand) for (const p of s.pts) if (p[2]) whisker(g, s, xs(p[0]), ys(Math.max(p[1] - p[2], 1e-9)), xs(p[0]), ys(p[1] + p[2]));
    }
    for (const s of series) for (const p of s.pts) marker(g, s, xs(p[0]), ys(p[1]));

    if (labels) {
      const ends = series.map((s) => ({ s, y: ys(s.pts.at(-1)[1]), x: xs(s.pts.at(-1)[0]) })).sort((a, b) => a.y - b.y);
      for (let i = 1; i < ends.length; i++) ends[i].ly = Math.max((ends[i - 1].ly ?? ends[i - 1].y) + 13, ends[i].y);
      ends[0].ly ??= ends[0].y;
      const over = (ends.at(-1).ly ?? 0) - (H - m.b);
      if (over > 0) for (const e of ends) e.ly -= over;
      for (const e of ends) {
        if (Math.abs(e.ly - e.y) > 3) svg("line", { class: "lead", x1: e.x + 7, y1: e.y, x2: e.x + 13, y2: e.ly }, g);
        text(g, e.x + 15, e.ly + 4, e.s.name, { class: `lbl${e.s.tool === "zsasa" ? " strong" : ""}` });
      }
    }

    const cross = svg("line", { class: "cross", y1: m.t, y2: H - m.b, visibility: "hidden" }, c.svg);
    let at = -1;
    const show = (i) => {
      at = i;
      const x = xv[i], px = xs(x);
      cross.setAttribute("x1", px);
      cross.setAttribute("x2", px);
      cross.setAttribute("visibility", "visible");
      const rows = [];
      let label = null, top = H;
      for (const s of series) {
        const p = s.pts.find((q) => q[0] === x);
        if (!p) continue;
        label ??= p[isBand ? 4 : 3];
        top = Math.min(top, ys(p[1]));
        rows.push({ s, value: isBand ? `${num(p[1], yDef.digits)} (${num(p[2], yDef.digits)} to ${num(p[3], yDef.digits)})` : val(p[1], yDef) + (p[2] ? ` ± ${num(p[2], yDef.digits)}` : ""), label: s.name });
      }
      showTip(c, px, top, `${label ? label + " · " : ""}${xDef.label} ${num(x, 0)}`, rows, isBand ? `Median (5th to 95th percentile), ${yDef.unit}` : null);
    };
    const hide = () => { at = -1; cross.setAttribute("visibility", "hidden"); hideTip(c); };
    c.svg.addEventListener("pointermove", (e) => {
      const px = e.clientX - c.svg.getBoundingClientRect().left;
      let best = 0;
      xv.forEach((v, i) => { if (Math.abs(xs(v) - px) < Math.abs(xs(xv[best]) - px)) best = i; });
      if (best !== at) show(best);
    });
    c.svg.addEventListener("pointerleave", hide);
    c.svg.setAttribute("tabindex", 0);
    c.svg.addEventListener("keydown", (e) => {
      if (e.key === "ArrowRight") show(Math.min(xv.length - 1, at + 1));
      else if (e.key === "ArrowLeft") show(Math.max(0, at < 0 ? 0 : at - 1));
      else if (e.key === "Escape") hide();
      else return;
      e.preventDefault();
    });
    c.svg.addEventListener("blur", hide);

    for (const s of series) {
      const item = h("span", "chart__legend-i");
      item.append(key(s, true), h("span", null, s.name));
      c.legend.append(item);
    }
    const cell = (s, x) => {
      const p = s.pts.find((q) => q[0] === x);
      return !p ? "–" : isBand ? `${num(p[1], yDef.digits)} (${num(p[2], yDef.digits)} to ${num(p[3], yDef.digits)})` : num(p[1], yDef.digits);
    };
    return { head: [xDef.label, ...series.map((s) => s.name)], rows: xv.map((x) => [num(x, 0), ...series.map((s) => cell(s, x))]) };
  }

  // ----- dense scatter -----

  function dense(c) {
    const v = c.view, xDef = c.spec.x, yDef = c.spec.y, H = 400, ids = c.spec.ids;
    c.svg.setAttribute("height", H);
    const views = Object.values(c.spec.views);
    const yAll = views.flatMap((w) => [Math.min(...w.y), Math.max(...w.y)]);
    const m = { t: 14, r: 18, b: 44, l: 0 };
    const ys = linear(Math.min(...yAll, 0), Math.max(...yAll, 0), H - m.b, m.t, 6);
    m.l = leftMargin(ys);
    const xAll = views.flatMap((w) => [Math.min(...w.x), Math.max(...w.x)]);
    const xs = log(Math.min(...xAll), Math.max(...xAll), m.l, c.W - m.r);
    frame(c, H, m, xs, ys, xDef, yDef, { zero: true });
    const px = v.x.map(xs), py = v.y.map(ys);
    let d = "";
    for (let i = 0; i < px.length; i++) d += `M${px[i].toFixed(1)},${py[i].toFixed(1)}h0`;
    svg("path", { class: `dots t-${v.tool}`, d }, c.svg);
    const ring = svg("circle", { class: `ring t-${v.tool}`, r: 5, visibility: "hidden" }, c.svg);
    c.svg.addEventListener("pointermove", (e) => {
      const r = c.svg.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
      let best = -1, bd = 30 * 30;
      for (let i = 0; i < px.length; i++) {
        const q = (px[i] - mx) ** 2 + (py[i] - my) ** 2;
        if (q < bd) { bd = q; best = i; }
      }
      if (best < 0) { ring.setAttribute("visibility", "hidden"); return hideTip(c); }
      ring.setAttribute("cx", px[best]);
      ring.setAttribute("cy", py[best]);
      ring.setAttribute("visibility", "visible");
      showTip(c, px[best], py[best], ids[best], [
        { value: (v.y[best] > 0 ? "+" : "") + num(v.y[best], 3) + " %", label: "difference" },
        { value: val(v.x[best], xDef), label: "FreeSASA" },
      ]);
    });
    c.svg.addEventListener("pointerleave", () => { ring.setAttribute("visibility", "hidden"); hideTip(c); });
    c.stats.textContent = v.stats;
    return null;
  }

  const RENDER = { bar, scatter, line: (c) => line(c, false), band: (c) => line(c, true), dense };

  // ----- figure -----

  function init(fig) {
    const spec = JSON.parse(fig.querySelector('script[type="application/json"]').textContent);
    const body = fig.querySelector(".chart__body"), legend = fig.querySelector(".chart__legend");
    const stats = fig.querySelector(".chart__stats"), data = fig.querySelector(".chart__data");
    const tip = h("div", "chart__tip");
    tip.hidden = true;
    const state = (spec.controls || []).map(() => 0);
    let width = 0;

    function render() {
      const W = Math.floor(body.clientWidth);
      if (!W) return;
      width = W;
      const view = spec.views[state.join("|") || "0"];
      const el = svg("svg", { class: "chart__svg", width: W, role: "img", "aria-label": spec.title });
      legend.replaceChildren();
      stats.textContent = "";
      body.replaceChildren(el, tip);
      tip.hidden = true;
      let table = null;
      if (!view || (view.series && !view.series.length)) {
        el.setAttribute("height", 120);
        text(el, W / 2, 64, "Not measured for this combination.", { class: "empty", "text-anchor": "middle" });
      } else table = RENDER[spec.type]({ spec, view, W, svg: el, tip, legend, stats });
      data.hidden = !table;
      if (table) {
        const t = h("table"), head = t.createTHead().insertRow();
        for (const cell of table.head) head.append(h("th", null, cell));
        const tb = t.createTBody();
        for (const r of table.rows) {
          const tr = tb.insertRow();
          for (const cell of r) tr.insertCell().textContent = cell;
        }
        data.querySelector(".table").replaceChildren(t);
      }
    }

    const controls = fig.querySelector(".chart__controls");
    (spec.controls || []).forEach((ctl, i) => {
      const group = h("div", "chart__ctl");
      group.setAttribute("role", "group");
      group.setAttribute("aria-label", ctl.label);
      group.append(h("span", "label", ctl.label));
      const seg = h("div", "seg");
      ctl.options.forEach((name, j) => {
        const b = h("button", null, name);
        b.type = "button";
        b.setAttribute("aria-pressed", String(j === 0));
        b.addEventListener("click", () => {
          state[i] = j;
          [...seg.children].forEach((x, k) => x.setAttribute("aria-pressed", String(k === j)));
          render();
        });
        seg.append(b);
      });
      group.append(seg);
      controls.append(group);
    });

    new ResizeObserver(() => {
      if (Math.floor(body.clientWidth) !== width) render();
    }).observe(body);
    document.fonts?.ready.then(render);
    render();
  }

  document.querySelectorAll("figure.chart").forEach(init);
})();
