// Hero figure: the solvent accessible surface of ubiquitin (1UBQ) as a rotating dot cloud.
// This is Shrake-Rupley run in the browser: test points are placed on every
// probe-inflated atom, and only those not covered by a neighbouring atom are kept.
(() => {
  const cv = document.getElementById("surface");
  const D = window.ZSASA_UBQ;
  if (!cv || !D) return;

  const N = D.length / 4;
  const VDW = [1.7, 1.55, 1.52, 1.8]; // C, N, O, S (Bondi)
  const PROBE = 1.4;
  const PTS = 320; // test points per atom
  const TILT = -0.3;
  const BINS = 5; // depth buckets, back to front

  const X = new Float32Array(N), Y = new Float32Array(N), Z = new Float32Array(N), R = new Float32Array(N);
  let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9, z0 = 1e9, z1 = -1e9;
  for (let i = 0; i < N; i++) {
    X[i] = D[i * 4] / 10;
    Y[i] = D[i * 4 + 1] / 10;
    Z[i] = D[i * 4 + 2] / 10;
    R[i] = VDW[D[i * 4 + 3]] + PROBE;
    x0 = Math.min(x0, X[i]); x1 = Math.max(x1, X[i]);
    y0 = Math.min(y0, Y[i]); y1 = Math.max(y1, Y[i]);
    z0 = Math.min(z0, Z[i]); z1 = Math.max(z1, Z[i]);
  }
  const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2, cz = (z0 + z1) / 2;
  const half = (x1 - x0) / 2 + 3.2;

  // Unit sphere test points (golden spiral).
  const ux = new Float32Array(PTS), uy = new Float32Array(PTS), uz = new Float32Array(PTS);
  for (let k = 0; k < PTS; k++) {
    const y = 1 - (2 * (k + 0.5)) / PTS, r = Math.sqrt(1 - y * y), a = k * 2.399963;
    ux[k] = r * Math.cos(a);
    uy[k] = y;
    uz[k] = r * Math.sin(a);
  }

  const px = [], py = [], pz = [];
  for (let i = 0; i < N; i++) {
    const near = [];
    for (let j = 0; j < N; j++) {
      if (j === i) continue;
      const d = R[i] + R[j], dx = X[i] - X[j], dy = Y[i] - Y[j], dz = Z[i] - Z[j];
      if (dx * dx + dy * dy + dz * dz < d * d) near.push(j);
    }
    for (let k = 0; k < PTS; k++) {
      const x = X[i] + R[i] * ux[k], y = Y[i] + R[i] * uy[k], z = Z[i] + R[i] * uz[k];
      let buried = false;
      for (const j of near) {
        const dx = x - X[j], dy = y - Y[j], dz = z - Z[j];
        if (dx * dx + dy * dy + dz * dz < R[j] * R[j]) { buried = true; break; }
      }
      if (!buried) { px.push(x - cx); py.push(y - cy); pz.push(z - cz); }
    }
  }
  const M = px.length;

  const set = (sel, text) => {
    const el = document.querySelector(sel);
    if (el) el.textContent = text;
  };
  set("[data-exposed]", M.toLocaleString("en"));
  set("[data-total]", (N * PTS).toLocaleString("en"));

  const ctx = cv.getContext("2d");
  const ct = Math.cos(TILT), st = Math.sin(TILT);
  let W = 0, H = 0, light = false, accent = "#f7a41d";

  function readTheme() {
    light = document.documentElement.dataset.theme === "light";
    accent = getComputedStyle(cv).getPropertyValue("--accent").trim() || accent;
  }

  function resize() {
    const r = cv.getBoundingClientRect(), dpr = Math.min(window.devicePixelRatio || 1, 2);
    W = r.width;
    H = r.height;
    cv.width = Math.round(W * dpr);
    cv.height = Math.round(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  const paths = [];
  function draw(angle, t) {
    const s = Math.min(W / (2 * (half + 2)), H / 42);
    const ca = Math.cos(angle), sa = Math.sin(angle);
    const scan = half * 0.92 * Math.sin(t * 0.00042);
    ctx.clearRect(0, 0, W, H);
    for (let b = 0; b <= BINS; b++) paths[b] = new Path2D();
    for (let i = 0; i < M; i++) {
      // Spin about the molecule's long (x) axis, then tilt in the screen plane.
      const x = px[i], y = py[i] * ca - pz[i] * sa, z = py[i] * sa + pz[i] * ca;
      const sx = W / 2 + (x * ct - y * st) * s, sy = H / 2 - (x * st + y * ct) * s;
      if (Math.abs(x - scan) < 0.55) {
        paths[BINS].rect(sx - 1.3, sy - 1.3, 2.6, 2.6);
        continue;
      }
      const b = Math.max(0, Math.min(BINS - 1, Math.floor(((z + 16) / 32) * BINS)));
      const d = 0.8 + b * 0.22;
      paths[b].rect(sx - d / 2, sy - d / 2, d, d);
    }
    ctx.globalCompositeOperation = light ? "source-over" : "lighter";
    ctx.fillStyle = light ? "#a35f00" : accent;
    for (let b = 0; b < BINS; b++) {
      ctx.globalAlpha = 0.2 + b * 0.2;
      ctx.fill(paths[b]);
    }
    ctx.globalAlpha = 1;
    ctx.fillStyle = light ? "#0e0e0e" : "#fff6e0";
    ctx.shadowColor = accent;
    ctx.shadowBlur = light ? 0 : 12;
    ctx.fill(paths[BINS]);
    ctx.shadowBlur = 0;

    // Scan plane marker.
    const mx = W / 2 + scan * ct * s, my = H / 2 - scan * st * s, len = 19.5 * s;
    ctx.globalCompositeOperation = "source-over";
    ctx.strokeStyle = accent;
    ctx.globalAlpha = 0.45;
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 5]);
    ctx.beginPath();
    ctx.moveTo(mx - st * len, my - ct * len);
    ctx.lineTo(mx + st * len, my + ct * len);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.globalAlpha = 1;
    set("[data-scan]", "x = " + (scan >= 0 ? "+" : "−") + Math.abs(scan).toFixed(1) + " Å");
  }

  const still = matchMedia("(prefers-reduced-motion: reduce)").matches;
  let angle = 0.9, last = 0, visible = true, raf = 0, drag = null;

  function tick(t) {
    raf = 0;
    const dt = last ? Math.min(t - last, 64) : 0;
    last = t;
    if (!drag && !still) angle += dt * 0.00028;
    draw(angle, still ? 0 : t);
    if (visible && !document.hidden && (!still || drag)) raf = requestAnimationFrame(tick);
  }

  const kick = () => {
    if (!raf) { last = 0; raf = requestAnimationFrame(tick); }
  };

  cv.addEventListener("pointerdown", (e) => {
    drag = { y: e.clientY, angle };
    cv.setPointerCapture(e.pointerId);
    kick();
  });
  cv.addEventListener("pointermove", (e) => {
    if (drag) angle = drag.angle + (e.clientY - drag.y) * 0.012;
  });
  const release = () => (drag = null);
  cv.addEventListener("pointerup", release);
  cv.addEventListener("pointercancel", release);

  new IntersectionObserver((es) => {
    visible = es[0].isIntersecting;
    if (visible) kick();
  }).observe(cv);
  document.addEventListener("visibilitychange", kick);
  new ResizeObserver(() => { resize(); kick(); }).observe(cv);
  window.addEventListener("themechange", () => { readTheme(); kick(); });

  readTheme();
  resize();
  kick();
})();
