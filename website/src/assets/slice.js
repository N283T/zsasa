// Hero figure: a moving cross-section through ubiquitin (1UBQ).
// Each atom contributes a circle of radius sqrt((r + probe)^2 - dz^2) to the plane.
// Test points on those circles that fall outside every other circle are solvent
// accessible: the 2D analogue of the Shrake-Rupley algorithm.
(() => {
  const cv = document.getElementById("slice");
  const D = window.ZSASA_UBQ;
  if (!cv || !D) return;

  const ctx = cv.getContext("2d");
  const N = D.length / 4;
  const VDW = [1.7, 1.55, 1.52, 1.8]; // C, N, O, S (Bondi)
  const PROBE = 1.4;
  const STEP = 0.42; // test-point spacing along each circle, in Å
  const PERIOD = 26000; // ms per full sweep

  const ax = new Float32Array(N);
  const ay = new Float32Array(N);
  const az = new Float32Array(N);
  const ar = new Float32Array(N);
  let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9, z0 = 1e9, z1 = -1e9;
  for (let i = 0; i < N; i++) {
    ax[i] = D[i * 4] / 10;
    ay[i] = D[i * 4 + 1] / 10;
    az[i] = D[i * 4 + 2] / 10;
    ar[i] = VDW[D[i * 4 + 3]];
    x0 = Math.min(x0, ax[i]); x1 = Math.max(x1, ax[i]);
    y0 = Math.min(y0, ay[i]); y1 = Math.max(y1, ay[i]);
    z0 = Math.min(z0, az[i]); z1 = Math.max(z1, az[i]);
  }
  const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
  const hx = (x1 - x0) / 2 + 4.6, hy = (y1 - y0) / 2 + 4.6;
  const zMid = (z0 + z1) / 2, zAmp = ((z1 - z0) / 2) * 0.86;

  const out = {
    z: document.querySelector("[data-z]"),
    atoms: document.querySelector("[data-atoms]"),
    pts: document.querySelector("[data-pts]"),
    len: document.querySelector("[data-len]"),
    scrub: document.querySelector("[data-scrub]"),
  };

  let W = 0, H = 0, scale = 1, ink = "#000", accent = "#f7a41d";

  function readTheme() {
    const s = getComputedStyle(cv);
    ink = s.getPropertyValue("--ink").trim();
    accent = s.getPropertyValue("--accent").trim();
  }

  function resize() {
    const r = cv.getBoundingClientRect();
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    W = r.width;
    H = r.height;
    cv.width = Math.round(W * dpr);
    cv.height = Math.round(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    scale = Math.min(W / (2 * hx), H / (2 * hy));
  }

  // Scratch buffers for atoms intersecting the current plane.
  const sx = new Float32Array(N), sy = new Float32Array(N);
  const sR = new Float32Array(N), sR2 = new Float32Array(N), sV = new Float32Array(N);

  function draw(z) {
    let n = 0;
    for (let i = 0; i < N; i++) {
      const dz = az[i] - z;
      const R = ar[i] + PROBE;
      if (dz > -R && dz < R) {
        sx[n] = W / 2 + (ax[i] - cx) * scale;
        sy[n] = H / 2 - (ay[i] - cy) * scale;
        const rs = Math.sqrt(R * R - dz * dz) * scale;
        sR[n] = rs;
        sR2[n] = rs * rs - 0.01;
        const v = ar[i] * ar[i] - dz * dz;
        sV[n] = v > 0 ? Math.sqrt(v) * scale : 0;
        n++;
      }
    }

    ctx.clearRect(0, 0, W, H);

    // Probe-inflated region (where the probe centre cannot go).
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      ctx.moveTo(sx[i] + sR[i], sy[i]);
      ctx.arc(sx[i], sy[i], sR[i], 0, 6.2832);
    }
    ctx.globalAlpha = 0.1;
    ctx.fillStyle = accent;
    ctx.fill();

    // van der Waals sections.
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      if (sV[i] <= 0) continue;
      ctx.moveTo(sx[i] + sV[i], sy[i]);
      ctx.arc(sx[i], sy[i], sV[i], 0, 6.2832);
    }
    ctx.globalAlpha = 0.13;
    ctx.fillStyle = ink;
    ctx.fill();
    ctx.globalAlpha = 0.34;
    ctx.strokeStyle = ink;
    ctx.lineWidth = 0.75;
    ctx.stroke();

    // Test points.
    const step = STEP * scale;
    let total = 0, exposed = 0, length = 0;
    const lit = new Path2D();
    ctx.globalAlpha = 0.22;
    ctx.fillStyle = ink;
    for (let i = 0; i < n; i++) {
      const m = Math.max(8, Math.round((6.2832 * sR[i]) / step));
      const da = 6.2832 / m;
      let e = 0;
      for (let k = 0; k < m; k++) {
        const px = sx[i] + sR[i] * Math.cos(k * da);
        const py = sy[i] + sR[i] * Math.sin(k * da);
        let buried = false;
        for (let j = 0; j < n; j++) {
          if (j === i) continue;
          const dx = px - sx[j], dy = py - sy[j];
          if (dx * dx + dy * dy < sR2[j]) { buried = true; break; }
        }
        if (buried) {
          ctx.fillRect(px - 0.5, py - 0.5, 1, 1);
        } else {
          lit.moveTo(px + 1.5, py);
          lit.arc(px, py, 1.5, 0, 6.2832);
          e++;
        }
      }
      total += m;
      exposed += e;
      length += (e / m) * 6.2832 * (sR[i] / scale);
    }
    ctx.globalAlpha = 1;
    ctx.fillStyle = accent;
    ctx.fill(lit);

    // 10 Å scale bar.
    const bar = 10 * scale, bx = 16, by = H - 16;
    ctx.strokeStyle = ink;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(bx, by - 4); ctx.lineTo(bx, by); ctx.lineTo(bx + bar, by); ctx.lineTo(bx + bar, by - 4);
    ctx.stroke();
    ctx.fillStyle = ink;
    ctx.font = '500 11px "JetBrains Mono", ui-monospace, monospace';
    ctx.fillText("10 Å", bx, by - 9);

    // Probe, drawn to scale.
    const pr = PROBE * scale, pxr = W - 16 - pr;
    ctx.strokeStyle = accent;
    ctx.beginPath();
    ctx.arc(pxr, by - pr, pr, 0, 6.2832);
    ctx.stroke();
    ctx.textAlign = "right";
    ctx.fillText("probe 1.4 Å", pxr - pr - 8, by - pr + 4);
    ctx.textAlign = "left";

    if (out.z) out.z.textContent = (z >= 0 ? "+" : "−") + Math.abs(z).toFixed(1) + " Å";
    if (out.atoms) out.atoms.textContent = n;
    if (out.pts) out.pts.innerHTML = exposed.toLocaleString("en") + " <i>/ " + total.toLocaleString("en") + "</i>";
    if (out.len) out.len.textContent = length.toFixed(1) + " Å";
  }

  const still = matchMedia("(prefers-reduced-motion: reduce)").matches;
  let phase = 0.6, last = 0, visible = true, held = 0, raf = 0;

  const zAt = (p) => zMid + zAmp * Math.sin(p * 6.2832);

  function tick(t) {
    raf = 0;
    const dt = last ? Math.min(t - last, 64) : 0;
    last = t;
    if (t > held) {
      phase = (phase + dt / PERIOD) % 1;
      if (out.scrub) out.scrub.value = String(Math.round(((zAt(phase) - zMid) / zAmp) * 500 + 500));
    }
    draw(zAt(phase));
    if (visible && !still && !document.hidden) raf = requestAnimationFrame(tick);
  }

  const kick = () => {
    if (!raf) { last = 0; raf = requestAnimationFrame(tick); }
  };

  out.scrub?.addEventListener("input", () => {
    const s = (out.scrub.value - 500) / 500; // -1..1
    const a = Math.asin(Math.max(-1, Math.min(1, s))) / 6.2832;
    // Stay on the branch of the sine the sweep is currently on.
    const rising = Math.cos(phase * 6.2832) >= 0;
    phase = ((rising ? a : 0.5 - a) + 1) % 1;
    held = performance.now() + 3500;
    kick();
  });

  new IntersectionObserver((es) => {
    visible = es[0].isIntersecting;
    if (visible) kick();
  }).observe(cv);
  document.addEventListener("visibilitychange", kick);
  new ResizeObserver(() => { resize(); kick(); }).observe(cv);
  window.addEventListener("themechange", () => { readTheme(); kick(); });
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { readTheme(); kick(); });
  document.fonts?.ready.then(kick);

  readTheme();
  resize();
  kick();
})();
