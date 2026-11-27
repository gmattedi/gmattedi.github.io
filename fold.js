// Monochrome reverse-SDE denoising animation of Villin headpiece HP35 (PDB 1YRF, CA trace).
(function () {
  const canvas = document.getElementById('fold-canvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  const CSS_W = 220;
  const CSS_H = 78;
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  canvas.width = Math.round(CSS_W * dpr);
  canvas.height = Math.round(CSS_H * dpr);
  ctx.scale(dpr, dpr);

  // Centered, PCA-aligned right-handed CA coordinates of Villin headpiece (1YRF, 35 residues).
  const TARGET = [
    [-7.42, -5.58, 4.16], [-10.31, -5.31, 1.7], [-9.5, -3.6, -1.56], [-11.64, -0.7, -0.4],
    [-9.74, -0.37, 2.87], [-6.42, -0.61, 1.06], [-7.34, 2.24, -1.22], [-8.45, 4.32, 1.77],
    [-5.22, 3.66, 3.72], [-2.72, 4.04, 0.88], [-4.5, 6.46, -1.44], [-4.32, 4.27, -4.52],
    [-5.38, 0.82, -5.61], [-3.36, -2.35, -4.99], [-2.61, -2.46, -8.67], [-1.15, 1.05, -8.56],
    [0.84, 0.37, -5.43], [2.32, -2.73, -6.96], [4.08, -0.51, -9.55], [6.12, 1.27, -6.9],
    [9.68, 0.06, -6.06], [9.58, -2.93, -3.78], [11.48, -1.17, -0.99], [8.99, 1.69, -1.21],
    [5.99, -0.64, -0.94], [7.61, -2.09, 2.18], [8.28, 1.38, 3.6], [4.81, 2.69, 2.95],
    [3.21, -0.33, 4.54], [5.47, -0.2, 7.62], [4.82, 3.52, 8.03], [1.09, 2.68, 8.26],
    [1.74, -0.18, 10.71], [0.99, -2.92, 8.1], [2.93, -5.85, 6.61],
  ];

  const N = TARGET.length;
  const STEPS = 110;
  const HOLD_MS = 240;
  const DENOISE_MS = 2400;
  const CAM = 58;
  const FOCAL = 3.55;
  const Z_RANGE = 11.6;

  let seedCounter = 17;

  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a |= 0;
      a = (a + 0x6d2b79f5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function randn(rng) {
    const u = Math.max(1e-7, rng());
    const v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  function clampSoft(val, limit) {
    return limit * Math.tanh(val / limit);
  }

  // Precompute a discrete Euler-Maruyama reverse-SDE trajectory so timing is refresh-rate agnostic.
  function buildTrajectory() {
    const rng = mulberry32(seedCounter++);
    const frames = new Array(STEPS + 1);

    // Bounded initial Gaussian cloud strictly inside the 220x78 viewport.
    let state = new Array(N);
    for (let i = 0; i < N; i++) {
      state[i] = [
        clampSoft(randn(rng) * 14.5, 22.0),
        clampSoft(randn(rng) * 4.4, 7.2),
        clampSoft(randn(rng) * 7.5, 11.5),
      ];
    }
    frames[0] = state;

    for (let k = 1; k <= STEPS; k++) {
      const u = k / STEPS;
      const remaining = 1 - u;
      // Stronger score pull as noise anneals, paired with cooling Langevin churn.
      const drift = 0.045 + 0.14 * u * u;
      const churn = 1.45 * Math.pow(remaining, 1.35);
      const next = new Array(N);

      for (let i = 0; i < N; i++) {
        const cur = state[i];
        const tgt = TARGET[i];
        let nx = cur[0] + drift * (tgt[0] - cur[0]) + randn(rng) * churn;
        let ny = cur[1] + drift * (tgt[1] - cur[1]) + randn(rng) * churn * 0.55;
        let nz = cur[2] + drift * (tgt[2] - cur[2]) + randn(rng) * churn * 0.75;

        // Blend smoothly into exact crystal coordinates over the final 20% of steps.
        if (u > 0.8) {
          const snap = (u - 0.8) / 0.2;
          const w = snap * snap;
          nx = nx * (1 - w) + tgt[0] * w;
          ny = ny * (1 - w) + tgt[1] * w;
          nz = nz * (1 - w) + tgt[2] * w;
        }
        next[i] = [clampSoft(nx, 23.0), clampSoft(ny, 7.6), clampSoft(nz, 12.0)];
      }
      state = next;
      frames[k] = state;
    }

    frames[STEPS] = TARGET;
    return frames;
  }

  function project(x, y, z, rotY, rotX) {
    const cy = Math.cos(rotY);
    const sy = Math.sin(rotY);
    const x1 = x * cy + z * sy;
    const z1 = -x * sy + z * cy;

    const cx = Math.cos(rotX);
    const sx = Math.sin(rotX);
    const y2 = y * cx - z1 * sx;
    const z2 = y * sx + z1 * cx;

    const persp = CAM / Math.max(20, CAM - z2);
    return {
      sx: CSS_W / 2 + x1 * persp * FOCAL,
      sy: CSS_H / 2 + y2 * persp * FOCAL,
      z: z2,
      persp: persp,
    };
  }

  let trajectory = buildTrajectory();
  const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  let elapsedMs = reducedMotion ? HOLD_MS + DENOISE_MS : 0;
  let lastTick = performance.now();
  let rafId = 0;
  let inView = true;
  let tabVisible = !document.hidden;
  let inkRgb = '17, 17, 17';
  let bgRgb = '255, 255, 255';

  const coords = new Array(N);
  const projected = new Array(N);
  const drawList = new Array(N + N - 1);

  function parseRgb(cssColor) {
    const match =
      cssColor && cssColor.match(/rgba?\(\s*(\d+)[\s,]+(\d+)[\s,]+(\d+)(?:[\s,/]+([\d.]+))?\s*\)/);
    if (!match) return null;
    if (match[4] !== undefined && parseFloat(match[4]) === 0) return null;
    return match[1] + ', ' + match[2] + ', ' + match[3];
  }

  function syncColors(redraw) {
    const bodyStyle = getComputedStyle(document.body);
    const htmlStyle = getComputedStyle(document.documentElement);
    const nextInk = parseRgb(bodyStyle.color) || '17, 17, 17';
    const nextBg =
      parseRgb(bodyStyle.backgroundColor) ||
      parseRgb(htmlStyle.backgroundColor) ||
      '255, 255, 255';
    if (nextInk !== inkRgb || nextBg !== bgRgb) {
      inkRgb = nextInk;
      bgRgb = nextBg;
      if (redraw) drawFrame(elapsedMs);
    }
  }

  function drawFrame(tMs) {
    ctx.clearRect(0, 0, CSS_W, CSS_H);

    const denoiseU = Math.max(0, Math.min(1, (tMs - HOLD_MS) / DENOISE_MS));
    const exactStep = denoiseU * STEPS;
    const s0 = Math.floor(exactStep);
    const s1 = Math.min(STEPS, s0 + 1);
    const frac = exactStep - s0;
    const f0 = trajectory[s0];
    const f1 = trajectory[s1];

    const rotY = reducedMotion ? 0.3 : 0.3 + tMs * 0.00024;
    const rotX = reducedMotion ? 0.08 : Math.sin(tMs * 0.00018) * 0.12;

    for (let i = 0; i < N; i++) {
      const x = f0[i][0] + (f1[i][0] - f0[i][0]) * frac;
      const y = f0[i][1] + (f1[i][1] - f0[i][1]) * frac;
      const z = f0[i][2] + (f1[i][2] - f0[i][2]) * frac;
      coords[i] = [x, y, z];
      projected[i] = project(x, y, z, rotY, rotX);
    }

    let itemCount = 0;

    // Add peptide bonds, fading each in strictly by its live 3D CA-CA distance.
    for (let i = 0; i < N - 1; i++) {
      const dx = coords[i + 1][0] - coords[i][0];
      const dy = coords[i + 1][1] - coords[i][1];
      const dz = coords[i + 1][2] - coords[i][2];
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const bondFrac = Math.max(0, Math.min(1, (9.0 - dist) / 4.2));
      if (bondFrac > 0.01) {
        const a = projected[i];
        const b = projected[i + 1];
        drawList[itemCount++] = {
          kind: 0,
          z: 0.5 * (a.z + b.z) - 0.05,
          a: a,
          b: b,
          bondFrac: bondFrac,
        };
      }
    }

    // Add CA atoms into the same depth-sorted pass.
    for (let i = 0; i < N; i++) {
      const p = projected[i];
      drawList[itemCount++] = { kind: 1, z: p.z, p: p };
    }

    const active = drawList.slice(0, itemCount).sort((u, v) => u.z - v.z);

    for (let k = 0; k < active.length; k++) {
      const item = active[k];
      const depth = Math.max(0, Math.min(1, (item.z + Z_RANGE) / (2 * Z_RANGE)));

      if (item.kind === 0) {
        const width = 0.85 + depth * 0.9;
        const alpha = (0.14 + depth * 0.52) * item.bondFrac;

        // Subtle background-matched casing behind foreground segments keeps crossing helices legible.
        if (depth > 0.45 && item.bondFrac > 0.5) {
          ctx.lineWidth = width + 1.5;
          ctx.strokeStyle = 'rgba(' + bgRgb + ', ' + (0.65 * item.bondFrac).toFixed(2) + ')';
          ctx.beginPath();
          ctx.moveTo(item.a.sx, item.a.sy);
          ctx.lineTo(item.b.sx, item.b.sy);
          ctx.stroke();
        }

        ctx.lineWidth = width;
        ctx.strokeStyle = 'rgba(' + inkRgb + ', ' + alpha.toFixed(3) + ')';
        ctx.beginPath();
        ctx.moveTo(item.a.sx, item.a.sy);
        ctx.lineTo(item.b.sx, item.b.sy);
        ctx.stroke();
      } else {
        const radius = (1.3 + depth * 1.25) * (0.88 + 0.12 * item.p.persp);
        const alpha = 0.22 + depth * 0.72;
        ctx.fillStyle = 'rgba(' + inkRgb + ', ' + alpha.toFixed(3) + ')';
        ctx.beginPath();
        ctx.arc(item.p.sx, item.p.sy, radius, 0, Math.PI * 2);
        ctx.fill();
      }
    }
  }

  function tick(now) {
    rafId = 0;
    if (!inView || !tabVisible || reducedMotion) return;
    elapsedMs += Math.min(50, now - lastTick);
    lastTick = now;
    drawFrame(elapsedMs);
    rafId = requestAnimationFrame(tick);
  }

  function startLoop() {
    if (reducedMotion || !inView || !tabVisible || rafId) return;
    lastTick = performance.now();
    rafId = requestAnimationFrame(tick);
  }

  function stopLoop() {
    if (rafId) {
      cancelAnimationFrame(rafId);
      rafId = 0;
    }
  }

  const replayBtn = document.getElementById('btn-replay-fold');
  if (replayBtn) {
    replayBtn.addEventListener('click', () => {
      if (reducedMotion) return;
      trajectory = buildTrajectory();
      elapsedMs = 0;
      lastTick = performance.now();
      drawFrame(0);
      startLoop();
    });
  }

  if ('IntersectionObserver' in window) {
    const observer = new IntersectionObserver(
      (entries) => {
        inView = entries.some((e) => e.isIntersecting);
        if (inView) {
          startLoop();
        } else {
          stopLoop();
        }
      },
      { threshold: 0.05 }
    );
    observer.observe(canvas);
  }

  document.addEventListener('visibilitychange', () => {
    tabVisible = !document.hidden;
    if (tabVisible) {
      startLoop();
    } else {
      stopLoop();
    }
  });

  if ('MutationObserver' in window) {
    const colorObserver = new MutationObserver(() => syncColors(true));
    colorObserver.observe(document.documentElement, { attributes: true });
    colorObserver.observe(document.body, { attributes: true });
    if (document.head) {
      colorObserver.observe(document.head, { childList: true, subtree: true, characterData: true });
    }
  }

  window.addEventListener('load', () => syncColors(true));

  syncColors(false);
  drawFrame(elapsedMs);
  startLoop();
})();
