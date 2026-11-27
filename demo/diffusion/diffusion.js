// Exact-score EDM reverse-SDE conformer sampler on Adenylate Kinase (4AKE open vs. 1AKE closed).
(function () {
  // Kabsch-aligned 214 CA atoms of 4AKE (open) and 1AKE (closed), oriented with 87% of LID/NMP swing in XY.
  const X_OPEN = [[13.3,-9.4,-7.1],[12.9,-6.4,-4.9],[9.9,-4.2,-5.6],[8.5,-1.1,-3.9],[4.7,-0.8,-4.7],[3.3,2.8,-4.5],[-0.3,3.8,-4.8],[-2.9,6.2,-3.5],[-5.3,5.1,-0.8],[-8.3,3.6,-2.6],[-6.5,2.9,-5.9],[-6.3,-0.9,-5.3],[-2.6,-1.7,-5.8],[-3.1,-4.7,-3.6],[-5.0,-6.3,-6.5],[-1.7,-7.2,-8.2],[0.2,-8.4,-5.1],[-0.7,-12.1,-5.1],[0.5,-12.2,-8.7],[3.8,-10.7,-7.9],[4.4,-13.0,-4.9],[3.3,-16.0,-7.1],[5.5,-15.0,-10.0],[8.7,-14.3,-8.1],[8.3,-16.3,-4.9],[8.7,-13.3,-2.6],[6.7,-12.1,0.4],[4.7,-8.9,1.0],[6.2,-6.7,3.7],[4.3,-3.8,5.0],[5.5,-1.7,8.0],[2.3,-2.3,10.0],[2.8,-6.0,9.7],[6.5,-5.8,10.4],[6.0,-3.7,13.5],[3.4,-6.0,14.9],[5.8,-8.9,14.3],[8.9,-7.3,16.0],[7.1,-6.7,19.3],[5.9,-10.4,19.5],[9.3,-11.8,18.8],[10.9,-9.9,21.6],[14.1,-9.0,19.8],[16.1,-5.9,20.7],[14.3,-3.6,18.2],[10.9,-5.3,18.8],[11.1,-4.4,22.4],[12.2,-0.8,21.5],[9.1,-0.6,19.2],[6.7,-2.3,21.6],[7.3,0.5,24.2],[6.9,3.3,21.6],[3.6,1.9,20.4],[2.4,1.1,24.0],[3.3,4.7,25.0],[1.2,5.9,22.1],[4.4,7.4,20.6],[5.2,7.3,16.9],[7.9,5.1,15.5],[10.8,7.0,14.2],[12.2,6.5,10.6],[15.3,5.3,12.3],[13.6,2.4,14.1],[11.5,1.7,11.0],[14.6,1.3,8.7],[16.1,-0.9,11.4],[13.1,-3.3,11.3],[13.2,-3.4,7.5],[16.9,-4.3,7.6],[16.2,-7.0,10.3],[13.1,-8.3,8.4],[14.5,-8.6,5.0],[17.4,-10.5,6.5],[14.9,-13.3,7.2],[15.4,-16.3,4.9],[11.9,-15.9,3.4],[13.2,-12.8,1.7],[16.0,-14.8,0.1],[14.3,-14.2,-3.4],[13.6,-10.5,-2.8],[10.5,-8.9,-1.4],[7.4,-6.7,-2.1],[7.5,-3.7,0.1],[4.0,-2.3,0.0],[3.6,1.5,0.6],[7.1,2.2,2.1],[9.0,4.4,1.8],[6.6,7.3,2.3],[9.1,10.2,2.1],[12.5,11.1,0.5],[14.1,11.0,4.0],[12.9,7.4,4.3],[13.9,6.5,0.8],[17.3,8.0,1.2],[17.9,6.1,4.4],[17.1,2.8,2.8],[19.8,3.6,0.2],[22.3,4.5,2.9],[21.7,1.2,4.7],[22.5,-0.7,1.5],[18.9,-1.7,0.8],[18.6,-1.4,-3.0],[15.4,-2.1,-5.0],[15.1,-3.3,-8.5],[11.6,-1.8,-9.1],[9.5,1.0,-7.9],[6.0,0.7,-9.4],[3.4,3.4,-9.1],[-0.1,2.1,-9.4],[-1.9,5.0,-11.0],[-5.7,5.6,-10.6],[-7.1,9.2,-11.0],[-8.3,10.9,-7.9],[-11.8,11.3,-9.3],[-12.1,7.7,-10.1],[-10.9,6.8,-6.6],[-13.4,9.3,-5.1],[-16.2,7.7,-7.1],[-15.1,4.2,-6.1],[-15.5,5.2,-2.5],[-17.9,7.9,-1.6],[-20.7,6.0,-3.4],[-20.7,3.0,-1.1],[-23.3,1.9,1.4],[-22.6,-0.7,4.1],[-24.7,-2.9,6.3],[-23.1,-2.0,9.7],[-23.9,-5.3,11.7],[-22.5,-7.7,9.1],[-19.9,-5.7,7.2],[-21.6,-6.2,3.8],[-20.7,-3.5,1.3],[-22.5,-2.1,-1.7],[-21.9,0.6,-4.3],[-24.8,2.6,-5.6],[-23.4,2.0,-9.1],[-21.1,-1.0,-8.9],[-23.0,-3.4,-6.7],[-26.3,-1.8,-5.6],[-28.5,-3.7,-3.2],[-31.6,-5.5,-4.7],[-33.8,-3.2,-2.6],[-31.7,0.0,-3.0],[-30.8,1.5,0.4],[-31.7,-1.9,2.0],[-29.4,-4.8,3.0],[-30.3,-7.8,0.9],[-29.9,-10.4,3.7],[-31.7,-8.7,6.6],[-33.9,-6.1,4.9],[-32.1,-3.7,7.2],[-31.0,-0.2,6.2],[-27.5,0.4,5.0],[-25.1,3.2,5.9],[-22.1,5.1,4.8],[-18.8,4.7,6.6],[-17.2,7.9,8.0],[-14.1,7.2,5.9],[-16.2,7.0,2.8],[-17.2,10.6,3.2],[-16.3,12.7,0.2],[-14.7,15.4,2.4],[-12.6,12.9,4.2],[-11.6,11.1,1.0],[-10.6,14.2,-0.9],[-8.5,15.3,2.0],[-6.7,11.9,2.3],[-5.7,12.2,-1.4],[-4.3,15.7,-1.3],[-2.4,14.8,1.8],[-1.1,11.7,-0.1],[0.0,13.8,-3.1],[1.8,16.6,-1.4],[3.8,13.9,0.5],[4.4,11.5,-2.4],[5.1,13.9,-5.3],[8.8,14.5,-4.4],[9.3,10.8,-3.7],[8.1,9.8,-7.1],[10.5,12.5,-8.4],[13.3,11.1,-6.3],[12.7,7.6,-7.7],[12.7,8.9,-11.3],[16.0,10.6,-10.8],[17.3,7.4,-9.3],[15.9,5.4,-12.3],[17.7,7.8,-14.7],[20.8,7.4,-12.7],[20.5,3.6,-13.3],[20.1,2.8,-9.6],[16.6,1.3,-9.9],[13.8,1.0,-12.5],[10.7,3.1,-12.0],[7.4,2.4,-13.9],[3.9,3.9,-13.8],[1.0,1.4,-14.1],[-2.5,2.4,-15.3],[-4.6,1.0,-12.5],[-7.7,1.9,-14.5],[-7.3,-1.2,-16.8],[-8.8,-4.6,-16.5],[-7.0,-6.7,-13.9],[-6.0,-9.0,-16.7],[-4.5,-6.2,-18.7],[-2.6,-4.7,-15.7],[-1.3,-8.2,-15.0],[-0.2,-8.4,-18.6],[1.5,-5.1,-18.3],[3.3,-6.1,-15.1],[4.9,-9.1,-16.8],[6.2,-6.9,-19.6],[7.8,-4.6,-17.0],[9.3,-7.2,-14.7],[10.1,-10.3,-16.7]];
  const X_CLOSED = [[10.4,-10.7,-7.9],[10.5,-7.7,-5.6],[8.0,-4.9,-6.3],[7.3,-1.4,-4.7],[3.8,0.1,-5.2],[3.9,3.9,-5.1],[0.8,6.2,-5.2],[-1.4,8.8,-3.4],[-3.9,7.6,-0.7],[-6.7,5.6,-2.4],[-4.9,5.2,-5.7],[-5.8,1.4,-5.8],[-2.4,0.0,-4.5],[-3.9,-2.4,-2.0],[-6.6,-3.6,-4.5],[-3.9,-4.4,-7.1],[-1.5,-5.9,-4.5],[-4.1,-8.8,-3.9],[-4.0,-9.6,-7.7],[-0.2,-9.6,-8.0],[0.0,-11.7,-4.9],[-2.5,-14.3,-6.1],[-1.1,-14.4,-9.7],[2.6,-14.3,-8.9],[2.7,-16.0,-5.6],[4.7,-13.2,-3.8],[3.7,-12.0,-0.4],[2.2,-8.6,0.1],[4.0,-6.7,3.0],[1.6,-3.8,3.7],[2.6,-1.5,6.4],[-0.9,-0.1,6.9],[-2.1,-3.7,7.5],[0.7,-4.5,9.8],[0.2,-1.3,11.9],[-3.7,-1.6,12.2],[-3.0,-5.2,13.3],[-0.4,-4.2,15.9],[-2.7,-1.5,17.2],[-5.5,-4.0,17.5],[-3.3,-6.7,19.2],[-1.5,-4.2,21.5],[1.8,-5.6,20.0],[4.8,-3.8,21.5],[6.2,-3.1,18.0],[3.5,-1.0,16.4],[0.7,-0.3,19.0],[2.2,3.2,19.3],[0.8,4.1,15.8],[-2.8,4.3,17.1],[-3.0,8.1,17.5],[-0.8,9.1,14.7],[-2.9,7.1,12.1],[-6.1,8.6,13.6],[-4.7,12.0,13.1],[-3.8,11.2,9.4],[-0.2,11.5,10.0],[2.9,9.4,9.0],[4.2,7.0,11.8],[8.0,7.6,12.6],[10.4,6.1,10.0],[12.3,4.4,12.8],[9.2,2.4,13.8],[8.0,1.4,10.4],[11.5,0.4,9.2],[11.8,-2.0,12.2],[8.5,-3.7,11.3],[9.5,-4.2,7.6],[12.9,-5.6,8.4],[11.3,-8.1,10.7],[8.8,-9.2,8.1],[11.3,-9.5,5.3],[13.2,-11.9,7.6],[10.3,-14.5,7.6],[10.9,-17.4,5.3],[8.2,-16.4,2.8],[10.2,-13.3,1.8],[13.4,-15.2,0.6],[11.8,-15.5,-2.9],[10.3,-12.0,-3.3],[7.7,-9.8,-1.9],[5.6,-6.9,-2.8],[5.9,-3.6,-0.7],[2.6,-1.9,-0.4],[2.8,1.4,1.5],[6.6,1.3,2.1],[8.8,3.2,1.5],[6.7,6.2,2.3],[9.5,8.7,2.8],[13.1,9.3,1.6],[14.6,8.0,5.0],[12.7,4.8,4.5],[14.0,4.5,0.9],[17.5,5.1,2.0],[17.2,2.5,4.8],[16.3,-0.1,2.1],[19.3,1.0,0.1],[21.7,0.7,3.1],[20.3,-2.7,3.9],[20.8,-4.4,0.4],[17.1,-4.5,-0.5],[17.0,-3.7,-4.1],[13.7,-4.2,-6.1],[12.8,-5.2,-9.6],[9.9,-2.9,-10.3],[8.5,0.3,-8.9],[4.9,0.9,-10.1],[3.4,4.3,-9.7],[-0.4,4.4,-9.6],[-1.4,7.9,-10.7],[-4.7,9.3,-9.7],[-5.9,13.0,-9.6],[-6.3,14.4,-6.0],[-9.9,15.3,-6.7],[-11.0,11.6,-7.0],[-9.6,10.5,-3.6],[-12.1,11.6,-1.1],[-15.0,9.6,-2.8],[-12.9,6.5,-3.3],[-12.2,6.5,0.3],[-15.7,7.1,1.6],[-17.0,4.0,-0.2],[-14.1,1.9,1.1],[-15.0,-1.2,3.2],[-12.4,-3.7,4.4],[-12.1,-7.2,5.8],[-9.6,-6.4,8.5],[-8.2,-9.8,9.3],[-7.1,-10.3,5.8],[-6.7,-6.7,4.6],[-8.9,-7.4,1.6],[-10.5,-4.2,0.4],[-13.8,-3.4,-1.4],[-15.7,-0.3,-2.6],[-19.5,-0.1,-2.3],[-19.7,1.3,-5.9],[-16.5,0.6,-7.7],[-15.7,-2.9,-6.6],[-18.4,-4.2,-4.4],[-18.4,-7.7,-2.5],[-20.5,-10.5,-4.2],[-22.8,-10.3,-1.2],[-23.2,-6.8,0.2],[-21.4,-6.2,3.4],[-19.6,-9.5,3.6],[-16.1,-10.6,2.6],[-16.1,-13.0,-0.4],[-13.9,-15.5,1.1],[-14.9,-15.8,4.7],[-18.4,-14.4,4.7],[-17.4,-12.1,7.5],[-18.9,-8.6,7.9],[-16.9,-5.8,6.3],[-16.1,-2.6,8.3],[-14.9,0.9,7.6],[-11.5,2.2,8.6],[-11.5,4.7,11.4],[-9.2,7.0,9.5],[-11.4,7.0,6.3],[-13.3,10.2,6.9],[-13.0,13.2,4.7],[-10.8,15.3,7.0],[-8.3,12.6,7.9],[-8.0,11.8,4.1],[-7.3,15.5,3.4],[-4.5,15.5,6.0],[-3.1,12.3,4.5],[-3.0,13.9,1.0],[-0.9,16.9,2.3],[1.5,14.5,4.0],[1.7,12.5,0.8],[2.6,15.5,-1.3],[5.1,16.7,1.2],[6.9,13.5,2.1],[6.0,10.7,-0.2],[5.6,12.3,-3.6],[9.4,12.9,-3.7],[9.9,9.1,-4.0],[8.7,9.4,-7.7],[11.9,11.4,-8.5],[14.0,9.2,-6.2],[13.1,5.9,-8.0],[12.9,7.5,-11.4],[16.6,8.3,-10.8],[17.1,4.7,-9.7],[15.6,3.5,-13.0],[17.8,5.8,-15.0],[20.8,4.3,-13.4],[19.8,0.7,-14.3],[19.0,-0.3,-10.6],[15.2,-1.0,-11.1],[12.4,-0.8,-13.7],[9.9,2.2,-13.1],[6.3,1.9,-14.7],[3.5,4.5,-14.3],[-0.1,3.2,-14.6],[-3.4,5.1,-14.9],[-5.1,4.5,-11.7],[-8.4,6.1,-13.0],[-9.2,3.1,-15.2],[-11.6,0.5,-14.0],[-10.1,-2.0,-11.5],[-10.0,-4.8,-14.2],[-8.4,-2.5,-16.7],[-5.7,-1.5,-14.2],[-5.3,-5.2,-13.6],[-4.7,-5.8,-17.3],[-2.0,-3.0,-17.7],[-0.2,-4.6,-14.8],[-0.1,-8.0,-16.5],[1.5,-6.3,-19.5],[4.2,-4.9,-17.4],[4.9,-8.1,-15.2],[4.2,-11.2,-17.4]];

  const N = X_OPEN.length;
  const SIGMA_MAX = 240.0;
  const SIGMA_MIN = 0.05;
  const RHO = 3.5;

  // Mobile wing mask: NMP domain (res 30..59 -> idx 29..58) and LID domain (res 122..159 -> idx 121..158).
  function isWing(idx) {
    return (idx >= 29 && idx <= 58) || (idx >= 121 && idx <= 158);
  }

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

  function l2Sq(a, b) {
    let s = 0;
    for (let i = 0; i < N; i++) {
      const dx = a[i][0] - b[i][0];
      const dy = a[i][1] - b[i][1];
      const dz = a[i][2] - b[i][2];
      s += dx * dx + dy * dy + dz * dz;
    }
    return s;
  }

  function exactDenoiser(x, sigma, bias) {
    const inv2s2 = 0.5 / (sigma * sigma);
    const dOpen = l2Sq(x, X_OPEN);
    const dClosed = l2Sq(x, X_CLOSED);
    const logit = bias - (dOpen - dClosed) * inv2s2;
    let wOpen;
    if (logit > 40) wOpen = 1;
    else if (logit < -40) wOpen = 0;
    else wOpen = 1 / (1 + Math.exp(-logit));
    const wClosed = 1 - wOpen;
    const x0 = new Array(N);
    for (let i = 0; i < N; i++) {
      x0[i] = [
        wOpen * X_OPEN[i][0] + wClosed * X_CLOSED[i][0],
        wOpen * X_OPEN[i][1] + wClosed * X_CLOSED[i][1],
        wOpen * X_OPEN[i][2] + wClosed * X_CLOSED[i][2],
      ];
    }
    return { x0, wOpen };
  }

  function simulateEDM(steps, churn, bias, seed) {
    const initRng = mulberry32(seed);
    const churnRng = mulberry32((seed * 2654435761) ^ 0x9e3779b9);

    const sigmas = new Array(steps + 1);
    const sMaxInv = Math.pow(SIGMA_MAX, 1 / RHO);
    const sMinInv = Math.pow(SIGMA_MIN, 1 / RHO);
    for (let i = 0; i < steps; i++) {
      const frac = i / (steps - 1);
      sigmas[i] = Math.pow(sMaxInv + frac * (sMinInv - sMaxInv), RHO);
    }
    sigmas[steps] = 0;

    let x = new Array(N);
    for (let i = 0; i < N; i++) {
      x[i] = [randn(initRng) * SIGMA_MAX, randn(initRng) * SIGMA_MAX, randn(initRng) * SIGMA_MAX];
    }

    const frames = [];
    for (let i = 0; i < steps; i++) {
      const s = sigmas[i];
      const sNext = sigmas[i + 1];

      // Record the monotonically cooling Markov state x at noise level s.
      const { x0: x0State, wOpen } = exactDenoiser(x, s, bias);
      const rmsdOpen = Math.sqrt(l2Sq(x0State, X_OPEN) / N);
      const rmsdClosed = Math.sqrt(l2Sq(x0State, X_CLOSED) / N);
      frames.push({
        step: i,
        sigma: s,
        wOpen,
        rmsdOpen,
        rmsdClosed,
        xt: x,
      });

      const gamma =
        churn > 0 && s >= 8.0 && s <= 220.0 ? Math.min(churn / steps, Math.SQRT2 - 1) : 0;
      const sHat = s * (1 + gamma);

      let xIn = x;
      let x0 = x0State;
      if (gamma > 0) {
        const std = Math.sqrt(Math.max(0, sHat * sHat - s * s));
        xIn = new Array(N);
        for (let j = 0; j < N; j++) {
          xIn[j] = [
            x[j][0] + randn(churnRng) * std,
            x[j][1] + randn(churnRng) * std,
            x[j][2] + randn(churnRng) * std,
          ];
        }
        x0 = exactDenoiser(xIn, sHat, bias).x0;
      }

      const d1 = new Array(N);
      const xEuler = new Array(N);
      const dt = sNext - sHat;
      for (let j = 0; j < N; j++) {
        const d0 = (xIn[j][0] - x0[j][0]) / sHat;
        const d1y = (xIn[j][1] - x0[j][1]) / sHat;
        const d1z = (xIn[j][2] - x0[j][2]) / sHat;
        d1[j] = [d0, d1y, d1z];
        xEuler[j] = [xIn[j][0] + dt * d0, xIn[j][1] + dt * d1y, xIn[j][2] + dt * d1z];
      }

      if (sNext > 0) {
        const { x0: x0Next } = exactDenoiser(xEuler, sNext, bias);
        const xNext = new Array(N);
        for (let j = 0; j < N; j++) {
          const d2x = (xEuler[j][0] - x0Next[j][0]) / sNext;
          const d2y = (xEuler[j][1] - x0Next[j][1]) / sNext;
          const d2z = (xEuler[j][2] - x0Next[j][2]) / sNext;
          xNext[j] = [
            xIn[j][0] + 0.5 * dt * (d1[j][0] + d2x),
            xIn[j][1] + 0.5 * dt * (d1[j][1] + d2y),
            xIn[j][2] + 0.5 * dt * (d1[j][2] + d2z),
          ];
        }
        x = xNext;
      } else {
        x = xEuler;
      }
    }

    const rmsdOpenF = Math.sqrt(l2Sq(x, X_OPEN) / N);
    const rmsdClosedF = Math.sqrt(l2Sq(x, X_CLOSED) / N);
    const wOpenF = rmsdOpenF <= rmsdClosedF ? 1 : 0;
    frames.push({
      step: steps,
      sigma: SIGMA_MIN,
      wOpen: wOpenF,
      rmsdOpen: rmsdOpenF,
      rmsdClosed: rmsdClosedF,
      xt: x,
    });

    return frames;
  }

  function setupCanvas(el, cssW, cssH) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    el.width = Math.round(cssW * dpr);
    el.height = Math.round(cssH * dpr);
    const ctx = el.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w: cssW, h: cssH };
  }

  const el3d = document.getElementById('view3d');
  const elRmsd = document.getElementById('plot-rmsd');
  if (!el3d || !elRmsd) return;

  const view3d = setupCanvas(el3d, 380, 290);
  const plotRmsd = setupCanvas(elRmsd, 320, 244);

  const status3d = document.getElementById('status-3d');
  const statusRmsd = document.getElementById('status-rmsd');

  const ctrlChurn = document.getElementById('ctrl-churn');
  const ctrlSteps = document.getElementById('ctrl-steps');
  const ctrlBias = document.getElementById('ctrl-bias');
  const valChurn = document.getElementById('val-churn');
  const valSteps = document.getElementById('val-steps');
  const valBias = document.getElementById('val-bias');

  const reducedMotion =
    window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  let inkRgb = '17, 17, 17';
  let bgRgb = '255, 255, 255';
  let mutedColor = '#595959';
  let ruleColor = '#b8b8b8';

  function parseRgb(cssColor) {
    if (!cssColor) return null;
    const m = cssColor.match(
      /rgba?\(\s*(\d+)[\s,]+(\d+)[\s,]+(\d+)(?:[\s,/]+([\d.]+))?\s*\)/
    );
    if (m) {
      if (m[4] !== undefined && parseFloat(m[4]) === 0) return null;
      return m[1] + ', ' + m[2] + ', ' + m[3];
    }
    if (document.createElement) {
      try {
        const c = document.createElement('canvas');
        c.width = 1;
        c.height = 1;
        const g = c.getContext('2d');
        g.fillStyle = cssColor;
        g.fillRect(0, 0, 1, 1);
        const d = g.getImageData(0, 0, 1, 1).data;
        if (d[3] === 0) return null;
        return d[0] + ', ' + d[1] + ', ' + d[2];
      } catch (_) {}
    }
    return null;
  }

  function syncThemeColors(redraw) {
    const bodyStyle = getComputedStyle(document.body);
    const htmlStyle = getComputedStyle(document.documentElement);
    const nextInk = parseRgb(bodyStyle.color) || '17, 17, 17';
    const nextBg =
      parseRgb(bodyStyle.backgroundColor) ||
      parseRgb(htmlStyle.backgroundColor) ||
      '255, 255, 255';
    const statusStyle = status3d ? getComputedStyle(status3d) : null;
    const nextMuted = (statusStyle && statusStyle.color) || '#595959';
    const nextRule = 'rgba(' + nextInk + ', 0.28)';
    const changed = nextInk !== inkRgb || nextBg !== bgRgb || nextMuted !== mutedColor;
    inkRgb = nextInk;
    bgRgb = nextBg;
    mutedColor = nextMuted;
    ruleColor = nextRule;
    if (redraw && changed && frames.length > 0) {
      renderAt(currentPos);
    }
  }

  let seed = 103;
  let rotY = 0.18;
  let rotX = -0.14;
  let dragging = false;
  let lastPtrX = 0;
  let lastPtrY = 0;

  let frames = [];
  let playStart = null;
  let currentPos = 0;
  let animMs = 4600;
  let rafId = 0;

  const coords = new Array(N);
  const projected = new Array(N);
  const drawList = new Array(N + N - 1);

  function project3D(x, y, z) {
    const cy = Math.cos(rotY);
    const sy = Math.sin(rotY);
    const x1 = x * cy + z * sy;
    const z1 = -x * sy + z * cy;
    const cx = Math.cos(rotX);
    const sx = Math.sin(rotX);
    const y2 = y * cx - z1 * sx;
    const z2 = y * sx + z1 * cx;
    const cam = 150;
    const persp = cam / Math.max(45, cam - z2);
    return {
      sx: view3d.w / 2 + x1 * persp * 4.4,
      sy: view3d.h / 2 - y2 * persp * 4.4,
      z: z2,
      persp,
    };
  }

  function draw3D(f0, f1, frac) {
    const { ctx, w, h } = view3d;
    ctx.clearRect(0, 0, w, h);

    const sigma = f0.sigma + (f1.sigma - f0.sigma) * frac;
    const normScale = 1 / Math.sqrt(1 + Math.pow(sigma / 24, 2));
    const settle = Math.max(0, Math.min(1, (14 - sigma) / 12));

    for (let i = 0; i < N; i++) {
      const x = (f0.xt[i][0] + (f1.xt[i][0] - f0.xt[i][0]) * frac) * normScale;
      const y = (f0.xt[i][1] + (f1.xt[i][1] - f0.xt[i][1]) * frac) * normScale;
      const z = (f0.xt[i][2] + (f1.xt[i][2] - f0.xt[i][2]) * frac) * normScale;
      coords[i] = [x, y, z];
      projected[i] = project3D(x, y, z);
    }

    let itemCount = 0;

    // Fade in backbone bonds only once adjacent atoms condense to realistic CA-CA spacing.
    if (sigma < 16) {
      const sigmaGate = Math.max(0, Math.min(1, (16 - sigma) / 10));
      for (let i = 0; i < N - 1; i++) {
        const dx = coords[i + 1][0] - coords[i][0];
        const dy = coords[i + 1][1] - coords[i][1];
        const dz = coords[i + 1][2] - coords[i][2];
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        const bondFrac = Math.max(0, Math.min(1, (9.5 - dist) / 4.5)) * sigmaGate;
        if (bondFrac > 0.01) {
          const a = projected[i];
          const b = projected[i + 1];
          drawList[itemCount++] = {
            kind: 0,
            z: 0.5 * (a.z + b.z) - 0.05,
            a,
            b,
            bondFrac,
            wing: isWing(i) || isWing(i + 1),
          };
        }
      }
    }

    for (let i = 0; i < N; i++) {
      const p = projected[i];
      drawList[itemCount++] = {
        kind: 1,
        z: p.z,
        p,
        wing: isWing(i),
      };
    }

    const active = drawList.slice(0, itemCount).sort((u, v) => u.z - v.z);
    ctx.lineCap = 'round';

    for (let k = 0; k < active.length; k++) {
      const item = active[k];
      const depth = Math.max(0.15, Math.min(1, (item.z + 28) / 56));

      if (item.kind === 0) {
        const width = (item.wing ? 1.4 : 0.85) + depth * (item.wing ? 1.05 : 0.55);
        if (depth > 0.45 && item.bondFrac > 0.5) {
          ctx.lineWidth = width + 1.5;
          ctx.strokeStyle = 'rgba(' + bgRgb + ', ' + (0.68 * item.bondFrac).toFixed(2) + ')';
          ctx.beginPath();
          ctx.moveTo(item.a.sx, item.a.sy);
          ctx.lineTo(item.b.sx, item.b.sy);
          ctx.stroke();
        }
        const alpha = (item.wing ? 0.34 + depth * 0.62 : 0.14 + depth * 0.34) * item.bondFrac;
        ctx.lineWidth = width;
        ctx.strokeStyle = 'rgba(' + inkRgb + ', ' + alpha.toFixed(3) + ')';
        ctx.beginPath();
        ctx.moveTo(item.a.sx, item.a.sy);
        ctx.lineTo(item.b.sx, item.b.sy);
        ctx.stroke();
      } else {
        const cloudAlpha = 0.22 + depth * 0.56;
        const settledAlpha = item.wing ? 0.36 + depth * 0.62 : 0.16 + depth * 0.36;
        const alpha = cloudAlpha * (1 - settle) + settledAlpha * settle;
        const r = (item.wing ? 1.25 : 1.0) + depth * (item.wing ? 0.95 : 0.7);
        ctx.fillStyle = 'rgba(' + inkRgb + ', ' + alpha.toFixed(3) + ')';
        ctx.beginPath();
        ctx.arc(item.p.sx, item.p.sy, r, 0, Math.PI * 2);
        ctx.fill();
      }
    }
  }

  function drawRmsdPlot(idx) {
    const { ctx, w, h } = plotRmsd;
    ctx.clearRect(0, 0, w, h);

    const padL = 34;
    const padR = 12;
    const padT = 14;
    const padB = 26;
    const pw = w - padL - padR;
    const ph = h - padT - padB;
    const total = Math.max(1, frames.length - 1);
    const maxR = 8.0;

    function xOfStep(i) {
      return padL + (i / total) * pw;
    }
    function yOfR(r) {
      return padT + (1 - Math.min(maxR, r) / maxR) * ph;
    }

    ctx.strokeStyle = ruleColor;
    ctx.lineWidth = 0.85;
    ctx.beginPath();
    ctx.moveTo(padL, padT);
    ctx.lineTo(padL, padT + ph);
    ctx.lineTo(padL + pw, padT + ph);
    ctx.stroke();

    ctx.setLineDash([2, 3]);
    ctx.strokeStyle = 'rgba(' + inkRgb + ', 0.12)';
    ctx.beginPath();
    ctx.moveTo(padL, yOfR(4));
    ctx.lineTo(padL + pw, yOfR(4));
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = mutedColor;
    ctx.font = '10px ui-monospace, Menlo, monospace';
    ctx.textAlign = 'right';
    ctx.fillText('8 Å', padL - 5, padT + 4);
    ctx.fillText('4 Å', padL - 5, yOfR(4) + 3);
    ctx.fillText('0 Å', padL - 5, padT + ph + 3);

    ctx.textAlign = 'center';
    ctx.fillText('step 0', padL + 18, h - 7);
    ctx.fillText('step ' + Math.round(total / 2), padL + pw / 2, h - 7);
    ctx.fillText('step ' + total, padL + pw - 18, h - 7);

    const closedStroke = 'rgba(' + inkRgb + ', 0.52)';

    ctx.setLineDash([4, 3]);
    ctx.strokeStyle = closedStroke;
    ctx.lineWidth = 1.4;
    ctx.beginPath();
    for (let i = 0; i <= idx; i++) {
      const px = xOfStep(i);
      const py = yOfR(frames[i].rmsdClosed);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    }
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.strokeStyle = 'rgb(' + inkRgb + ')';
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    for (let i = 0; i <= idx; i++) {
      const px = xOfStep(i);
      const py = yOfR(frames[i].rmsdOpen);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    }
    ctx.stroke();

    const cur = frames[idx];
    const curX = xOfStep(idx);
    ctx.fillStyle = closedStroke;
    ctx.beginPath();
    ctx.arc(curX, yOfR(cur.rmsdClosed), 2.7, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = 'rgb(' + inkRgb + ')';
    ctx.beginPath();
    ctx.arc(curX, yOfR(cur.rmsdOpen), 3.0, 0, Math.PI * 2);
    ctx.fill();
  }

  function renderAt(pos) {
    const maxIdx = frames.length - 1;
    if (maxIdx < 0) return;
    const clampedPos = Math.max(0, Math.min(maxIdx, pos));
    const s0 = Math.floor(clampedPos);
    const s1 = Math.min(maxIdx, s0 + 1);
    const frac = clampedPos - s0;
    const f0 = frames[s0];
    const f1 = frames[s1];
    const idx = Math.round(clampedPos);
    const f = frames[idx];

    const basinSuffix =
      f.wOpen > 0.85 ? ' · Open (4AKE)' : f.wOpen < 0.15 ? ' · Closed (1AKE)' : '';
    status3d.textContent =
      'step ' + f.step + ' / ' + maxIdx + ' · σ = ' + (f.sigma >= 10 ? f.sigma.toFixed(0) : f.sigma.toFixed(2)) + ' Å' + basinSuffix;
    statusRmsd.textContent =
      'Open ' + f.rmsdOpen.toFixed(2) + ' · Closed ' + f.rmsdClosed.toFixed(2);

    draw3D(f0, f1, frac);
    drawRmsdPlot(idx);
  }

  function tick(now) {
    rafId = 0;
    if (document.hidden) return;
    if (playStart === null) {
      playStart = now - (currentPos / Math.max(1, frames.length - 1)) * animMs;
    }
    const u = Math.max(0, Math.min(1, (now - playStart) / animMs));
    currentPos = u * (frames.length - 1);
    if (!dragging) {
      rotY += 0.0025;
    }
    renderAt(currentPos);
    if (u < 1) {
      rafId = requestAnimationFrame(tick);
    }
  }

  function runSimulation(newSeed) {
    try {
      if (rafId) {
        cancelAnimationFrame(rafId);
        rafId = 0;
      }
      if (newSeed) seed += 1;
      const steps = parseInt(ctrlSteps.value, 10);
      const churn = parseFloat(ctrlChurn.value);
      const bias = parseFloat(ctrlBias.value);
      valChurn.textContent = churn === 0 ? '0 (ODE)' : String(churn);
      valSteps.textContent = String(steps);
      valBias.textContent = (bias > 0 ? '+' : '') + bias.toFixed(1);

      frames = simulateEDM(steps, churn, bias, seed);
      playStart = null;
      animMs = Math.max(2600, Math.min(5400, steps * 105));

      if (reducedMotion) {
        currentPos = frames.length - 1;
        renderAt(currentPos);
        return;
      }

      currentPos = 0;
      renderAt(0);
      rafId = requestAnimationFrame(tick);
    } catch (err) {
      if (status3d) {
        status3d.textContent = 'Error: ' + (err && err.message ? err.message : String(err));
      }
    }
  }

  el3d.addEventListener('pointerdown', (e) => {
    dragging = true;
    lastPtrX = e.clientX;
    lastPtrY = e.clientY;
    try {
      el3d.setPointerCapture(e.pointerId);
    } catch (_) {}
  });
  el3d.addEventListener('pointermove', (e) => {
    if (!dragging) return;
    rotY += (e.clientX - lastPtrX) * 0.01;
    rotX = Math.max(-1.2, Math.min(1.2, rotX + (e.clientY - lastPtrY) * 0.01));
    lastPtrX = e.clientX;
    lastPtrY = e.clientY;
    renderAt(currentPos);
  });
  el3d.addEventListener('pointerup', () => {
    dragging = false;
  });
  el3d.addEventListener('pointercancel', () => {
    dragging = false;
  });

  document.getElementById('btn-resample').addEventListener('click', () => runSimulation(true));
  ctrlChurn.addEventListener('input', () => runSimulation(false));
  ctrlSteps.addEventListener('input', () => runSimulation(false));
  ctrlBias.addEventListener('input', () => runSimulation(false));

  document.addEventListener('visibilitychange', () => {
    if (!document.hidden && !reducedMotion && currentPos < frames.length - 1 && !rafId) {
      playStart = null;
      rafId = requestAnimationFrame(tick);
    }
  });

  syncThemeColors(false);
  window.addEventListener('themechange', () => syncThemeColors(true));
  if (window.matchMedia) {
    const darkMq = window.matchMedia('(prefers-color-scheme: dark)');
    if (darkMq.addEventListener) {
      darkMq.addEventListener('change', () => syncThemeColors(true));
    } else if (darkMq.addListener) {
      darkMq.addListener(() => syncThemeColors(true));
    }
  }
  if ('MutationObserver' in window) {
    const themeObs = new MutationObserver(() => syncThemeColors(true));
    themeObs.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ['style', 'class', 'data-darkreader-scheme'],
    });
    themeObs.observe(document.body, {
      attributes: true,
      attributeFilter: ['style', 'class'],
    });
  }

  runSimulation(false);
})();
