// Client-side quantized multi-scale 1D CNN for peptide solubility and 20 x L mutational scanning.
(function () {
  const MODEL_SPEC = window.SOLUBILITY_MODEL;
  if (!MODEL_SPEC) return;

  const PRESETS = {
    insulin_b: 'FVNQHLCGSHLVEALYLVCGERGFFYTPKT',
    melittin: 'GIGAVLKVLTTGLPALISWIKRKRQQ',
    glp1: 'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR',
    villin: 'LSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF',
    ubiquitin: 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
    walp23: 'GWWLALALALALALALALALWWA',
  };

  const AA = MODEL_SPEC.alphabet;
  const AA_IDX = {};
  for (let i = 0; i < AA.length; i++) AA_IDX[AA[i]] = i;

  // Rows ordered biochemically: charged -> polar -> special -> aliphatic -> aromatic.
  const ROW_AA = 'EDKRQNHSTPGCAVMILFYW';
  const ROW_TOK = new Int32Array(20);
  for (let r = 0; r < 20; r++) ROW_TOK[r] = AA_IDX[ROW_AA[r]];

  function decodeB64Int8(b64) {
    const bin = atob(b64);
    const bytes = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    return new Int8Array(bytes.buffer);
  }

  function decodeLayer(layer) {
    const q = decodeB64Int8(layer.q_b64);
    const outCh = layer.shape[0];
    const inCh = layer.shape[1];
    const k = layer.shape[2];
    const w = new Float32Array(outCh * inCh * k);
    const stride = inCh * k;
    for (let o = 0; o < outCh; o++) {
      const s = layer.scales[o];
      const base = o * stride;
      for (let j = 0; j < stride; j++) w[base + j] = q[base + j] * s;
    }
    return { outCh, inCh, k, w, bias: Float32Array.from(layer.bias) };
  }

  const conv3 = decodeLayer(MODEL_SPEC.layers.conv3);
  const conv5 = decodeLayer(MODEL_SPEC.layers.conv5);
  const conv7 = decodeLayer(MODEL_SPEC.layers.conv7);
  const convMid = decodeLayer(MODEL_SPEC.layers.conv_mid);
  const head = decodeLayer(MODEL_SPEC.layers.head);
  const C1 = conv3.outCh + conv5.outCh + conv7.outCh;
  const INV_T = 1 / MODEL_SPEC.temperature;

  // Reorder convMid weights [outCh, inCh, k] -> [outCh, k, inCh] for contiguous unit-stride inner loops.
  const convMidContig = new Float32Array(convMid.w.length);
  for (let o = 0; o < convMid.outCh; o++) {
    for (let c = 0; c < convMid.inCh; c++) {
      for (let d = 0; d < convMid.k; d++) {
        convMidContig[(o * convMid.k + d) * convMid.inCh + c] =
          convMid.w[(o * convMid.inCh + c) * convMid.k + d];
      }
    }
  }

  function clamp(v, lo, hi) {
    return v < lo ? lo : v > hi ? hi : v;
  }

  function conv1dOneHotRep(tokens, layer, outArr, outOffset) {
    const L = tokens.length;
    const { outCh, inCh, k, w, bias } = layer;
    const pad = (k - 1) >> 1;
    const strideO = inCh * k;
    for (let pos = 0; pos < L; pos++) {
      const baseOut = pos * C1 + outOffset;
      for (let o = 0; o < outCh; o++) {
        let acc = bias[o];
        const wO = o * strideO;
        for (let d = 0; d < k; d++) {
          const src = clamp(pos + d - pad, 0, L - 1);
          acc += w[wO + tokens[src] * k + d];
        }
        outArr[baseOut + o] = acc > 0 ? acc : 0;
      }
    }
  }

  function conv1dDenseReluRep(inp, L) {
    const outCh = convMid.outCh;
    const inCh = convMid.inCh;
    const k = convMid.k;
    const bias = convMid.bias;
    const pad = (k - 1) >> 1;
    const out = new Float32Array(L * outCh);
    for (let pos = 0; pos < L; pos++) {
      const baseOut = pos * outCh;
      for (let o = 0; o < outCh; o++) {
        let acc = bias[o];
        const wBaseO = o * k * inCh;
        for (let d = 0; d < k; d++) {
          const src = clamp(pos + d - pad, 0, L - 1);
          const baseIn = src * inCh;
          const wOD = wBaseO + d * inCh;
          for (let c = 0; c < inCh; c++) {
            acc += inp[baseIn + c] * convMidContig[wOD + c];
          }
        }
        out[baseOut + o] = acc > 0 ? acc : 0;
      }
    }
    return out;
  }

  function predictSequence(seq) {
    const L = seq.length;
    const tokens = new Int32Array(L);
    for (let i = 0; i < L; i++) {
      tokens[i] = AA_IDX[seq[i]];
    }
    const h1 = new Float32Array(L * C1);
    conv1dOneHotRep(tokens, conv3, h1, 0);
    conv1dOneHotRep(tokens, conv5, h1, conv3.outCh);
    conv1dOneHotRep(tokens, conv7, h1, conv3.outCh + conv5.outCh);
    const h2 = conv1dDenseReluRep(h1, L);

    const inCh = head.inCh;
    const scores = new Float32Array(L);
    let sum = 0;
    for (let pos = 0; pos < L; pos++) {
      const baseIn = pos * inCh;
      let acc = head.bias[0];
      for (let c = 0; c < inCh; c++) acc += h2[baseIn + c] * head.w[c];
      scores[pos] = acc;
      sum += acc;
    }
    const logit = sum / L;
    const prob = 1 / (1 + Math.exp(-(logit * INV_T)));
    return { seq, tokens, h1, scores, logit, prob };
  }

  // Exact local single-mutant scan: only recomputes receptive field [p-5, p+5].
  function computeMutantScan(pred, topLimit) {
    const L = pred.seq.length;
    const tokens = new Int32Array(pred.tokens);
    const h1Base = pred.h1;
    const h1Tmp = new Float32Array(h1Base.length);
    h1Tmp.set(h1Base);
    const baseSum = pred.logit * L;
    const probsMatrix = new Float32Array(L * 20);
    const layers = [
      [conv3, 0],
      [conv5, conv3.outCh],
      [conv7, conv3.outCh + conv5.outCh],
    ];
    const midOut = convMid.outCh;
    const midIn = convMid.inCh;
    const midK = convMid.k;
    const midPad = (midK - 1) >> 1;
    const midBias = convMid.bias;
    const headW = head.w;
    const headB = head.bias[0];

    function recomputeH1Range(lo, hi) {
      for (let li = 0; li < layers.length; li++) {
        const layer = layers[li][0];
        const off = layers[li][1];
        const { outCh, inCh, k, w, bias } = layer;
        const pad = (k - 1) >> 1;
        const strideO = inCh * k;
        for (let pos = lo; pos <= hi; pos++) {
          const baseOut = pos * C1 + off;
          for (let o = 0; o < outCh; o++) {
            let acc = bias[o];
            const wO = o * strideO;
            for (let d = 0; d < k; d++) {
              const src = clamp(pos + d - pad, 0, L - 1);
              acc += w[wO + tokens[src] * k + d];
            }
            h1Tmp[baseOut + o] = acc > 0 ? acc : 0;
          }
        }
      }
    }

    function scoreRange(lo, hi) {
      let total = 0;
      for (let pos = lo; pos <= hi; pos++) {
        let s = headB;
        for (let o = 0; o < midOut; o++) {
          let acc = midBias[o];
          const wBaseO = o * midK * midIn;
          for (let d = 0; d < midK; d++) {
            const src = clamp(pos + d - midPad, 0, L - 1);
            const baseIn = src * midIn;
            const wOD = wBaseO + d * midIn;
            for (let c = 0; c < midIn; c++) {
              acc += h1Tmp[baseIn + c] * convMidContig[wOD + c];
            }
          }
          if (acc > 0) s += acc * headW[o];
        }
        total += s;
      }
      return total;
    }

    const top = [];
    let maxAbsDelta = 0.08;
    for (let p = 0; p < L; p++) {
      const orig = tokens[p];
      const h1Lo = Math.max(0, p - 3);
      const h1Hi = Math.min(L - 1, p + 3);
      const sLo = Math.max(0, p - 5);
      const sHi = Math.min(L - 1, p + 5);
      let oldSub = 0;
      for (let pos = sLo; pos <= sHi; pos++) oldSub += pred.scores[pos];

      for (let r = 0; r < 20; r++) {
        const aa = ROW_TOK[r];
        if (aa === orig) {
          probsMatrix[p * 20 + r] = pred.prob;
          continue;
        }
        tokens[p] = aa;
        recomputeH1Range(h1Lo, h1Hi);
        const newSub = scoreRange(sLo, sHi);
        const newLogit = (baseSum - oldSub + newSub) / L;
        const newProb = 1 / (1 + Math.exp(-(newLogit * INV_T)));
        probsMatrix[p * 20 + r] = newProb;
        const delta = newProb - pred.prob;
        const absD = Math.abs(delta);
        if (absD > maxAbsDelta) maxAbsDelta = absD;
        if (delta > 0.0005) {
          top.push({ pos: p, from: AA[orig], to: AA[aa], prob: newProb, delta });
        }
      }
      tokens[p] = orig;
      for (let pos = h1Lo; pos <= h1Hi; pos++) {
        const start = pos * C1;
        for (let c = 0; c < C1; c++) h1Tmp[start + c] = h1Base[start + c];
      }
    }
    top.sort((a, b) => b.delta - a.delta);
    return { probsMatrix, maxAbsDelta, top: top.slice(0, topLimit) };
  }

  const seqInput = document.getElementById('seq-input');
  const inputNotice = document.getElementById('input-notice');
  const btnReset = document.getElementById('btn-reset');
  const presetBtns = document.querySelectorAll('.preset-btn');
  const plotEl = document.getElementById('plot-profile');
  const heatmapEl = document.getElementById('plot-heatmap');
  const heatmapWrap = document.getElementById('heatmap-wrap');
  const heatmapPrompt = document.getElementById('heatmap-prompt');
  const btnRunScan = document.getElementById('btn-run-scan');
  const statusProfile = document.getElementById('status-profile');
  const statusHover = document.getElementById('status-hover');
  const statusLen = document.getElementById('status-len');
  const statusHeatmap = document.getElementById('status-heatmap');
  const scoreVal = document.getElementById('score-val');
  const scoreBand = document.getElementById('score-band');
  const gaugeFill = document.getElementById('gauge-fill');
  const mutsList = document.getElementById('muts-list');
  const mutsHint = document.getElementById('muts-hint');
  const mutTrail = document.getElementById('mut-trail');
  const cbarMin = document.getElementById('cbar-min');
  const cbarMax = document.getElementById('cbar-max');

  if (!plotEl || !heatmapEl || !seqInput) return;

  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const PROF_W = 420;
  const PROF_H = 190;
  plotEl.width = Math.round(PROF_W * dpr);
  plotEl.height = Math.round(PROF_H * dpr);
  const ctx = plotEl.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  const HEAT_H = 268;
  let heatW = 780;
  const hCtx = heatmapEl.getContext('2d');

  function resizeHeatmapCanvas(cssW) {
    heatW = cssW;
    heatmapEl.width = Math.round(cssW * dpr);
    heatmapEl.height = Math.round(HEAT_H * dpr);
    heatmapEl.style.width = cssW > 780 ? cssW + 'px' : '100%';
    hCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }
  resizeHeatmapCanvas(780);

  let inkRgb = '17, 17, 17';
  let posRgb = '45, 106, 138';
  let negRgb = '181, 86, 63';
  let mutedColor = '#595959';
  let ruleColor = '#b8b8b8';

  function parseRgbArr(cssColor) {
    if (!cssColor) return null;
    const m = cssColor.match(
      /rgba?\(\s*(\d+)[\s,]+(\d+)[\s,]+(\d+)(?:[\s,/]+([\d.]+))?\s*\)/
    );
    if (m) {
      if (m[4] !== undefined && parseFloat(m[4]) === 0) return null;
      return [parseInt(m[1], 10), parseInt(m[2], 10), parseInt(m[3], 10)];
    }
    return null;
  }

  function updateGaugeColor(p) {
    gaugeFill.style.backgroundColor =
      p >= 0.6 ? 'rgb(' + posRgb + ')' : p < 0.4 ? 'rgb(' + negRgb + ')' : 'rgb(' + inkRgb + ')';
  }

  function syncThemeColors(redraw) {
    const bodyStyle = getComputedStyle(document.body);
    const htmlStyle = getComputedStyle(document.documentElement);
    const inkArr = parseRgbArr(bodyStyle.color) || [17, 17, 17];
    const bgArr =
      parseRgbArr(bodyStyle.backgroundColor) ||
      parseRgbArr(htmlStyle.backgroundColor) || [255, 255, 255];
    const bgLum = 0.299 * bgArr[0] + 0.587 * bgArr[1] + 0.114 * bgArr[2];
    const isDark = bgLum < 128;

    const nextInk = inkArr[0] + ', ' + inkArr[1] + ', ' + inkArr[2];
    const nextPos = isDark ? '112, 178, 214' : '45, 106, 138';
    const nextNeg = isDark ? '226, 134, 112' : '181, 86, 63';
    const stStyle = statusLen ? getComputedStyle(statusLen) : null;
    const nextMuted = (stStyle && stStyle.color) || '#595959';
    const changed =
      nextInk !== inkRgb || nextPos !== posRgb || nextNeg !== negRgb || nextMuted !== mutedColor;
    inkRgb = nextInk;
    posRgb = nextPos;
    negRgb = nextNeg;
    mutedColor = nextMuted;
    ruleColor = 'rgba(' + nextInk + ', 0.28)';
    if (document.documentElement && document.documentElement.style && document.documentElement.style.setProperty) {
      document.documentElement.style.setProperty('--sol-pos', 'rgb(' + posRgb + ')');
      document.documentElement.style.setProperty('--sol-neg', 'rgb(' + negRgb + ')');
    }
    if (redraw && changed && currentPred) {
      updateGaugeColor(currentPred.prob);
      drawProfile();
      drawHeatmap();
    }
  }

  let activePresetKey = 'insulin_b';
  let currentPred = null;
  let currentScan = null;
  let mutHistory = [];
  let hoverPos = -1;
  let hoverRow = -1;
  let barGeom = { padL: 34, pw: 374 };
  let heatGeom = { padL: 24, padT: 20, cellW: 10, cellH: 11, cols: 0 };

  function cleanSequence(raw) {
    const lines = raw.split(/\r?\n/);
    let joined = '';
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i].trim();
      if (!line || line.startsWith('>')) continue;
      joined += line;
    }
    const upper = joined.toUpperCase();
    let clean = '';
    let ignored = 0;
    for (let i = 0; i < upper.length; i++) {
      const ch = upper[i];
      if (AA_IDX[ch] !== undefined) clean += ch;
      else if (ch >= 'A' && ch <= 'Z') ignored++;
    }
    return { clean, ignored };
  }

  function drawProfile() {
    ctx.clearRect(0, 0, PROF_W, PROF_H);
    if (!currentPred) return;
    const scores = currentPred.scores;
    const seq = currentPred.seq;
    const L = seq.length;

    const padL = 34;
    const padR = 10;
    const padT = 14;
    const padB = L <= 46 ? 24 : 20;
    const pw = PROF_W - padL - padR;
    const ph = PROF_H - padT - padB;
    barGeom = { padL, pw };

    let maxAbs = 2.5;
    for (let i = 0; i < L; i++) {
      const a = Math.abs(scores[i]);
      if (a > maxAbs) maxAbs = a;
    }
    maxAbs = Math.ceil(maxAbs * 2) / 2;
    const midY = padT + ph / 2;

    ctx.strokeStyle = ruleColor;
    ctx.lineWidth = 0.85;
    ctx.beginPath();
    ctx.moveTo(padL, padT);
    ctx.lineTo(padL, padT + ph);
    ctx.lineTo(padL + pw, padT + ph);
    ctx.stroke();

    ctx.setLineDash([2, 3]);
    ctx.strokeStyle = 'rgba(' + inkRgb + ', 0.25)';
    ctx.beginPath();
    ctx.moveTo(padL, midY);
    ctx.lineTo(padL + pw, midY);
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = mutedColor;
    ctx.font = '9.5px ui-monospace, Menlo, monospace';
    ctx.textAlign = 'right';
    ctx.fillText('+' + maxAbs.toFixed(1), padL - 4, padT + 4);
    ctx.fillText('0', padL - 4, midY + 3);
    ctx.fillText('-' + maxAbs.toFixed(1), padL - 4, padT + ph + 3);

    const step = pw / L;
    const barW = Math.max(1.2, step * (L <= 46 ? 0.68 : 0.78));

    if (hoverPos >= 0 && hoverPos < L) {
      const cx = padL + (hoverPos + 0.5) * step;
      ctx.fillStyle = 'rgba(' + inkRgb + ', 0.08)';
      ctx.fillRect(cx - step / 2, padT, step, ph);
    }

    for (let i = 0; i < L; i++) {
      const val = scores[i];
      const cx = padL + (i + 0.5) * step;
      const x = cx - barW / 2;
      const h = (Math.abs(val) / maxAbs) * (ph / 2);
      const isFocus = i === hoverPos;
      if (val >= 0) {
        ctx.fillStyle = isFocus ? 'rgb(' + posRgb + ')' : 'rgba(' + posRgb + ', 0.78)';
        ctx.fillRect(x, midY - h, barW, h);
      } else {
        ctx.fillStyle = isFocus ? 'rgb(' + negRgb + ')' : 'rgba(' + negRgb + ', 0.78)';
        ctx.fillRect(x, midY, barW, h);
      }
    }

    ctx.textAlign = 'center';
    if (L <= 46) {
      ctx.font = (L <= 32 ? '9.5px' : '8.5px') + ' ui-monospace, Menlo, monospace';
      for (let i = 0; i < L; i++) {
        const cx = padL + (i + 0.5) * step;
        ctx.fillStyle = i === hoverPos ? 'rgb(' + inkRgb + ')' : mutedColor;
        ctx.fillText(seq[i], cx, padT + ph + 13);
      }
    } else {
      ctx.font = '9.5px ui-monospace, Menlo, monospace';
      ctx.fillStyle = mutedColor;
      for (let t = 0; t <= 5; t++) {
        const pos = Math.round((t / 5) * (L - 1));
        const cx = padL + (pos + 0.5) * step;
        ctx.fillText(String(pos + 1), cx, padT + ph + 13);
      }
    }
  }

  function drawHeatmap() {
    hCtx.clearRect(0, 0, heatW, HEAT_H);
    if (!currentPred || !currentScan) return;

    const seq = currentPred.seq;
    const L = seq.length;
    const padL = 24;
    const padR = 10;
    const padT = 22;
    const padB = 18;
    const desiredW = Math.max(780, padL + padR + L * 11);
    if (desiredW !== heatW) resizeHeatmapCanvas(desiredW);

    const gw = heatW - padL - padR;
    const gh = HEAT_H - padT - padB;
    const cellW = gw / L;
    const cellH = gh / 20;
    heatGeom = { padL, padT, cellW, cellH, cols: L };

    const probs = currentScan.probsMatrix;
    const baseP = currentPred.prob;
    const maxD = currentScan.maxAbsDelta;
    if (cbarMin && cbarMax) {
      cbarMin.textContent = '−' + maxD.toFixed(2);
      cbarMax.textContent = '+' + maxD.toFixed(2);
    }

    hCtx.strokeStyle = ruleColor;
    hCtx.lineWidth = 0.8;
    hCtx.strokeRect(padL, padT, gw, gh);

    for (let p = 0; p < L; p++) {
      const origAA = seq[p];
      const x = padL + p * cellW;
      for (let r = 0; r < 20; r++) {
        const y = padT + r * cellH;
        const aa = ROW_AA[r];
        if (aa === origAA) {
          hCtx.fillStyle = 'rgba(' + inkRgb + ', 0.04)';
          hCtx.fillRect(x, y, cellW, cellH);
          const dotR = Math.max(1.6, Math.min(2.6, Math.min(cellW, cellH) * 0.24));
          hCtx.fillStyle = 'rgb(' + inkRgb + ')';
          hCtx.beginPath();
          hCtx.arc(x + cellW / 2, y + cellH / 2, dotR, 0, Math.PI * 2);
          hCtx.fill();
          continue;
        }
        const delta = probs[p * 20 + r] - baseP;
        if (delta > 0.002) {
          const norm = Math.min(1, delta / maxD);
          hCtx.fillStyle = 'rgba(' + posRgb + ', ' + (0.08 + 0.86 * norm).toFixed(3) + ')';
          hCtx.fillRect(x + 0.3, y + 0.3, Math.max(0.5, cellW - 0.6), Math.max(0.5, cellH - 0.6));
        } else if (delta < -0.002) {
          const norm = Math.min(1, -delta / maxD);
          hCtx.fillStyle = 'rgba(' + negRgb + ', ' + (0.08 + 0.86 * norm).toFixed(3) + ')';
          hCtx.fillRect(x + 0.3, y + 0.3, Math.max(0.5, cellW - 0.6), Math.max(0.5, cellH - 0.6));
        }
      }
    }

    hCtx.strokeStyle = 'rgba(' + inkRgb + ', 0.16)';
    hCtx.lineWidth = 0.7;
    const breaks = [4, 9, 12];
    for (let b = 0; b < breaks.length; b++) {
      const y = padT + breaks[b] * cellH;
      hCtx.beginPath();
      hCtx.moveTo(padL, y);
      hCtx.lineTo(padL + gw, y);
      hCtx.stroke();
    }

    hCtx.font = '9.5px ui-monospace, Menlo, monospace';
    hCtx.textAlign = 'right';
    for (let r = 0; r < 20; r++) {
      const y = padT + (r + 0.5) * cellH + 3;
      hCtx.fillStyle = r === hoverRow ? 'rgb(' + inkRgb + ')' : mutedColor;
      hCtx.fillText(ROW_AA[r], padL - 5, y);
    }

    hCtx.textAlign = 'center';
    hCtx.font = (cellW >= 14 ? '9.5px' : '8.5px') + ' ui-monospace, Menlo, monospace';
    for (let p = 0; p < L; p++) {
      const cx = padL + (p + 0.5) * cellW;
      if (cellW >= 9.5 || p === hoverPos) {
        hCtx.fillStyle = p === hoverPos ? 'rgb(' + inkRgb + ')' : mutedColor;
        hCtx.fillText(seq[p], cx, padT - 6);
      }
      if (p === 0 || (p + 1) % 10 === 0 || p === L - 1) {
        hCtx.fillStyle = mutedColor;
        hCtx.fillText(String(p + 1), cx, padT + gh + 12);
      }
    }

    if (hoverPos >= 0 && hoverPos < L) {
      const hx = padL + hoverPos * cellW;
      hCtx.strokeStyle = 'rgba(' + inkRgb + ', 0.25)';
      hCtx.lineWidth = 0.8;
      hCtx.strokeRect(hx, padT, cellW, gh);
      if (hoverRow >= 0 && hoverRow < 20) {
        const hy = padT + hoverRow * cellH;
        hCtx.strokeStyle = 'rgb(' + inkRgb + ')';
        hCtx.lineWidth = 1.4;
        hCtx.strokeRect(hx + 0.5, hy + 0.5, Math.max(1, cellW - 1), Math.max(1, cellH - 1));
      }
    }
  }

  function applyMutation(pos, newAA) {
    if (!currentPred || pos < 0 || pos >= currentPred.seq.length) return;
    const oldAA = currentPred.seq[pos];
    if (oldAA === newAA) return;
    const nextSeq = currentPred.seq.slice(0, pos) + newAA + currentPred.seq.slice(pos + 1);
    mutHistory.push(oldAA + (pos + 1) + newAA);
    if (mutHistory.length > 5) mutHistory = mutHistory.slice(-5);
    seqInput.value = nextSeq;
    updateAll(true);
  }

  function renderMutants(mutants) {
    while (mutsList.firstChild) mutsList.removeChild(mutsList.firstChild);
    if (!mutants.length) {
      const p = document.createElement('span');
      p.className = 'muts-hint';
      p.textContent = 'No single-point mutation raises predicted solubility further.';
      mutsList.appendChild(p);
      return;
    }
    for (let i = 0; i < mutants.length; i++) {
      const m = mutants[i];
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'mut-row';

      const left = document.createElement('span');
      left.className = 'mut-name';
      left.textContent = m.from + (m.pos + 1) + m.to;

      const right = document.createElement('span');
      right.className = 'mut-delta';
      right.textContent =
        currentPred.prob.toFixed(2) + ' → ' + m.prob.toFixed(2) + ' (+' + m.delta.toFixed(2) + ')';

      btn.appendChild(left);
      btn.appendChild(right);
      btn.addEventListener('click', () => applyMutation(m.pos, m.to));
      mutsList.appendChild(btn);
    }
  }

  function updateHoverText() {
    if (!currentPred || hoverPos < 0 || hoverPos >= currentPred.seq.length) {
      statusProfile.textContent = 'Hover to inspect';
      statusHover.textContent = 'Additive residue scores';
      statusHeatmap.textContent = 'Click any cell to mutate';
      return;
    }
    const aa = currentPred.seq[hoverPos];
    const sc = currentPred.scores[hoverPos];
    const sign = sc >= 0 ? '+' : '';
    statusProfile.textContent = aa + (hoverPos + 1) + ' · ' + sign + sc.toFixed(2);
    statusHover.textContent = 'Residue ' + aa + (hoverPos + 1) + ' (' + sign + sc.toFixed(2) + ')';

    if (currentScan && hoverRow >= 0 && hoverRow < 20) {
      const toAA = ROW_AA[hoverRow];
      const mutP = currentScan.probsMatrix[hoverPos * 20 + hoverRow];
      if (toAA === aa) {
        statusHeatmap.textContent = aa + (hoverPos + 1) + ' (current) · P = ' + mutP.toFixed(2);
      } else {
        const d = mutP - currentPred.prob;
        const ds = d >= 0 ? '+' : '';
        statusHeatmap.textContent =
          aa + (hoverPos + 1) + toAA + ' · ' + currentPred.prob.toFixed(2) + ' → ' + mutP.toFixed(2) + ' (' + ds + d.toFixed(2) + ')';
      }
    } else {
      statusHeatmap.textContent = 'Column ' + aa + (hoverPos + 1) + ' · click any row to mutate';
    }
  }

  function runAndShowScan() {
    if (!currentPred) return;
    currentScan = computeMutantScan(currentPred, 5);
    heatmapWrap.hidden = false;
    heatmapPrompt.hidden = true;
    mutsHint.textContent = 'Click to apply';
    renderMutants(currentScan.top);
    drawHeatmap();
    updateHoverText();
  }

  function updateAll(keepHistory) {
    try {
      if (!keepHistory) mutHistory = [];
      const { clean, ignored } = cleanSequence(seqInput.value);
      btnReset.hidden = !(activePresetKey && PRESETS[activePresetKey] !== clean);
      mutTrail.textContent = mutHistory.length ? 'WT → ' + mutHistory.join(' → ') : '';

      if (!clean.length) {
        currentPred = null;
        currentScan = null;
        hoverPos = -1;
        hoverRow = -1;
        scoreVal.textContent = '--';
        scoreBand.textContent = 'Enter A–Y sequence';
        gaugeFill.style.width = '0%';
        statusLen.textContent = '0 aa';
        inputNotice.textContent = '';
        while (mutsList.firstChild) mutsList.removeChild(mutsList.firstChild);
        ctx.clearRect(0, 0, PROF_W, PROF_H);
        hCtx.clearRect(0, 0, heatW, HEAT_H);
        return;
      }

      let seq = clean;
      const notes = [];
      if (ignored > 0) notes.push('Ignored ' + ignored + ' non-canonical character(s).');
      if (seq.length > 1500) {
        seq = seq.slice(0, 1500);
        notes.push('Truncated to first 1,500 residues.');
      }
      inputNotice.textContent = notes.join(' ');

      currentPred = predictSequence(seq);
      if (hoverPos >= seq.length) hoverPos = -1;

      const p = currentPred.prob;
      scoreVal.textContent = p.toFixed(2);
      scoreBand.textContent =
        p >= 0.65 ? 'High solubility' : p >= 0.35 ? 'Moderate' : 'Low solubility';
      gaugeFill.style.width = Math.round(p * 100) + '%';
      updateGaugeColor(p);
      statusLen.textContent = seq.length + ' aa';

      drawProfile();

      if (seq.length <= 180) {
        runAndShowScan();
      } else {
        currentScan = null;
        heatmapWrap.hidden = true;
        heatmapPrompt.hidden = false;
        btnRunScan.textContent = 'Compute 20 × ' + seq.length + ' mutational scan';
        mutsHint.textContent = 'Long sequence';
        while (mutsList.firstChild) mutsList.removeChild(mutsList.firstChild);
        const hint = document.createElement('span');
        hint.className = 'muts-hint';
        hint.textContent = 'Click below to run the 20 × ' + seq.length + ' scan.';
        mutsList.appendChild(hint);
        updateHoverText();
      }
    } catch (err) {
      scoreVal.textContent = 'Err';
      inputNotice.textContent = err && err.message ? err.message : String(err);
    }
  }

  function posFromProfileEvent(e) {
    if (!currentPred) return -1;
    const rect = plotEl.getBoundingClientRect();
    if (!rect.width) return -1;
    const cssX = ((e.clientX - rect.left) / rect.width) * PROF_W;
    const rel = (cssX - barGeom.padL) / barGeom.pw;
    if (rel < 0 || rel > 1) return -1;
    return Math.min(currentPred.seq.length - 1, Math.floor(rel * currentPred.seq.length));
  }

  function cellFromHeatmapEvent(e) {
    if (!currentPred || !currentScan) return { col: -1, row: -1 };
    const rect = heatmapEl.getBoundingClientRect();
    if (!rect.width || !rect.height) return { col: -1, row: -1 };
    const cssX = ((e.clientX - rect.left) / rect.width) * heatW;
    const cssY = ((e.clientY - rect.top) / rect.height) * HEAT_H;
    const col = Math.floor((cssX - heatGeom.padL) / heatGeom.cellW);
    const row = Math.floor((cssY - heatGeom.padT) / heatGeom.cellH);
    if (col < 0 || col >= heatGeom.cols || row < 0 || row >= 20) return { col: -1, row: -1 };
    return { col, row };
  }

  plotEl.addEventListener('pointermove', (e) => {
    const next = posFromProfileEvent(e);
    if (next !== hoverPos || hoverRow !== -1) {
      hoverPos = next;
      hoverRow = -1;
      drawProfile();
      drawHeatmap();
      updateHoverText();
    }
  });
  plotEl.addEventListener('pointerleave', () => {
    if (hoverPos !== -1) {
      hoverPos = -1;
      hoverRow = -1;
      drawProfile();
      drawHeatmap();
      updateHoverText();
    }
  });

  heatmapEl.addEventListener('pointermove', (e) => {
    const cell = cellFromHeatmapEvent(e);
    if (cell.col !== hoverPos || cell.row !== hoverRow) {
      hoverPos = cell.col;
      hoverRow = cell.row;
      drawProfile();
      drawHeatmap();
      updateHoverText();
    }
  });
  heatmapEl.addEventListener('pointerleave', () => {
    if (hoverPos !== -1 || hoverRow !== -1) {
      hoverPos = -1;
      hoverRow = -1;
      drawProfile();
      drawHeatmap();
      updateHoverText();
    }
  });
  heatmapEl.addEventListener('click', (e) => {
    const cell = cellFromHeatmapEvent(e);
    if (cell.col >= 0 && cell.row >= 0) {
      applyMutation(cell.col, ROW_AA[cell.row]);
    }
  });

  btnRunScan.addEventListener('click', () => runAndShowScan());

  seqInput.addEventListener('input', () => {
    presetBtns.forEach((b) => b.setAttribute('aria-pressed', 'false'));
    updateAll(false);
  });

  presetBtns.forEach((btn) => {
    btn.addEventListener('click', () => {
      const key = btn.getAttribute('data-preset');
      if (!key || !PRESETS[key]) return;
      activePresetKey = key;
      hoverPos = -1;
      hoverRow = -1;
      seqInput.value = PRESETS[key];
      presetBtns.forEach((b) => b.setAttribute('aria-pressed', String(b === btn)));
      updateAll(false);
    });
  });

  btnReset.addEventListener('click', () => {
    if (activePresetKey && PRESETS[activePresetKey]) {
      seqInput.value = PRESETS[activePresetKey];
      hoverPos = -1;
      hoverRow = -1;
      presetBtns.forEach((b) =>
        b.setAttribute('aria-pressed', String(b.getAttribute('data-preset') === activePresetKey))
      );
      updateAll(false);
    }
  });

  syncThemeColors(false);
  window.addEventListener('themechange', () => syncThemeColors(true));
  if (window.matchMedia) {
    const mq = window.matchMedia('(prefers-color-scheme: dark)');
    if (mq.addEventListener) mq.addEventListener('change', () => syncThemeColors(true));
  }
  if ('MutationObserver' in window) {
    const obs = new MutationObserver(() => syncThemeColors(true));
    obs.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ['style', 'class', 'data-darkreader-scheme'],
    });
    obs.observe(document.body, {
      attributes: true,
      attributeFilter: ['style', 'class'],
    });
  }

  updateAll(false);
})();
