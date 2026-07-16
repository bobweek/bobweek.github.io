// =====================================================================
//  Stochastic PDE for a 1-D trait density  ν(z, t)
//  (stochastic generalization of multivariate DAGA; Eq. 63 of the
//   Supplementary Material — the SPDE that the measure-valued
//   martingale characterization reduces to in one trait dimension)
//
//     ∂ν/∂t = m(ν,z) ν + ½ M ∂²ν/∂z²  +  √(v ν) · ξ(z,t)
//     m(ν,z) = r − ½ ψ (θ − z)²  − c n            (b = 0, no directional selection)
//     n = ∫ ν(z) dz
//
//  ξ(z,t) is space–time white noise (the supplement's Ẇ_t(z)); the
//  √(v ν) prefactor is the demographic-stochasticity / drift term.
//
//  Solved by Euler–Maruyama on a fixed grid.  Space–time white noise
//  discretizes so the per-cell noise s.d. is √(v ν_i / Δz)·√(Δt), which
//  reproduces the abundance SDE dn = m̄ n dt + √(v n) dB (Eq. 18).
//
//  Fixed parameters (below) put the deterministic equilibrium at the
//  Gaussian mutation–selection–competition balance:
//     variance  P̂ = √(M/ψ) = 1 ,   mean  ẑ̄ = θ = 0 ,   abundance n̂ = 40.
//
//  The three check boxes overlay test functions f(z) = 1, z, z² onto the
//  density and open a trace of the associated inner product ⟨f, ν⟩:
//     ⟨1, ν⟩ = n ,   ⟨z, ν⟩/n = z̄ ,   ⟨(z−z̄)², ν⟩/n = P.
//  This is the Hilbert-space view — moments are projections of the
//  population measure onto test functions.
//
//  Single file, p5.js 1.x global mode.  Starts PAUSED.
// =====================================================================

// ---------- canvas ----------
const cw = 960, ch = 540;

// ---------- fixed model parameters ----------
const psi   = 0.5;     // ψ  strength of stabilizing selection
const Mmut  = 0.5;     // M  mutational variance (diffusion coefficient)
const rGrow = 1.25;    // r  intrinsic growth
const cComp = 0.025;   // c  competition strength
const theta = 0.0;     // θ  phenotypic optimum
const vRep  = 0.6;     // v  reproductive variance (drives drift)

// ---------- numerics ----------
const zmin = -4.5, zmax = 4.5;
const N    = 200;                       // grid cells
const dz   = (zmax - zmin) / N;
const invdz2 = 1 / (dz * dz);
const DT   = 0.0015;                    // Euler–Maruyama sub-step
const SUBSTEPS = 35;                    // sub-steps advanced per drawn frame
const noiseCoef = Math.sqrt(vRep * DT / dz);
const densMax = 26;                     // fixed y-axis for the density

// deterministic equilibrium (for reference lines)
const P_hat = Math.sqrt(Mmut / psi);    // = 1
const n_hat = (rGrow - 0.5 * Math.sqrt(psi * Mmut)) / cComp; // = 40

// ---------- state ----------
let z = new Float64Array(N);
let selTerm = new Float64Array(N);      // −½ ψ (θ − z)²   (fixed)
let nu = new Float64Array(N);
let nuNew = new Float64Array(N);
let simTime = 0;
let curN = 0, curZ = 0, curP = 0;

// ---------- moment histories (fixed-length rings; newest at end) ----------
const W = 280;
let histN, histZ, histP;

// ---------- UI ----------
let cbs = [];                           // three check boxes
let show = [false, false, false];
let playBtn, resetBtn;
let playing = false;
let host = null;

// ---------- palette (deck "Daytime" tokens) ----------
const INK  = '#1b1f24';
const GRY  = '#8a97a6';
const LGRY = '#d5dde6';
const CY   = '#007c91';   // n   / f(z)=1
const VIO  = '#5a4cc7';   // z̄  / f(z)=z
const PK   = '#c2185b';   // P   / f(z)=z²
const DENS_STROKE = '#3a4350';
const DENS_FILL   = 'rgba(58,67,80,0.10)';
const THETA_GUIDE = '#c3ccd6';
const PANEL_BG    = '#fafbfc';
const FONT = "'Cambria Math','STIX Two Math','Segoe UI',system-ui,sans-serif";

// ---------- main-plot geometry ----------
const px0 = 56, px1 = 600, py0 = 48, py1 = 468;

// ---------- right-column geometry ----------
const RX0 = 620, RX1 = 944, slotTop = 48, slotBot = 532, GAP = 12;
let SLOTS = [];

// =====================================================================
function setup() {
  host = (typeof document !== 'undefined') ? document.getElementById('sim') : null;
  const cnv = createCanvas(cw, ch);
  if (host) cnv.parent(host);

  for (let i = 0; i < N; i++) {
    z[i] = zmin + (i + 0.5) * dz;
    selTerm[i] = -0.5 * psi * (theta - z[i]) * (theta - z[i]);
  }
  layoutSlots();
  buildControls();
  seed();

  if (typeof document !== 'undefined') {
    document.addEventListener('visibilitychange', () => {
      if (document.hidden) noLoop();
      else { if (playing) loop(); else redraw(); }
    });
  }
  setPlaying(false);   // start paused (no background animation in a slide deck)
  redraw();            // render the initial frame once
}

function layoutSlots() {
  const slotH = (slotBot - slotTop - 2 * GAP) / 3;
  SLOTS = [];
  for (let k = 0; k < 3; k++) {
    const y0 = slotTop + k * (slotH + GAP);
    SLOTS.push({ y0, y1: y0 + slotH });
  }
}

// ---------- controls ----------
function buildControls() {
  const labels = ['f(z) = 1', 'f(z) = z', 'f(z) = z²'];
  const cols   = [CY, VIO, PK];
  for (let k = 0; k < 3; k++) {
    const cb = createCheckbox(labels[k], false);
    cb.addClass('chk');
    cb.style('color', cols[k]);
    cb.position(RX0 + 2, SLOTS[k].y0 + 2);
    if (host) cb.parent(host);
    // force a redraw when toggled while paused (physics is off, but overlays must update)
    cb.changed(() => { if (!playing) redraw(); });
    cbs.push(cb);
  }

  playBtn = createButton('\u25B6 play');
  playBtn.position(cw - 156, 12);
  playBtn.mousePressed(() => setPlaying(!playing));
  if (host) playBtn.parent(host);

  resetBtn = createButton('\u21BA reset');
  resetBtn.position(cw - 78, 12);
  resetBtn.mousePressed(() => { seed(); redraw(); });
  if (host) resetBtn.parent(host);
}

function setPlaying(p) {
  playing = p;
  if (playBtn) playBtn.html(p ? '\u23F8 pause' : '\u25B6 play');
  if (p) loop(); else noLoop();
}

// ---------- initial condition ----------
function seed() {
  const c0 = -0.9, sd0 = 0.4, mass0 = 10;     // narrow, off-centre, under-populated
  const A = mass0 / (Math.sqrt(2 * Math.PI) * sd0);
  for (let i = 0; i < N; i++) {
    const e = (z[i] - c0) / sd0;
    nu[i] = A * Math.exp(-0.5 * e * e);
  }
  histN = new Array(W).fill(NaN);
  histZ = new Array(W).fill(NaN);
  histP = new Array(W).fill(NaN);
  simTime = 0;
  const st = computeStats();
  curN = st.n; curZ = st.zbar; curP = st.P;
}

// ---------- one Euler–Maruyama sub-step ----------
function step() {
  let n = 0;
  for (let i = 0; i < N; i++) n += nu[i];
  n *= dz;

  for (let i = 0; i < N; i++) {
    const m = rGrow + selTerm[i] - cComp * n;
    let lap;
    if (i === 0)          lap = (nu[1] - nu[0]) * invdz2;          // Neumann (zero-flux)
    else if (i === N - 1) lap = (nu[N - 2] - nu[N - 1]) * invdz2;
    else                  lap = (nu[i + 1] - 2 * nu[i] + nu[i - 1]) * invdz2;
    const drift = m * nu[i] + 0.5 * Mmut * lap;
    const g = (nu[i] > 0) ? noiseCoef * Math.sqrt(nu[i]) * randomGaussian() : 0;
    const val = nu[i] + drift * DT + g;
    nuNew[i] = val > 0 ? val : 0;                                   // density stays ≥ 0
  }
  const tmp = nu; nu = nuNew; nuNew = tmp;
}

function computeStats() {
  let n = 0;
  for (let i = 0; i < N; i++) n += nu[i];
  n *= dz;
  if (n <= 1e-9) return { n: n, zbar: NaN, P: NaN };
  let zb = 0;
  for (let i = 0; i < N; i++) zb += z[i] * nu[i];
  zb = zb * dz / n;
  let P = 0;
  for (let i = 0; i < N; i++) { const e = z[i] - zb; P += e * e * nu[i]; }
  P = P * dz / n;
  return { n: n, zbar: zb, P: P };
}

function pushHist(h, v) { h.push(v); h.shift(); }

// =====================================================================
function draw() {
  if (playing) {
    for (let s = 0; s < SUBSTEPS; s++) step();
    const st = computeStats();
    pushHist(histN, st.n); pushHist(histZ, st.zbar); pushHist(histP, st.P);
    simTime += SUBSTEPS * DT;
    curN = st.n; curZ = st.zbar; curP = st.P;
  }
  for (let k = 0; k < 3; k++) show[k] = cbs[k].checked();
  render();
}

// =====================================================================
//  Rendering
// =====================================================================
function render() {
  background(255);
  drawHeader();
  drawMainPlot();
  drawSlots();
}

function drawHeader() {
  T('\u03BD(z) = abundance density',
     16, 30, 14.5, { c: INK });
  // T('\u2202\u03BD/\u2202t  =  m(\u03BD,z)\u2009\u03BD  +  \u00BD M \u2202\u00B2\u03BD/\u2202z\u00B2  +  \u221A(v\u2009\u03BD) \u00B7 \u03BE(z,t)',
  //   16, 30, 14.5, { c: INK });
}

function drawMainPlot() {
  const ctx = drawingContext;

  // panel
  ctx.fillStyle = PANEL_BG;
  ctx.fillRect(px0, py0, px1 - px0, py1 - py0);

  // horizontal gridlines + y ticks (density)
  for (const dv of [0, 10, 20]) {
    const y = yDens(dv);
    dline(px0, y, px1, y, LGRY, 1, [3, 4]);
    T(String(dv), px0 - 8, y + 4, 11, { c: GRY, align: 'right' });
  }
  // x ticks
  for (const zt of [-4, -2, 0, 2, 4]) {
    const x = xZ(zt);
    ctx.strokeStyle = LGRY; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x, py1); ctx.lineTo(x, py1 + 5); ctx.stroke();
    T(String(zt), x, py1 + 18, 11, { c: GRY, align: 'center' });
  }

  // optimum θ guide
  const xt = xZ(theta);
  dline(xt, py0, xt, py1, THETA_GUIDE, 1.2, [4, 4]);
  T('\u03B8', xt, py1 + 18, 12, { c: GRY, align: 'center' });

  // frame
  ctx.strokeStyle = LGRY; ctx.lineWidth = 1;
  ctx.strokeRect(px0, py0, px1 - px0, py1 - py0);

  // density: filled area
  ctx.beginPath();
  ctx.moveTo(px0, yDens(0));
  for (let i = 0; i < N; i++) ctx.lineTo(xZ(z[i]), clampY(yDens(nu[i])));
  ctx.lineTo(px1, yDens(0));
  ctx.closePath();
  ctx.fillStyle = DENS_FILL; ctx.fill();
  dline(px0, yDens(0), px1, yDens(0), DENS_STROKE, 1, []);   // solid ν = 0 boundary
  // density: top curve
  ctx.beginPath();
  ctx.moveTo(xZ(z[0]), clampY(yDens(nu[0])));
  for (let i = 1; i < N; i++) ctx.lineTo(xZ(z[i]), clampY(yDens(nu[i])));
  ctx.strokeStyle = DENS_STROKE; ctx.lineWidth = 1.7; ctx.stroke();

  // test-function overlays
  drawOverlays();

  // axis labels
  T('trait value   z', (px0 + px1) / 2, py1 + 36, 13, { c: INK, align: 'center' });
  // y-axis label (rotated)
  ctx.save();
  ctx.translate(16, (py0 + py1) / 2);
  ctx.rotate(-Math.PI / 2);
  T('abundance density   \u03BD(z)', 0, 0, 13, { c: INK, align: 'center' });
  ctx.restore();

  // fixed-parameter caption
  // T('m(\u03BD,z) = r \u2212 \u00BD\u2009\u03C8\u2009(\u03B8\u2212z)\u00B2 \u2212 c\u2009n      fixed:  r=1.25  \u03C8=0.5  M=0.5  c=0.025  \u03B8=0  v=0.6',
  //   px0, ch - 12, 11.5, { c: GRY });

  // onboarding hint
  if (!playing && simTime === 0) {
    T('press  \u25B6 play  \u00B7  then tick a box \u2192', (px0 + px1) / 2, py0 + 26, 13,
      { c: GRY, align: 'center' });
  }
}

// ---- dashed test functions, drawn on the TRUE density y-axis ----
//  x = z via xZ(),  y = f(z) via yDens()  — the same value→pixel map as the density
function drawOverlays() {
  const ctx = drawingContext;
  ctx.save();
  ctx.beginPath();
  ctx.rect(px0, py0, px1 - px0, py1 - py0);   // clip to the plot frame
  ctx.clip();

  // f(z) = 1 : horizontal line at value 1
  if (show[0]) {
    const y = yDens(1);
    dline(px0, y, px1, y, CY, 1.6, [7, 5]);
    T('f(z) = 1', px1 - 6, y - 6, 12, { c: CY, align: 'right', w: 'bold' });
  }

  // f(z) = z : straight line through the origin (0,0)
  if (show[1]) {
    dline(px0, yDens(0), px1, yDens(0), VIO, 1, [2, 4]);   // value = 0 line
    ctx.save(); ctx.setLineDash([7, 5]); ctx.strokeStyle = VIO; ctx.lineWidth = 1.6;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = xZ(z[i]), y = yDens(z[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke(); ctx.restore();
    T('f(z) = z', px1 - 6, yDens(zmax) - 6, 12, { c: VIO, align: 'right', w: 'bold' });
  }

  // f(z) = z² : parabola touching the origin at z = 0
  if (show[2]) {
    ctx.save(); ctx.setLineDash([7, 5]); ctx.strokeStyle = PK; ctx.lineWidth = 1.6;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = xZ(z[i]), y = yDens(z[i] * z[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke(); ctx.restore();
    T('f(z) = z\u00B2', px1 - 6, yDens(zmax * zmax) - 6, 12, { c: PK, align: 'right', w: 'bold' });
  }

  ctx.restore();
}
// ---- right-hand moment windows ----
function drawSlots() {
  drawMomentSlot(0, {
    color: CY, active: show[0],
    runs: [ {s: 'n', c: CY, w: 'bold'},
            {s: '  =  \u222B \u03BD(z) dz', c: INK} ],
    hist: histN, cur: curN, vlo: 0, vhi: 56, ref: n_hat,
    ticks: [0, 20, 40], fmt: (x) => x.toFixed(1)
  });
  drawMomentSlot(1, {
    color: VIO, active: show[1],
    runs: [ {s: 'z', c: VIO, w: 'bold', over: true},
            {s: '  =  \u222B z \u00B7 \u03BD(z)/n dz', c: INK} ],
    hist: histZ, cur: curZ, vlo: -1.5, vhi: 1.5, ref: theta,
    ticks: [-1, 0, 1], fmt: (x) => x.toFixed(2)
  });
  drawMomentSlot(2, {
    color: PK, active: show[2],
    runs: [ {s: 'P', c: PK, w: 'bold'},
            {s: '  =  \u222B (z \u2212 ', c: INK},
            {s: 'z', c: INK, over: true},
            {s: ' )\u00B2 \u00B7 \u03BD(z)/n dz', c: INK} ],
    hist: histP, cur: curP, vlo: 0, vhi: 2, ref: P_hat,
    ticks: [0, 1, 2], fmt: (x) => x.toFixed(2)
  });
}

function drawMomentSlot(k, o) {
  const s = SLOTS[k];
  const fy = s.y0 + 40;                          // formula baseline (below the check box)

  if (!o.active) {
    // greyed formula so the layout & meaning are visible before ticking
    const grey = o.runs.map(r => ({ s: r.s, over: r.over, c: GRY, w: 'normal' }));
    drawRuns(grey, RX0 + 6, fy, 14.5);
    return;
  }

  drawRuns(o.runs, RX0 + 6, fy, 14.5);

  // strip chart of the statistic
  const box = { x0: RX0 + 38, y0: s.y0 + 52, x1: RX1 - 8, y1: s.y1 - 12 };
  stripChart(box, o);
}

function stripChart(box, o) {
  const ctx = drawingContext;
  const bw = box.x1 - box.x0, bh = box.y1 - box.y0;
  const yVal = (v) => box.y1 - (v - o.vlo) / (o.vhi - o.vlo) * bh;
  const xIdx = (i) => box.x0 + (i / (W - 1)) * bw;

  // panel + y ticks
  ctx.fillStyle = PANEL_BG; ctx.fillRect(box.x0, box.y0, bw, bh);
  for (const tv of o.ticks) {
    const y = yVal(tv);
    dline(box.x0, y, box.x1, y, LGRY, 1, [3, 4]);
    T(o.fmt(tv), box.x0 - 6, y + 4, 10, { c: GRY, align: 'right' });
  }
  ctx.strokeStyle = LGRY; ctx.lineWidth = 1; ctx.strokeRect(box.x0, box.y0, bw, bh);

  // deterministic-equilibrium reference line
  const yr = yVal(o.ref);
  dline(box.x0, yr, box.x1, yr, o.color, 1.1, [5, 4]);

  // data trace
  ctx.save();
  ctx.beginPath(); ctx.rect(box.x0, box.y0, bw, bh); ctx.clip();
  ctx.strokeStyle = o.color; ctx.lineWidth = 1.6;
  ctx.beginPath();
  let started = false;
  for (let i = 0; i < W; i++) {
    const v = o.hist[i];
    if (v == null || Number.isNaN(v)) { started = false; continue; }
    const x = xIdx(i), y = yVal(v);
    if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
  }
  ctx.stroke();
  ctx.restore();

  // latest value: dot at right edge + label
  if (!Number.isNaN(o.cur)) {
    const yc = Math.max(box.y0 + 2, Math.min(box.y1 - 2, yVal(o.cur)));
    ctx.fillStyle = o.color;
    ctx.beginPath(); ctx.arc(box.x1, yc, 2.6, 0, 2 * Math.PI); ctx.fill();
    T(o.fmt(o.cur), box.x1 - 5, yc - 6, 12, { c: o.color, align: 'right', w: 'bold' });
  }

  // time axis hint
  T('time \u2192', box.x1, box.y1 + 11, 10, { c: GRY, align: 'right' });
}

// =====================================================================
//  Text / drawing helpers (canvas 2D for full glyph + overline control)
// =====================================================================
function T(str, x, y, px, opt) {
  opt = opt || {};
  const ctx = drawingContext;
  ctx.font = `${opt.w || 'normal'} ${px}px ${FONT}`;
  ctx.fillStyle = opt.c || INK;
  ctx.textAlign = opt.align || 'left';
  ctx.textBaseline = opt.baseline || 'alphabetic';
  ctx.fillText(str, x, y);
  ctx.textAlign = 'left';
}

// draw a sequence of coloured/weighted runs; runs with {over:true} get an overline
function drawRuns(runs, x, y, px) {
  const ctx = drawingContext;
  ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
  for (const r of runs) {
    ctx.font = `${r.w || 'normal'} ${px}px ${FONT}`;
    ctx.fillStyle = r.c || INK;
    ctx.fillText(r.s, x, y);
    const w = ctx.measureText(r.s).width;
    if (r.over) {
      ctx.strokeStyle = r.c || INK;
      ctx.lineWidth = Math.max(1, px * 0.06);
      const oy = y - px * 0.74;
      ctx.beginPath();
      ctx.moveTo(x + px * 0.05, oy);
      ctx.lineTo(x + w - px * 0.02, oy);
      ctx.stroke();
    }
    x += w;
  }
  return x;
}

function dline(x0, y0, x1, y1, colorS, wt, dash) {
  const ctx = drawingContext;
  ctx.save();
  ctx.setLineDash(dash || [5, 4]);
  ctx.strokeStyle = colorS; ctx.lineWidth = wt || 1;
  ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
  ctx.restore();
}

// ---------- coordinate transforms (main plot) ----------
function xZ(zz) { return px0 + (zz - zmin) / (zmax - zmin) * (px1 - px0); }
// function yDens(d) { return py1 - (d / densMax) * (py1 - py0); }
function yDens(d) { const Y_MIN = -zmax, Y_MAX = densMax; return py1 - (d - Y_MIN) / (Y_MAX - Y_MIN) * (py1 - py0); }
function clampY(y) { return y < py0 ? py0 : (y > py1 ? py1 : y); }

// =====================================================================
function keyPressed() {
  if (key === ' ')                { setPlaying(!playing); return false; }
  if (key === 'r' || key === 'R') { seed(); redraw(); return false; }
  if (key === '1') { toggleCb(0); return false; }
  if (key === '2') { toggleCb(1); return false; }
  if (key === '3') { toggleCb(2); return false; }
}
function toggleCb(k) {
  cbs[k].checked(!cbs[k].checked());
  show[k] = cbs[k].checked();
  if (!playing) redraw();
}
