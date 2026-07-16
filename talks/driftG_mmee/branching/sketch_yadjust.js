// =====================================================================
//  Measure-valued branching diffusion  —  individual-based approximation
//  Trait-through-time "waterfall".
//
//  Light theme: WHITE background, population drawn as accumulating GREY->BLACK ink.
//  Starts PAUSED (so it never runs in the background of a slide deck until asked).
//  Single file; uses p5.js 1.x global mode. No other libraries.
//
//  RENDERING: individuals are stored as points (per time-column) and re-drawn
//  every frame as fixed-SCREEN-size dots through the current view transform.
//  So zooming changes where individuals sit and how far apart they are, but NOT
//  how thick their marks are -- it is a real zoom, not a vertical stretch.
//
//  AUTO-ZOOM (y-axis) with a DEADBAND: the view stays completely still while the
//  bulk of the population (central 95%) sits comfortably inside it. It only
//  re-frames when the bulk spills toward an edge or shrinks enough to leave too
//  much empty space -- and even then it eases over with a critically-damped
//  spring. A few outliers shooting off-screen don't move it (the trigger is a
//  quantile, i.e. count-weighted).
// ---------------------------------------------------------------------
//  Scaling factor n refines the particle system toward its measure-valued
//  diffusion limit (Fournier-Meleard 2004; Champagnat-Ferriere-Meleard;
//  Etheridge's logistic super-Brownian motion; Dawson-Watanabe).
//    birth rate b(N) = sqrt(n) * Birth_Rate * exp(-Competition * N / n)
//    death rate d(y) = sqrt(n) * Death_Rate + Abiotic_Selection * y^2   (optimum at 0)
//    mutation    y' ~ Normal(parent.y, sqrt(mu / b))
// =====================================================================

// ---------- canvas ----------
let cw = 800;
let ch = 300;

// ---------- plot layout (margins carve out labelled axes around the waterfall) ----------
const PAD_L = 60;                       // left  margin: trait (y) axis
const PAD_R = 16;                       // right margin
const PAD_T = 44;                       // top   margin: control strip (scaling slider + play)
const PAD_B = 38;                       // bottom margin: time (x) axis
const PLOT_X = PAD_L;
const PLOT_Y = PAD_T;
const PLOT_W = cw - PAD_L - PAD_R;      // waterfall width  == ring-buffer length (columns)
const PLOT_H = ch - PAD_T - PAD_B;      // waterfall height (pixels on screen)

// ---------- trait axis ----------
const TRAIT_UNIT_PX = 50;               // trait-pixels per labelled unit (optimum = trait 0)

// ---------- auto-zoom view: centred on the optimum (trait 0); only the ZOOM adapts ----------
//  The view holds still until the population's spread crosses a deadband, then eases
//  to a new zoom with a spring. The centre is pinned, so the frame never chases the
//  wandering sample mean -- the population drifts freely within a stable frame.
const VIEW_QLO = 0.025, VIEW_QHI = 0.975;   // "bulk" = central 95% of individuals
const VIEW_ALPHA = 0.025;                   // smoothing of the needed half-extent (per frame)
const VIEW_EDGE  = 0.13;                    // zoom OUT if bulk reaches (1-EDGE) of the half-view
const VIEW_SHRINK = 0.40;                   // zoom IN  if bulk fills < SHRINK of the half-view
const VIEW_PAD = 1.50;                      // padding when (re)framing
const VIEW_OMEGA = 0.10;                    // spring rate (per frame); smaller = floatier
const VIEW_MIN_HALF = 30;                   // trait-px: tightest zoom (+/- 0.6 units)
const VIEW_MAX_HALF = 400;                  // trait-px: widest zoom
let hV = 90, velH = 0, emaH = 90;           // half-view height, spring velocity, smoothed extent
let viewLo = -hV, viewHi = hV;              // window is always centred on the optimum (0)
function screenYofRow(y) { return PLOT_Y + (viewHi - y) / (viewHi - viewLo) * PLOT_H; }

// ---------- parameters ----------
// Only the scaling factor has a slider; the rest are fixed at these values.
let Scaling_Factor = 1.0;
let Birth_Rate = 0.05;
let Death_Rate = 0.02;
let Abiotic_Selection = 1.5e-7;
let Competition = 0.03;
let mu = 2.0;

// ---------- fixed params ----------
let N0 = 25.0;
let inds = [];

// ---------- point storage (preallocated ring buffers; column i lives at offset i*K) ----------
const K_STORE = 600;                    // max individuals stored/drawn per column (subsample above)
const CON_MAX = 80;                     // max birth connectors stored per column
let ptsBuf, colCount, colInk;           // individual trait positions + per-column ink
let conBuf, conCount;                   // birth-connector endpoint pairs
let screenImg = null;                   // screen-space plot raster (rebuilt each frame)
let cursor = 0, newest = 0, simTime = 0;

// ---------- display toggles ----------
let bpshow = true;     // lifelines (individual marks)
let blshow = true;     // dashed birth connectors
let zshow  = false;    // mean trait marker
let vshow  = false;    // +/- std marker
let tshow  = false;    // optimum tick

// ---------- GUI / playback ----------
let ui = [];
let btns = [];
let scalingSlider, scalingUI, playBtn;
let guiVisible = true;
let host = null;
let playing = false;   // starts paused

// =====================================================================
function setup() {
  pixelDensity(1);
  host = (typeof document !== 'undefined') ? document.getElementById('sim') : null;
  let cnv = createCanvas(cw, ch);
  if (host) cnv.parent(host);
  noSmooth();
  background(255);

  ptsBuf   = new Float32Array(PLOT_W * K_STORE);
  colCount = new Int32Array(PLOT_W);
  colInk   = new Float32Array(PLOT_W);
  conBuf   = new Float32Array(PLOT_W * CON_MAX * 2);
  conCount = new Int32Array(PLOT_W);
  screenImg = drawingContext.createImageData(PLOT_W, PLOT_H);

  buildGUI();
  seed();

  if (typeof document !== 'undefined') {
    document.addEventListener('visibilitychange', () => {
      if (document.hidden) noLoop();
      else if (playing) loop();
    });
  }

  setPlaying(false);   // initialize PAUSED
}

// ---------- playback ----------
function setPlaying(p) {
  playing = p;
  if (p) loop(); else noLoop();
  if (playBtn) playBtn.html(p ? '\u23F8 pause' : '\u25B6 play');
}

// ---------- control panel ----------
function addSlider(name, mn, mx, val, step, fmt, set) {
  const i = ui.length, x = 12, y = 8 + i * 38;
  const lab = createSpan('');
  lab.position(x, y);
  lab.style('color', '#111'); lab.style('font-family', 'monospace');
  lab.style('font-size', '11px'); lab.style('text-shadow', '0 0 2px #fff, 0 0 2px #fff');
  const s = createSlider(mn, mx, val, step);
  s.position(x, y + 15); s.style('width', '170px');

  const u = { slider: s, label: lab, name, fmt, set };
  const refresh = () => { const v = s.value(); set(v); lab.html(name + ' = ' + fmt(v)); };
  u.refresh = refresh;
  refresh();
  s.input(refresh);

  if (host) { lab.parent(host); s.parent(host); }
  ui.push(u);
  return u;
}

function buildGUI() {
  scalingUI = addSlider('scaling k', 1, 50, Scaling_Factor, 1, v => v.toFixed(0), v => Scaling_Factor = v);
  scalingSlider = scalingUI.slider;

  const by = 8 + (ui.length - 1) * 38 + 15;
  const mk = (lbl, x, fn) => {
    const btn = createButton(lbl);
    btn.position(x, by);
    btn.style('font', '12px monospace'); btn.style('padding', '0 6px');
    btn.mousePressed(fn);
    if (host) btn.parent(host);
    btns.push(btn);
  };
  mk('\u25C0', 190, () => nudgeN(-1));
  mk('\u25B6', 222, () => nudgeN(+1));

  playBtn = createButton('\u25B6 play');
  playBtn.position(cw - 86, 10);
  playBtn.style('font', '12px monospace'); playBtn.style('padding', '2px 8px');
  playBtn.mousePressed(() => setPlaying(!playing));
  if (host) playBtn.parent(host);
  btns.push(playBtn);
}

function nudgeN(d) {
  let v = constrain(round(scalingSlider.value()) + d, 1, 50);
  scalingSlider.value(v);
  scalingUI.refresh();
}

function readGUI() {
  for (const u of ui) { const v = u.slider.value(); u.set(v); u.label.html(u.name + ' = ' + u.fmt(v)); }
}

function seed() {
  inds = [];
  let N = floor(N0 * Scaling_Factor);
  for (let i = 0; i < N; i++) inds.push(new Individual(0));   // optimum = trait 0
  cursor = 0; newest = 0; simTime = 0;
  if (colCount) colCount.fill(0);
  if (conCount) conCount.fill(0);
  hV = 90; velH = 0; emaH = 90;                               // reset the view spring
  viewLo = -hV; viewHi = hV;
}

// ---------- n-dependent rates (optimum at 0) ----------
function birthRate(N, n) { return Math.sqrt(n) * Birth_Rate * Math.exp(-Competition * N / n); }
function deathRate(y, n) { return Math.sqrt(n) * Death_Rate + Abiotic_Selection * y * y; }

class Individual { constructor(y) { this.y = y; this.dead = false; } }

// ---------- store one column's snapshot (reservoir-subsample) + its connectors ----------
function storeColumn(ci, arr, ink, cons) {
  const pop = arr.length;
  const cnt = Math.min(pop, K_STORE);
  const base = ci * K_STORE;
  if (pop <= K_STORE) {
    for (let i = 0; i < pop; i++) ptsBuf[base + i] = arr[i].y;
  } else {                                        // uniform random subsample (preserves shape)
    for (let i = 0; i < K_STORE; i++) ptsBuf[base + i] = arr[i].y;
    for (let i = K_STORE; i < pop; i++) {
      const j = (Math.random() * (i + 1)) | 0;
      if (j < K_STORE) ptsBuf[base + j] = arr[i].y;
    }
  }
  colCount[ci] = cnt;
  colInk[ci] = ink * (pop / cnt);                 // boost ink so darkness ~ true density
  const nc = Math.min((cons.length / 2) | 0, CON_MAX);
  const cbase = ci * CON_MAX * 2;
  for (let i = 0; i < nc * 2; i++) conBuf[cbase + i] = cons[i];
  conCount[ci] = nc;
}

// ---------- auto-zoom: deadband target + critically-damped spring ----------
function quantileSorted(a, p) {
  const n = a.length;
  if (n === 1) return a[0];
  const idx = (n - 1) * p, lo = Math.floor(idx), hi = Math.ceil(idx), frac = idx - lo;
  return a[lo] * (1 - frac) + a[hi] * frac;
}

function updateView() {
  const N = inds.length;
  let reqH = emaH;                                    // half-extent needed to hold the bulk
  if (N > 0) {
    const rows = new Float64Array(N);
    for (let i = 0; i < N; i++) rows[i] = inds[i].y;
    rows.sort();
    const ql = quantileSorted(rows, VIEW_QLO), qh = quantileSorted(rows, VIEW_QHI);
    reqH = Math.max(Math.abs(ql), Math.abs(qh));      // distance from optimum to the 95% edge
  }
  emaH += VIEW_ALPHA * (reqH - emaH);                 // smooth it (ignore brief excursions)

  let tH = hV;                                        // default: HOLD STILL (deadband)
  if (emaH > hV * (1 - VIEW_EDGE) || emaH < hV * VIEW_SHRINK) {
    tH = constrain(emaH * VIEW_PAD, VIEW_MIN_HALF, VIEW_MAX_HALF);
  }
  // critically-damped zoom spring (semi-implicit Euler, dt = 1 frame): smooth, no overshoot
  velH += VIEW_OMEGA * VIEW_OMEGA * (tH - hV) - 2 * VIEW_OMEGA * velH;
  hV += velH;
  hV = constrain(hV, VIEW_MIN_HALF, VIEW_MAX_HALF);
  viewLo = -hV; viewHi = hV;                          // centre pinned at the optimum
}

// =====================================================================
function draw() {
  readGUI();

  let n = Scaling_Factor;
  let rootN = Math.sqrt(n);
  let sub = max(1, round(rootN));

  colAccAdvance(rootN);
  let dt = 1.0 / sub;
  let addPt = 180 / rootN;             // per-individual ink at full population
  let addBr = 150 / n;                 // connector ink (fades fast in n)

  for (let c = 0; c < columnsPerFrame; c++) {
    let cons = [];
    for (let s = 0; s < sub; s++) {
      let N = inds.length;
      let b = birthRate(N, n);
      let pBirth = 1 - Math.exp(-b * dt);
      let sigma = Math.sqrt(mu / Math.max(b, 1e-9));

      let born = [];
      for (let i = 0; i < N; i++) {
        let ind = inds[i];
        if (Math.random() < 1 - Math.exp(-deathRate(ind.y, n) * dt)) { ind.dead = true; continue; }
        if (Math.random() < pBirth) {
          let yo = randomGaussian(ind.y, sigma);
          born.push(yo);
          if (blshow && addBr > 1.5 && cons.length < CON_MAX * 2) { cons.push(ind.y, yo); }
        }
      }
      let alive = [];
      for (let i = 0; i < N; i++) if (!inds[i].dead) alive.push(inds[i]);
      for (let i = 0; i < born.length; i++) alive.push(new Individual(born[i]));
      inds = alive;
      if (inds.length === 0) inds.push(new Individual(0));
    }
    storeColumn(cursor, inds, addPt, cons);
    newest = cursor;
    cursor = (cursor + 1) % PLOT_W;
  }

  simTime += columnsPerFrame;
  updateView();                        // pan/zoom the y-axis (deadband + spring)

  background(255);
  renderWaterfall(rootN);              // re-draw individuals through the current view
  drawAxes();
  drawOverlays();
}

// smooth-scroll accumulator (advance ~sqrt(n) columns/frame on average)
let colAcc = 0, columnsPerFrame = 1;
function colAccAdvance(rootN) {
  colAcc += rootN;
  columnsPerFrame = Math.floor(colAcc);
  colAcc -= columnsPerFrame;
  if (columnsPerFrame < 1) columnsPerFrame = 1;
}

// ---------- waterfall: fixed-size dots at view-mapped positions (true zoom) ----------
function renderWaterfall(rootN) {
  const sd = screenImg.data;
  sd.fill(255);
  const rDot = 1.3 / rootN;                      // FIXED screen half-height of a mark
  const scale = PLOT_H / (viewHi - viewLo);
  const showCon = (150 / Scaling_Factor) > 1.5;
  const addBrNow = 150 / Scaling_Factor;

  for (let i = 0; i < PLOT_W; i++) {
    const cnt = colCount[i];
    if (cnt === 0) continue;
    const sx = (i - newest - 1 + PLOT_W) % PLOT_W;   // column -> screen x (newest at right)

    if (blshow && showCon) {
      const nc = conCount[i], cbase = i * CON_MAX * 2;
      for (let k = 0; k < nc; k++) {
        const p = conBuf[cbase + 2 * k], q = conBuf[cbase + 2 * k + 1];
        paintConnImg(sd, sx, (viewHi - p) * scale, (viewHi - q) * scale, addBrNow);
      }
    }
    if (bpshow) {
      const ink = colInk[i], base = i * K_STORE;
      for (let k = 0; k < cnt; k++) {
        paintDotImg(sd, sx, (viewHi - ptsBuf[base + k]) * scale, rDot, ink);
      }
    }
  }
  drawingContext.putImageData(screenImg, PLOT_X, PLOT_Y);
}

function paintDotImg(sd, sx, syf, r, ink) {
  const yc = Math.round(syf), R = Math.ceil(r + 0.5);
  for (let dy = -R; dy <= R; dy++) {
    const y = yc + dy; if (y < 0 || y >= PLOT_H) continue;
    const cov = r + 0.5 - Math.abs(dy); if (cov <= 0) continue;
    const idx = (y * PLOT_W + sx) * 4, a = ink * Math.min(1, cov);
    sd[idx] = Math.max(0, sd[idx] - a);
    sd[idx + 1] = Math.max(0, sd[idx + 1] - a);
    sd[idx + 2] = Math.max(0, sd[idx + 2] - a);
  }
}
function paintConnImg(sd, sx, y0f, y1f, ink) {
  let y0 = Math.round(y0f), y1 = Math.round(y1f);
  if (y0 > y1) { const t = y0; y0 = y1; y1 = t; }
  y0 = Math.max(0, y0); y1 = Math.min(PLOT_H - 1, y1);
  for (let y = y0; y <= y1; y++) {
    if (((y - y0) & 3) < 2) {
      const idx = (y * PLOT_W + sx) * 4;
      sd[idx] = Math.max(0, sd[idx] - ink);
      sd[idx + 1] = Math.max(0, sd[idx + 1] - ink);
      sd[idx + 2] = Math.max(0, sd[idx + 2] - ink);
    }
  }
}

// ---------- nice tick helpers ----------
function niceStep(raw) {
  if (!(raw > 0)) return 1;
  const m = Math.pow(10, Math.floor(Math.log10(raw)));
  const nrm = raw / m;
  return (nrm < 1.5 ? 1 : nrm < 3.5 ? 2 : nrm < 7.5 ? 5 : 10) * m;
}
function fmtTick(u, step) {
  const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return u.toFixed(d);
}

// ---------- axes ----------
function drawAxes() {
  push();
  textFont('monospace');
  noFill(); stroke(150); strokeWeight(1);
  rect(PLOT_X - 0.5, PLOT_Y - 0.5, PLOT_W, PLOT_H);

  // y-axis: trait value (0 at optimum), ticks follow the auto-zoom
  let uTop = viewHi / TRAIT_UNIT_PX, uBot = viewLo / TRAIT_UNIT_PX;
  let step = niceStep((uTop - uBot) / 6);
  textAlign(RIGHT, CENTER); textSize(11);
  for (let u = Math.ceil(uBot / step) * step; u <= uTop + 1e-9; u += step) {
    let ys = screenYofRow(u * TRAIT_UNIT_PX);
    if (ys < PLOT_Y - 0.5 || ys > PLOT_Y + PLOT_H + 0.5) continue;
    let isZero = Math.abs(u) < 1e-9;
    stroke(150); strokeWeight(1); line(PLOT_X - 4, ys, PLOT_X, ys);
    noStroke(); fill(isZero ? 30 : 95);
    text(fmtTick(u, step), PLOT_X - 7, ys);
  }
  push();
  translate(15, PLOT_Y + PLOT_H / 2); rotate(-HALF_PI);
  noStroke(); fill(40); textAlign(CENTER, CENTER); textSize(12);
  text('trait value  z', 0, 0);
  pop();

  // x-axis: time (moving; ticks scroll left with the data)
  let axisY = PLOT_Y + PLOT_H;
  const dTime = 100;
  textAlign(CENTER, TOP); textSize(11);
  let firstTick = Math.ceil((simTime - (PLOT_W - 1)) / dTime) * dTime;
  for (let T = firstTick; T <= simTime; T += dTime) {
    if (T < 0) continue;
    let xs = PLOT_X + (PLOT_W - 1) - (simTime - T);
    if (xs < PLOT_X - 0.5 || xs > PLOT_X + PLOT_W + 0.5) continue;
    stroke(150); strokeWeight(1); line(xs, axisY, xs, axisY + 4);
    noStroke(); fill(95);
    text(T, xs, axisY + 6);
  }
  noStroke(); fill(40); textAlign(CENTER, TOP); textSize(12);
  text('time', PLOT_X + PLOT_W / 2, axisY + 20);
  pop();
}

// ---------- overlays on the newest (right-edge) column, mapped through the view ----------
function drawOverlays() {
  let xr = PLOT_X + PLOT_W - 1;
  if (tshow) {
    let yo = screenYofRow(0);
    stroke(220, 0, 0); strokeWeight(4); line(xr - 6, yo, xr - 2, yo); strokeWeight(1);
  }
  if (zshow && inds.length > 0) {
    let m = meanTrait();
    if (vshow) {
      let sd = stdTrait(m);
      stroke(120, 120, 120, 180);
      line(xr - 6, screenYofRow(m + sd), xr - 6, screenYofRow(m - sd));
    }
    stroke(0); line(xr - 8, screenYofRow(m), xr - 2, screenYofRow(m));
  }
}
function meanTrait() { let m = 0; for (let i = 0; i < inds.length; i++) m += inds[i].y; return m / inds.length; }
function stdTrait(m) { let v = 0; for (let i = 0; i < inds.length; i++) { let e = inds[i].y - m; v += e * e; } return Math.sqrt(v / inds.length); }

// ---------- keyboard ----------
function keyPressed() {
  if (key === ' ')             { setPlaying(!playing); return false; }
  if (keyCode === LEFT_ARROW)  { nudgeN(-1); return false; }
  if (keyCode === RIGHT_ARROW) { nudgeN(+1); return false; }
  switch (key) {
    case 'p': setPlaying(!playing); break;
    case 'g': guiVisible = !guiVisible;
      for (const u of ui) { if (guiVisible) { u.slider.show(); u.label.show(); } else { u.slider.hide(); u.label.hide(); } }
      for (const btn of btns) { if (guiVisible) btn.show(); else btn.hide(); } break;
    case 'b': bpshow = !bpshow; break;
    case 'l': blshow = !blshow; break;
    case 'z': zshow  = !zshow;  break;
    case 'v': vshow  = !vshow;  break;
    case 't': tshow  = !tshow;  break;
    case 'r': seed();           break;
  }
}
