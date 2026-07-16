// =====================================================================
//  Measure-valued branching diffusion  —  individual-based approximation
//  Trait-through-time "waterfall", rendered with an O(width) ring buffer.
//
//  Light theme: WHITE background, population drawn as accumulating GREY->BLACK ink.
//  Starts PAUSED (so it never runs in the background of a slide deck until asked).
//  Single file; uses p5.js 1.x global mode. No other libraries.
//
//  UI: only the scaling factor is exposed as a slider — every other parameter
//  is fixed. The waterfall sits inside a plot rectangle with a labelled TRAIT
//  axis (y) and a moving TIME axis (x, one pixel column = one time unit).
//
//  AUTO-ZOOM: the dynamics live in a tall fixed "trait buffer"; the y-axis is a
//  window into that buffer that pans/zooms to frame the population. The target
//  window is the central 95% of individuals (a quantile, so a few outliers
//  shooting off-screen exert almost no pull — the pull is count-weighted), and
//  the view is driven toward it by a critically-damped spring (smooth, no jerk).
// ---------------------------------------------------------------------
//  Scaling factor n refines the particle system toward its measure-valued
//  diffusion limit (Fournier-Meleard 2004; Champagnat-Ferriere-Meleard;
//  Etheridge's logistic super-Brownian motion; Dawson-Watanabe).
//    initial count  N0 * n
//    birth rate b(N) = sqrt(n) * Birth_Rate * exp(-Competition * N / n)
//    death rate d(y) = sqrt(n) * Death_Rate + Abiotic_Selection*(theta - y)^2
//    mutation    y' ~ Normal(parent.y, sqrt(mu / b))     (variance ~ 1/sqrt(n))
//  sigma^2 = mu/b keeps (rate)*(variance) = mu, so the trait -> Brownian motion
//  with diffusion coefficient mu in the limit.
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
const PLOT_W = cw - PAD_L - PAD_R;      // waterfall width  == ring-buffer length
const PLOT_H = ch - PAD_T - PAD_B;      // on-screen plot height (a WINDOW into the buffer)

// ---------- trait buffer (dynamics live here; the view is a sub-window of it) ----------
const TH = 1200;                        // buffer height in "trait pixels" (tall, fixed)
const Abiotic_Optimum = TH / 2;         // optimum trait at the buffer centre; FIXED
const TRAIT_UNIT_PX = 50;               // trait pixels per labelled unit
function rowToTrait(row) { return (Abiotic_Optimum - row) / TRAIT_UNIT_PX; }

// ---------- auto-zoom view (a window [viewLo, viewHi] of buffer rows) ----------
const VIEW_QLO = 0.025, VIEW_QHI = 0.975;   // frame the central 95% of individuals
const VIEW_PAD = 1.30;                      // padding around that band
const VIEW_MIN_HALF = 35;                   // rows: don't zoom in tighter than +/-0.7 units
const VIEW_MAX_HALF = 500;                  // rows: don't zoom out past this
const VIEW_OMEGA_C = 0.18;                  // spring rate for the centre (per frame)
const VIEW_OMEGA_H = 0.14;                  // spring rate for the half-height (calmer)
let viewC = Abiotic_Optimum, viewH = 120;   // spring position (centre, half-height)
let viewVC = 0, viewVH = 0;                 // spring velocities
let viewLo = viewC - viewH, viewHi = viewC + viewH;   // current displayed window
function screenYofRow(row) { return PLOT_Y + (row - viewLo) / (viewHi - viewLo) * PLOT_H; }

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

// ---------- ring-buffer rendering ----------
let pg, col, cursor = 0;
let colAcc = 0;
let simTime = 0;                        // total columns rendered == elapsed time (1 column = 1 unit)

// ---------- display toggles ----------
let bpshow = true;     // lifelines
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

  pg = createGraphics(PLOT_W, TH);       // tall trait buffer
  pg.pixelDensity(1);
  pg.noSmooth();
  pg.background(255);
  col = pg.drawingContext.createImageData(1, TH);

  buildGUI();
  seed();

  // pause whenever the browser tab is hidden; resume only if the user had it playing
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

// ---------- control panel (native p5 sliders; setter closure writes the binding) ----------
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
  for (let i = 0; i < N; i++) inds.push(new Individual(Abiotic_Optimum));
  cursor = 0;
  colAcc = 0;
  simTime = 0;
  viewC = Abiotic_Optimum; viewH = 120; viewVC = 0; viewVH = 0;   // reset the view spring
  viewLo = viewC - viewH; viewHi = viewC + viewH;
  if (pg) pg.background(255);
}

// ---------- n-dependent rates ----------
function birthRate(N, n) { return Math.sqrt(n) * Birth_Rate * Math.exp(-Competition * N / n); }
function deathRate(y, n) { let e = Abiotic_Optimum - y; return Math.sqrt(n) * Death_Rate + Abiotic_Selection * e * e; }

class Individual { constructor(y) { this.y = y; this.dead = false; } }

// ---------- auto-zoom: robust target window + critically-damped spring ----------
function quantileSorted(a, p) {
  const n = a.length;
  if (n === 1) return a[0];
  const idx = (n - 1) * p, lo = Math.floor(idx), hi = Math.ceil(idx), frac = idx - lo;
  return a[lo] * (1 - frac) + a[hi] * frac;
}

function updateView() {
  let targetC = viewC, targetH = viewH;
  const N = inds.length;
  if (N > 0) {
    const rows = new Float64Array(N);
    for (let i = 0; i < N; i++) rows[i] = inds[i].y;
    rows.sort();                                    // ascending
    const ql = quantileSorted(rows, VIEW_QLO);
    const qh = quantileSorted(rows, VIEW_QHI);
    targetC = 0.5 * (ql + qh);                      // centre of the bulk
    targetH = constrain(0.5 * (qh - ql) * VIEW_PAD, VIEW_MIN_HALF, VIEW_MAX_HALF);
  }
  // critically-damped springs (semi-implicit Euler, dt = 1 frame): smooth, no overshoot
  viewVC += VIEW_OMEGA_C * VIEW_OMEGA_C * (targetC - viewC) - 2 * VIEW_OMEGA_C * viewVC;
  viewC  += viewVC;
  viewVH += VIEW_OMEGA_H * VIEW_OMEGA_H * (targetH - viewH) - 2 * VIEW_OMEGA_H * viewVH;
  viewH  += viewVH;

  viewH = constrain(viewH, VIEW_MIN_HALF, VIEW_MAX_HALF);
  viewC = constrain(viewC, 0, TH);
  viewLo = viewC - viewH; viewHi = viewC + viewH;
  if (viewLo < 0)  { viewLo = 0; }
  if (viewHi > TH) { viewHi = TH; }
  if (viewHi - viewLo < 2) viewHi = viewLo + 2;
}

// =====================================================================
function draw() {
  readGUI();

  let n = Scaling_Factor;
  let rootN = Math.sqrt(n);
  let sub = max(1, round(rootN));

  colAcc += rootN;
  let columnsPerFrame = Math.floor(colAcc);
  colAcc -= columnsPerFrame;
  if (columnsPerFrame < 1) columnsPerFrame = 1;
  let dt = 1.0 / sub;

  let radius = 1.3 / rootN;
  let addPt = 180 / rootN;
  let addBr = 150 / n;

  let newest = cursor;
  let data = col.data;

  for (let c = 0; c < columnsPerFrame; c++) {
    for (let k = 0; k < data.length; k += 4) { data[k]=255; data[k+1]=255; data[k+2]=255; data[k+3]=255; } // white column

    for (let s = 0; s < sub; s++) {
      let N = inds.length;
      let b = birthRate(N, n);
      let pBirth = 1 - Math.exp(-b * dt);
      let sigma = Math.sqrt(mu / Math.max(b, 1e-9));
      let isLast = (s === sub - 1);

      let born = [];
      for (let i = 0; i < N; i++) {
        let ind = inds[i];
        if (Math.random() < 1 - Math.exp(-deathRate(ind.y, n) * dt)) { ind.dead = true; continue; }
        if (Math.random() < pBirth) {
          let yo = randomGaussian(ind.y, sigma);
          born.push(yo);
          if (blshow && addBr > 1.5) paintDashed(data, ind.y, yo, addBr);
        }
      }

      let alive = [];
      for (let i = 0; i < N; i++) if (!inds[i].dead) alive.push(inds[i]);
      for (let i = 0; i < born.length; i++) alive.push(new Individual(born[i]));
      inds = alive;
      if (inds.length === 0) inds.push(new Individual(Abiotic_Optimum));

      if (isLast && bpshow) {
        for (let i = 0; i < inds.length; i++) paintMark(data, inds[i].y, radius, addPt);
      }
    }

    pg.drawingContext.putImageData(col, cursor, 0);
    newest = cursor;
    cursor = (cursor + 1) % PLOT_W;
  }

  simTime += columnsPerFrame;
  updateView();                        // pan/zoom the y-axis toward the population

  background(255);
  // blit the current vertical WINDOW of the ring buffer, scaled to the plot rect.
  // horizontal is 1:1 (crisp time columns); only the vertical trait axis is scaled.
  const vy = viewLo, vh = viewHi - viewLo;
  const dctx = drawingContext;
  dctx.imageSmoothingEnabled = true;   // smooth vertical zoom
  let aW = PLOT_W - 1 - newest;
  if (aW > 0) image(pg, PLOT_X, PLOT_Y, aW, PLOT_H, newest + 1, vy, aW, vh);
  image(pg, PLOT_X + aW, PLOT_Y, newest + 1, PLOT_H, 0, vy, newest + 1, vh);
  dctx.imageSmoothingEnabled = false;

  drawAxes();
  drawOverlays();
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

// ---------- axes (drawn in the margins, on top of the white background) ----------
function drawAxes() {
  push();
  textFont('monospace');

  // plot frame
  noFill(); stroke(150); strokeWeight(1);
  rect(PLOT_X - 0.5, PLOT_Y - 0.5, PLOT_W, PLOT_H);

  // ---- y-axis: trait value (0 at the optimum), ticks follow the auto-zoom ----
  let uTop = rowToTrait(viewLo), uBot = rowToTrait(viewHi);   // uTop > uBot
  let step = niceStep((uTop - uBot) / 6);
  textAlign(RIGHT, CENTER); textSize(11);
  for (let u = Math.ceil(uBot / step) * step; u <= uTop + 1e-9; u += step) {
    let ys = screenYofRow(Abiotic_Optimum - u * TRAIT_UNIT_PX);
    if (ys < PLOT_Y - 0.5 || ys > PLOT_Y + PLOT_H + 0.5) continue;
    let isZero = Math.abs(u) < 1e-9;
    stroke(150); strokeWeight(1); line(PLOT_X - 4, ys, PLOT_X, ys);
    noStroke(); fill(isZero ? 30 : 95);
    text(fmtTick(u, step), PLOT_X - 7, ys);
  }
  // y-axis title (rotated)
  push();
  translate(15, PLOT_Y + PLOT_H / 2); rotate(-HALF_PI);
  noStroke(); fill(40); textAlign(CENTER, CENTER); textSize(12);
  text('trait value  z', 0, 0);
  pop();

  // ---- x-axis: time (moving; ticks scroll left with the data) ----
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

// ---------- column painters (subtractive: darken white toward black) ----------
function ink(data, y, amt) {
  let idx = y * 4;
  data[idx]   = Math.max(0, data[idx]   - amt);
  data[idx+1] = Math.max(0, data[idx+1] - amt);
  data[idx+2] = Math.max(0, data[idx+2] - amt);
}
function paintMark(data, yf, r, add) {
  let yc = Math.round(yf);
  let R = Math.ceil(r + 0.5);
  for (let dy = -R; dy <= R; dy++) {
    let y = yc + dy; if (y < 0 || y >= TH) continue;
    let cov = r + 0.5 - Math.abs(dy);
    if (cov > 0) ink(data, y, add * Math.min(1, cov));
  }
}
function paintDashed(data, y0f, y1f, add) {
  let y0 = y0f | 0, y1 = y1f | 0; if (y0 > y1) { let t = y0; y0 = y1; y1 = t; }
  y0 = Math.max(0, y0); y1 = Math.min(TH - 1, y1);
  for (let y = y0; y <= y1; y++) if (((y - y0) & 3) < 2) ink(data, y, add);
}

// ---------- overlays on the newest (right-edge) column, mapped through the view ----------
function drawOverlays() {
  let xr = PLOT_X + PLOT_W - 1;
  if (tshow) {
    let yo = screenYofRow(Abiotic_Optimum);
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
