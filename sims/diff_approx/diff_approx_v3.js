// =====================================================================
//  Measure-valued branching diffusion  —  individual-based approximation
//  Trait-through-time "waterfall", rendered with an O(width) ring buffer.
//
//  SINGLE-FILE: paste straight into the p5.js web editor and press play.
//  No external libraries (the GUI uses p5's built-in createSlider).
// ---------------------------------------------------------------------
//  Model (per unit time, per individual):
//    birth rate  b(N) = Birth_Rate * exp(-Competition * N / n)        <- competition LOWERS birth
//    death rate  d(y) = Death_Rate + Abiotic_Selection*(theta - y)^2  <- selection RAISES death
//    At a birth event: the parent persists and exactly ONE offspring
//    buds off with trait ~ Normal(parent.y, sqrt(mu / b)).  Because births
//    happen at rate b, a lineage accrues trait-variance at rate b*(mu/b)=mu,
//    so mu is the (scale-invariant) trait diffusion coefficient.
//
//  n = Scaling_Factor refines the IBM toward its diffusion limit:
//    - density terms carry /n so a carrying capacity N* ~ n emerges;
//    - more simulation time is advanced per frame (~sqrt(n)) -> zoom out in TIME;
//    - the trait axis is compressed (~n^{1/4}) -> zoom out in SPACE.
// =====================================================================

// ---------- canvas ----------
let cw = 800;
let ch = 500;

// ---------- live parameters (driven by the sliders below) ----------
let Scaling_Factor = 1.0;     // n : refinement toward the diffusion limit
let Birth_Rate = 0.05;        // intrinsic (uncrowded) birth rate b0
let Death_Rate = 0.02;        // baseline / intrinsic mortality d0
let Abiotic_Selection = 1e-6; // strength of stabilizing selection (on death)
let Competition = 6e-3;       // strength of competition (suppresses birth)
let Abiotic_Optimum = 5 * ch / 8; // optimum trait value (trait units ~ pixels)

// ---------- fixed params ----------
let N0 = 100.0;        // initial abundance per unit of scaling
let mu = 1.0;          // trait diffusion coefficient (variance / unit time)
let baseDt = 1.0;      // sim-time advanced per column when Scaling_Factor == 1
let dtCap = 1.0;       // max integration sub-step (keeps event probabilities small)
let inds = [];         // the population (array of Individual)

// ---------- ring-buffer rendering ----------
let pg;                // offscreen waterfall buffer (cw x ch)
let col;               // ImageData for the single freshly-written column (1 x ch)
let cursor = 0;        // write position (newest column) in the ring buffer

// ---------- display toggles ----------
let bpshow = true;     // population trace
let blshow = true;     // branch (budding) segments
let zshow  = false;    // mean trait marker
let vshow  = false;    // +/- std marker
let tshow  = false;    // optimum tick

// ---------- GUI (built-in p5 sliders) ----------
let ui = [];           // [{key, slider, label, name, fmt}]
let guiVisible = true;

// =====================================================================
function setup() {
  pixelDensity(1);
  createCanvas(cw, ch);
  noSmooth();
  background(0);

  // offscreen ring buffer + reusable 1-pixel-wide column
  pg = createGraphics(cw, ch);
  pg.pixelDensity(1);
  pg.noSmooth();
  pg.background(0);
  col = pg.drawingContext.createImageData(1, ch);

  buildGUI();
  seed();
}

// ---------- minimal control panel using native p5 sliders ----------
function addSlider(key, name, mn, mx, val, step, fmt) {
  const i = ui.length;
  const x = 12, y = 12 + i * 40;
  const lab = createSpan('');
  lab.position(x, y);
  lab.style('color', '#fff');
  lab.style('font-family', 'monospace');
  lab.style('font-size', '11px');
  lab.style('text-shadow', '0 0 3px #000');
  const s = createSlider(mn, mx, val, step);
  s.position(x, y + 16);
  s.style('width', '170px');
  ui.push({ key, slider: s, label: lab, name, fmt });
}

function buildGUI() {
  addSlider('Birth_Rate',        'birth b0',      1e-3, 0.2,   Birth_Rate,        1e-3, v => v.toFixed(3));
  addSlider('Death_Rate',        'death d0',      1e-3, 0.2,   Death_Rate,        1e-3, v => v.toFixed(3));
  addSlider('Abiotic_Selection', 'selection',     0,    2e-6,  Abiotic_Selection, 1e-7, v => v.toExponential(1));
  addSlider('Competition',       'competition',   0,    2e-2,  Competition,       1e-4, v => v.toFixed(4));
  addSlider('Abiotic_Optimum',   'optimum theta', 0,    ch,    Math.round(Abiotic_Optimum), 1, v => v.toFixed(0));
  addSlider('Scaling_Factor',    'scaling n',     1,    50,    Scaling_Factor,    1,    v => v.toFixed(0));
}

function readGUI() {
  for (const u of ui) {
    const v = u.slider.value();
    window[u.key] = v;                       // write into the live global
    u.label.html(u.name + ' = ' + u.fmt(v));
  }
}

function seed() {
  inds = [];
  let N = floor(N0 * Scaling_Factor);
  for (let i = 0; i < N; i++) inds.push(new Individual(Abiotic_Optimum));
  cursor = 0;
}

// ---------- rates ----------
function birthRate(N, n) {
  return Birth_Rate * Math.exp(-Competition * N / n);   // competition lowers birth
}
function deathRate(y) {
  let e = Abiotic_Optimum - y;
  return Death_Rate + Abiotic_Selection * e * e;        // selection raises death
}

// ---------- trait <-> screen mapping (spatial zoom-out with scaling) ----------
//  Continuous in n: time per column ~ sqrt(n)  ->  space ~ sqrt(time) ~ n^{1/4}.
function spaceZoom(n) {
  return Math.pow(n, 0.25);
}
function traitToY(y, zoom) {
  return Abiotic_Optimum + (y - Abiotic_Optimum) / zoom; // identity at zoom == 1
}

// ---------- individual ----------
class Individual {
  constructor(y) { this.y = y; this.dead = false; }
}

// =====================================================================
function draw() {
  readGUI();

  let n = Scaling_Factor;
  let timePerColumn   = baseDt * Math.sqrt(n);          // TIME zoom-out (continuous in n)
  let columnsPerFrame = max(1, round(Math.sqrt(n)));    // literal scroll speed (px / frame)
  let zoom            = spaceZoom(n);                    // SPACE zoom-out (continuous in n)

  // integrate each column in sub-steps small enough to keep event probs low
  let nSub = Math.max(1, Math.ceil(timePerColumn / dtCap));
  let dt   = timePerColumn / nSub;

  // additive brightness, normalised by sqrt(n)
  let addPt = 60 / Math.sqrt(n);   // population trace
  let addBr = 22 / Math.sqrt(n);   // branch segments

  let newest = cursor;
  let data = col.data;

  for (let c = 0; c < columnsPerFrame; c++) {
    // clear the reusable column to opaque black
    for (let k = 0; k < data.length; k += 4) {
      data[k] = 0; data[k + 1] = 0; data[k + 2] = 0; data[k + 3] = 255;
    }

    // advance the process by timePerColumn (one displayed column = one sampled time)
    for (let s = 0; s < nSub; s++) {
      let N = inds.length;
      let b = birthRate(N, n);
      let pBirth = 1 - Math.exp(-b * dt);
      let sigma = Math.sqrt(mu / Math.max(b, 1e-9));   // keeps trait diffusion = mu
      let isLast = (s === nSub - 1);

      let born = [];
      for (let i = 0; i < N; i++) {
        let ind = inds[i];

        // death (stabilizing selection + baseline)
        let pDeath = 1 - Math.exp(-deathRate(ind.y) * dt);
        if (Math.random() < pDeath) { ind.dead = true; continue; }

        // birth: parent stays, exactly one Gaussian offspring buds off
        if (Math.random() < pBirth) {
          let yo = randomGaussian(ind.y, sigma);
          born.push(yo);
          if (blshow) paintSegment(data, traitToY(ind.y, zoom), traitToY(yo, zoom), addBr);
        }
      }

      // compact out the dead, then append the newborns
      let alive = [];
      for (let i = 0; i < N; i++) if (!inds[i].dead) alive.push(inds[i]);
      for (let i = 0; i < born.length; i++) alive.push(new Individual(born[i]));
      inds = alive;

      if (inds.length === 0) inds.push(new Individual(Abiotic_Optimum)); // anti-extinction

      // paint the living population once, on the final sub-step of this column
      if (isLast && bpshow) {
        for (let i = 0; i < inds.length; i++) {
          paintPoint(data, traitToY(inds[i].y, zoom), addPt);
        }
      }
    }

    // commit the column to the ring buffer and advance the write cursor
    pg.drawingContext.putImageData(col, cursor, 0);
    newest = cursor;
    cursor = (cursor + 1) % cw;
  }

  // --- display via two GPU blits so the buffer appears to scroll left ---
  //     newest column sits at the right edge; older data wraps to the left.
  background(0);
  let aW = cw - 1 - newest;
  if (aW > 0) image(pg, 0, 0, aW, ch, newest + 1, 0, aW, ch);
  image(pg, aW, 0, newest + 1, ch, 0, 0, newest + 1, ch);

  drawOverlays(zoom);
}

// ---------- column painters (write directly into the ImageData) ----------
function paintPoint(data, yf, add) {
  let y = yf | 0;
  if (y < 0 || y >= ch) return;
  let idx = y * 4;
  data[idx]     = Math.min(255, data[idx]     + add);
  data[idx + 1] = Math.min(255, data[idx + 1] + add);
  data[idx + 2] = Math.min(255, data[idx + 2] + add);
}
function paintSegment(data, y0f, y1f, add) {
  let y0 = y0f | 0, y1 = y1f | 0;
  if (y0 > y1) { let t = y0; y0 = y1; y1 = t; }
  y0 = Math.max(0, y0); y1 = Math.min(ch - 1, y1);
  for (let y = y0; y <= y1; y++) {
    let idx = y * 4;
    data[idx]     = Math.min(255, data[idx]     + add);
    data[idx + 1] = Math.min(255, data[idx + 1] + add);
    data[idx + 2] = Math.min(255, data[idx + 2] + add);
  }
}

// ---------- overlays on the newest (right-edge) column ----------
function drawOverlays(zoom) {
  let xr = cw - 1;
  if (tshow) {
    let oy = traitToY(Abiotic_Optimum, zoom);
    stroke(255, 0, 0); strokeWeight(4);
    line(xr - 6, oy, xr - 2, oy);
    strokeWeight(1);
  }
  if (zshow && inds.length > 0) {
    let m = meanTrait();
    if (vshow) {
      let sd = stdTrait(m);
      stroke(120, 120, 120, 180);
      line(xr - 6, traitToY(m + sd, zoom), xr - 6, traitToY(m - sd, zoom));
    }
    stroke(255);
    let my = traitToY(m, zoom);
    line(xr - 8, my, xr - 2, my);
  }
}

function meanTrait() {
  let m = 0;
  for (let i = 0; i < inds.length; i++) m += inds[i].y;
  return m / inds.length;
}
function stdTrait(m) {
  let v = 0;
  for (let i = 0; i < inds.length; i++) { let e = inds[i].y - m; v += e * e; }
  return Math.sqrt(v / inds.length);
}

// ---------- keyboard ----------
function keyPressed() {
  switch (key) {
    case 'g':
      guiVisible = !guiVisible;
      for (const u of ui) { if (guiVisible) { u.slider.show(); u.label.show(); } else { u.slider.hide(); u.label.hide(); } }
      break;
    case 'b': bpshow = !bpshow; break;   // population trace
    case 'l': blshow = !blshow; break;   // branch segments
    case 'z': zshow  = !zshow;  break;   // mean
    case 'v': vshow  = !vshow;  break;   // +/- std
    case 't': tshow  = !tshow;  break;   // optimum tick
    case 'r': seed();           break;   // reseed
  }
}
