// =====================================================================
//  Measure-valued branching diffusion  —  individual-based approximation
//  Trait-through-time "waterfall", rendered with an O(width) ring buffer.
// ---------------------------------------------------------------------
//  Model (per unit time, per individual):
//    birth rate  b(N) = Birth_Rate * exp(-Competition * N / n)   <- competition LOWERS birth
//    death rate  d(y) = Death_Rate + Abiotic_Selection*(theta - y)^2  <- selection RAISES death
//    At a birth event: the parent persists and exactly ONE offspring
//    buds off with trait ~ Normal(parent.y, sqrt(mu / b)).  Because births
//    happen at rate b, a lineage accrues trait-variance at rate b*(mu/b)=mu,
//    so mu is the (scale-invariant) trait diffusion coefficient.
//
//  n = Scaling_Factor refines the IBM toward its diffusion limit:
//    - density terms carry /n so a carrying capacity N* ~ n emerges
//      (each individual effectively has mass 1/n);
//    - more simulation time is advanced per frame (~sqrt(n)) -> zoom out in TIME;
//    - the trait axis is compressed (~n^{1/4}) so a Brownian cloud keeps a
//      constant on-screen width -> zoom out in SPACE.
// =====================================================================

// ---------- canvas ----------
let cw = 800;
let ch = 500;

// ---------- GUI params (p5.gui convention: X / Xmin / Xmax / Xstep) ----------
var Scaling_Factor = 1.0;            // n : refinement toward the diffusion limit
var Scaling_FactorMin = 1.0;
var Scaling_FactorMax = 50.0;
var Scaling_FactorStep = 1.0;

var Birth_Rate = 0.05;               // intrinsic (uncrowded) birth rate b0
var Birth_RateMin = 1e-3;
var Birth_RateMax = 0.2;
var Birth_RateStep = 1e-3;

var Death_Rate = 0.02;               // baseline / intrinsic mortality d0
var Death_RateMin = 1e-3;
var Death_RateMax = 0.2;
var Death_RateStep = 1e-3;

var Abiotic_Selection = 1e-6;        // strength of stabilizing selection (on death)
var Abiotic_SelectionMin = 0.0;
var Abiotic_SelectionMax = 2e-6;
var Abiotic_SelectionStep = 1e-7;

var Competition = 6e-3;              // strength of competition (suppresses birth)
var CompetitionMin = 0.0;
var CompetitionMax = 2e-2;
var CompetitionStep = 1e-4;

var Abiotic_Optimum = 5 * ch / 8;    // optimum trait value (trait units ~ pixels)
var Abiotic_OptimumMin = 0;
var Abiotic_OptimumMax = ch;
var Abiotic_OptimumStep = 1;

// ---------- fixed params ----------
let N0 = 100.0;        // initial abundance per unit of scaling
let mu = 1.0;          // trait diffusion coefficient (variance / unit time)
let dt = 1.0;          // integration sub-step (per-unit-time rates * dt = event prob)
let baseSteps = 1;     // sub-steps per frame when Scaling_Factor == 1
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
let visible = true;
let gui;

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

  gui = createGui('Parameters');
  gui.addGlobals(
    'Birth_Rate', 'Death_Rate', 'Abiotic_Selection',
    'Competition', 'Abiotic_Optimum', 'Scaling_Factor'
  );

  seed();
  gui.show();
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
function spaceZoom(stepsPerFrame) {
  return Math.sqrt(stepsPerFrame / baseSteps);          // ~ n^{1/4}
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
  let n = Scaling_Factor;
  let stepsPerFrame = max(1, round(baseSteps * Math.sqrt(n))); // more sim-time / frame
  let zoom = spaceZoom(stepsPerFrame);

  // --- clear the reusable column to opaque black ---
  let data = col.data;
  for (let k = 0; k < data.length; k += 4) {
    data[k] = 0; data[k + 1] = 0; data[k + 2] = 0; data[k + 3] = 255;
  }

  // additive brightness, normalised by sqrt(n) so total ink ~ constant as N ~ n
  let addPt = 60 / Math.sqrt(n);   // population trace
  let addBr = 22 / Math.sqrt(n);   // branch segments

  // --- advance the process; each frame writes ONE column (a sampled snapshot) ---
  for (let s = 0; s < stepsPerFrame; s++) {
    let N = inds.length;
    let b = birthRate(N, n);
    let pBirth = 1 - Math.exp(-b * dt);
    let sigma = Math.sqrt(mu / Math.max(b, 1e-9));   // keeps trait diffusion = mu
    let isLast = (s === stepsPerFrame - 1);

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

    // paint the living population once, on the final sub-step
    if (isLast && bpshow) {
      for (let i = 0; i < inds.length; i++) {
        paintPoint(data, traitToY(inds[i].y, zoom), addPt);
      }
    }
  }

  // --- commit the column to the ring buffer ---
  pg.drawingContext.putImageData(col, cursor, 0);

  // --- display via two GPU blits so the buffer appears to scroll left ---
  //     newest column (cursor) sits at the right edge; older data wraps to the left.
  background(0);
  let aW = cw - 1 - cursor;                            // width of the "older" left block
  if (aW > 0) image(pg, 0, 0, aW, ch, cursor + 1, 0, aW, ch);
  image(pg, aW, 0, cursor + 1, ch, 0, 0, cursor + 1, ch);

  drawOverlays(zoom);

  cursor = (cursor + 1) % cw;
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
    case 'g': visible = !visible; if (visible) gui.show(); else gui.hide(); break;
    case 'b': bpshow = !bpshow; break;   // population trace
    case 'l': blshow = !blshow; break;   // branch segments
    case 'z': zshow  = !zshow;  break;   // mean
    case 'v': vshow  = !vshow;  break;   // +/- std
    case 't': tshow  = !tshow;  break;   // optimum tick
    case 'r': seed();           break;   // reseed
  }
}
