// =====================================================================
//  Measure-valued branching diffusion  —  individual-based approximation
//  Trait-through-time "waterfall", rendered with an O(width) ring buffer.
//
//  SINGLE FILE: paste into the p5.js web editor and press play. No libraries.
// ---------------------------------------------------------------------
//  Scaling factor n refines the particle system toward its measure-valued
//  diffusion limit (Fournier-Meleard 2004; Champagnat-Ferriere-Meleard;
//  Etheridge's logistic super-Brownian motion; Dawson-Watanabe).  n enters
//  EVERYWHERE:
//    initial count  N0 * n
//    birth rate b(N) = sqrt(n) * Birth_Rate * exp(-Competition * N / n)
//    death rate d(y) = sqrt(n) * Death_Rate + Abiotic_Selection*(theta - y)^2
//    mutation    y' ~ Normal(parent.y, sqrt(mu / b))     (variance ~ 1/sqrt(n))
//  sigma^2 = mu/b keeps (rate)*(variance) = mu, so the trait -> Brownian motion
//  with diffusion coefficient mu in the limit.  Equilibrium N* = (n/c)ln(b/d0) ~ n.
//
//  At n = 1 the picture is a sparse, clearly DISCRETE branching process:
//  thick horizontal lifelines (each individual lives at a fixed trait) with
//  dashed vertical connectors at births (parent -> budded offspring).  As n
//  rises the lifelines thin, the connectors fade, and the empirical measure
//  becomes the smooth diffusion cloud -- that morph IS the diffusion approximation.
// =====================================================================

// ---------- canvas ----------
let cw = 800;
let ch = 500;

// ---------- live parameters (driven by the sliders) ----------
let Scaling_Factor = 1.0;     // n : refinement toward the diffusion limit
let Birth_Rate = 0.05;        // intrinsic per-capita birth rate (before sqrt(n) & competition)
let Death_Rate = 0.02;        // intrinsic baseline mortality (before sqrt(n))
let Abiotic_Selection = 1.5e-7; // stabilizing selection on death (kept O(1) in n)
let Competition = 0.03;       // competition strength (acts on birth, per-capita ~ 1/n)
let mu = 5.0;                 // trait diffusion coefficient (variance accrued / unit time)
let Abiotic_Optimum = 5 * ch / 8; // optimum trait value (trait units == pixels)

// ---------- fixed params ----------
let N0 = 25.0;         // initial abundance per unit of scaling
let inds = [];         // the population

// ---------- ring-buffer rendering ----------
let pg, col, cursor = 0;

// ---------- display toggles ----------
let bpshow = true;     // lifelines
let blshow = true;     // dashed birth connectors
let zshow  = false;    // mean trait marker
let vshow  = false;    // +/- std marker
let tshow  = false;    // optimum tick

// ---------- GUI ----------
let ui = [];
let guiVisible = true;

// =====================================================================
function setup() {
  pixelDensity(1);
  createCanvas(cw, ch);
  noSmooth();
  background(0);

  pg = createGraphics(cw, ch);
  pg.pixelDensity(1);
  pg.noSmooth();
  pg.background(0);
  col = pg.drawingContext.createImageData(1, ch);

  buildGUI();
  seed();
}

// ---------- control panel (native p5 sliders; setter closure writes the binding) ----------
function addSlider(name, mn, mx, val, step, fmt, set) {
  const i = ui.length, x = 12, y = 12 + i * 38;
  const lab = createSpan('');
  lab.position(x, y);
  lab.style('color', '#fff'); lab.style('font-family', 'monospace');
  lab.style('font-size', '11px'); lab.style('text-shadow', '0 0 3px #000');
  const s = createSlider(mn, mx, val, step);
  s.position(x, y + 15); s.style('width', '170px');
  ui.push({ slider: s, label: lab, name, fmt, set });
}

function buildGUI() {
  addSlider('birth b0',      1e-3, 0.2,  Birth_Rate,        1e-3, v => v.toFixed(3),       v => Birth_Rate = v);
  addSlider('death d0',      1e-3, 0.2,  Death_Rate,        1e-3, v => v.toFixed(3),       v => Death_Rate = v);
  addSlider('selection',     0,    1e-6, Abiotic_Selection, 1e-8, v => v.toExponential(1), v => Abiotic_Selection = v);
  addSlider('competition',   1e-3, 0.06, Competition,       1e-3, v => v.toFixed(3),       v => Competition = v);
  addSlider('diffusion mu',  0.5,  12,   mu,                0.5,  v => v.toFixed(1),        v => mu = v);
  addSlider('optimum theta', 0,    ch,   Math.round(Abiotic_Optimum), 1, v => v.toFixed(0), v => Abiotic_Optimum = v);
  addSlider('scaling n',     1,    50,   Scaling_Factor,    1,    v => v.toFixed(0),        v => Scaling_Factor = v);
}

function readGUI() {
  for (const u of ui) { const v = u.slider.value(); u.set(v); u.label.html(u.name + ' = ' + u.fmt(v)); }
}

function seed() {
  inds = [];
  let N = floor(N0 * Scaling_Factor);
  for (let i = 0; i < N; i++) inds.push(new Individual(Abiotic_Optimum));
  cursor = 0;
}

// ---------- n-dependent rates ----------
function birthRate(N, n) { return Math.sqrt(n) * Birth_Rate * Math.exp(-Competition * N / n); }
function deathRate(y, n) { let e = Abiotic_Optimum - y; return Math.sqrt(n) * Death_Rate + Abiotic_Selection * e * e; }

class Individual { constructor(y) { this.y = y; this.dead = false; } }

// =====================================================================
function draw() {
  readGUI();

  let n = Scaling_Factor;
  let sub = max(1, round(Math.sqrt(n)));      // sub-steps keep event probs bounded
  let dt = 1.0 / sub;                          // one column == one unit of sim-time
  let columnsPerFrame = max(1, round(Math.sqrt(n)));

  // --- n-adaptive aesthetics: thick bright lifelines at n=1 -> fine faint cloud ---
  let thickness = max(1, round(3 / Math.sqrt(n)));   // lifeline half-height below
  let half = thickness >> 1;
  let addPt = 180 / Math.sqrt(n);                    // lifeline brightness
  let addBr = 150 / n;                               // connector brightness (fades fast in n)

  let newest = cursor;
  let data = col.data;

  for (let c = 0; c < columnsPerFrame; c++) {
    for (let k = 0; k < data.length; k += 4) { data[k]=0; data[k+1]=0; data[k+2]=0; data[k+3]=255; }

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
          if (blshow && addBr > 1.5) paintDashed(data, ind.y, yo, addBr); // dashed parent->offspring
        }
      }

      let alive = [];
      for (let i = 0; i < N; i++) if (!inds[i].dead) alive.push(inds[i]);
      for (let i = 0; i < born.length; i++) alive.push(new Individual(born[i]));
      inds = alive;
      if (inds.length === 0) inds.push(new Individual(Abiotic_Optimum));

      if (isLast && bpshow) {
        for (let i = 0; i < inds.length; i++) paintMark(data, inds[i].y, half, addPt);
      }
    }

    pg.drawingContext.putImageData(col, cursor, 0);
    newest = cursor;
    cursor = (cursor + 1) % cw;
  }

  background(0);
  let aW = cw - 1 - newest;
  if (aW > 0) image(pg, 0, 0, aW, ch, newest + 1, 0, aW, ch);
  image(pg, aW, 0, newest + 1, ch, 0, 0, newest + 1, ch);

  drawOverlays();
}

// ---------- column painters ----------
function addRGB(data, y, add) {
  let idx = y * 4;
  data[idx]   = Math.min(255, data[idx]   + add);
  data[idx+1] = Math.min(255, data[idx+1] + add);
  data[idx+2] = Math.min(255, data[idx+2] + add);
}
function paintMark(data, yf, half, add) {           // thick lifeline mark, height 2*half+1
  let yc = yf | 0;
  for (let dy = -half; dy <= half; dy++) { let y = yc + dy; if (y >= 0 && y < ch) addRGB(data, y, add); }
}
function paintDashed(data, y0f, y1f, add) {          // dashed vertical connector (2 on / 2 off)
  let y0 = y0f | 0, y1 = y1f | 0; if (y0 > y1) { let t = y0; y0 = y1; y1 = t; }
  y0 = Math.max(0, y0); y1 = Math.min(ch - 1, y1);
  for (let y = y0; y <= y1; y++) if (((y - y0) & 3) < 2) addRGB(data, y, add);
}

// ---------- overlays on the newest (right-edge) column ----------
function drawOverlays() {
  let xr = cw - 1;
  if (tshow) { stroke(255, 0, 0); strokeWeight(4); line(xr - 6, Abiotic_Optimum, xr - 2, Abiotic_Optimum); strokeWeight(1); }
  if (zshow && inds.length > 0) {
    let m = meanTrait();
    if (vshow) { let sd = stdTrait(m); stroke(120,120,120,180); line(xr - 6, m + sd, xr - 6, m - sd); }
    stroke(255); line(xr - 8, m, xr - 2, m);
  }
}
function meanTrait() { let m = 0; for (let i = 0; i < inds.length; i++) m += inds[i].y; return m / inds.length; }
function stdTrait(m) { let v = 0; for (let i = 0; i < inds.length; i++) { let e = inds[i].y - m; v += e * e; } return Math.sqrt(v / inds.length); }

// ---------- keyboard ----------
function keyPressed() {
  switch (key) {
    case 'g': guiVisible = !guiVisible;
      for (const u of ui) { if (guiVisible) { u.slider.show(); u.label.show(); } else { u.slider.hide(); u.label.hide(); } } break;
    case 'b': bpshow = !bpshow; break;   // lifelines
    case 'l': blshow = !blshow; break;   // connectors
    case 'z': zshow  = !zshow;  break;   // mean
    case 'v': vshow  = !vshow;  break;   // +/- std
    case 't': tshow  = !tshow;  break;   // optimum tick
    case 'r': seed();           break;   // reseed at current scale
  }
}
