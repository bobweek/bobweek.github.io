// =====================================================================
//  Measure-valued branching diffusion  —  individual-based approximation
//  Trait-through-time "waterfall", rendered with an O(width) ring buffer.
//
//  Light theme: WHITE background, population drawn as accumulating GREY->BLACK ink.
//  Starts PAUSED (so it never runs in the background of a slide deck until asked).
//  Single file; uses p5.js 1.x global mode. No other libraries.
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
//  At n = 1 the picture is a sparse, clearly DISCRETE branching process: thick
//  lifelines (each individual lives at a fixed trait) with dashed birth
//  connectors.  As n rises the lifelines thin, the connectors fade, and the
//  empirical measure becomes the smooth diffusion cloud -- the approximation.
// =====================================================================

// ---------- canvas ----------
let cw = 800;
let ch = 500;

// ---------- live parameters (driven by the sliders) ----------
let Scaling_Factor = 1.0;
let Birth_Rate = 0.05;
let Death_Rate = 0.02;
let Abiotic_Selection = 1.5e-7;
let Competition = 0.03;
let mu = 5.0;
let Abiotic_Optimum = 5 * ch / 8;

// ---------- fixed params ----------
let N0 = 25.0;
let inds = [];

// ---------- ring-buffer rendering ----------
let pg, col, cursor = 0;
let colAcc = 0;

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

  pg = createGraphics(cw, ch);
  pg.pixelDensity(1);
  pg.noSmooth();
  pg.background(255);
  col = pg.drawingContext.createImageData(1, ch);

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
  const i = ui.length, x = 12, y = 12 + i * 38;
  const lab = createSpan('');
  lab.position(x, y);
  lab.style('color', '#111'); lab.style('font-family', 'monospace');
  lab.style('font-size', '11px'); lab.style('text-shadow', '0 0 2px #fff, 0 0 2px #fff');
  const s = createSlider(mn, mx, val, step);
  s.position(x, y + 15); s.style('width', '170px');

  const u = { slider: s, label: lab, name, fmt, set };
  const refresh = () => { const v = s.value(); set(v); lab.html(name + ' = ' + fmt(v)); };
  u.refresh = refresh;
  refresh();              // initialize label + binding now (so it shows while paused)
  s.input(refresh);       // live update on drag, even when paused

  if (host) { lab.parent(host); s.parent(host); }
  ui.push(u);
  return u;
}

function buildGUI() {
  addSlider('birth b0',      1e-3, 0.2,  Birth_Rate,        1e-3, v => v.toFixed(3),       v => Birth_Rate = v);
  addSlider('death d0',      1e-3, 0.2,  Death_Rate,        1e-3, v => v.toFixed(3),       v => Death_Rate = v);
  addSlider('selection',     0,    1e-6, Abiotic_Selection, 1e-8, v => v.toExponential(1), v => Abiotic_Selection = v);
  addSlider('competition',   1e-3, 0.06, Competition,       1e-3, v => v.toFixed(3),       v => Competition = v);
  addSlider('diffusion mu',  0.5,  12,   mu,                0.5,  v => v.toFixed(1),        v => mu = v);
  addSlider('optimum theta', 0,    ch,   Math.round(Abiotic_Optimum), 1, v => v.toFixed(0), v => Abiotic_Optimum = v);
  scalingUI = addSlider('scaling n', 1, 50, Scaling_Factor, 1, v => v.toFixed(0), v => Scaling_Factor = v);
  scalingSlider = scalingUI.slider;

  // step-n-by-one buttons beside the scaling slider
  const by = 12 + (ui.length - 1) * 38 + 15;
  const mk = (lbl, x, fn) => {
    const btn = createButton(lbl);
    btn.position(x, by);
    btn.style('font', '12px monospace'); btn.style('padding', '0 6px');
    btn.mousePressed(fn);
    if (host) btn.parent(host);
    btns.push(btn);
  };
  mk('\u25C0', 190, () => nudgeN(-1));   // decrease n
  mk('\u25B6', 222, () => nudgeN(+1));   // increase n

  // play / pause button, top-right of the window
  playBtn = createButton('\u25B6 play');
  playBtn.position(cw - 86, 12);
  playBtn.style('font', '12px monospace'); playBtn.style('padding', '2px 8px');
  playBtn.mousePressed(() => setPlaying(!playing));
  if (host) playBtn.parent(host);
  btns.push(playBtn);
}

function nudgeN(d) {
  let v = constrain(round(scalingSlider.value()) + d, 1, 50);
  scalingSlider.value(v);
  scalingUI.refresh();   // value() set programmatically doesn't fire input(); refresh manually
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
  let rootN = Math.sqrt(n);
  let sub = max(1, round(rootN));

  // smooth scroll: advance sqrt(n) columns/frame ON AVERAGE (fractional accumulator)
  colAcc += rootN;
  let columnsPerFrame = Math.floor(colAcc);
  colAcc -= columnsPerFrame;
  if (columnsPerFrame < 1) columnsPerFrame = 1;
  let dt = 1.0 / sub;

  // n-adaptive aesthetics: thick dark lifelines at n=1 -> fine faint cloud
  let radius = 1.3 / rootN;       // CONTINUOUS half-width with soft (anti-aliased) edges
  let addPt = 180 / rootN;        // lifeline ink (amount of darkening)
  let addBr = 150 / n;            // connector ink (fades fast in n)

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
    cursor = (cursor + 1) % cw;
  }

  background(255);
  let aW = cw - 1 - newest;
  if (aW > 0) image(pg, 0, 0, aW, ch, newest + 1, 0, aW, ch);
  image(pg, aW, 0, newest + 1, ch, 0, 0, newest + 1, ch);

  drawOverlays();
}

// ---------- column painters (subtractive: darken white toward black) ----------
function ink(data, y, amt) {
  let idx = y * 4;
  data[idx]   = Math.max(0, data[idx]   - amt);
  data[idx+1] = Math.max(0, data[idx+1] - amt);
  data[idx+2] = Math.max(0, data[idx+2] - amt);
}
function paintMark(data, yf, r, add) {              // soft lifeline; coverage fades at the edge
  let yc = Math.round(yf);
  let R = Math.ceil(r + 0.5);
  for (let dy = -R; dy <= R; dy++) {
    let y = yc + dy; if (y < 0 || y >= ch) continue;
    let cov = r + 0.5 - Math.abs(dy);
    if (cov > 0) ink(data, y, add * Math.min(1, cov));
  }
}
function paintDashed(data, y0f, y1f, add) {          // dashed vertical connector (2 on / 2 off)
  let y0 = y0f | 0, y1 = y1f | 0; if (y0 > y1) { let t = y0; y0 = y1; y1 = t; }
  y0 = Math.max(0, y0); y1 = Math.min(ch - 1, y1);
  for (let y = y0; y <= y1; y++) if (((y - y0) & 3) < 2) ink(data, y, add);
}

// ---------- overlays on the newest (right-edge) column ----------
function drawOverlays() {
  let xr = cw - 1;
  if (tshow) { stroke(220, 0, 0); strokeWeight(4); line(xr - 6, Abiotic_Optimum, xr - 2, Abiotic_Optimum); strokeWeight(1); }
  if (zshow && inds.length > 0) {
    let m = meanTrait();
    if (vshow) { stroke(120, 120, 120, 180); let sd = stdTrait(m); line(xr - 6, m + sd, xr - 6, m - sd); }
    stroke(0); line(xr - 8, m, xr - 2, m);
  }
}
function meanTrait() { let m = 0; for (let i = 0; i < inds.length; i++) m += inds[i].y; return m / inds.length; }
function stdTrait(m) { let v = 0; for (let i = 0; i < inds.length; i++) { let e = inds[i].y - m; v += e * e; } return Math.sqrt(v / inds.length); }

// ---------- keyboard ----------
function keyPressed() {
  if (key === ' ')             { setPlaying(!playing); return false; }   // play/pause
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
