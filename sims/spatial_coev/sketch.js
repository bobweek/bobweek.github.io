/* =====================================================================
   sketch.js  —  p5 glue only.  The model lives in field.js.

   Everything the user can touch is a real DOM control, so it keeps
   working inside an iframe, on a touchscreen, and with a keyboard.
   ===================================================================== */

let F;                 // the lattice
let img;               // p5.Image at lattice resolution, blown up to fill
let resIdx;            // index into RESOLUTIONS, read from the markup
let playing = true;

/* Set this to true to start with the control panel tucked away. The same
   thing can be done per-embed with ?panel=hidden on the iframe src. */
const PANEL_HIDDEN_BY_DEFAULT = false;
let stepOnce = false;
let generation = 0;
let rate = 0, rateAt = 0, rateShown = 0;
let genBudget = 0;     // carries the fractional part of generations/frame

const $ = (id) => document.getElementById(id);

/* Defaults live in ONE place: the value= / selected / checked attributes on
   the controls in index.html.  They are snapshotted here at startup so the
   "Restore defaults" button has something to go back to.  Add a control to
   the markup and it is picked up automatically. */
const STATEFUL = '#panel input[type=range], #panel input[type=checkbox], #panel select';
let DEFAULTS = {};

function captureDefaults() {
  DEFAULTS = {};
  document.querySelectorAll(STATEFUL).forEach((el) => {
    if (el.id) DEFAULTS[el.id] = el.type === 'checkbox' ? el.checked : el.value;
  });
}

/* A slider may declare data-map="sq", in which case the value it feeds the
   model is the square of its position.  That buys fine control at the low end
   of a wide range: drift runs 0..25 but the interesting part is under 2. */
function mapped(el) {
  const v = +el.value;
  return el.dataset.map === 'sq' ? v * v : v;
}

// how many decimals to show for a slider, taken from its own step attribute
function decimals(step) {
  const s = String(step);
  const dot = s.indexOf('.');
  return dot < 0 ? 0 : s.length - dot - 1;
}

/* ------------------------------------------------------------- p5 setup */

function setup() {
  pixelDensity(1);
  const c = createCanvas(400, 700);
  c.parent('holder');
  c.elt.style.touchAction = 'manipulation';

  // the starting grid is a data attribute on the stepper, so it lives in the
  // markup with every other default
  resIdx = constrain(+document.querySelector('.stepper').dataset.level || 0,
                     0, RESOLUTIONS.length - 1);
  buildField(RESOLUTIONS[resIdx], null);
  F.randomize();

  wireControls();
  fitCanvas();

  new ResizeObserver(() => fitCanvas()).observe($('stage'));

  const q = new URLSearchParams(location.search);
  if (PANEL_HIDDEN_BY_DEFAULT || q.get('panel') === 'hidden' || q.get('panel') === '0') tuck(true);

  // the simulation *is* the motion, so honour the system setting
  if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) setPlaying(false);
}

function buildField(dims, old) {
  F = new Field(dims[0], dims[1]);
  if (old) F.resampleFrom(old); else F.randomize();
  img = createImage(F.wid, F.hi);
  img.loadPixels();
  $('resv').value = `${F.wid} \u00d7 ${F.hi}`;
  rate = 0; rateAt = 0;
}

/* -------------------------------------------------------------- params */

const COEF = ['t12', 'f12', 't13', 'f13', 't23', 'f23'];

/* The 2016 sketch's slider positions, worked back through its own
   normalisation.  "Even cycle" is the same structure with one asymmetry
   throughout; match and avoid flip the sign of the push. */
const PRESETS = {
  original: { t12: 0.095, f12: 0.095, t13: 0.095, f13: 0.19, t23: 0, f23: 0 },
  cycle:    { t12: 0.02,  f12: 0.04,  t13: 0.02,  f13: 0.04, t23: 0.02, f23: 0.04 },
  match:    { t12: 0.05,  f12: -0.05, t13: 0.05,  f13: -0.05, t23: 0.05, f23: -0.05 },
  avoid:    { t12: -0.05, f12: 0.05,  t13: -0.05, f13: 0.05,  t23: -0.05, f23: 0.05 },
};

function readParams() {
  const c = {};
  COEF.forEach((id) => { c[id] = +$(id).value; });
  return {
    m: [+$('m1').value, +$('m2').value, +$('m3').value],
    M: buildM(c),
    gamma: +$('gamma').value,
    sigma: mapped($('sigma')),
    localNe: $('localNe').checked,
    mode: $('bounds').value,
    sequential: $('order').value === 'sequential',
  };
}

/* ---------------------------------------------------------------- draw */

function draw() {
  if (playing || stepOnce) {
    const p = readParams();

    /* Generations per frame can be fractional. The budget carries the
       remainder between frames, so 0.75 steps on three frames out of four
       rather than rounding up to one step every frame. That is the only way
       to run slower than the display refresh. */
    let n;
    if (stepOnce) {
      n = 1;
      stepOnce = false;
    } else {
      genBudget += +$('spf').value;
      n = Math.floor(genBudget);
      genBudget -= n;
    }

    for (let q = 0; q < n; q++) F.step(p);
    if (n) {
      generation += n;
      $('gen').textContent = generation.toLocaleString();
    }

    /* Smoothed generations per second. This is averaged over every frame,
       including the ones that take no step, or skipped frames would bias it
       upward at fractional speeds. */
    const now = performance.now();
    if (rateAt) {
      rate = rate * 0.93 + ((n * 1000) / Math.max(1, now - rateAt)) * 0.07;
      if (now - rateShown > 400) {
        $('rate').textContent = rate < 10 ? rate.toFixed(1) : rate.toFixed(0);
        rateShown = now;
      }
    }
    rateAt = now;
  } else {
    rateAt = 0;
    genBudget = 0;
    $('rate').textContent = '';
  }

  F.writePixels(img.pixels, $('view').value);
  img.updatePixels();
  image(img, 0, 0, width, height);
}

/* ------------------------------------------------------------- sizing */

function fitCanvas() {
  const stage = $('stage');
  const cs = getComputedStyle(stage);
  const availW = stage.clientWidth - parseFloat(cs.paddingLeft) - parseFloat(cs.paddingRight);
  const availH = stage.clientHeight - parseFloat(cs.paddingTop) - parseFloat(cs.paddingBottom);
  if (availW < 20 || availH < 20) return;

  const aspect = F.wid / F.hi;
  let w = availW, h = availW / aspect;
  if (h > availH) { h = availH; w = availH * aspect; }
  w = Math.max(20, Math.round(w));
  h = Math.max(20, Math.round(h));
  // resizing the canvas resets the 2D context state, and the lattice can
  // change without the canvas changing, so smoothing is set either way
  if (w !== width || h !== height) resizeCanvas(w, h);
  setSmoothing();
}

/* A lattice coarser than the canvas is magnified and should stay as crisp
   blocks.  A lattice finer than the canvas has to be shrunk, and sampling it
   nearest-neighbour would alias the fine structure into noise, so hand the
   downsampling to the browser instead. */
function setSmoothing() {
  if (F.wid > width) smooth(); else noSmooth();
}

function windowResized() { fitCanvas(); }

/* ------------------------------------------------------------ controls */

function setPlaying(on) {
  playing = on;
  $('play').textContent = on ? 'Pause' : 'Play';
}

function setResolution(idx) {
  const next = constrain(idx, 0, RESOLUTIONS.length - 1);
  if (next === resIdx) return;
  resIdx = next;
  buildField(RESOLUTIONS[resIdx], F);
  fitCanvas();
}

function showValues() {
  document.querySelectorAll('#panel input[type=range]').forEach((el) => {
    const out = $(el.id + 'v');
    if (out) out.value = mapped(el).toFixed(el.dataset.dp ? +el.dataset.dp : decimals(el.step));
  });
}

// mark the preset select as Custom once the coefficients no longer match it
function syncPreset() {
  const now = {};
  COEF.forEach((id) => { now[id] = +$(id).value; });
  const hit = Object.keys(PRESETS).find((name) =>
    COEF.every((id) => Math.abs(PRESETS[name][id] - now[id]) < 1e-9));
  $('preset').value = hit || 'custom';
}

function applyPreset(name) {
  const p = PRESETS[name];
  if (!p) return;
  COEF.forEach((id) => { $(id).value = p[id]; });
  showValues();
}

function reseed(how) {
  if (how === 'random') F.randomize();
  else if (how === 'optima') F.seedOptima();
  else F.seedFlat();
  generation = 0;
  $('gen').textContent = '0';
}

function wireControls() {
  captureDefaults();
  document.querySelectorAll('#panel input[type=range]')
    .forEach((el) => el.addEventListener('input', showValues));
  COEF.forEach((id) => $(id).addEventListener('input', syncPreset));
  $('preset').addEventListener('change', () => applyPreset($('preset').value));
  showValues();

  $('play').onclick = () => setPlaying(!playing);
  $('stepBtn').onclick = () => { setPlaying(false); stepOnce = true; };

  $('rand').onclick = () => reseed('random');
  $('seedOpt').onclick = () => reseed('optima');
  $('flat').onclick = () => reseed('flat');

  $('coarser').onclick = () => setResolution(resIdx - 1);
  $('finer').onclick = () => setResolution(resIdx + 1);

  $('defaults').onclick = () => {
    for (const [id, v] of Object.entries(DEFAULTS)) {
      const el = $(id);
      if (el.type === 'checkbox') el.checked = v; else el.value = v;
    }
    showValues();
  };

  $('hide').onclick = () => tuck(true);
  $('reveal').onclick = () => tuck(false);

  // pausing belongs to the canvas, not to every click on the page
  document.querySelector('#holder canvas')
    .addEventListener('pointerdown', () => setPlaying(!playing));
}

function tuck(on) {
  document.body.classList.toggle('tucked', on);
  $('reveal').hidden = !on;
  requestAnimationFrame(fitCanvas);
}

/* ----------------------------------------------------------- keyboard */

window.addEventListener('keydown', (e) => {
  // let arrow keys nudge a focused slider instead of hijacking them
  const a = document.activeElement;
  if (a && /^(INPUT|SELECT|TEXTAREA|BUTTON)$/.test(a.tagName) && e.key !== 'Escape') {
    if (e.key !== ' ') return;
    if (a.tagName === 'BUTTON') return;
  }
  switch (e.key) {
    case ' ': setPlaying(!playing); e.preventDefault(); break;
    case 's': setPlaying(false); stepOnce = true; break;
    case 'r': reseed('random'); break;
    case 'h': tuck(!document.body.classList.contains('tucked')); break;
    case '[': setResolution(resIdx - 1); break;
    case ']': setResolution(resIdx + 1); break;
  }
});
