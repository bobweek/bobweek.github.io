/* =====================================================================
   sketch.js  —  p5 glue only.  The model lives in field.js.

   Everything the user can touch is a real DOM control, so it keeps
   working inside an iframe, on a touchscreen, and with a keyboard.
   ===================================================================== */

let F;                 // the lattice
let img;               // p5.Image at lattice resolution, blown up to fill
let resIdx = 3;        // index into RESOLUTIONS
let playing = true;
let stepOnce = false;
let generation = 0;

const $ = (id) => document.getElementById(id);

const DEFAULTS = {
  m1: 0.23, m2: 0.43, m3: 0.30,
  s12: 0.02, s13: 0.064, s23: 0.03,
  ratio: 0.5, mode: 'chase',
  gamma: 0, sigma: 2.26, localNe: false,
  view: 'traits', bounds: 'clamp', spf: 1,
};

/* ------------------------------------------------------------- p5 setup */

function setup() {
  pixelDensity(1);
  const c = createCanvas(400, 700);
  c.parent('holder');
  c.elt.style.touchAction = 'manipulation';

  buildField(RESOLUTIONS[resIdx], null);
  F.randomize();

  wireControls();
  fitCanvas();

  new ResizeObserver(() => fitCanvas()).observe($('stage'));

  tuck(true);
}

function buildField(dims, old) {
  F = new Field(dims[0], dims[1]);
  if (old) F.resampleFrom(old); else F.randomize();
  img = createImage(F.wid, F.hi);
  img.loadPixels();
  $('resv').value = `${F.wid} \u00d7 ${F.hi}`;
}

/* -------------------------------------------------------------- params */

function readParams() {
  const mode = $('mode').value;
  return {
    m: [+$('m1').value, +$('m2').value, +$('m3').value],
    M: buildM(+$('s12').value, +$('s13').value, +$('s23').value, +$('ratio').value, mode),
    gamma: +$('gamma').value,
    sigma: +$('sigma').value,
    localNe: $('localNe').checked,
    mode: $('bounds').value,
  };
}

/* ---------------------------------------------------------------- draw */

function draw() {
  if (playing || stepOnce) {
    const p = readParams();
    const n = stepOnce ? 1 : +$('spf').value;
    for (let q = 0; q < n; q++) F.step(p);
    generation += n;
    stepOnce = false;
    $('gen').textContent = generation.toLocaleString();
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
  if (w === width && h === height) return;

  resizeCanvas(w, h);
  noSmooth();          // resizing the canvas resets the 2D context state
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
  const fmt = (id, d) => { $(id + 'v').value = (+$(id).value).toFixed(d); };
  fmt('m1', 2); fmt('m2', 2); fmt('m3', 2);
  fmt('s12', 3); fmt('s13', 3); fmt('s23', 3);
  fmt('ratio', 2); fmt('gamma', 3); fmt('sigma', 2);
  $('spfv').value = $('spf').value;
  // the flee:track ratio only means anything for the chase interaction
  $('ratioWrap').classList.toggle('off', $('mode').value !== 'chase');
}

function reseed(how) {
  if (how === 'random') F.randomize();
  else if (how === 'optima') F.seedOptima();
  else F.seedFlat();
  generation = 0;
  $('gen').textContent = '0';
}

function wireControls() {
  const sliders = ['m1', 'm2', 'm3', 's12', 's13', 's23', 'ratio', 'gamma', 'sigma', 'spf'];
  sliders.forEach((id) => $(id).addEventListener('input', showValues));
  $('mode').addEventListener('change', showValues);
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
