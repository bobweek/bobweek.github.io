// =====================================================================
//  Drift builds genetic correlations  —  Wright-Fisher, no recombination
//  A population of haploid sequences on the left; the dynamics of every
//  pairwise locus correlation on the right.
//
//  SINGLE FILE: paste into the p5.js web editor and press play. No libraries.
// ---------------------------------------------------------------------
//  Model (haploid, discrete non-overlapping Wright-Fisher, N constant,
//  NO mutation, NO recombination):
//    - genome = L biallelic loci, ALL polymorphic.  Each allele is drawn as a
//      bright / dark shade of that locus' hue.
//    - start in linkage equilibrium: every locus is an independent
//      Bernoulli(1/2) across individuals, so all pairwise correlations ~ 0.
//    - every generation the next population is N whole-haplotype copies of
//      randomly chosen parents (clonal resampling).  With no recombination the
//      entire genome rides ONE shared genealogy, so as lineages coalesce the
//      loci drift into correlation purely from finite-population sampling --
//      the build-up of linkage disequilibrium in the fully linked limit.
//
//  The right panel tracks, for every locus pair (i,j), the signed allelic
//  correlation
//        r_ij = D_ij / sqrt( p_i(1-p_i) p_j(1-p_j) ),   D_ij = p_ij - p_i p_j,
//  the usual LD correlation coefficient (r^2 is its square).  There are
//  L-choose-2 of them, so the plot is a fan of overlapping time series.  A
//  line breaks when one of its loci fixes (r undefined).  Note r does not sit
//  at exactly 0 even in equilibrium: the finite-population baseline is
//  E[r^2] ~ 1/(N-1), i.e. E|r| ~ sqrt(2/pi)/sqrt(N-1).
//
//  This is transient: with no mutation the only true endpoint is fixation on a
//  single haplotype (every locus monomorphic).  The interesting regime is the
//  intermediate one, where the population sits on ~2 complementary haplotypes
//  and the surviving polymorphic loci are perfectly correlated.
//
//  Controls (left): population N, number of loci L, speed; play / step /
//  reset (same seed) / new seed.  Click a locus column to spotlight every pair
//  that involves it.  Keys: space play/pause, s step, r reset, n new seed.
// =====================================================================

// ---------- canvas ----------
let cw = 980;
let ch = 600;

// ---------- live parameters (driven by the sliders) ----------
let N = 120;          // population size
let L = 10;           // number of (all polymorphic) loci
let speed = 12;       // generations per second

// ---------- Rosé Pine Dawn palette (light) ----------
const COL = {
  base:"#faf4ed", surface:"#fffaf3", overlay:"#f2e9e1", hlMed:"#dfdad9",
  muted:"#9893a5", subtle:"#797593", text:"#575279"
};
// accent hues: love, gold, rose, pine, foam, iris
const ACCENTS = ["#b4637a","#ea9d34","#d7827e","#286983","#56949f","#907aa9"];
let shadeDeep = [], shadeBright = [], pairList = [], pairHue = [];

// ---------- model state ----------
let pop = [];          // pop[i][loc] in {0,1}
let gen = 0;
let playing = true;
let seedValue = 20250609;
let rngState = 1;

let genHist = [];      // recorded generation indices
let rHist = [];        // rHist[k] = array of r (or null) for pair k
const HIST_MAX = 1400;
let stepAcc = 0;

let highlightLocus = -1;

// ---------- layout (computed in setup) ----------
let sideX = 12, sideW = 196;
let popX0, popW, plotX0, plotW, vizTop, vizBot;

// ---------- GUI ----------
let ui = [], btns = [];

// =====================================================================
function setup() {
  pixelDensity(1);
  const cnv = createCanvas(cw, ch);
  cnv.style('display', 'block');
  textFont('monospace');

  // precompute shade pairs for the maximum number of loci we might show
  rebuildPalette();

  // viz geometry
  popX0  = sideX + sideW + 26;
  popW   = 300;
  plotX0 = popX0 + popW + 40;
  plotW  = cw - plotX0 - 18;
  vizTop = 74;
  vizBot = ch - 30;

  buildGUI();
  initPop();
}

function rebuildPalette() {
  shadeDeep = []; shadeBright = [];
  for (let k = 0; k < ACCENTS.length; k++) {
    const h = color(ACCENTS[k]);
    shadeDeep.push(lerpColor(h, color(COL.text), 0.18));
    shadeBright.push(lerpColor(h, color(255), 0.62));
  }
}
function locusDeep(loc)   { return shadeDeep[loc % ACCENTS.length]; }
function locusBright(loc) { return shadeBright[loc % ACCENTS.length]; }
function locusHue(loc)    { return ACCENTS[loc % ACCENTS.length]; }

// ---------- seeded RNG (mulberry32) ----------
function srand(s){ rngState = s >>> 0; }
function rnd(){
  rngState |= 0; rngState = (rngState + 0x6D2B79F5) | 0;
  let t = Math.imul(rngState ^ (rngState >>> 15), 1 | rngState);
  t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
  return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
}
const ri = n => Math.floor(rnd() * n);

// ---------- model ----------
function buildPairs(){
  pairList = []; pairHue = [];
  for (let i = 0; i < L; i++)
    for (let j = i + 1; j < L; j++){ pairList.push([i, j]); pairHue.push(locusHue(i)); }
}

function initPop(){
  srand(seedValue);
  buildPairs();
  pop = [];
  for (let i = 0; i < N; i++){
    const row = new Uint8Array(L);
    for (let l = 0; l < L; l++) row[l] = rnd() < 0.5 ? 1 : 0;   // linkage equilibrium
    pop.push(row);
  }
  gen = 0; stepAcc = 0; highlightLocus = -1;
  genHist = []; rHist = pairList.map(() => []);
  record();
}

function stepGen(){
  const next = [];
  for (let i = 0; i < N; i++) next.push(pop[ri(N)].slice());  // clonal WF copy
  pop = next;
  sortByHaplotype();
  gen++;
  record();
}

function sortByHaplotype(){
  pop.sort((a, b) => { for (let l = 0; l < L; l++){ if (a[l] !== b[l]) return a[l] - b[l]; } return 0; });
}

// signed allelic correlation between loci i,j; null if either is fixed
function corr(i, j){
  let cA = 0, cB = 0, cAB = 0;
  for (let k = 0; k < N; k++){ const a = pop[k][i], b = pop[k][j]; cA += a; cB += b; cAB += a & b; }
  const pa = cA / N, pb = cB / N;
  if (pa === 0 || pa === 1 || pb === 0 || pb === 1) return null;
  const D = cAB / N - pa * pb;
  return D / Math.sqrt(pa * (1 - pa) * pb * (1 - pb));
}

function record(){
  genHist.push(gen);
  for (let k = 0; k < pairList.length; k++)
    rHist[k].push(corr(pairList[k][0], pairList[k][1]));
  if (genHist.length > HIST_MAX){          // rebin (drop every other point)
    genHist = genHist.filter((_, idx) => idx % 2 === 0);
    rHist = rHist.map(arr => arr.filter((_, idx) => idx % 2 === 0));
  }
}

function nHaplotypes(){
  const set = new Set();
  for (let k = 0; k < N; k++) set.add(pop[k].join(''));
  return set.size;
}

// =====================================================================
function draw(){
  readGUI();
  background(COL.base);

  // advance generations on a wall-clock accumulator (speed independent of fps)
  if (playing){
    stepAcc += (deltaTime / 1000) * speed;
    let budget = 0;
    while (stepAcc >= 1 && budget < 40){ stepGen(); stepAcc -= 1; budget++; }
  }

  drawSidebar();
  drawHeader();
  drawPopulation();
  drawCorrelations();
}

// ---------- header ----------
function drawHeader(){
  noStroke(); textAlign(LEFT, BASELINE);
  fill(COL.text); textSize(20); textStyle(BOLD);
  text("drift builds genetic correlations", popX0, 34);
  textStyle(NORMAL); fill(COL.subtle); textSize(11.5);
  text("Wright-Fisher resampling of whole haplotypes  ·  no mutation, no recombination",
       popX0, 52);

  fill(COL.text); textSize(12.5); textStyle(BOLD); textAlign(LEFT, BASELINE);
  text("POPULATION", popX0, vizTop - 8);
  textStyle(NORMAL); fill(COL.subtle); textSize(11);
  textAlign(RIGHT, BASELINE);
  text("gen " + gen, popX0 + popW, vizTop - 8);

  fill(COL.text); textSize(12.5); textStyle(BOLD); textAlign(LEFT, BASELINE);
  text("CORRELATIONS", plotX0, vizTop - 8);
  textStyle(NORMAL); fill(COL.subtle); textSize(11);
  text("one line per locus pair", plotX0 + 104, vizTop - 8);
}

// ---------- population grid ----------
function drawPopulation(){
  const H = vizBot - vizTop;
  const rowH = H / N;
  const cellW = popW / L;
  const segIdxHi = highlightLocus;

  // surface backing
  noStroke(); fill(COL.surface);
  rect(popX0, vizTop, popW, H, 4);

  noStroke();
  for (let r = 0; r < N; r++){
    const by = vizTop + r * rowH;
    for (let l = 0; l < L; l++){
      const allele = pop[r][l];
      const c = allele ? locusBright(l) : locusDeep(l);
      fill(c);
      rect(popX0 + l * cellW, by, cellW, Math.ceil(rowH) + 0.5);
    }
  }
  // column gutters (base colour) -> norns-style segmentation
  stroke(COL.base); strokeWeight(1.4);
  for (let l = 0; l <= L; l++){ const x = popX0 + l * cellW; line(x, vizTop, x, vizBot); }

  // spotlight a clicked locus column
  if (segIdxHi >= 0 && segIdxHi < L){
    noFill(); stroke(COL.text); strokeWeight(1.6);
    rect(popX0 + segIdxHi * cellW, vizTop, cellW, H);
    // caret above
    const cx = popX0 + (segIdxHi + 0.5) * cellW;
    line(cx - cellW * 0.18, vizTop - 5, cx, vizTop - 11);
    line(cx, vizTop - 11, cx + cellW * 0.18, vizTop - 5);
  }
  noStroke();
}

// ---------- correlation dynamics ----------
function drawCorrelations(){
  const H = vizBot - vizTop;
  const x0 = plotX0, y0 = vizTop, w = plotW, h = H;
  const yMid = y0 + h / 2;

  noStroke(); fill(COL.surface); rect(x0, y0, w, h, 4);

  // guides: r = +1 / 0 / -1
  stroke(COL.hlMed); strokeWeight(1);
  drawingContext.setLineDash([2, 4]);
  line(x0, y0, x0 + w, y0); line(x0, y0 + h, x0 + w, y0 + h);
  drawingContext.setLineDash([]);
  stroke(COL.muted); strokeWeight(1); line(x0, yMid, x0 + w, yMid);
  noStroke(); fill(COL.subtle); textSize(9.5); textAlign(RIGHT, CENTER);
  text("+1", x0 - 4, y0); text("0", x0 - 4, yMid); text("−1", x0 - 4, y0 + h);
  textAlign(CENTER, TOP); text("generation", x0 + w / 2, y0 + h + 8);
  push(); translate(x0 - 26, yMid); rotate(-HALF_PI);
  textAlign(CENTER, BOTTOM); text("allelic correlation  r", 0, 0); pop();

  const n = genHist.length;
  if (n >= 2){
    const gmax = Math.max(1, genHist[n - 1]);
    const xx = g => x0 + (g / gmax) * w;
    const yy = v => yMid - v * (h / 2);

    for (let k = 0; k < pairList.length; k++){
      const involved = highlightLocus >= 0 &&
                       (pairList[k][0] === highlightLocus || pairList[k][1] === highlightLocus);
      if (highlightLocus >= 0 && !involved){
        stroke(red(color(pairHue[k])), green(color(pairHue[k])), blue(color(pairHue[k])), 28);
        strokeWeight(1);
      } else {
        const a = highlightLocus >= 0 ? 230 : 130;
        const c = color(pairHue[k]);
        stroke(red(c), green(c), blue(c), a);
        strokeWeight(involved ? 2.0 : 1.3);
      }
      noFill();
      const arr = rHist[k];
      let open = false;
      for (let t = 0; t < n; t++){
        const v = arr[t];
        if (v === null || v === undefined){ if (open){ endShape(); open = false; } continue; }
        if (!open){ beginShape(); open = true; }
        vertex(xx(genHist[t]), yy(v));
      }
      if (open) endShape();
    }
  }

  // readouts
  let mAbs = 0, cnt = 0, committed = 0;
  for (let k = 0; k < pairList.length; k++){
    const v = rHist[k][rHist[k].length - 1];
    if (v !== null && v !== undefined){ mAbs += Math.abs(v); cnt++; if (Math.abs(v) >= 0.9) committed++; }
  }
  mAbs = cnt ? mAbs / cnt : 0;
  noStroke(); textAlign(LEFT, BASELINE); fill(COL.text); textSize(11);
  text("mean |r| = " + nf(mAbs, 1, 2), x0, y0 - 22 + 0);   // (kept inside header band)

  const nh = nHaplotypes();
  textAlign(RIGHT, BASELINE); fill(COL.subtle); textSize(10.5);
  let note = committed + " pairs at |r| ≥ 0.9";
  if (nh === 1) note = "fixed — one haplotype · press r / n";
  text(note, x0 + w, y0 - 22 + 0);

  // tiny haplotype-count tag, bottom-right inside plot
  textAlign(RIGHT, BASELINE); fill(COL.muted); textSize(10);
  text(nh + " haplotype" + (nh === 1 ? "" : "s"), x0 + w - 6, y0 + h - 6);
}

// ---------- sidebar panel behind the DOM controls ----------
function drawSidebar(){
  noStroke(); fill(COL.surface);
  rect(sideX - 6, 6, sideW + 16, ch - 12, 8);
  stroke(COL.hlMed); strokeWeight(1); noFill();
  rect(sideX - 6, 6, sideW + 16, ch - 12, 8);
  noStroke();
}

// =====================================================================
// GUI : native p5 sliders + buttons stacked in the left sidebar
// =====================================================================
function addSlider(name, mn, mx, val, step, fmt, set){
  const i = ui.length, x = sideX, y = 18 + i * 46;
  const lab = createSpan('');
  lab.position(x, y);
  lab.style('color', COL.text); lab.style('font-family', 'monospace');
  lab.style('font-size', '11px');
  const s = createSlider(mn, mx, val, step);
  s.position(x, y + 16); s.style('width', sideW + 'px');
  ui.push({ slider: s, label: lab, name, fmt, set, last: val });
}

function buildGUI(){
  addSlider('population N', 20, 240, N, 2, v => v.toFixed(0),
            v => { if (v !== N){ N = v; initPop(); } });
  addSlider('loci L',        4,  14, L, 1, v => v.toFixed(0),
            v => { if (v !== L){ L = v; initPop(); } });
  addSlider('speed',         1,  30, speed, 1, v => v.toFixed(0) + ' gen/s',
            v => { speed = v; });

  const by = 18 + ui.length * 46 + 4;
  const mk = (lbl, x, w, fn) => {
    const b = createButton(lbl);
    b.position(sideX + x, by);
    b.style('font', '12px monospace'); b.style('padding', '5px 0');
    b.style('width', w + 'px'); b.style('cursor', 'pointer');
    b.style('color', COL.text); b.style('background', COL.base);
    b.style('border', '1px solid ' + COL.hlMed); b.style('border-radius', '6px');
    b.mousePressed(fn); btns.push(b);
  };
  mk('play/pause', 0,   sideW, togglePlay);
  const by2 = by + 30;
  const mk2 = (lbl, x, w, fn) => {
    const b = createButton(lbl);
    b.position(sideX + x, by2);
    b.style('font', '12px monospace'); b.style('padding', '5px 0');
    b.style('width', w + 'px'); b.style('cursor', 'pointer');
    b.style('color', COL.text); b.style('background', COL.base);
    b.style('border', '1px solid ' + COL.hlMed); b.style('border-radius', '6px');
    b.mousePressed(fn); btns.push(b);
  };
  const third = Math.floor((sideW - 12) / 3);
  mk2('step',  0,             third, () => { if (!playing) stepGen(); });
  mk2('reset', third + 6,     third, () => initPop());
  mk2('seed',  2 * third + 12, third, newSeed);

  // hint line under the buttons
  const hint = createSpan('space play · s step · r reset · n seed<br>click a locus to spotlight its pairs');
  hint.position(sideX, by2 + 34);
  hint.style('color', COL.subtle); hint.style('font-family', 'monospace');
  hint.style('font-size', '10px'); hint.style('width', sideW + 'px');
  hint.style('line-height', '1.5');
}

function readGUI(){
  for (const u of ui){
    const v = u.slider.value();
    u.label.html(u.name + ' = ' + u.fmt(v));
    if (v !== u.last){ u.last = v; u.set(v); }
  }
}

function togglePlay(){ playing = !playing; }
function newSeed(){ seedValue = (Math.random() * 4294967296) >>> 0; initPop(); }

// ---------- interaction ----------
function mousePressed(){
  if (mouseX < popX0 || mouseX > popX0 + popW || mouseY < vizTop || mouseY > vizBot) return;
  const cellW = popW / L;
  const loc = Math.floor((mouseX - popX0) / cellW);
  highlightLocus = (loc === highlightLocus) ? -1 : loc;   // toggle
}

function keyPressed(){
  if (key === ' ') { togglePlay(); return false; }
  if (key === 's' || key === 'S') { if (!playing) stepGen(); }
  if (key === 'r' || key === 'R') initPop();
  if (key === 'n' || key === 'N') newSeed();
}
