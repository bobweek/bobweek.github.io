const { Field, buildM, SQUASH, ZMAX } = require('./field.js');

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const ZERO = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
let fails = 0;
function ok(name, cond, extra = '') {
  console.log(`${cond ? 'PASS' : 'FAIL'}  ${name}${extra ? '   ' + extra : ''}`);
  if (!cond) fails++;
}

// ---------------------------------------------------------------- 1. squash
{
  const big = [-1e4, -300, -1, 0, 128, 255, 256, 900, 1e4];
  let inRange = true;
  for (const name of ['clamp', 'wrap', 'fold', 'soft']) {
    for (const v of big) {
      const w = SQUASH[name](v);
      // wrap lives on the circle [0,256); the others on [0,255]
      const hi = name === 'wrap' ? 256 : ZMAX + 1e-9;
      if (!(w >= 0 && w < hi && Number.isFinite(w))) { inRange = false; console.log('  bad', name, v, w); }
    }
  }
  ok('every squash maps R -> [0,255]', inRange);
  // original wrap:  (255 + v) % 255
  const oldWrap = (v) => (255 + v) % 255;
  ok('original wrap leaks negative for v=-300', oldWrap(-300) < 0, `old=${oldWrap(-300)}  new=${SQUASH.wrap(-300).toFixed(1)}`);
  const oldFold = (v) => (v > 255 ? 510 - v : v < 0 ? -v : v);
  ok('original fold leaks for v=-600', oldFold(-600) > 255, `old=${oldFold(-600)}  new=${SQUASH.fold(-600).toFixed(1)}`);
  const eps = 1e-3, mid = ZMAX / 2;
  const slope = (SQUASH.soft(mid + eps) - SQUASH.soft(mid - eps)) / (2 * eps);
  ok('soft squash has unit slope at mid-range', Math.abs(slope - 1) < 1e-4, `slope=${slope.toFixed(6)}`);
}

// ------------------------------------------------- 2. dispersal only: mixing
{
  const F = new Field(40, 70, mulberry32(1));
  F.randomize();
  const sd = () => {
    let s = 0, ss = 0;
    for (let k = 0; k < F.cur.length; k += 3) { s += F.cur[k]; ss += F.cur[k] * F.cur[k]; }
    const n = F.n; const mu = s / n; return { mu, sd: Math.sqrt(ss / n - mu * mu) };
  };
  const a = sd();
  const p = { m: [0.5, 0.5, 0.5], M: ZERO, gamma: 0, sigma: 0, mode: 'clamp', localNe: false };
  let mono = true, prev = a.sd;
  for (let t = 0; t < 400; t++) { F.step(p); const c = sd().sd; if (c > prev + 1e-6) mono = false; prev = c; }
  const b = sd();
  ok('dispersal alone reduces spatial sd monotonically', mono && b.sd < 0.05 * a.sd, `sd ${a.sd.toFixed(2)} -> ${b.sd.toFixed(3)}`);
  ok('dispersal alone conserves the spatial mean', Math.abs(b.mu - a.mu) < 1e-2, `mean ${a.mu.toFixed(4)} -> ${b.mu.toFixed(4)}`);
}

// -------------------------------------- 3. stabilising selection only: rate
{
  const gam = 0.2;
  const F = new Field(24, 42, mulberry32(2));
  F.seedFlat(0);
  const p = { m: [0, 0, 0], M: ZERO, gamma: gam, sigma: 0, mode: 'clamp', localNe: false };
  const err = () => { let e = 0; for (let k = 0; k < F.cur.length; k++) e = Math.max(e, Math.abs(F.cur[k] - F.theta[k])); return e; };
  const e0 = err();
  const T = 30;
  for (let t = 0; t < T; t++) F.step(p);
  const eT = err();
  const pred = e0 * Math.pow(1 - gam, T);
  ok('z -> theta at rate (1-gamma)^t', Math.abs(eT - pred) / pred < 1e-4, `obs=${eT.toExponential(3)} pred=${pred.toExponential(3)}`);
}

// --------------------------------------------- 4. drift variance and scaling
{
  const sig = 1.5;
  const F = new Field(56, 98, mulberry32(3));
  F.seedFlat(128);
  const p = { m: [0, 0, 0], M: ZERO, gamma: 0, sigma: sig, mode: 'clamp', localNe: false };
  const T = 100;
  for (let t = 0; t < T; t++) F.step(p);
  let ss = 0, n = 0;
  for (let k = 0; k < F.cur.length; k++) { const d = F.cur[k] - 128; ss += d * d; n++; }
  const obs = ss / n, pred = sig * sig * T;
  ok('drift variance grows as sigma^2 * t', Math.abs(obs - pred) / pred < 0.05, `obs=${obs.toFixed(1)} pred=${pred.toFixed(1)}`);

  // local Ne: sigma(x) = sigma0 * sqrt(Ncore/Ne(x)), stronger at the margins
  const G = new Field(56, 98, mulberry32(4));
  G.seedFlat(128);
  const q2 = { ...p, localNe: true, sigma: 0.4 };
  for (let t = 0; t < 60; t++) G.step(q2);
  const cq = (Math.floor(G.wid / 2) + G.wid * Math.floor(G.hi / 2));
  const mq = 0; // corner
  const varAt = (q) => { let s = 0; for (let d = 0; d < 3; d++) { const e = G.cur[3 * q + d] - 128; s += e * e; } return s / 3; };
  ok('drift is stronger at the range margin than the core',
    varAt(mq) > varAt(cq),
    `scale core=${G.driftScale[cq].toFixed(2)} margin=${G.driftScale[mq].toFixed(2)}`);
}

// --------------------------------- 5. which pair each coevolution slider hits
{
  const probe = (s12, s13, s23) => {
    const M = buildM(s12, s13, s23, 2, 'chase');
    const touched = new Set();
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) if (M[a][b] !== 0) { touched.add(a); touched.add(b); }
    return [...touched].sort().join('');
  };
  ok('s12 slider drives species 1 and 2 only', probe(0.1, 0, 0) === '01', probe(0.1, 0, 0));
  ok('s13 slider drives species 1 and 3 only', probe(0, 0.1, 0) === '02', probe(0, 0.1, 0));
  ok('s23 slider drives species 2 and 3 only', probe(0, 0, 0.1) === '12', probe(0, 0, 0.1));

  // the original mapping, for comparison
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  const Msel13 = 0.1, Msel12 = 0, Msel23 = 0, Mratio = 2;
  M[1][0] = Msel12; M[2][1] = Msel13; M[0][2] = Msel23;
  M[0][1] = -Msel12; M[1][2] = -Mratio * Msel13; M[2][0] = -Mratio * Msel23;
  const t = new Set();
  for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) if (M[a][b] !== 0) { t.add(a); t.add(b); }
  ok('original "spp 1 & 3" slider actually drove species 2 and 3', [...t].sort().join('') === '12', [...t].sort().join(''));
}

// ------------------------------------------------ 6. synchronous vs in-place
{
  // one step of the original in-place rule vs the synchronous rule,
  // on a single cell with no dispersal
  const M = buildM(0.08, 0.08, 0.08, 2, 'chase');
  const z = [40, 150, 220];
  // in place (original loop order)
  const a = z.slice();
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) a[i] += M[j][i] * (a[j] - a[i]);
  // synchronous
  const b = z.slice();
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) b[i] += M[j][i] * (z[j] - z[i]);
  const d = Math.max(...a.map((v, i) => Math.abs(v - b[i])));
  ok('in-place update differs from synchronous update', d > 0.1,
    `max|diff| = ${d.toFixed(3)} trait units per generation`);
}

// ------------------------------------------------------- 7. stability bounds
{
  // worst case checkerboard mode: neighbour mean = -z, so dispersal
  // contributes -2m on the diagonal.  |1-2m| < 1 requires m < 1.
  const mMax = 0.75;
  ok('dispersal cap keeps the checkerboard mode contracting', Math.abs(1 - 2 * mMax) < 1, `|1-2m| = ${Math.abs(1 - 2 * mMax).toFixed(2)}`);

  // full run at the most aggressive settings: must stay finite and in range
  const F = new Field(40, 70, mulberry32(9));
  F.randomize();
  for (const mode of ['clamp', 'wrap', 'fold', 'soft']) {
    const p = {
      m: [mMax, mMax, mMax],
      M: buildM(0.1, 0.1, 0.1, 3, 'chase'),
      gamma: 0.4, sigma: 3, mode, localNe: true,
    };
    for (let t = 0; t < 300; t++) F.step(p);
    let bad = 0;
    const hi = mode === 'wrap' ? 256 : 255 + 1e-9;
    for (let k = 0; k < F.cur.length; k++) if (!Number.isFinite(F.cur[k]) || F.cur[k] < 0 || F.cur[k] >= hi) bad++;
    ok(`sliders pinned, mode=${mode}: state stays finite and on scale`, bad === 0, `${bad} bad cells`);
    F.randomize();
  }
}

// ------------------------------------------------------------- 8. throughput
{
  const F = new Field(80, 140, mulberry32(11));
  F.randomize();
  const p = { m: [0.3, 0.3, 0.3], M: buildM(0.05, 0.02, 0.05, 2, 'chase'), gamma: 0.05, sigma: 0.5, mode: 'clamp', localNe: true };
  const t0 = process.hrtime.bigint();
  for (let t = 0; t < 200; t++) F.step(p);
  const ms = Number(process.hrtime.bigint() - t0) / 1e6 / 200;
  ok('80x140 generation under 4 ms', ms < 4, `${ms.toFixed(2)} ms/generation`);
}

console.log(fails === 0 ? '\nall checks passed' : `\n${fails} FAILED`);
process.exit(fails ? 1 : 0);
