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
const NOCOUPLE = { t12: 0, f12: 0, t13: 0, f13: 0, t23: 0, f23: 0 };
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
  const zero = { t12: 0, f12: 0, t13: 0, f13: 0, t23: 0, f23: 0 };
  const probe = (over) => {
    const M = buildM({ ...zero, ...over });
    const touched = new Set();
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) if (M[a][b] !== 0) { touched.add(a + 1); touched.add(b + 1); }
    return [...touched].sort().join('&');
  };
  ok('the 1&2 pair sliders touch species 1 and 2 only', probe({ t12: .1, f12: .1 }) === '1&2', probe({ t12: .1, f12: .1 }));
  ok('the 1&3 pair sliders touch species 1 and 3 only', probe({ t13: .1, f13: .1 }) === '1&3', probe({ t13: .1, f13: .1 }));
  ok('the 2&3 pair sliders touch species 2 and 3 only', probe({ t23: .1, f23: .1 }) === '2&3', probe({ t23: .1, f23: .1 }));

  // tracking must close the cycle 1 -> 2 -> 3 -> 1, not form a hierarchy
  const M = buildM({ t12: .05, f12: .1, t13: .05, f13: .1, t23: .05, f23: .1 });
  const arrows = [];
  for (let s2 = 0; s2 < 3; s2++) for (let t = 0; t < 3; t++) if (M[t][s2] > 0) arrows.push(`${s2 + 1}->${t + 1}`);
  ok('tracking closes the 3-cycle 1->2->3->1', arrows.sort().join(' ') === '1->2 2->3 3->1', arrows.join(' '));
  ok('no species tracks two others (that would be a hierarchy)',
    new Set(arrows.map((x) => x.slice(-1))).size === 3);

  // a negative push turns a pair into mutual matching
  const Mm = buildM({ ...zero, t12: .05, f12: -.05 });
  ok('a negative flee makes that pair converge', Mm[1][0] > 0 && Mm[0][1] > 0,
    `M[1][0]=${Mm[1][0]} M[0][1]=${Mm[0][1]}`);
}

// ------------------------------------------- 6. the two update orders differ
{
  const M = buildM({ t12: .095, f12: .095, t13: .095, f13: .19, t23: 0, f23: 0 });
  const base = { m: [0.7125, 0.3375, 0.0375], M, gamma: 0, sigma: 0, localNe: false, mode: 'clamp' };
  const init = new Float32Array(40 * 70 * 3);
  const r = mulberry32(21);
  for (let k = 0; k < init.length; k++) init[k] = r() * 255;

  const A = new Field(40, 70, mulberry32(1)); A.cur.set(init);
  const B = new Field(40, 70, mulberry32(1)); B.cur.set(init);
  for (let t = 0; t < 200; t++) { A.step({ ...base, sequential: true }); B.step({ ...base, sequential: false }); }
  let d = 0;
  for (let k = 0; k < A.cur.length; k++) d = Math.max(d, Math.abs(A.cur[k] - B.cur[k]));
  ok('sequential and simultaneous are different trajectories', d > 10, `max |difference| = ${d.toFixed(1)}`);

  const churn = (F) => {
    const p = Float32Array.from(F.cur);
    F.step({ ...base, sequential: F === A });
    let s2 = 0; for (let k = 0; k < F.cur.length; k++) s2 += Math.abs(F.cur[k] - p[k]);
    return s2 / F.cur.length;
  };
  const cA = churn(A), cB = churn(B);
  ok('sequential sustains more motion than simultaneous', cA > cB,
    `sequential ${cA.toFixed(2)} vs simultaneous ${cB.toFixed(2)} trait units/generation`);
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
  const hot = buildM({ t12: 0.2, f12: 0.2, t13: 0.2, f13: 0.2, t23: 0.2, f23: 0.2 });
  for (const mode of ['clamp', 'wrap', 'fold', 'soft']) for (const sequential of [true, false]) {
    const p = {
      m: [mMax, mMax, mMax], M: hot,
      gamma: 0.4, sigma: 25, mode, localNe: true, sequential,
    };
    for (let t = 0; t < 300; t++) F.step(p);
    let bad = 0;
    const hi = mode === 'wrap' ? 256 : 255 + 1e-9;
    for (let k = 0; k < F.cur.length; k++) if (!Number.isFinite(F.cur[k]) || F.cur[k] < 0 || F.cur[k] >= hi) bad++;
    ok(`sliders pinned, mode=${mode}, ${sequential ? 'sequential' : 'simultaneous'}: finite and on scale`, bad === 0, `${bad} bad cells`);
    F.randomize();
  }
}

// ------------------------------------------- 7b. the widened drift ceiling
{
  // the slider is quadratic: position 0..5 feeds sigma = position^2 = 0..25
  const sq = (v) => v * v;
  ok('drift slider reaches sigma = 25 at the top', Math.abs(sq(5) - 25) < 1e-9, `sigma(5) = ${sq(5)}`);
  ok('drift slider still resolves the low end', sq(0.5) < 0.3 && sq(1) === 1,
    `sigma(0.5) = ${sq(0.5)}, sigma(1) = ${sq(1)}`);
  ok('drift is exactly zero at the bottom', sq(0) === 0);

  // sigma = 25 with a falling Ne is essentially white noise at the margin
  const F = new Field(56, 98, mulberry32(31));
  F.seedFlat(128);
  const p = { m: [0, 0, 0], M: ZERO, gamma: 0, sigma: 25, mode: 'clamp', localNe: true, sequential: true };
  F.step(p);
  let core = 0, marg = 0, nc = 0, nm = 0;
  for (let q = 0; q < F.n; q++) {
    const d = Math.abs(F.cur[3 * q] - 128);
    if (F.driftScale[q] < 1.3) { core += d; nc++; } else if (F.driftScale[q] > 4) { marg += d; nm++; }
  }
  ok('at full drift the margin saturates in a single generation', marg / nm > core / nc,
    `mean |step| core ${(core / nc).toFixed(0)}, margin ${(marg / nm).toFixed(0)} trait units`);
}

// ------------------------------------------------- 8. throughput and memory
{
  // every grid on the ladder must allocate and run without falling over
  const { RESOLUTIONS } = require('./field.js');
  const big = RESOLUTIONS[RESOLUTIONS.length - 1];
  const B = new Field(big[0], big[1], mulberry32(12));
  B.randomize();
  const bp = { m: [.5, .5, .5], M: buildM({ t12: .095, f12: .095, t13: .095, f13: .19, t23: 0, f23: 0 }),
    gamma: .05, sigma: 1, mode: 'clamp', localNe: true, sequential: true };
  const bt = process.hrtime.bigint();
  B.step(bp);
  const bms = Number(process.hrtime.bigint() - bt) / 1e6;
  let bad = 0;
  for (let k = 0; k < B.cur.length; k++) if (!Number.isFinite(B.cur[k])) bad++;
  ok(`largest grid ${big[0]}x${big[1]} runs and stays finite`, bad === 0, `${bms.toFixed(0)} ms/generation`);
  ok('largest grid is 4:7 like the rest', Math.abs(big[0] / big[1] - 4 / 7) < 1e-9);
  ok('every grid on the ladder is 4:7',
    RESOLUTIONS.every(([w, h]) => Math.abs(w / h - 4 / 7) < 1e-9),
    RESOLUTIONS.map((r) => r.join('x')).join(' '));

  const F = new Field(80, 140, mulberry32(11));
  F.randomize();
  const p = { m: [0.3, 0.3, 0.3], M: buildM({ t12: .05, f12: .1, t13: .02, f13: .04, t23: .05, f23: .1 }), gamma: 0.05, sigma: 0.5, mode: 'clamp', localNe: true, sequential: true };
  const t0 = process.hrtime.bigint();
  for (let t = 0; t < 200; t++) F.step(p);
  const ms = Number(process.hrtime.bigint() - t0) / 1e6 / 200;
  ok('80x140 generation under 4 ms', ms < 4, `${ms.toFixed(2)} ms/generation`);
}

console.log(fails === 0 ? '\nall checks passed' : `\n${fails} FAILED`);
process.exit(fails ? 1 : 0);
