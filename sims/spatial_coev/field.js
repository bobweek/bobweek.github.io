/* =====================================================================
   field.js  —  lattice state + one-generation update kernel
   (replaces cell.js)

   No p5 dependency, on purpose: the kernel can be require()'d and
   checked numerically under node (see test_field.js).

   STATE
     Three species, each with a mean phenotype z_s(x) on a wid x hi torus,
     stored interleaved in a Float32Array:  z0,z1,z2, z0,z1,z2, ...
     Traits live on 0..255 so a cell's state *is* its colour.

   ONE GENERATION  (synchronous: read cur, write next, then swap)

     z_s  <-  z_s
            + m_s ( <z_s>_nbrs - z_s )        dispersal (4-neighbour mean)
            + gamma ( theta_s(x) - z_s )      stabilising selection
            + sum_t M[t][s] ( z_t - z_s )     coevolutionary selection
            + sqrt( G / Ne(x) ) * xi          drift,  xi ~ N(0,1)
     then
       z_s  <-  squash(z_s)                   keep it on the 0..255 scale

     M[t][s] is the pressure species t exerts on species s (the convention
     in the comment of the original sketch.js).  Every term on the right is
     evaluated at the *previous* state, so the three species are updated
     simultaneously rather than in index order.
   ===================================================================== */

const NSP = 3;      // species
const ZMAX = 255;   // trait scale == colour channel

// Lattice sizes; all exactly 4:7 so cells stay square as resolution changes.
// Cost is roughly linear in cell count up to ~180k, then goes superlinear as
// the state arrays stop fitting in cache.  Measured step+blit per generation:
// 80x140 1.5ms, 160x280 6ms, 320x560 23ms, 448x784 60ms, 640x1120 91ms.  The
// last two run below 30fps; the patterns move slowly enough to still read.
const RESOLUTIONS = [
  [24, 42], [40, 70], [56, 98], [80, 140], [112, 196], [160, 280],
  [224, 392], [320, 560], [448, 784], [640, 1120],
];

// --- standard normal, Box-Muller with a cached spare -------------------
let _spare = null;
function gauss(rand) {
  if (_spare !== null) { const s = _spare; _spare = null; return s; }
  let u, v, q;
  do {
    u = 2 * rand() - 1;
    v = 2 * rand() - 1;
    q = u * u + v * v;
  } while (q === 0 || q >= 1);
  const f = Math.sqrt((-2 * Math.log(q)) / q);
  _spare = v * f;
  return u * f;
}

// --- ways of keeping the trait on 0..255 -------------------------------
// "soft" is calibrated so its slope at mid-range is exactly 1: it is the
// identity to first order in the middle and only bites near the bounds.
const SOFT_K = 4 / ZMAX;
const SQUASH = {
  clamp: (v) => (v < 0 ? 0 : v > ZMAX ? ZMAX : v),
  wrap: (v) => { const p = ZMAX + 1; return ((v % p) + p) % p; },
  fold: (v) => { const p = 2 * ZMAX; const x = ((v % p) + p) % p; return x > ZMAX ? p - x : x; },
  soft: (v) => ZMAX / (1 + Math.exp(-SOFT_K * (v - ZMAX / 2))),
};
const SQUASH_NAMES = ['clamp', 'wrap', 'fold', 'soft'];

// --- coevolution matrix ------------------------------------------------
//
// M[t][s] is the pressure species t exerts on species s, entering the update
// for z_s as M[t][s] * (z_t - z_s).
//
// Each pair carries two signed coefficients: how hard the tracker pulls
// toward its target, and how hard the target pushes away.  Orientations are
// fixed so the tracking relations close the cycle 1 -> 2 -> 3 -> 1, which is
// the rock-paper-scissors structure:
//
//     pair 1&2:   1 tracks 2,   2 flees 1
//     pair 1&3:   3 tracks 1,   1 flees 3
//     pair 2&3:   2 tracks 3,   3 flees 2
//
// Because both numbers are signed and independent, a pair can run at any
// asymmetry (the 2016 defaults used 1:1 on pair 1&2 and 1:2 on pair 1&3, and
// that mismatch is a large part of how the original looked).  A negative flee
// turns a pair into mutual matching; a negative track into mutual avoidance.
const PAIRS = [
  // [tracker, target, track key, flee key]
  [0, 1, 't12', 'f12'],
  [2, 0, 't13', 'f13'],
  [1, 2, 't23', 'f23'],
];

function buildM(c) {
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const [a, b, tk, fk] of PAIRS) {
    M[b][a] = +c[tk];    // a moves toward b
    M[a][b] = -c[fk];    // b moves away from a
  }
  return M;
}

class Field {
  constructor(wid, hi, rand) {
    this.rand = rand || Math.random;
    this.alloc(wid, hi);
  }

  alloc(wid, hi) {
    this.wid = wid;
    this.hi = hi;
    this.n = wid * hi;
    this.cur = new Float32Array(this.n * NSP);
    this.next = new Float32Array(this.n * NSP);
    this.theta = new Float32Array(this.n * NSP);   // local optima
    this.driftScale = new Float32Array(this.n);    // sqrt(Ncore / Ne(x))
    this.buildLandscape();
  }

  // Optima and abundance are defined in normalised lattice coordinates,
  // so the landscape does not change shape when the canvas or the
  // resolution changes.
  buildLandscape() {
    const peaks = [[3 / 8, 2 / 5], [5 / 8, 2 / 5], [1 / 2, 3 / 5]];
    const sTh = 0.16;                 // width of each optimum
    const sN = 0.14;                  // width of the abundance hump
    const Nmin = 10, Nmax = 500, Ncore = Nmin + Nmax;
    for (let j = 0; j < this.hi; j++) {
      const y = (j + 0.5) / this.hi;
      for (let i = 0; i < this.wid; i++) {
        const x = (i + 0.5) / this.wid;
        const q = i + this.wid * j;
        const k = NSP * q;
        for (let s = 0; s < NSP; s++) {
          const dx = x - peaks[s][0], dy = y - peaks[s][1];
          this.theta[k + s] = 55 + 200 * Math.exp(-(dx * dx + dy * dy) / (2 * sTh * sTh));
        }
        const dx = x - 0.5, dy = y - 0.5;
        const Ne = Nmin + Nmax * Math.exp(-(dx * dx + dy * dy) / (2 * sN * sN));
        // sigma = sqrt(G/Ne); the slider sets sigma at the range core, and
        // it grows like sqrt(Ncore/Ne) out toward the margins
        this.driftScale[q] = Math.sqrt(Ncore / Ne);
      }
    }
  }

  randomize() {
    for (let k = 0; k < this.cur.length; k++) this.cur[k] = this.rand() * ZMAX;
  }

  seedOptima() { this.cur.set(this.theta); }

  seedFlat(v) {
    const z = v === undefined ? ZMAX / 2 : v;
    this.cur.fill(z);
  }

  // nearest-neighbour rescale from another field (used by the resolution
  // buttons) — a true rescale, not a tiling
  resampleFrom(old) {
    for (let j = 0; j < this.hi; j++) {
      const oj = Math.min(old.hi - 1, Math.floor((j * old.hi) / this.hi));
      for (let i = 0; i < this.wid; i++) {
        const oi = Math.min(old.wid - 1, Math.floor((i * old.wid) / this.wid));
        const k = NSP * (i + this.wid * j);
        const ok = NSP * (oi + old.wid * oj);
        for (let s = 0; s < NSP; s++) this.cur[k + s] = old.cur[ok + s];
      }
    }
  }

  step(p) {
    const wid = this.wid, hi = this.hi;
    const cur = this.cur, next = this.next, theta = this.theta;
    const m = p.m, M = p.M, gam = p.gamma, sig = p.sigma;
    const squash = SQUASH[p.mode] || SQUASH.clamp;
    const localNe = !!p.localNe;
    // 'sequential' reproduces the original sketch: species are updated in
    // index order, so species s sees the already-updated z_t for t < s, and
    // z_s is re-read after each term of the coupling sum.  'simultaneous'
    // evaluates every term at the previous state.  The two give visibly
    // different dynamics; sequential sustains more motion.
    const seq = p.sequential !== false;
    const rand = this.rand;
    const z = [0, 0, 0], a = [0, 0, 0];

    for (let j = 0; j < hi; j++) {
      const jm = (j - 1 + hi) % hi, jp = (j + 1) % hi;
      for (let i = 0; i < wid; i++) {
        const im = (i - 1 + wid) % wid, ip = (i + 1) % wid;
        const q = i + wid * j;
        const k = NSP * q;
        const kL = NSP * (im + wid * j), kR = NSP * (ip + wid * j);
        const kU = NSP * (i + wid * jm), kD = NSP * (i + wid * jp);

        for (let s = 0; s < NSP; s++) {
          z[s] = cur[k + s];
          a[s] = 0.25 * (cur[kL + s] + cur[kR + s] + cur[kU + s] + cur[kD + s]);
        }

        const sc = localNe ? this.driftScale[q] : 1;

        if (seq) {
          for (let s = 0; s < NSP; s++) {
            z[s] += m[s] * (a[s] - z[s]);
            if (gam !== 0) z[s] += gam * (theta[k + s] - z[s]);
            if (sig > 0) z[s] += sig * sc * gauss(rand);
            for (let t = 0; t < NSP; t++) z[s] += M[t][s] * (z[t] - z[s]);
          }
          for (let s = 0; s < NSP; s++) next[k + s] = squash(z[s]);
        } else {
          for (let s = 0; s < NSP; s++) {
            let v = z[s];
            v += m[s] * (a[s] - z[s]);
            if (gam !== 0) v += gam * (theta[k + s] - z[s]);
            for (let t = 0; t < NSP; t++) v += M[t][s] * (z[t] - z[s]);
            if (sig > 0) v += sig * sc * gauss(rand);
            next[k + s] = squash(v);
          }
        }
      }
    }

    const tmp = this.cur; this.cur = this.next; this.next = tmp;
  }

  // fill an RGBA byte array (a p5.Image pixel buffer) at lattice resolution
  writePixels(px, view) {
    const src = view === 'optima' ? this.theta : this.cur;
    for (let q = 0, k = 0; q < this.n; q++, k += NSP) {
      const o = q << 2;
      px[o] = src[k];
      px[o + 1] = src[k + 1];
      px[o + 2] = src[k + 2];
      px[o + 3] = 255;
    }
  }
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { Field, buildM, PAIRS, SQUASH, SQUASH_NAMES, RESOLUTIONS, NSP, ZMAX, gauss };
}
