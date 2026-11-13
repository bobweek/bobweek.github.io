const fs = require('fs');
const html = fs.readFileSync('index.html', 'utf8');
const js = fs.readFileSync('sketch.js', 'utf8');

const htmlIds = new Set([...html.matchAll(/\bid="([^"]+)"/g)].map((m) => m[1]));
let fails = 0;
const ok = (n, c, x = '') => { console.log(`${c ? 'PASS' : 'FAIL'}  ${n}${x ? '   ' + x : ''}`); if (!c) fails++; };

// 1. every $('id') in the sketch resolves
const direct = [...js.matchAll(/\$\('([A-Za-z0-9_]+)'\)/g)].map((m) => m[1]);
const missing = [...new Set(direct)].filter((id) => !htmlIds.has(id));
ok('every $(id) in sketch.js exists in index.html', missing.length === 0, missing.join(', '));

// 2. every range slider has a matching <output id="Xv"> next to it
const ranges = [...html.matchAll(/<input id="([^"]+)" type="range"([^>]*)>/g)]
  .map((m) => ({ id: m[1], attrs: m[2] }));
const missingOut = ranges.filter((r) => !htmlIds.has(r.id + 'v')).map((r) => r.id);
ok('every slider has a value readout', missingOut.length === 0, missingOut.join(', '));
ok('there is at least one slider', ranges.length > 0, `${ranges.length} sliders`);

// 3. defaults are read from the markup, not duplicated in the sketch
ok('sketch.js has no hardcoded copy of the defaults',
  !/const DEFAULTS = \{[^}]*[0-9]/.test(js) && js.includes('captureDefaults'),
  'markup is the single source of truth');
ok('sketch.js has no hardcoded slider list',
  !/const sliders = \[/.test(js) && js.includes("querySelectorAll('#panel input[type=range]')"));

// 4. no slider starts outside its own range, and step divides the value
const outOfRange = [];
for (const r of ranges) {
  const lo = +r.attrs.match(/min="([^"]+)"/)[1];
  const hi = +r.attrs.match(/max="([^"]+)"/)[1];
  const step = +r.attrs.match(/step="([^"]+)"/)[1];
  const v = +r.attrs.match(/value="([^"]+)"/)[1];
  if (!(v >= lo && v <= hi)) outOfRange.push(`${r.id}: ${v} outside [${lo},${hi}]`);
  const n = (v - lo) / step;
  if (Math.abs(n - Math.round(n)) > 1e-6) outOfRange.push(`${r.id}: ${v} is not a multiple of step ${step}`);
}
ok('no slider starts outside its range or off its step', outOfRange.length === 0, outOfRange.join(' | '));

// 5. every stateful control sits inside #panel so captureDefaults() sees it
const panelStart = html.indexOf('<aside id="panel">');
const panelEnd = html.indexOf('</aside>');
const outside = [...html.matchAll(/<(?:input|select)\b[^>]*\bid="([^"]+)"/g)]
  .filter((m) => m.index < panelStart || m.index > panelEnd).map((m) => m[1]);
ok('every stateful control lives inside #panel', outside.length === 0, outside.join(', '));

// 6. no control in the markup is left unwired
const controlIds = [...html.matchAll(/<(?:input|select|button)\b[^>]*\bid="([^"]+)"/g)].map((m) => m[1]);
const orphans = controlIds.filter((id) => !js.includes(`'${id}'`));
ok('no control in index.html is left unwired', orphans.length === 0, orphans.join(', '));

// 7. every <label for=> points at something real
const fors = [...html.matchAll(/\bfor="([^"]+)"/g)].map((m) => m[1]);
const danglingFor = fors.filter((id) => !htmlIds.has(id));
ok('every label points at a real control', danglingFor.length === 0, danglingFor.join(', '));

// 8. the select options match what the model actually accepts
const { SQUASH_NAMES } = require('./field.js');
const boundsOpts = [...html.matchAll(/<select id="bounds">([\s\S]*?)<\/select>/g)][0][1]
  .match(/value="([^"]+)"/g).map((s) => s.slice(7, -1));
ok('bounds options match the squash modes',
  boundsOpts.every((v) => SQUASH_NAMES.includes(v)) && boundsOpts.length === SQUASH_NAMES.length,
  boundsOpts.join(', '));

const presetOpts = [...html.matchAll(/<select id="preset">([\s\S]*?)<\/select>/g)][0][1]
  .match(/value="([^"]+)"/g).map((x) => x.slice(7, -1));
const presetKeys = [...js.matchAll(/^  (\w+):\s*\{ t12:/gm)].map((m) => m[1]);
ok('every preset option except Custom has a definition',
  presetOpts.filter((v) => v !== 'custom').every((v) => presetKeys.includes(v))
  && presetKeys.every((k) => presetOpts.includes(k)),
  `options ${presetOpts.join(',')} | defined ${presetKeys.join(',')}`);

const orderOpts = [...html.matchAll(/<select id="order">([\s\S]*?)<\/select>/g)][0][1]
  .match(/value="([^"]+)"/g).map((x) => x.slice(7, -1));
ok('update-order options are the two the kernel knows',
  orderOpts.join(',') === 'sequential,simultaneous', orderOpts.join(', '));

// every preset value must be reachable on the slider it writes to
const coefTags = {};
for (const m of html.matchAll(/<input id="([tf]\d\d)"([^>]*)>/g)) coefTags[m[1]] = m[2];
const unreachable = [];
for (const m of js.matchAll(/^  (\w+):\s*\{([^}]*)\}/gm)) {
  for (const kv of m[2].matchAll(/(\w+):\s*(-?[\d.]+)/g)) {
    const tag = coefTags[kv[1]];
    if (!tag) { unreachable.push(`${m[1]}.${kv[1]}: no such slider`); continue; }
    const lo = +tag.match(/min="([^"]+)"/)[1], hi2 = +tag.match(/max="([^"]+)"/)[1];
    const st = +tag.match(/step="([^"]+)"/)[1], v = +kv[2];
    if (v < lo || v > hi2) unreachable.push(`${m[1]}.${kv[1]}=${v} outside [${lo},${hi2}]`);
    const n = (v - lo) / st;
    if (Math.abs(n - Math.round(n)) > 1e-6) unreachable.push(`${m[1]}.${kv[1]}=${v} off step ${st}`);
  }
}
ok('every preset value is reachable on its slider', unreachable.length === 0, unreachable.join(' | '));

// the markup defaults must equal the "original" preset, or the sim opens
// somewhere the preset dropdown does not admit to
const orig = {};
for (const kv of js.match(/original:\s*\{([^}]*)\}/)[1].matchAll(/(\w+):\s*(-?[\d.]+)/g)) orig[kv[1]] = +kv[2];
const mismatch = Object.entries(orig)
  .filter(([id, v]) => Math.abs(+coefTags[id].match(/value="([^"]+)"/)[1] - v) > 1e-9)
  .map(([id, v]) => `${id}: markup ${coefTags[id].match(/value="([^"]+)"/)[1]} vs preset ${v}`);
ok('markup opens on the Original preset', mismatch.length === 0, mismatch.join(' | '));

// 8b. the starting grid level, and the fractional-speed plumbing
const level = +html.match(/class="stepper" data-level="(\d+)"/)[1];
const { RESOLUTIONS } = require('./field.js');
ok('the starting grid level names a real resolution',
  Number.isInteger(level) && level >= 0 && level < RESOLUTIONS.length,
  `level ${level} = ${RESOLUTIONS[level] ? RESOLUTIONS[level].join('x') : 'out of range'}`);
ok('sketch.js reads the grid level from the markup',
  /dataset\.level/.test(js) && !/let resIdx = \d/.test(js));

const spf = html.match(/<input id="spf"([^>]*)>/)[1];
const spfMin = +spf.match(/min="([^"]+)"/)[1];
const spfVal = +spf.match(/value="([^"]+)"/)[1];
ok('generations per frame can go below 1', spfMin < 1, `min ${spfMin}`);
ok('the sketch carries a fractional budget rather than rounding',
  /genBudget/.test(js) && /Math\.floor\(genBudget\)/.test(js),
  `opens at ${spfVal} generations/frame`);

// the drift slider is quadratic, so check the markup lands where intended
const sig = html.match(/<input id="sigma"([^>]*)>/)[1];
const sigPos = +sig.match(/value="([^"]+)"/)[1];
const sigma = sigPos * sigPos;
const shown = +sigma.toFixed(2);
ok('drift opens close enough to the value the readout shows',
  Math.abs(sigma - shown) / Math.max(shown, 1e-9) < 0.005,
  `position ${sigPos} -> sigma ${sigma.toFixed(4)}, shown as ${shown.toFixed(2)}`);

// 9. stylesheet and scripts are actually linked
ok('style.css is linked', /href="style\.css"/.test(html));
ok('field.js and sketch.js are both loaded, field first',
  html.indexOf('field.js') < html.indexOf('sketch.js') && html.includes('sketch.js'));
ok('cell.js is no longer loaded', !html.includes('cell.js'));

console.log(fails === 0 ? '\nwiring ok' : `\n${fails} FAILED`);
process.exit(fails ? 1 : 0);
