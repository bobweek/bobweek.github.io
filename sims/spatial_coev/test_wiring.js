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

// 2. the value readouts, which are built as id + 'v'
const fmt = [...js.matchAll(/fmt\('([A-Za-z0-9_]+)'/g)].map((m) => m[1] + 'v');
const sliderList = js.match(/const sliders = \[([^\]]+)\]/)[1]
  .split(',').map((s) => s.trim().replace(/'/g, ''));
const readouts = [...new Set([...fmt, ...sliderList.map((s) => s + 'v')])];
const missingOut = readouts.filter((id) => !htmlIds.has(id) && id !== 'spfvv');
ok('every value readout <output> exists', missingOut.length === 0, missingOut.join(', '));

// 3. every slider named in the sketch is a real range input
const badRange = sliderList.filter((id) => {
  const tag = html.match(new RegExp(`<input id="${id}"[^>]*>`));
  return !tag || !/type="range"/.test(tag[0]);
});
ok('every wired slider is a range input', badRange.length === 0, badRange.join(', '));

// 4. DEFAULTS keys all address real controls, with in-range values
const block = js.match(/const DEFAULTS = \{([\s\S]*?)\n\};/)[1];
const keys = [...block.matchAll(/([A-Za-z0-9_]+):/g)].map((m) => m[1]);
const badKeys = keys.filter((k) => !htmlIds.has(k));
ok('every DEFAULTS key names a control', badKeys.length === 0, badKeys.join(', '));

const defaults = {};
for (const m of block.matchAll(/([A-Za-z0-9_]+):\s*('?[A-Za-z0-9_.]+'?)/g)) {
  defaults[m[1]] = m[2].replace(/'/g, '');
}
const outOfRange = [];
for (const id of sliderList) {
  const tag = html.match(new RegExp(`<input id="${id}"[^>]*>`))[0];
  const lo = +tag.match(/min="([^"]+)"/)[1];
  const hi = +tag.match(/max="([^"]+)"/)[1];
  const markup = +tag.match(/value="([^"]+)"/)[1];
  if (markup < lo || markup > hi) outOfRange.push(`${id}: markup ${markup} outside [${lo},${hi}]`);
  if (id in defaults) {
    const d = +defaults[id];
    if (d < lo || d > hi) outOfRange.push(`${id}: default ${d} outside [${lo},${hi}]`);
    if (d !== markup) outOfRange.push(`${id}: markup ${markup} != default ${d}`);
  }
}
ok('no slider starts outside its own range', outOfRange.length === 0, outOfRange.join(' | '));

// 5. no control in the markup is left unwired
const controlIds = [...html.matchAll(/<(?:input|select|button)\b[^>]*\bid="([^"]+)"/g)].map((m) => m[1]);
const orphans = controlIds.filter((id) => !js.includes(`'${id}'`));
ok('no control in index.html is left unwired', orphans.length === 0, orphans.join(', '));

// 6. every <label for=> points at something real
const fors = [...html.matchAll(/\bfor="([^"]+)"/g)].map((m) => m[1]);
const danglingFor = fors.filter((id) => !htmlIds.has(id));
ok('every label points at a real control', danglingFor.length === 0, danglingFor.join(', '));

// 7. the select options match what the model actually accepts
const { SQUASH_NAMES } = require('./field.js');
const boundsOpts = [...html.matchAll(/<select id="bounds">([\s\S]*?)<\/select>/g)][0][1]
  .match(/value="([^"]+)"/g).map((s) => s.slice(7, -1));
ok('bounds options match the squash modes',
  boundsOpts.every((v) => SQUASH_NAMES.includes(v)) && boundsOpts.length === SQUASH_NAMES.length,
  boundsOpts.join(', '));

const modeOpts = [...html.matchAll(/<select id="mode">([\s\S]*?)<\/select>/g)][0][1]
  .match(/value="([^"]+)"/g).map((s) => s.slice(7, -1));
ok('interaction options match buildM', modeOpts.join(',') === 'chase,match,avoid', modeOpts.join(', '));

// 8. stylesheet and scripts are actually linked
ok('style.css is linked', /href="style\.css"/.test(html));
ok('field.js and sketch.js are both loaded, field first',
  html.indexOf('field.js') < html.indexOf('sketch.js') && html.includes('sketch.js'));
ok('cell.js is no longer loaded', !html.includes('cell.js'));

console.log(fails === 0 ? '\nwiring ok' : `\n${fails} FAILED`);
process.exit(fails ? 1 : 0);
