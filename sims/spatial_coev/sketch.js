// --- globals (ported from Processing) ---
let grid, dummy;

let w, h,
  wid = 80,
  hi = 140,
  old_wid,
  old_hi,
  neighborhood = 0,
  hud = 1,
  picCount,
  k = 0,
  e = 1,
  count = 0;

let C = 0.01,
  mix = 0,
  loge = 0.01,
  Mratio = 2,
  Msel12,
  Msel13,
  Msel23;

let m = [0.1, 0.1, 0.1];
let drift = [0, 0, 0];
let G = [1, 1, 1],
  Gamma = [0, 0, 0];

// matching selection matrix: M[i][j] is selection pressure from species i onto species j
let M = [
  [0, 0, 0],
  [0, 0, 0],
  [0, 0, 0],
];

let m_SL = "",
  mig = "",
  M_SL = "",
  Mratio_SL = "",
  A_SL = "",
  Image = "",
  Ran = "",
  Size = "",
  contrast = "",
  Mix = "",
  log_e,
  filename;

let number;
// let size = 10;
let c;

let pic = false,
  rando = false,
  vid = false,
  rev = false,
  capture = false,
  pause = false,
  pause1 = false;

// UI
let mig1, mig2, mig3, m12, m13, m23;

function setup() {
  document.body.style.margin = "0";
  document.body.style.padding = "0";
  document.body.style.overflow = "hidden";

  pixelDensity(1);
  
  createCanvas(400, 700);
  noStroke();
  background(0);

  // sliders
  mig1 = new HScrollbar(25, 50, 120, 20, 5, 120);
  mig2 = new HScrollbar(25, 100, 120, 20, 5, 70);
  mig3 = new HScrollbar(25, 150, 120, 20, 5, 30);
  m12 = new HScrollbar(25, 200, 120, 20, 5, 120);
  m13 = new HScrollbar(25, 250, 120, 20, 5, 20);
  m23 = new HScrollbar(25, 300, 120, 20, 5, 120);

  // grid init
  grid = new Array(wid);
  for (let i = 0; i < wid; i++) {
    grid[i] = new Array(hi);
    for (let j = 0; j < hi; j++) {
      grid[i][j] = new Cell(random(0, 255), random(0, 255), random(0, 255), i, j);
    }
  }
  
  
  // initGridFromInstructions();
}

function draw() {
  if (!pause) {
    // render
    for (let i = 0; i < wid; i++) {
      for (let j = 0; j < hi; j++) grid[i][j].render();
    }

    // read slider values -> parameters
    m[0] = 0.75 * (mig1.spos - mig1.sposMin) / (mig1.sposMax - mig1.sposMin);
    m[1] = 0.75 * (mig2.spos - mig2.sposMin) / (mig2.sposMax - mig2.sposMin);
    m[2] = 0.75 * (mig3.spos - mig3.sposMin) / (mig3.sposMax - mig3.sposMin);

    Mratio = 2;
    Msel12 = 0.1 * (m12.spos - m12.sposMin) / (m12.sposMax - m12.sposMin);
    Msel13 = 0.1 * (m13.spos - m13.sposMin) / (m13.sposMax - m13.sposMin);
    Msel23 = 0.1 * (m23.spos - m23.sposMin) / (m23.sposMax - m23.sposMin);

    // fill selection matrix
    M[1][0] = Msel12;
    M[2][1] = Msel13;
    M[0][2] = Msel23;
    M[0][1] = -Msel12;
    M[1][2] = -Mratio * Msel13;
    M[2][0] = -Mratio * Msel23;

    // dummy copy
    dummy = new Array(wid);
    for (let i = 0; i < wid; i++) {
      dummy[i] = new Array(hi);
      for (let j = 0; j < hi; j++) {
        dummy[i][j] = new Cell(
          grid[i][j].species[0],
          grid[i][j].species[1],
          grid[i][j].species[2],
          i,
          j
        );
      }
    }

    // update (reverse loops kept)
    for (let i = wid - 1; i >= 0; i--) {
      for (let j = hi - 1; j >= 0; j--) grid[i][j].update(dummy);
    }

    // HUD
    if (hud === 0) {
      mig1.update(); mig1.display();
      mig2.update(); mig2.display();
      mig3.update(); mig3.display();
      m12.update();  m12.display();
      m13.update();  m13.display();
      m23.update();  m23.display();

      fill(0);
      rect(mig1.xpos - 2, mig1.ypos - 20, 90, 12);
      rect(mig2.xpos - 2, mig2.ypos - 20, 90, 12);
      rect(mig3.xpos - 2, mig3.ypos - 20, 90, 12);
      rect(m12.xpos - 2, m12.ypos - 20, 125, 12);
      rect(m13.xpos - 2, m13.ypos - 20, 125, 12);
      rect(m23.xpos - 2, m23.ypos - 20, 125, 12);

      fill(255);
      text("spp1 migration", mig1.xpos, mig1.ypos - 10);
      text("spp2 migration", mig2.xpos, mig2.ypos - 10);
      text("spp3 migration", mig3.xpos, mig3.ypos - 10);
      text("spp 1 & 2 coevolution", m12.xpos, m12.ypos - 10);
      text("spp 1 & 3 coevolution", m13.xpos, m13.ypos - 10);
      text("spp 2 & 3 coevolution", m23.xpos, m23.ypos - 10);
    }
  } else {
    if (pause1) {
      // one-step redraw/update while paused (mirrors your Processing logic)
      for (let i = 0; i < wid; i++) {
        for (let j = 0; j < hi; j++) grid[i][j].render();
      }

      m[0] = 0.75 * (mig1.spos - mig1.sposMin) / (mig1.sposMax - mig1.sposMin);
      m[1] = 0.75 * (mig2.spos - mig2.sposMin) / (mig2.sposMax - mig2.sposMin);
      m[2] = 0.75 * (mig3.spos - mig3.sposMin) / (mig3.sposMax - mig3.sposMin);

      Mratio = 2;
      Msel12 = 0.1 * (m12.spos - m12.sposMin) / (m12.sposMax - m12.sposMin);
      Msel13 = 0.1 * (m13.spos - m13.sposMin) / (m13.sposMax - m13.sposMin);
      Msel23 = 0.1 * (m23.spos - m23.sposMin) / (m23.sposMax - m23.sposMin);

      M[1][0] = Msel12;
      M[2][1] = Msel13;
      M[0][2] = Msel23;
      M[0][1] = -Msel12;
      M[1][2] = -Mratio * Msel13;
      M[2][0] = -Mratio * Msel23;

      dummy = new Array(wid);
      for (let i = 0; i < wid; i++) {
        dummy[i] = new Array(hi);
        for (let j = 0; j < hi; j++) {
          dummy[i][j] = new Cell(
            grid[i][j].species[0],
            grid[i][j].species[1],
            grid[i][j].species[2],
            i,
            j
          );
        }
      }
      for (let i = wid - 1; i >= 0; i--) {
        for (let j = hi - 1; j >= 0; j--) grid[i][j].update(dummy);
      }

      mig1.update(); mig1.display();
      mig2.update(); mig2.display();
      mig3.update(); mig3.display();
      m12.update();  m12.display();
      m13.update();  m13.display();
      m23.update();  m23.display();

      fill(0);
      rect(mig1.xpos - 2, mig1.ypos - 20, 90, 12);
      rect(mig2.xpos - 2, mig2.ypos - 20, 90, 12);
      rect(mig3.xpos - 2, mig3.ypos - 20, 90, 12);
      rect(m12.xpos - 2, m12.ypos - 20, 125, 12);
      rect(m13.xpos - 2, m13.ypos - 20, 125, 12);
      rect(m23.xpos - 2, m23.ypos - 20, 125, 12);

      fill(255);
      text("spp1 migration", mig1.xpos, mig1.ypos - 10);
      text("spp2 migration", mig2.xpos, mig2.ypos - 10);
      text("spp3 migration", mig3.xpos, mig3.ypos - 10);
      text("spp 1 & 2 coevolution", m12.xpos, m12.ypos - 10);
      text("spp 1 & 3 coevolution", m13.xpos, m13.ypos - 10);
      text("spp 2 & 3 coevolution", m23.xpos, m23.ypos - 10);
    }

    fill(0);
    rect(140, 360, 220, 90);
    fill(255);
    text("press h to toggle controls", 150, 380);
    text("press p to pause/play", 150, 400);
    text("press +/- to incr/decr resolution", 150, 420);
    text("press r to reset with random values", 150, 440);
    pause1 = false;
  }

}

function keyReleased() {
  switch (key) {
    case '+': {
      // copy current grid into dummy
      dummy = new Array(wid);
      for (let i = 0; i < wid; i++) {
        dummy[i] = new Array(hi);
        for (let j = 0; j < hi; j++) {
          dummy[i][j] = new Cell(
            grid[i][j].species[0],
            grid[i][j].species[1],
            grid[i][j].species[2],
            i,
            j
          );
        }
      }

      old_wid = wid;
      old_hi = hi;

      hi += 3;
      wid += 4;

      grid = new Array(wid);
      for (let i = 0; i < wid; i++) {
        grid[i] = new Array(hi);
        for (let j = 0; j < hi; j++) {
          grid[i][j] = new Cell(
            dummy[i % (old_wid - 1)][j % (old_hi - 1)].species[0],
            dummy[i % (old_wid - 1)][j % (old_hi - 1)].species[1],
            dummy[i % (old_wid - 1)][j % (old_hi - 1)].species[2],
            i,
            j
          );
        }
      }
      dummy = null;
      break;
    }

    case '-': {
      dummy = new Array(wid);
      for (let i = 0; i < wid; i++) {
        dummy[i] = new Array(hi);
        for (let j = 0; j < hi; j++) {
          dummy[i][j] = new Cell(
            grid[i][j].species[0],
            grid[i][j].species[1],
            grid[i][j].species[2],
            i,
            j
          );
        }
      }

      old_wid = wid;
      old_hi = hi;

      hi -= 3;
      if (hi < 3) hi = 3;

      wid -= 4;
      if (wid < 4) wid = 4;

      grid = new Array(wid);
      for (let i = 0; i < wid; i++) {
        grid[i] = new Array(hi);
        for (let j = 0; j < hi; j++) {
          grid[i][j] = new Cell(
            dummy[i % (old_wid - 1)][j % (old_hi - 1)].species[0],
            dummy[i % (old_wid - 1)][j % (old_hi - 1)].species[1],
            dummy[i % (old_wid - 1)][j % (old_hi - 1)].species[2],
            i,
            j
          );
        }
      }
      dummy = null;
      break;
    }

    case 'h':
      hud = (1 + hud) % 2;
      break;

    case 'p':
      pause = !pause;
      pause1 = true;
      break;

    case 'r':
      for (let i = 0; i < wid; i++) {
        for (let j = 0; j < hi; j++) grid[i][j].set(random(0, 255), random(0, 255), random(0, 255));
      }
      break;
  }
}

// --- UI class (ported) ---
class HScrollbar {
  constructor(xp, yp, sw, sh, l, ispos) {
    this.swidth = sw;
    this.sheight = sh;
    const widthtoheight = sw - sh;
    this.ratio = sw / widthtoheight;

    this.xpos = xp;
    this.ypos = yp - sh / 2;

    this.spos = ispos;
    this.newspos = this.spos;

    this.sposMin = this.xpos;
    this.sposMax = this.xpos + this.swidth - this.sheight;

    this.loose = l;

    this.over = false;
    this.locked = false;
  }

  update() {
    this.over = this.overEvent();

    if (mouseIsPressed && this.over) this.locked = true;
    if (!mouseIsPressed) this.locked = false;

    if (this.locked) {
      this.newspos = constrain(mouseX - this.sheight / 2, this.sposMin, this.sposMax);
    }

    if (abs(this.newspos - this.spos) > 1) {
      this.spos = this.spos + (this.newspos - this.spos) / this.loose;
    }
  }

  overEvent() {
    return (
      mouseX > this.xpos &&
      mouseX < this.xpos + this.swidth &&
      mouseY > this.ypos &&
      mouseY < this.ypos + this.sheight
    );
  }

  display() {
    fill(204);
    rect(this.xpos, this.ypos, this.swidth, this.sheight);

    if (this.over || this.locked) fill(0);
    else fill(102);

    rect(this.spos, this.ypos, this.sheight, this.sheight);
  }

  getPos() {
    return this.spos * this.ratio;
  }
}

// --- spatial functions (ported) ---
function abiotic_selection(x, y, i) {
  return 0.1 * exp(-((y - height / 2) * (y - height / 2)) / (2 * 50000)) *
    exp(-((x - width / 2) * (x - width / 2)) / (2 * 50000));
}

function biotic_selection(x, y, i) {
  return 0.1 * exp(-((y - height / 2) * (y - height / 2)) / (2 * 50000)) *
    exp(-((x - width / 2) * (x - width / 2)) / (2 * 50000));
}

function theta(x, y, i) {
  let value = 0;
  if (i === 0)
    value =
      55 +
      200 *
        exp(-((y - (2 * height) / 5) * (y - (2 * height) / 5)) / (2 * 10000)) *
        exp(-((x - (3 * width) / 8) * (x - (3 * width) / 8)) / (2 * 10000));
  if (i === 1)
    value =
      55 +
      200 *
        exp(-((y - (2 * height) / 5) * (y - (2 * height) / 5)) / (2 * 10000)) *
        exp(-((x - (5 * width) / 8) * (x - (5 * width) / 8)) / (2 * 10000));
  if (i === 2)
    value =
      55 +
      200 *
        exp(-((y - (3 * height) / 5) * (y - (3 * height) / 5)) / (2 * 10000)) *
        exp(-((x - width / 2) * (x - width / 2)) / (2 * 10000));
  return value;
}

function N(x, y, i) {
  return 10 + 500 * exp(-((y - height / 2) * (y - height / 2)) / (2 * 5000)) *
    exp(-((x - width / 2) * (x - width / 2)) / (2 * 5000));
}


function initGridFromInstructions() {
  const msg =
    "press p to pause/play\n" +
    "press +/- to incr/decr resolution\n" +
    "press r to reset with random values";

  // offscreen buffer at lattice resolution
  const pg = createGraphics(wid, hi);
  pg.pixelDensity(1);
  pg.background(0);
  pg.fill(255);
  pg.noStroke();
  pg.textAlign(CENTER, CENTER);

  // pick a readable size in lattice-pixels
  pg.textSize(max(8, floor(min(wid, hi) * 0.12)));

  // draw text centered
  pg.text(msg, wid / 2, hi / 2);

  pg.loadPixels();

  // map buffer pixels -> Cell RGB
  for (let i = 0; i < wid; i++) {
    for (let j = 0; j < hi; j++) {
      const idx = 4 * (i + wid * j);
      const r = pg.pixels[idx + 0];
      const g = pg.pixels[idx + 1];
      const b = pg.pixels[idx + 2];
      grid[i][j].set(r, g, b);
    }
  }
}
