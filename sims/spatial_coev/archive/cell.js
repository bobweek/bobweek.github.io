class Cell {
  constructor(R, G, B, i, j) {
    this.species = [0, 0, 0];
    this.picColor = [0, 0, 0];
    this.vidColor = [0, 0, 0];
    this.av_species = [0, 0, 0];
    this.prev_species = [0, 0, 0];
    this.position = [0, 0];
    this.bore_count = 0;
    this.boring = 10;

    this.species[0] = R;
    this.species[1] = G;
    this.species[2] = B;

    this.prev_species[0] = this.species[0];
    this.prev_species[1] = this.species[1];
    this.prev_species[2] = this.species[2];

    this.position[0] = i;
    this.position[1] = j;
  }

  render() {
    fill(this.species[0], this.species[1], this.species[2]);
    rect(
      (this.position[0] * width) / wid,
      (this.position[1] * height) / hi,
      width / wid + 2,
      height / hi + 2
    );
  }

  set(R, G, B) {
    this.species[0] = R;
    this.species[1] = G;
    this.species[2] = B;
  }

  update(gridRef) {
    c = 0;

    for (let i = 0; i < 3; i++) {
      this.picColor[i] = 0;
      this.vidColor[i] = 0;
      this.av_species[i] = 0;
      this.prev_species[i] = this.species[i];
    }

    // nearest-neighbor average (toroidal wrap)
    for (let h = 0; h < 3; h++) this.av_species[h] += gridRef[(wid + this.position[0] - 1) % wid][(hi + this.position[1]) % hi].species[h];
    for (let h = 0; h < 3; h++) this.av_species[h] += gridRef[(wid + this.position[0] + 1) % wid][(hi + this.position[1]) % hi].species[h];
    for (let h = 0; h < 3; h++) this.av_species[h] += gridRef[(wid + this.position[0]) % wid][(hi + this.position[1] + 1) % hi].species[h];
    for (let h = 0; h < 3; h++) this.av_species[h] += gridRef[(wid + this.position[0]) % wid][(hi + this.position[1] - 1) % hi].species[h];
    for (let h = 0; h < 3; h++) this.av_species[h] /= 4;

    // dynamics
    for (let i = 0; i < 3; i++) {
      this.species[i] +=
        m[i] * (this.av_species[i] - this.species[i]) +
        Gamma[i] * (theta(((this.position[0] + 1) * width) / wid, ((this.position[1] + 1) * height) / hi, i) - this.species[i]);

      // drift if you re-enable it:
      // this.species[i] += drift[i] * randomGaussian();

      for (let j = 0; j < 3; j++) this.species[i] += M[j][i] * (this.species[j] - this.species[i]);
    }

    // "flippy flops"
    switch (e) {
      case 1:
        for (let h = 0; h < 3; h++) {
          if (this.species[h] > 255) this.species[h] = 255;
          if (this.species[h] < 0) this.species[h] = 0;
        }
        break;

      case 2:
        for (let h = 0; h < 3; h++) this.species[h] = (255 + this.species[h]) % 255;
        break;

      case 3:
        for (let h = 0; h < 3; h++) {
          this.species[h] = 255 / (1 + exp(-1 * loge * (this.species[h] - 127)));
        }
        break;

      case 4:
        for (let h = 0; h < 3; h++) {
          if (this.species[h] > 255) this.species[h] = 510 - this.species[h];
          if (this.species[h] < 0) this.species[h] = -this.species[h];
        }
        break;
    }
  }
}
