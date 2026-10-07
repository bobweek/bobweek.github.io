#!/usr/bin/env python3
"""Build every seminar figure made so far (batches 0 to 4).

    python3 scripts/make_all.py
    SEMINAR_MONO=1   python3 scripts/make_all.py   # pure greyscale
    SEMINAR_STRICT=1 python3 scripts/make_all.py   # fail on any text collision
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import b0_glyph_kit
import s03_map
import s04_18_pipeline
import s05_drift_G
import s05_roadmap
import s06_interface
import s07_interfaces
import s08_arms_race
import s09_11_zoom
import s14_17_space
import s19_likelihood
import s20_abc
import s22_hinge
import s24_ancestral_concordance
import s25_27_rescue_qme
import s_misc_batch4

MODULES = (b0_glyph_kit, s03_map, s04_18_pipeline, s05_drift_G, s05_roadmap, s06_interface,
           s07_interfaces, s08_arms_race, s09_11_zoom, s14_17_space,
           s19_likelihood, s20_abc,
           s22_hinge, s24_ancestral_concordance, s25_27_rescue_qme,
           s_misc_batch4)

if __name__ == "__main__":
    for m in MODULES:
        getattr(m, "main_all", m.main)()
