Kristen same-deck paired control
================================

Contents
--------
- `bO3_ion_on.csv`
  Clean rerun of Kristen's manuscript-native `Model 1.zip` deck.
- `bO3_ion_off.csv`
  Same-deck control with only the charged-product photoprocess `delta` values zeroed.
- `eval_points.csv`
  The six Gerakines `O3/O2_0` comparison points used in the manuscript-facing fit checks.

Authoritative numbers
---------------------
- Ion-on:
  `chi = 0.484166`, weighted RMSD `= 0.547233`, peak `O3/O2_0 = 23.547%`.
- Ion-off:
  `chi = 13.9327`, weighted RMSD `= 16.5804`, peak `O3/O2_0 = 3.243%`.

Notes
-----
- These paired files are the manuscript-facing source of truth for the same-model ion-on / ion-off comparison.
- They are distinct from the older reviewer-control paired runs in `comment_1.9_paired/`, which use a different `model.inp` and a different `photo_processes.dat`.
