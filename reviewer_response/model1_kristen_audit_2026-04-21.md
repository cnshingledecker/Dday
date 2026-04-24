Kristen Model 1 photoprocess audit
==================================

Scope
-----
This audit compares the corrected Kristen deck in
`Model_1_kristen_corrected/photo_processes.dat`
against the photoprocess table in
`Ion_Ice_Paper/revisedVersion.tex`.

Result
------
After the 2026-04-21 decimal update to the manuscript table, the 14 manuscript
photoprocess entries match the corrected Kristen deck in the quantities that
matter for the paper:
- column 1: branching fraction `f_br`
- column 2: average cross section `sigma`

Row mapping
-----------
- `P1`: `0001` / `1001` -> `0.50`, `0.00`
- `P2`: `0002` / `1002` -> `0.25`, `3.86e-20`
- `P3`: `0003` / `1003` -> `0.25`, `3.86e-20`
- `P4`: `0004` / `1004` -> `0.167`, `0.00`
- `P5`: `0005` / `1005` -> `0.167`, `0.00`
- `P6`: `0006` / `1006` -> `0.167`, `0.00`
- `P7`: `0007` / `1007` -> `0.50`, `0.00`
- `P8`: `0008` / `1008` -> `0.50`, `3.86e-20`
- `P9`: `0009` / `1009` -> `0.50`, `0.00`
- `P10`: `0010` / `1010` -> `1.00`, `0.00`
- `P11`: `0011` / `1011` -> `0.50`, `2.13e-18`
- `P12`: `0012` / `1012` -> `0.50`, `5.60e-18`
- `P13`: `0013` / `1013` -> `0.50`, `2.13e-18`
- `P14`: `0014` / `1014` -> `0.50`, `5.60e-18`

Branching-sum check
-------------------
- `O`: `R1 + R2 = 1.0`, `R3 + R4 = 1.0`
- `O2`: `R1 + R2 = 1.0`, `R3 + R4 = 1.0`
- `O3`: `R1 + R2 = 1.001`, `R3 + R4 = 1.0`

Interpretation
--------------
The `O3` ionization sum is `1.001` because an equal three-way split of `1/2`
cannot be represented exactly in the compact decimal notation used in the input
deck; each `R1` branch is stored as `1.67E-01`. This preserves the intended
symmetry of the three `R1` channels.
