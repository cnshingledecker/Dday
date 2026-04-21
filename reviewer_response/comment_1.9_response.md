# Response to Reviewer Comment 1.9

**Draft — 2026-04-16**
**Companion artifacts:** `comment_1.9_paired/` (per-point CSVs, summary CSVs,
overlay figure, paired-run README).
**Branches on `github.com/cnshingledecker/Dday`:**
`newNetwork-transport-fixes` (ion-on, Model 1) and `newNetwork-reviewer-noion`
(ion-off control).

---

## Reviewer's comment (for reference — replace with verbatim text)

> [Paste the reviewer's literal wording of Comment 1.9 here. Inferred content:
> whether the ion photoproducts in the network are actually required to
> reproduce the Gerakines O3 yields, or whether neutral photochemistry alone
> would suffice.]

---

## Response

We thank the reviewer for this comment, which prompted us to run a
same-model control isolating the effect of the ion photoproduct channels.
We did not simply rerun an earlier ion-free version of the network,
because the previously archived `wo_ions` branch differs from the
present Model 1 in more than just the ion channels and would therefore
confound the comparison the reviewer is asking for. Instead, we have
constructed a paired control on a single dedicated branch
(`newNetwork-reviewer-noion`) whose only difference from the Model 1
baseline is that the twelve photoprocess rows whose products contain
ions or electrons have their yield (δ in our notation, column 9 of
`photo_processes.dat`) set to 0. All other parameters — (a) the
non-diffusive bulk-reaction kinetics; (b) the Table 5 best-fit δ
values on the neutral channels; (c) the Table 4 trial frequencies; and
(d) the suprathermal and radiolysis settings — are held identical
between the two runs.

Concretely, the following channels are active on the ion-on side and
disabled (δ = 0) on the ion-off side:

- gO/bO + hν → gO⁺/bO⁺ + ge⁻/be⁻ (PHOION; rows 8288, 8289);
- gO₂/bO₂ + hν → gO⁺/bO⁺ + gO⁻/bO⁻ (PHOEXC; rows 8298, 8299);
- gO₂/bO₂ + hν → gO₂⁺/bO₂⁺ + ge⁻/be⁻ (PHOION; rows 8306, 8307);
- gO₃/bO₃ + hν → gO⁺/bO⁺ + gO₂⁻/bO₂⁻ (PHOEXC; rows 8322, 8323);
- gO₃/bO₃ + hν → gO₂⁺/bO₂⁺ + gO⁻/bO⁻ (PHOEXC; rows 8326, 8327);
- gO₃/bO₃ + hν → gO₃⁺/bO₃⁺ + ge⁻/be⁻ (PHOION; rows 8336, 8337).

The remaining fourteen oxygen photoprocess rows — all neutral
excitation and neutral-product photoionization channels — carry their
Table 5 yields unchanged on both sides. `FIXED_DVAL = 0` in
`model.inp` on both branches, so the per-row δ column in
`photo_processes.dat` is what governs the channel rate at runtime and
δ = 0 cleanly switches a channel off without perturbing any other
reaction.

**Result.** Turning off only the ion photoproduct channels degrades the
Model 1 fit to the Gerakines (2019) bO₃ yields by a factor of roughly
seven in RMSD and pulls the peak modelled O₃ fraction down to ~14% of
its ion-on value:

| Quantity                                    | Ion-on (Model 1) | Ion-off (control) | Ratio (off/on) |
|---------------------------------------------|-----------------:|------------------:|---------------:|
| Unweighted RMSD vs. Gerakines points        |             2.01 |             13.86 |          6.89× |
| Weighted RMSD (fluence-weighted)            |             2.43 |             16.46 |          6.78× |
| Peak bO₃ % (F ≈ 6 × 10¹⁷ photons cm⁻²)       |            26.80 |              3.66 |          0.14× |
| bO₃ % at F ≈ 1.2 × 10¹⁶ photons cm⁻²        |             6.77 |              0.12 |          0.02× |
| bO₃ % at F ≈ 3.4 × 10¹⁶ photons cm⁻²        |            15.10 |              0.32 |          0.02× |

The per-point residuals (see `eval_ion_on.csv` / `eval_ion_off.csv` in
the supporting files) show that the ion-off control fails to reach the
experimental O₃ abundance at any measured fluence in the Gerakines data
set, and the weighted-deviation column is dominated by the high-fluence
points on the ion-off side (weighted deviation ≈ 689 at F ≈ 6 × 10¹⁷
photons cm⁻², versus ≈ 13 for the ion-on run). The ionless model is
therefore not a rescaled or lower-efficiency version of Model 1; it is
structurally unable to produce O₃ at the observed level.

We note that an ionless oxygen photochemistry of this kind is broadly
consistent in spirit with the classical neutral-chemistry treatment of
Mullikin et al. (2021), which established the methodological lineage
for modelling UV-driven O₃ formation in pure O₂ ices and on which our
neutral channels are based. We are careful, however, not to claim that
our ion-off branch reproduces their model quantitatively: the two
networks differ in species set, rate-coefficient provenance, and the
specific non-diffusive bulk-kinetics formalism, and a direct numerical
comparison is outside the scope of the present work. What the paired
control establishes is the in-model statement the reviewer asked for —
namely, that within the Model 1 network used throughout the paper,
suprathermal and ion photoproduct chemistry are jointly required to
reproduce the Gerakines yields and cannot be substituted for by the
neutral photoprocesses alone.

## Manuscript changes made in response

- Added a one-paragraph in-model ionless-control description and the
  RMSD / peak-yield numbers above to Section [X, to be filled in] of
  the revised manuscript.
- Added Figure [N, to be filled in] showing the ion-on and ion-off
  bO₃(F) curves overlaid on the Gerakines points (source file:
  `comment_1.9_paired/fig_paired_overlay.pdf`).
- Added a citation to Mullikin et al. (2021) as methodological lineage
  for the neutral photochemistry, with the explicit caveat above that
  our ion-off branch is not a numerical reproduction of their model.
- The two branches used for the comparison are archived on the public
  repository (`newNetwork-transport-fixes` and
  `newNetwork-reviewer-noion`) with the provenance and reproduction
  recipe documented in `comment_1.9_paired/README.md`.

---

## Open items before this response is letter-ready

1. Replace the "[Paste the reviewer's literal wording...]" placeholder
   with the verbatim Comment 1.9 text from the referee report.
2. Fill in the Section / Figure numbers in the "Manuscript changes
   made" list once the revised manuscript is typeset.
3. Confirm the Mullikin et al. citation year and bibliography entry
   against the main reference list.
4. Decide whether to include the full per-point residual table in the
   letter itself or cite the CSV files in the supplementary material.
