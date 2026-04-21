# Resubmission Handoff

Date: 2026-04-20

Primary sources used for this handoff:
- Revised manuscript PDF: `/Users/cnshingledecker/Downloads/Ion_Ice_Paper-3.pdf`
- Comment 1.9 response draft: [comment_1.9_response.md](/Users/cnshingledecker/Research/Dday/reviewer_response/comment_1.9_response.md)
- Paired ion-on / ion-off artifacts: [comment_1.9_paired](/Users/cnshingledecker/Research/Dday/reviewer_response/comment_1.9_paired)
- Reproducibility sweep: [reproducibility_sweep_2026-04-16](/Users/cnshingledecker/Research/Dday/reviewer_response/reproducibility_sweep_2026-04-16)

## Executive Summary

The revised manuscript is scientifically stronger than the earlier drafts, especially in the route-analysis sections. The biggest remaining issue is that the body of the paper and the Conclusions are not yet saying exactly the same thing. In the body, Model 1 clearly shows a fluence-dependent transition in which ionic chemistry becomes mechanistically important. In the Conclusions, that gets softened too far into a claim that neutral networks are "adequate" for `O3` here.

The best current synthesis is:
- neutral channels can still dominate the gross terminal formation of `O3` at some fluences;
- ionic channels remain kinetically load-bearing for reproducing the observed `O3` trajectory in Model 1;
- therefore, "ions matter" does not require that most `O3` be formed directly in ion-tagged reactions.

This should be made explicit in the revised manuscript and in the reply to reviewer Comment 1.9.

## Highest-Priority Fixes Before Resubmission

### 1. Real internal inconsistency in Section 3.2 / Table 4

Section 3.2 currently says that Models 2-28 used low / medium / high trial frequencies of:
- `1.00 x 10^8 s^-1`
- `1.00 x 10^10 s^-1`
- `1.00 x 10^13 s^-1`

But Table 4 in the revised PDF shows:
- `1.00 x 10^11 s^-1`
- `1.00 x 10^13 s^-1`
- `1.00 x 10^15 s^-1`

This is not just a typo-level issue. It is a real internal inconsistency and should be corrected from the actual source-of-truth model set before resubmission.

Suggested repair:

```text
In Models 2-28, all delta values were set to unity, while nu_n-n, nu_i-n, and nu_i-i were tested using all possible combinations of three different frequency values: low (1.00 x 10^11 s^-1), medium (1.00 x 10^13 s^-1), and high (1.00 x 10^15 s^-1).
```

Only use those numbers if Table 4 is indeed the correct source of truth.

### 2. Conclusions currently understate the ionic result

The manuscript body says:
- `O3` chemistry follows the earlier neutral pattern at early and late fluence;
- at intermediate fluence, the dominant `O3` formation pathway becomes `O3+ + O2 -> O2+ + O3`;
- `O3` destruction also involves ionic channels.

But the Conclusions still say:

```text
Fortuitously, our results suggest that purely neutral grain-chemical networks are adequate for certain species, such as O3 here.
```

That is now too soft given the paired same-model ion-off control and also softer than the manuscript's own route analysis.

Suggested replacement:

```text
Our results show that ionic pathways are not uniformly dominant at all fluences, but they are kinetically important for reproducing the observed O3 growth in Model 1. Neutral channels dominate some terminal O3-forming steps, especially at early and late fluence, whereas ion-neutral and ion-ion processes regulate the mid-fluence chemistry and materially improve agreement with experiment. Thus, neutral-only oxygen chemistry captures the qualitative possibility of ozone formation, but the explicit inclusion of ionic pathways is required to recover the observed abundance trajectory in the present model family.
```

### 3. Make the route-dominance versus sensitivity point explicit

The students' observation is scientifically interesting, not contradictory:
- dominant gross `O3`-forming routes can still be neutral;
- strong dependence on `nu_i-n` and `nu_i-i` can coexist with that if ionic chemistry controls precursor supply, recycling, or destruction balance.

Recommended bridge sentence for either Section 3.3.1 or the Conclusions:

```text
Although the dominant gross O3-forming route is neutral at early and late fluence, the much stronger sensitivity of the model to nu_i-n and nu_i-i than to nu_n-n implies that ionic channels regulate the precursor and recycling chemistry that feeds those neutral terminal steps.
```

This is likely the cleanest mechanism-level interpretation of the current results.

## Suggested Section-Level Edits

### Abstract

Current phrasing is mostly fine, but two lines could be sharpened.

Suggested edits:
- `To this end, we here present` -> `To this end, we present`
- `ion-neutral and ion-molecule reactions` -> consider `ion-neutral and ion-ion reactions`, or just `ionic reactions`, for consistency with the rest of the paper

Possible revised ending:

```text
A comparison of our calculations with previous data suggests that ion-neutral and ion-ion reactions may play a critical role in the chemistry of irradiated interstellar ice.
```

### Section 2.2

The entropy sentence is survivable, but only as a heuristic.

Current line:

```text
The underlying causes for such variation are not, as far as we know, currently understood but could be related to effects such as the entropy of activation.
```

Suggested revision:

```text
The underlying causes of such variation are not yet clear, but may reflect the strong configurational constraints imposed by the condensed-phase environment. In transition-state-theory language, this could be viewed as an effective entropy-of-activation effect, although we do not treat that interpretation as unique here.
```

That keeps the idea without overclaiming.

### Section 3.1

The discussion of the fitted `nu` values is already useful. The paragraph would benefit from one sentence making clear that the ionic sector is kinetically important even though neutral chemistry still contributes strongly to the route network.

Suggested addition after the fitted `nu` discussion:

```text
The resulting parameter hierarchy suggests that the ionic sector is not merely a minor correction to an otherwise sufficient neutral network, but instead materially influences the kinetic pathways required to recover the measured O3 abundance profile.
```

### Section 3.3.1

This section is one of the strongest in the manuscript. It already states the key mechanistic transition well.

I would add one short synthesis sentence after the paragraph describing `N33`:

```text
This transition shows that the reactions carrying the largest gross O3-forming flux need not be the same processes that most strongly control the overall model sensitivity; in Model 1, ionic chemistry appears to regulate the state of the network even when neutral reactions remain important terminal formation steps.
```

### Conclusions

The conclusions should be brought into alignment with both the route analysis and the no-ion control.

Suggested replacement block:

```text
The model described here shows that inclusion of ionic species and related reactions is possible in slightly modified three-phase astrochemical codes, albeit with a substantial increase in network size, complexity, and uncertainty. Our results further indicate that ionic pathways are kinetically important for reproducing the observed O3 growth in Model 1. Although neutral channels remain important, and in some fluence regimes dominate the gross terminal formation of O3, the explicit inclusion of ion-neutral and ion-ion chemistry is required to recover the measured abundance trajectory in the present model family. In this sense, the ionic sector is not simply an optional refinement to the neutral chemistry, but a load-bearing part of the successful kinetics.
```

## Comment 1.9: Recommended Response-Letter Text

The current draft in [comment_1.9_response.md](/Users/cnshingledecker/Research/Dday/reviewer_response/comment_1.9_response.md) is already pretty good. The strongest version should do three things:
- make clear that the test is same-model;
- keep Mullikin as lineage, not as the actual control;
- state the quantitative collapse in the ion-off case directly.

Suggested text:

```text
We thank the reviewer for this suggestion. To address it, we performed a same-model ion-off control in which only the twelve charged-product photoprocess channels were suppressed, while all neutral channels, trial frequencies, and other Model 1 settings were held fixed. This is a cleaner test than reverting to an earlier ion-free network, because it isolates the ionic contribution within the same model framework.

In this control, the agreement with the Gerakines et al. O3 data degrades strongly: the unweighted RMSD increases from 2.01 in the full ion-enabled Model 1 calculation to 13.86 in the ion-off control, while the peak modeled O3 abundance falls from 26.80% to 3.66% of the initial O2 abundance. Thus, neutral chemistry alone remains capable of producing some O3, but is quantitatively insufficient within the present Model 1 network to reproduce the observed abundance trajectory. In this sense, the ionic channels are not a minor refinement to an otherwise sufficient neutral model; they provide kinetically important pathways required to recover both the magnitude and fluence dependence of the measured O3 growth.

We note that Mullikin et al. provides important methodological lineage for the neutral-core chemistry in the present model. However, it is not the operative control for the present reviewer question, because the current Model 1 network and parameterization were developed in an ion-capable framework. For that reason, the same-model ion-off calculation above is the appropriate comparison.
```

## Should the No-Ions Run Be In a Figure?

Recommendation: yes, include it in a figure in the manuscript.

I would not leave the no-ion control only in the response letter or describe it only in prose. This is exactly the kind of reviewer-targeted control that becomes much more persuasive when the curve is visible.

Best option:
- use the existing overlay figure at [fig_paired_overlay.pdf](/Users/cnshingledecker/Research/Dday/reviewer_response/comment_1.9_paired/fig_paired_overlay.pdf)
- place it either:
  - as a new panel attached to Figure 1, or
  - as a short standalone figure immediately after the Model 1 results section

Why a figure is worth it:
- the result is simple and decisive;
- the visual comparison communicates the effect immediately;
- it prevents the reviewer from feeling the control was buried in prose;
- it strengthens the manuscript itself, not just the response letter.

If space is tight, a panel is better than prose-only.

## Suggested Manuscript Paragraph Introducing the No-Ions Control

This could go near the end of the Model 1 discussion or in a dedicated response-driven insertion:

```text
To assess directly the role of ionic channels in Model 1, we also carried out a same-model ion-off control in which only the twelve photoprocess channels yielding charged products were suppressed, while all neutral channels and all other model parameters were held fixed. In this control, the agreement with the Gerakines et al. O3 data deteriorates markedly: the unweighted RMSD increases from 2.01 to 13.86, and the peak modeled O3 abundance decreases from 26.80% to 3.66% of the initial O2 abundance. This result shows that, although neutral chemistry alone can still produce some ozone, it is quantitatively insufficient within the present Model 1 framework to reproduce the observed abundance trajectory without the contribution of the ionic sector.
```

If a figure is added, end that paragraph with:

```text
The ion-on and ion-off O3 curves are compared directly in Figure X.
```

## Copyediting / Typo Pass

These should be fixed before resubmission:

1. `persuavively` -> `persuasively`
2. `yeilded` -> `yielded`
3. `A comparison ... show good agreement` -> `shows good agreement`
4. `the chemistry of neutrals and ions are analyzed` -> `is analyzed` or `are discussed`
5. `the results of the best-fit results of Model 1` -> `the results of the best-fit Model 1`
6. `The model described here show` -> `shows`
7. `codes ,` -> `codes,`
8. `To this end, we here present` -> `To this end, we present`

Potential style-level cleanup:

1. `Terrestrial-like` -> probably lowercase `terrestrial-like`
2. `gas-phase.` versus `gas phase` usage should be made consistent throughout
3. If desired, soften `the first computational study` to `to our knowledge, the first computational study`

## Scientific Takeaway to Preserve

The clean mechanistic reading to preserve across the manuscript and response letter is:
- the neutral Mullikin-like backbone is real and still chemically important;
- the ionic sector is not merely decorative or a small correction;
- dominant gross formation routes and strongest model sensitivities need not be the same thing;
- the most likely interpretation is that ionic chemistry controls precursor supply, recycling, and possibly survival, even when the terminal `O3`-forming step is neutral.

That is a stronger and more interesting result than either of the oversimplified alternatives:
- "ions directly form all the ozone"
- "neutral chemistry alone is sufficient"

## Recommended Next Steps

1. Correct the Section 3.2 / Table 4 frequency inconsistency from the actual source-of-truth files.
2. Revise the Conclusions to align with the route analysis and the paired no-ion control.
3. Add the no-ion control as a figure or figure panel in the manuscript.
4. Insert a short bridging sentence distinguishing route dominance from sensitivity.
5. Apply the typo / grammar fixes listed above.
6. Update the response to Comment 1.9 using the same-model language and quantitative results.
