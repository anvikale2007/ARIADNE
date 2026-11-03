# Does removing more NETWORK-abnormal tissue predict seizure freedom?

**Pre-registered** in PREREGISTRATION.md Amendment 4 (Analysis A), written before the analysis was run.

This is the test Amendment 3 could not do. Amendment 3 refuted the *morphometric* version (cortical
thickness, area, volume; n=405). This tests the **network** version, which is what Sinha (2021) used and
what ARIADNE actually computes: each patient's regional nodal strength is z-scored against **97 control
connectomes**, after the same normalisation ARIADNE applies internally.

**n = 182** (127 seizure-free, 55 not seizure-free).

| Test | Seizure-free | Not seizure-free | AUC | p |
|---|---|---|---|---|
| H4 network abnormality of resected tissue (PRIMARY) | 0.925 | 0.870 | **0.535** | 0.226 |
| H5 whole-brain network abnormality (higher = worse) | 0.825 | 0.821 | **0.523** | 0.31 |
| confound: resected volume alone | 3.677 | 4.014 | **0.437** | 0.911 |
| secondary: total network abnormality removed | 3.415 | 3.602 | **0.474** | 0.715 |
| H4 restricted to temporal surgery (n=142) | 0.909 | 0.911 | **0.499** | 0.51 |

Adjusted logistic regression (density, controlling for resected volume and operated side):
β = +0.512, p = 0.30. Significance threshold was pre-set at p < 0.025 (two analyses, Bonferroni).

## Verdict: not supported

The primary measure points in the hypothesised direction — seizure-free patients had slightly more
network-abnormal tissue removed (0.925 vs 0.870) — but AUC 0.535 at p = 0.23 is far from the threshold,
and restricting to temporal surgery gives AUC 0.499, which is exactly chance.

**Both versions of the abnormality hypothesis now fail at adequate power:** morphometric (n = 405,
AUC 0.476-0.499) and network (n = 182, AUC 0.535). Reported as pre-registered negative results.

## What this does not claim

It does not refute Sinha (2021). The abnormality measure used here is deliberately simple — a z-score of
nodal strength against controls — whereas that work used a more developed patient-specific metric. What
is refuted is the simple version, which is the one ARIADNE would have adopted.

The confound is also worth noting: resected volume alone gives AUC 0.437, i.e. patients who did worse had
*larger* resections, consistent with more extensive disease rather than with the surgery itself failing.
