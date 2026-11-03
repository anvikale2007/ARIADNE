# Virtual surgery vs real outcome (IDEAS II) - pre-registered analysis
Plan fixed in PREREGISTRATION.md before this analysis was run.

Patients: 6 (seizure-free 3, not seizure-free 3)
Per-patient calibration brought pre-surgery spread capacity to 0.54-0.69 (pilot was 0.93-1.00, saturated).

**PRIMARY: reduction in spread capacity (graded resection)**
- mean: seizure-free 0.273, not seizure-free 0.355
- Mann-Whitney U=4, one-sided p=0.65, **AUC=0.44**
- Pre-registered decision rule: H1 not supported.

**SENSITIVITY: same, original binary resection model**
- mean: seizure-free 0.081, not seizure-free 0.143
- Mann-Whitney U=1, one-sided p=0.95, **AUC=0.11**
- Pre-registered decision rule: H1 not supported.

**Power.** Smallest attainable one-sided p at n=3 vs 3 is 0.050. 80% power needs about 15 per group at AUC 0.80, 25 per group at AUC 0.70 (simulated, see PREREGISTRATION.md). This sample is a pilot and cannot confirm the hypothesis.

**Confounds (reported regardless of significance)**

| Variable | Seizure-free | Not seizure-free |
|---|---|---|
| Regions >=50% resected | 3.000 | 4.000 |
| Total resected fraction | 3.292 | 4.714 |
| Pre-surgery capacity | 0.590 | 0.631 |
| Calibrated excitability | -2.121 | -2.116 |

| Patient | ILAE | Surgery | Pre | Post | Reduction |
|---|---|---|---|---|---|
| 1 | 4 | T Lx | 0.663 | 0.290 | 0.372 |
| 2 | 1 | T Lx | 0.590 | 0.206 | 0.384 |
| 3 | 1 | T Lesx | 0.542 | 0.470 | 0.072 |
| 7 | 1 | T Lx | 0.638 | 0.276 | 0.363 |
| 143 | 4 | T Lx | 0.542 | 0.231 | 0.311 |
| 145 | 3 | T Lesx | 0.687 | 0.305 | 0.381 |

Exploratory beyond the primary test. Model parameters were never tuned to outcomes.
The earlier pilot (fixed excitability, binary resection, n=6) gave AUC 0.33 and is reported alongside.
## Exploratory (added after seeing the result; see PREREGISTRATION.md Amendment 1)

The confound table above shows that patients who did **not** become seizure-free had larger resections
(total resected fraction 4.71 vs 3.29).
The primary measure therefore partly rewards removing more tissue rather than removing the *right* tissue.
Spearman(resected volume, reduction) = 0.49 (p=0.33).

Dividing the reduction by the total resected fraction gives seizure-free 0.356 vs not
seizure-free 0.077, AUC 0.78, p=0.20 - the hypothesised direction.

**This is not a finding.** It is driven almost entirely by one patient (sub-3: resected fraction 0.079,
reduction 0.072, ratio 0.915). With three patients per group it is noise. It has been pre-specified as a
co-primary outcome for the next adequately powered run, at alpha 0.025, and is claimed as nothing more.
