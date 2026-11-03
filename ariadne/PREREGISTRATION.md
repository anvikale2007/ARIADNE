# ARIADNE: pre-registered analysis plan

**Written 31 October 2026, 02:40 IST (Asia/Kolkata), before the analysis described here was run.**
Author: Anvi Kale. The purpose of this document is to fix the hypothesis, the analysis and the decision
rule in advance, so that the result reported on 3 November cannot have been chosen after seeing the data.

Nothing in this document may be edited after the first run of `virtual_surgery2.py`. Changes go in a
dated amendment section at the end.

---

## 1. Background and what is already known

ARIADNE simulates how a focal seizure spreads through a patient's structural connectome, and asks how
much of each antiseizure drug is needed to stop that spread. The drug ranking itself **cannot currently
be validated**: no public dataset links a brain scan to which drug worked in that patient. This is a
property of the field, not of this project.

The one prediction that *can* be tested with public data is about surgery. IDEAS II (OpenNeuro
ds007401) provides pre-surgical connectomes, the tissue each patient actually had resected, and ILAE
outcome at one year. If ARIADNE's network model captures something real about seizure generation, then
virtually removing the resected tissue should weaken the network more in patients who actually became
seizure-free.

**Pilot already run (30 October 2026, n = 6, reported in `results/virtual_surgery.csv`):**
mean reduction in spread capacity 0.102 in seizure-free patients versus 0.153 in not-seizure-free
patients; AUC 0.33; one-sided p = 0.8. The pilot therefore pointed in the *opposite* direction to the
hypothesis.

**The pilot could not have succeeded.** With 3 patients per group, the smallest one-sided
Mann-Whitney p-value obtainable is 1/C(6,3) = 1/20 = 0.05, which is not below alpha = 0.05 for any
possible arrangement of the data. Statistical power at AUC 0.80 was therefore **0.00**. The pilot was
uninformative by construction, and its AUC of 0.33 is not evidence against the hypothesis.

Simulated power (Mann-Whitney, one-sided, alpha = 0.05, normal shift model, 1500 replicates):

| True AUC | n per group for 80% power |
|---|---|
| 0.70 | 25 |
| 0.75 | 20 |
| 0.80 | 15 |
| 0.85 | 10 |

---

## 2. Two known defects in the pilot, and the fixes

Both were identified from the pilot's internal diagnostics, **without reference to outcome labels**.

**Defect 1 — ceiling effect.** Every patient was simulated at the same fixed background excitability
(`X0H = -2.10`). At that setting, pre-surgery spread capacity was 0.93–1.00 for all six patients: the
networks were already saturated, so removing tissue had almost nowhere to move the measurement.

*Fix.* Calibrate the background excitability **per patient**, using only the pre-surgery connectome, so
that each patient's pre-surgery spread capacity sits near the middle of the measurable range (target
0.50). The resection is not used in the calibration, and the outcome label is not used anywhere, so
the test remains non-circular.

**Defect 2 — partial resections under-counted.** Resection was modelled by scaling a region's
connections by (1 − fraction resected) and silencing only regions that were at least 95% removed. A
region 60% removed therefore kept full seizure-generating capacity.

*Fix.* Graded resection: a region's excitability moves toward the silenced value in proportion to the
fraction resected, `x0 = x0_base + f·(x0_removed − x0_base)`, alongside the existing connection scaling.

---

## 3. Hypothesis

**H1 (directional).** In patients who became seizure-free (ILAE 1–2 at year 1), virtual resection of
the tissue actually removed produces a **larger** reduction in network spread capacity than in patients
who did not become seizure-free (ILAE 3–5).

**H0.** The reduction does not differ between outcome groups.

---

## 4. Analysis plan, fixed in advance

- **Data.** Every IDEAS II subject with a `36-Lausanne_Count` connectome, an entry in
  `table_resected.csv`, and a non-missing `ILAE_Year1`. No subject is excluded on any other basis. If a
  connectome fails to load or calibration does not converge, the subject is reported as excluded, with
  the reason and the count stated in the result.
- **Groups.** Seizure-free = ILAE 1–2. Not seizure-free = ILAE 3–5.
- **Primary outcome measure.** Reduction in spread capacity = capacity(pre) − capacity(post), where
  capacity is the mean fraction of the network recruited over seizures seeded in each region of the
  operated hemisphere in turn.
- **Primary test.** One-sided Mann-Whitney U (seizure-free > not seizure-free), alpha = 0.05, reported
  with AUC = U / (n₁·n₂) as the effect size.
- **Calibration target.** Per-patient bisection on background excitability to bring pre-surgery
  capacity closest to 0.50, evaluated on a fixed pseudo-random subsample of 12 onset regions
  (seed 0). Final pre- and post-surgery capacities are then measured on the full operated hemisphere.
- **Resection model.** Graded (defect 2 fix) as the primary analysis. The original binary model is
  reported alongside as a sensitivity analysis.
- **Confounds reported whether or not they are significant.** Number of regions at least 50% resected,
  operated hemisphere, surgery type, pathology, and calibrated excitability, each compared between
  outcome groups.

## 5. Decision rule

- **AUC > 0.5 with p < 0.05** → H1 supported at this sample size.
- **AUC > 0.5 with p ≥ 0.05** → direction consistent with H1, not statistically supported; reported as
  inconclusive, with the sample size needed stated.
- **AUC ≤ 0.5** → H1 not supported. Reported as a negative result.

## 6. What would falsify the underlying idea

If, at adequate sample size (n ≥ 15 per group by the table above), AUC is at or below 0.5, then
ARIADNE's spread-capacity measure carries no information about surgical outcome, and the network model
as currently formulated should not be used to reason about which tissue matters. That would be a real
finding and will be reported as such.

## 7. Commitments

1. The result will be reported whichever way it comes out, including a null or a reversed result.
2. The pilot's AUC of 0.33 will continue to be reported alongside the new analysis.
3. The sample is expected to remain under-powered (roughly 10 subjects are currently available locally).
   The analysis is therefore labelled a **pilot**, and no claim of confirmation will be made from it
   regardless of the p-value obtained.
4. The primary analysis is the one specified above. Any other analysis performed is exploratory and
   will be labelled exploratory.

---

## 8. Separate pre-registered prediction: the drug ranking

This cannot be tested with available data, and is recorded here so that it is falsifiable later.

**Prediction.** In a cohort of adults with drug-resistant focal epilepsy who have diffusion MRI, a
localised seizure-onset zone, and a documented response to a specific antiseizure drug, ARIADNE's
predicted drug order will correlate with observed response better than chance (pre-specified: AUC > 0.5
for distinguishing responders from non-responders to the drug ARIADNE ranks first).

**Dataset that would test it.** A retrospective epilepsy-surgery or drug-resistance cohort with
pre-treatment diffusion MRI and per-drug response recorded; n ≥ 15 responders and 15 non-responders by
the power table above. This requires ethics approval and does not exist publicly as of 31 October 2026.

**What would falsify it.** AUC ≤ 0.5 in such a cohort.

**Known reason it might fail even if the mechanism is right.** ARIADNE models efficacy only. Real drug
choice and real discontinuation are driven substantially by tolerability, which ARIADNE does not
represent at all.

---

## Amendment 1 - 31 October 2026, 04:35 IST, written AFTER seeing the pilot result

This amendment is dated and kept separate because it was written after the data were seen. Nothing in
sections 1-8 was changed.

**What the pre-registered analysis returned.** n = 6. Per-patient calibration worked: pre-surgery spread
capacity fell from the pilot's saturated 0.93-1.00 to 0.54-0.69, so the ceiling effect is fixed. Primary
outcome AUC 0.44 (one-sided p = 0.65); binary-resection sensitivity AUC 0.11. By the decision rule in
section 5, **H1 is not supported** at this sample size. This is reported as a negative result.

**What the confound table shows.** Patients who did *not* become seizure-free had larger resections
(4.0 vs 3.0 regions at least 50% removed; total resected fraction 4.71 vs 3.29). The primary measure
rewards removing more tissue, so it partly measures resection size rather than whether the *right*
tissue was removed. Spearman correlation between resected volume and reduction is 0.49 (p = 0.33) in
this sample - suggestive, not established.

**Exploratory, not a finding.** Dividing the reduction by the total resected fraction gives AUC 0.78
(p = 0.20) in the hypothesised direction. **This is driven almost entirely by one patient** (sub-3,
resected fraction 0.079, reduction 0.072, ratio 0.915) and with three patients per group it is noise.
It is recorded here as a hypothesis, not a result, and it is not claimed anywhere as support for H1.

**Pre-specified for the next, adequately powered run.** When n >= 15 per group is available:

- **Primary (unchanged):** reduction in spread capacity, graded resection, one-sided Mann-Whitney.
- **Co-primary, added here:** reduction divided by total resected fraction - "how much seizure-generating
  capacity was removed per unit of tissue removed". Same test, same alpha. Because two primaries are now
  specified, each is tested at alpha = 0.025 (Bonferroni).
- Resected volume will be reported as a covariate for both.

This amendment is itself falsifiable: if the normalised measure gives AUC <= 0.5 at adequate sample
size, the idea that ARIADNE identifies *which* tissue matters is not supported, and will be reported
as such.

---

## Amendment 2 — 31 October 2026, prior-art acknowledgement

**This analysis is not novel methodology, and should not be presented as such.**

Nishant Sinha's PhD thesis (*Predicting seizure spread and neurosurgical outcomes in epilepsy by combining
neuroimaging, machine learning, and computer modelling*, Newcastle University, July 2021, supervised by
Peter N. Taylor) addresses the same question — predicting post-surgical outcome from structural brain
networks — in substantially the same patient population. Taylor's group released the IDEAS dataset used
here, and the thesis co-authors (de Tisi, Vos, Duncan, Winston) are the UCL/Queen Square surgical cohort
that IDEAS comes from.

Key differences in method, and what ARIADNE should borrow:

- Sinha standardises each patient's structural network **against healthy controls** to derive a regional
  abnormality score, rather than using raw network measures. ARIADNE currently uses raw spread capacity.
  That is very likely part of why the measure here tracks resection volume (Amendment 1): a
  control-normalised measure would quantify how abnormal the removed tissue was, not simply how much of
  it there was. One IDEAS control (sub-4097) is already in `data/ideas/`; a proper normalisation needs a
  control group.
- Sinha reports that structural network abnormality **was** associated with post-surgical seizure
  recurrence. ARIADNE's n=6 negative result should therefore be read as a likely failure of this
  implementation at this sample size, not as evidence against the underlying premise.

**Positioning.** Surgical-outcome prediction from structural networks is an established field with
published positive results. ARIADNE's virtual surgery is a **validation test borrowed from that
literature**, used to check whether the network model carries real information — it is not ARIADNE's
contribution. The contribution is the drug layer: translating published drug potencies into predicted
effects on seizure spread. No part of the prior work above models drugs.

---

## Amendment 3 — 31 October 2026, ~13:00 IST. Written BEFORE running this analysis.

New data arrived: IDEAS z-scored FreeSurfer measures (figshare), giving regional cortical **thickness, surface
area and volume z-scores against 100 normative controls** for 428 patients. Combined with the resection table
and outcomes already held, **415 patients** have all three. This is the first adequately powered sample in the
project (292/135 split overall; 80% power needs ~15–25 per group).

**This analysis does not use the ARIADNE model.** No connectome, no Epileptor, no drugs. It tests the
*principle* behind the fix proposed in Amendment 1 — that what matters is how abnormal the removed tissue was,
not how much of it there was — on a properly powered sample, before that principle is built into ARIADNE's
virtual surgery. It is methodologically close to Sinha (2021) (Amendment 2) and is a replication, not a novel
method.

### Measures (fixed now)

For patient *i*, over the 68 cortical Desikan-Killiany regions (`ctx-lh-X` ↔ `lh_X`):

- `vol_resected` = Σ_r f_ir  (f = fraction of region r resected) — the existing confounded measure
- **`abnormality_removed` = Σ_r f_ir · |z_ir|** — total abnormality removed
- **`abnormality_density` = Σ_r f_ir · |z_ir| / Σ_r f_ir** — mean abnormality of the tissue removed (PRIMARY)
- `global_abnormality` = mean_r |z_ir| over all regions — whole-brain disease burden

Primary z-metric: **cortical thickness** (the most established structural marker in epilepsy). Surface area and
volume are secondary and reported whatever they show.

### Hypotheses

- **H1 (primary).** Seizure-free patients (ILAE 1–2) have higher `abnormality_density` than not-seizure-free
  (ILAE 3–5). One-sided Mann-Whitney, α = 0.05, effect size AUC.
- **H2 (control direction, from Sinha 2021).** Higher `global_abnormality` is associated with *not* being
  seizure-free. One-sided in that direction. H1 and H2 are not in conflict: one is about targeting, the other
  about overall disease burden.
- **H3 (confound).** `vol_resected` alone does not separate the groups (if it does, it must be adjusted for).

### Adjustment, fixed now

Logistic regression of seizure freedom on `abnormality_density`, with `vol_resected` and operated hemisphere as
covariates, reported alongside the unadjusted test. Most operations in this cohort are temporal, so a
temporal-only sensitivity analysis is also reported.

### Decision rule

- AUC > 0.5 and p < 0.05 for H1 → the abnormality-weighting principle is supported in this cohort, and it
  should be built into ARIADNE's virtual-surgery measure.
- AUC ≤ 0.5, or not significant at n = 415 → the principle is not supported, Amendment 1's proposed co-primary
  is withdrawn, and that is reported as a negative result.

Three hypotheses are tested, so each is read at α = 0.0167 (Bonferroni) for any claim of significance; the
uncorrected p-values are reported too.

### Amendment 3 — RESULT (run 31 Oct 2026, ~13:05 IST, immediately after the plan above)

n = 405 (275 seizure-free, 130 not), 68 cortical regions matched.

| Hypothesis | Thickness | Area | Volume |
|---|---|---|---|
| **H1** density of abnormality in resected tissue (higher → seizure-free) | AUC 0.476, p=0.79 | AUC 0.499, p=0.51 | AUC 0.484, p=0.70 |
| **H2** whole-brain abnormality (higher → not seizure-free) | AUC 0.559, p=0.027 | AUC 0.507, p=0.41 | AUC 0.531, p=0.16 |
| **H3** resected volume alone (confound check) | AUC 0.465, p=0.87 | — | — |

Adjusted logistic regression (density, controlling for resected volume and operated side): β = −0.277,
p = 0.15 for thickness; not significant for any metric. Temporal-surgery-only sensitivity (n = 328): AUC
0.512, p = 0.36.

**Verdict by the rule fixed above: H1 is NOT supported, at adequate power.** AUC is at or below chance for
all three z-metrics. **Amendment 1's proposed co-primary outcome (reduction normalised by resected fraction)
is therefore WITHDRAWN.** The exploratory AUC 0.78 seen in the n = 6 connectome pilot was noise, as was
stated at the time.

**H2** is directionally consistent with Sinha (2021) — more whole-brain cortical thinning, worse outcome —
but at p = 0.027 it does not clear the pre-specified Bonferroni threshold of 0.0167 and is reported as
suggestive only.

**H3 is informative in the other direction:** resected volume alone does *not* separate the outcome groups
here (AUC 0.465). The volume confound that appeared to dominate the n = 6 pilot does not hold at n = 405,
which means that confound was also a small-sample artefact.

**Scope of this negative result.** It tests *morphometric* abnormality — cortical thickness, surface area and
volume. It does **not** test *network* abnormality, which is what Sinha used and what ARIADNE would compute.
So the result rules out the simple morphometric weighting, not the network version. The network version
remains untested and still requires the connectome files.

---

## Amendment 4 — 31 October 2026, ~13:15 IST. Written BEFORE running either analysis below.

The full Lausanne-36 deterministic-tractography release arrived: **288 connectomes — 191 patients and 97
controls**. Of these, **183 patients have a connectome, a resection map and a year-1 outcome** (128
seizure-free, 55 not). This is 30x the pilot and is adequately powered (80% power needs 15–25 per group).

Two analyses are pre-specified. Neither has been run at the time of writing.

### Analysis A — network abnormality (the test Amendment 3 could not do)

Amendment 3 refuted the *morphometric* version of the abnormality hypothesis. It could not test the
*network* version, which is what Sinha (2021) used and what ARIADNE actually computes. The 97 control
connectomes now make that possible.

- For each region r, build a normative distribution of its **nodal strength** across the 97 controls
  (after the same `connectome.prepare` normalisation ARIADNE uses), and z-score each patient's region
  against it: `z_ir`.
- `net_abnormality_density` = Σ_r f_ir·|z_ir| / Σ_r f_ir  — mean network abnormality of resected tissue (PRIMARY)
- `net_global_abnormality` = mean_r |z_ir| — whole-brain network abnormality
- **H4:** seizure-free patients have higher `net_abnormality_density`. One-sided Mann-Whitney, AUC.
- **H5:** higher `net_global_abnormality` associates with *not* seizure-free (Sinha's direction).
- Adjusted for resected volume and operated side by logistic regression, as in Amendment 3.

### Analysis B — virtual surgery, adequately powered (the original H1)

The pre-registered test from sections 3–5, now at n = 183 instead of 6. Method unchanged: per-patient
calibration blind to outcome, graded resection, reduction in spread capacity, one-sided Mann-Whitney, AUC.

**One method change, made for compute cost and fixed now:** the onset sweep uses a fixed pseudo-random
subsample of **15** regions of the operated hemisphere (seed 0) rather than all of them, and **the same
subsample is used before and after resection** so the paired difference cancels most of the sampling noise.
Calibration uses 4 bisection steps on 10 onsets. The binary-resection sensitivity analysis is dropped.
This is a precision trade, not a change of hypothesis: it makes 183 patients feasible in ~3 hours where the
full sweep would take ~21.

### Decision rules

- A: AUC > 0.5 with p < 0.025 (two analyses, Bonferroni) → network abnormality weighting is supported and
  should enter ARIADNE's virtual-surgery measure. AUC ≤ 0.5 → not supported; reported as negative.
- B: AUC > 0.5 with p < 0.025 → H1 supported at adequate power, reversing the pilot's inconclusive result.
  AUC ≤ 0.5 → H1 refuted at adequate power, and that is the project's headline negative result.

Both results will be reported whichever way they come out. Given Amendment 3's outcome, a negative result
for A is the more likely of the two.
