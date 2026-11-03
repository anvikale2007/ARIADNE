# Does removing more abnormal tissue predict seizure freedom?

**Pre-registered** in PREREGISTRATION.md Amendment 3, written before the analysis was run.
**This analysis does not use the ARIADNE model** - no connectome, no Epileptor, no drugs. It tests the
principle behind the proposed fix to ARIADNE's virtual surgery, on the first adequately powered sample.

Data: IDEAS z-scored FreeSurfer measures (428 patients vs 100 normative controls) + resection percentages
+ ILAE year-1 outcome. **n = 405** with all three (275 seizure-free, 130 not).

## Result: the hypothesis is not supported

| Measure | z-metric | Seizure-free | Not seizure-free | AUC | p |
|---|---|---|---|---|---|
| H1 abnormality density of resected tissue | thickness | 1.028 | 1.112 | **0.476** | 0.787 |
| H2 whole-brain abnormality (higher = worse outcome) | thickness | 0.890 | 0.938 | **0.559** | 0.027 |
| H3 resected volume alone (confound check) | thickness | 2.843 | 2.970 | **0.465** | 0.87 |
| total abnormality removed (secondary) | thickness | 2.921 | 3.206 | **0.491** | 0.621 |
| H1 restricted to temporal surgery | thickness | 1.021 | 1.063 | **0.512** | 0.361 |
| H1 abnormality density of resected tissue | area | 0.995 | 1.012 | **0.499** | 0.515 |
| H2 whole-brain abnormality (higher = worse outcome) | area | 0.870 | 0.900 | **0.507** | 0.412 |
| H3 resected volume alone (confound check) | area | 2.843 | 2.970 | **0.465** | 0.87 |
| total abnormality removed (secondary) | area | 2.806 | 3.088 | **0.481** | 0.732 |
| H1 restricted to temporal surgery | area | 0.985 | 0.987 | **0.511** | 0.375 |
| H1 abnormality density of resected tissue | volume | 0.981 | 1.034 | **0.484** | 0.703 |
| H2 whole-brain abnormality (higher = worse outcome) | volume | 0.884 | 0.901 | **0.531** | 0.158 |
| H3 resected volume alone (confound check) | volume | 2.843 | 2.970 | **0.465** | 0.87 |
| total abnormality removed (secondary) | volume | 2.695 | 3.022 | **0.486** | 0.68 |
| H1 restricted to temporal surgery | volume | 0.931 | 0.966 | **0.518** | 0.306 |

AUC 0.5 is chance. Three hypotheses were pre-specified, so significance requires p < 0.0167.

## What this establishes

1. **Removing more abnormal tissue does not predict seizure freedom** in this cohort, on any of the three
   morphometric measures, at a sample size with ample power. Reported as a pre-registered negative result.
2. **The exploratory AUC 0.78 from the 6-patient pilot was noise**, exactly as stated when it was found.
   The co-primary outcome proposed on the strength of it has been formally withdrawn.
3. **The resected-volume confound was also a small-sample artefact.** At n = 405, resected volume alone does
   not separate the outcome groups (AUC 0.465).
4. Whole-brain cortical thinning is weakly associated with worse outcome (AUC 0.559, p = 0.027), the
   direction Sinha (2021) reports, but it does not clear the corrected threshold here.

## What it does not establish

This tests **morphometric** abnormality (thickness, area, volume). It does not test **network**
abnormality, which is what Sinha used and what ARIADNE computes. The network version is still untested
and still needs the connectome files.
