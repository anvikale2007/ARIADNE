# Clinical benchmarks for ARIADNE (read in full by Claude from user-supplied documents unless marked)

## Nevitt et al. 2022, Cochrane CD011412 (IPD network meta-analysis, monotherapy, newly diagnosed) - READ
- Focal seizures, time to treatment failure (any reason, includes side effects), HR of lamotrigine vs:
  levetiracetam 1.01 (0.88-1.20), zonisamide 1.18, lacosamide 1.19, carbamazepine 1.26 (1.10-1.44),
  oxcarbazepine 1.30, valproate 1.35 (1.09-1.69), phenytoin 1.44 (1.11-1.85), topiramate 1.50, gabapentin 1.53, phenobarbitone 1.97.
- Focal, remission (efficacy) outcomes: no differences between lamotrigine and any AED; carbamazepine better than valproate (6-month remission).
- Focal, time to first seizure (efficacy): phenytoin and phenobarbitone performed better than others.
- Generalised: valproate best profile; no difference vs lamotrigine or levetiracetam.

## WHO mhGAP evidence profile EPI 5 (2015), drug-resistant convulsive epilepsy (add-on) - READ
- Lamotrigine add-on vs placebo, >=50% seizure reduction: OR 2.51 (1.86-3.4), 11 RCTs (Ramaratnam 2010).
- No head-to-head trials between essential and newer ASMs in drug-resistant epilepsy.

## Hadady et al. 2022 (EpiPick external validation) - VERIFIED numbers
- Group-1 vs other: retention 79.4 vs 67.2%; seizure freedom 76.0 vs 61.6%; AE discontinuation 12.0 vs 29.2%.

## How ARIADNE compares (honest reading)
| Check | Evidence | ARIADNE (Stanford scan) | Verdict |
|---|---|---|---|
| Lamotrigine effective in drug-resistant focal epilepsy | WHO / Cochrane add-on OR 2.51 | stops spread near nominal dose (0.72) | CONSISTENT (qualitative) |
| Lamotrigine vs carbamazepine vs valproate (focal) | Treatment failure favours LTG, but driven largely by tolerability; on remission there are no LTG differences | LTG 0.72 < VPA 1.53 ~ CBZ 1.59 | PARTLY: ARIADNE shows a larger efficacy gap than efficacy endpoints support; it does not model side effects |
| Phenytoin best for time to first seizure (focal) | Nevitt 2022 | phenytoin not yet modelled | NOT TESTABLE yet |
| Cenobamate, pregabalin | not in these monotherapy reviews | cenobamate 1st, pregabalin last | NOT TESTABLE with these sources (see Zhang 2024 add-on NMA, not yet read) |

Key caveat for the poster: population mismatch. Trials are newly diagnosed monotherapy; ARIADNE simulates drug-resistant focal spread and models efficacy only (no tolerability).

## Marson et al. 2021, SANAD II (focal), Lancet 397:1363-74 - READ (abstract)
- 990 newly diagnosed focal epilepsy patients randomised to lamotrigine / levetiracetam / zonisamide (open-label).
- 12-month remission, ITT: levetiracetam NOT non-inferior to lamotrigine (HR 1.18, 97.5% CI 0.95-1.47); zonisamide non-inferior (1.03, 0.83-1.28).
- Per-protocol: lamotrigine superior to levetiracetam (HR 1.32, 1.05-1.66) and zonisamide (1.37, 1.08-1.73).
- Adverse reactions: lamotrigine 33%, levetiracetam 44%, zonisamide 45%.
- Conclusion: lamotrigine should remain first-line for focal epilepsy.
- Relevance to ARIADNE: supports lamotrigine ranking high among drugs that are modelled; LEV/ZNS not in ARIADNE's panel.

## Xie et al. 2025, Brain 148:4288-4298 (emulated trials from EHR, single US centre) - READ (abstract)
- 1596 patients (seizure freedom, 2 yr), 2945 (retention, 5 yr).
- Focal epilepsy: first-line ASMs achieved SIMILAR seizure freedom; levetiracetam better retention than lamotrigine/oxcarbazepine.
- Third-line (drug-resistant): cenobamate, brivaracetam, clobazam achieved similar retention and seizure freedom in focal epilepsy;
  across all epilepsy types cenobamate had longer retention than brivaracetam but LOWER seizure freedom than clobazam and brivaracetam.
- Relevance to ARIADNE: real-world data do NOT show cenobamate clearly superior for seizure freedom in this cohort.
  This tempers ARIADNE's cenobamate-first result (and Zhang 2024 NMA's favourable cenobamate estimate). Rankings differ by outcome.

## Zhang et al. 2024, Epilepsia Open 9:1550, doi 10.1002/epi4.12997 (NMA of add-on ASMs, focal epilepsy) - READ
RR vs placebo (95% CI), as reported in the text:
| Drug | >=50% responder | 100% responder (seizure freedom) |
|---|---|---|
| Cenobamate | 2.90 (2.04-4.13) | 29.11 (5.20-794.9) |
| Valproate | 3.00 (2.00-4.44) | not reported |
| Pregabalin | 2.82 (2.30-3.70) | 12.16 (1.92-330.2) |
| Lamotrigine | 2.2 (1.64-3.04) | 3.84 (1.71-10.40) |
| Carbamazepine | not reported (no add-on placebo data) | not reported |
- Authors: cenobamate had the greatest likelihood of seizure freedom. CIs are very wide (indirect comparisons).

### Direction check vs ARIADNE (Stanford scan: CNB > LTG > VPA ~ CBZ > PGB)
- Seizure freedom (3 drugs): Zhang CNB > PGB > LTG vs ARIADNE CNB > LTG > PGB. Spearman rho = 0.5 (n=3, meaningless statistically).
- >=50% responder (4 drugs): Zhang VPA > CNB > PGB > LTG vs ARIADNE CNB > LTG > VPA > PGB. Spearman rho = 0.0 (n=4).
- MATCH: cenobamate top for seizure freedom. MISMATCH: ARIADNE puts pregabalin last; trials put it at or above lamotrigine.
  Likely cause: pregabalin is modelled only through an ASSUMED release reduction (0.6), with no excitability effect -> probably underrated.
- Bottom line: no meaningful agreement can be claimed; overlap is 3-4 drugs and CIs overlap heavily.

## Löscher & Klein 2026, Epilepsia (doi 10.1002/epi.70210) - READ (abstract + Table 2)
- METHOD PRECEDENT: ranked 23 ASMs by potency in rodent seizure models and correlated the ranks with clinical seizure-freedom
  rates in add-on RCTs in drug-resistant focal epilepsy. Single models mostly failed; combinations of 2-4 focal models correlated
  significantly. ARIADNE can be evaluated the same way (model rank vs clinical rank).
- Table 2, placebo-subtracted seizure freedom at max dose (add-on, drug-resistant focal):
  cenobamate 20%, levetiracetam 7.3%, vigabatrin 6.7%, topiramate 5.4%, lacosamide 3.9%, lamotrigine 1.6%, pregabalin 1.5%,
  zonisamide 1.4%, gabapentin 0.8% (valproate, carbamazepine not included).
- Direction check vs ARIADNE (CNB > LTG > PGB): MATCHES on placebo-subtracted seizure freedom (20 >> 1.6 ~ 1.5).
  On the NMA risk ratio (Zhang 2024) pregabalin is above lamotrigine, so the LTG-vs-PGB order depends on the metric (effectively a tie).
