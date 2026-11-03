# ARIADNE value verification (30 Oct 2026)

Each number used by the model, where it was checked, and the result.
VERIFIED = read in the source (or an independent copy of it) this session.
CONSISTENT = the primary source could not be opened, but independent data agree.
NOT VERIFIED = could not be checked; treat as provisional.

## Drug potencies
| Value | Used | Check | Status |
|---|---|---|---|
| CBZ persistent INa IC50 | 77.5 uM | Kahlig 2022 Table 1 + Kahlig AAN 2023 poster (77,490 nM, slope 1.1) | VERIFIED |
| LTG persistent INa IC50 | 78.5 uM | same (78,480 nM, 1.0) | VERIFIED |
| Cenobamate persistent INa IC50 | 71.7 uM | same (71,690 nM, 1.1) | VERIFIED |
| CBZ / LTG / cenobamate inactivated-state (VDB) IC50 | 44.4 / 39.1 / 66.7 uM | same (44,370 / 39,090 / 66,710 nM, slopes 0.9) | VERIFIED |
| CBZ / LTG / cenobamate 10 Hz UDB IC50 | 1418 / 516 / 749 uM | same | VERIFIED |
| Valproate Na block (2% persistent, 18% VDB at 1 mM) | ~1-1.5% at 69 uM | Kahlig 2022; poster says column alignment for valproate is unclear | PARTLY VERIFIED (negligible either way) |
| Cenobamate tonic GABA-A EC50 | 36.6 uM | Sharma et al. Eur J Pharmacol 2020;879:173117 abstract (36.63 uM; IGABA 164 uM; human subtypes 42-194 uM) | VERIFIED |
| Cenobamate 53.1 uM (rat CA3) | not used | quoted by MDPI review; Nakamura 2019 abstract gives no number | SECONDHAND |

## Concentrations and binding
| Value | Used | Check | Status |
|---|---|---|---|
| CBZ range 4-12 mg/L, 75% bound | yes | Patsalos & Berry 2013; Patsalos 2018 | VERIFIED |
| LTG range 2.5-15 mg/L | yes | Patsalos & Berry 2013; Patsalos 2018 | VERIFIED |
| LTG 55% bound (fu 0.45) | yes | FDA lamotrigine label (12.3: ~55% at 1-10 mcg/mL) and Patsalos & Berry 2013 = 55%; Patsalos 2018 = 66% | VERIFIED (2 of 3 sources); 66% kept as sensitivity |
| VPA range 50-100 mg/L | yes | Patsalos 2013, 2018 | VERIFIED |
| VPA free fraction 10% @40, 18.5% @130 mg/L | yes (fu 0.133 at 75 mg/L) | valproic acid label | VERIFIED (formula is arithmetic on these) |
| Cenobamate 60% bound | yes | FDA XCOPRI label 12.3; drugs.com monograph | VERIFIED |
| Cenobamate free Cave 30.1 uM (200 mg/day) | yes | Vernillet 2020 / FDA review not openable | CONSISTENT: patient levels 17.8-20.9 mg/L at ~200-250 mg/day (Steinhoff 2024) give ~27-31 uM free |
| Cenobamate 400 mg/day 54.4 uM | sensitivity only | as above | NOT VERIFIED directly |
| Pregabalin unbound | yes | Lyrica labels (US, EU, AU) | VERIFIED |
| Pregabalin Cave 19.8 uM (300 mg/day) | yes | monograph table not openable; AU label: avg 1.5 and 6.0 ug/mL at 150 and 600 mg/day | CONSISTENT (dose-proportional 300 mg -> ~3.0 ug/mL = 18.8 uM) |
| Pregabalin KD 6.0 nM (a2d-1), 7.2 nM (a2d-2) | context only | Li 2011 table caption; Piechan 2004 abstract gives 5, 7, 9 nM | VERIFIED |
| Pregabalin range | not used | Reimers 2018: 2-6 mg/L; Patsalos 2018: 2-8 mg/L | sources differ |
| Molecular weights | yes | computed from formulas; cenobamate formula/mass from Wikipedia (267.67) | VERIFIED |

## Gene rules
| Rule | Check | Status |
|---|---|---|
| HLA-B*15:02 -> exclude CBZ | CPIC: Phillips et al. Clin Pharmacol Ther 2018;103(4):574-581, doi 10.1002/cpt.1004 | VERIFIED citation and association |
| Dravet / SCN1A LoF -> avoid Na blockers | Wirrell 2016 Can J Neurol Sci statement read; Wirrell 2022 consensus citation verified but its abstract does not mention Na blockers | VERIFIED via 2016 source |
| POLG -> exclude valproate | Depakote label Contraindications; Stewart et al. Hepatology 2010 (OR 23.6, 95% CI 8.4-65.8) | VERIFIED |

## Other
| Item | Check | Status |
|---|---|---|
| CHB-MIT seizure times (2996-3036, 1732-1772, 1015-1066, 1862-1963 s) | independent Zenodo timestamp table | VERIFIED |
| Epileptor-2D equations | identical to TVB code (tests/test_against_tvb.py) | VERIFIED |
| Single-node seizure threshold | simulated: between x0 -2.07 and -2.06 (earlier statement "~-2.05" was slightly off) | VERIFIED numerically |
| "Over 30% of patients are refractory" | stated in Nakamura et al. 2019 abstract | VERIFIED (cite a primary epidemiology paper for the poster) |
| Sanz-Leon 2013 Front Neuroinform 7:10; Garyfallidis 2014 Front Neuroinform 8:8 | DOIs 10.3389/fninf.2013.00010, 10.3389/fninf.2014.00008 | VERIFIED |
| Proix 2014 J Neurosci 34(45):15009 | journal URL path + title | VERIFIED |
| Jirsa 2014 Brain; Jirsa 2017 NeuroImage (VEP) | title, journal, year verified; volume/pages not | PARTLY VERIFIED |

## Modelling assumptions (not facts; no source exists by design)
GABA weight 0.5, NAV_REF 0.30, release reductions (LTG 0.3, VPA 0.2, PGB 0.6), valproate other excitability 0.4, coupling K, effect units,
calibration grid, Desikan-to-TVB label crosswalk. See results/robustness_report.md for which of these change the ranking.

## Added 31 Oct from FDA labels supplied by the user (DailyMed)
| Value | Source | Status |
|---|---|---|
| Phenytoin therapeutic range 10-20 mcg/mL; 'highly protein bound' (no % in this label) | Extended phenytoin sodium capsules label | VERIFIED range; binding % still needed |
| Valproate 50-100 mcg/mL; free fraction ~10% at 40 to 18.5% at 130 mcg/mL | Valproic acid capsules label | VERIFIED (second source) |
| Levetiracetam <10% protein bound; binds SV2A; label: 'molecular significance ... is not understood' | Levetiracetam tablets label 12.1/12.3 | VERIFIED; supports not modelling LEV yet |

## Levetiracetam (added 31 Oct 2026) — derived, not assumed

| Value | Source | Status |
|---|---|---|
| Therapeutic range 12–46 mg/L | Patsalos & Gulcebi 2018 (also reports 3–34 and 20–40 mg/L) | sourced, range disputed |
| Protein binding ~0 ("not plasma protein bound") | Patsalos & Gulcebi 2018 | sourced |
| Free concentration 170.4 µM (mid-range) | computed: 29 mg/L ÷ 170.2 g/mol × fu 1.0 | derived |
| N-type Ca current block ~18% (total HVA), ~37% (isolated N-type) | Lukyanetz, Shkryl & Kostyuk, Epilepsia 2002;43:9–18 (rat CA1, 1–200 µM) | sourced; **no IC50 published**, so the measured block is used directly rather than fitting an occupancy curve |
| Release cooperativity n = 2.6 | Reid, Bekkers & Clements, J Neurosci 1998;18:2849–2855 (Hill coefficient, N-type at rat hippocampal autapses; 2.5 for P/Q; 3.1–3.3 unblocked) | sourced |
| Coupling effect 0.403 | derived: 1 − (1 − 0.18)^2.6 | **derived from two sourced measurements — the only coupling value in ARIADNE that is not assumed** |
| SV2A effect | **none published** | NOT modelled. This makes the figure a LOWER BOUND on levetiracetam's true effect. |

Sensitivity: using the 37% isolated N-type block gives coupling 0.699 and a threshold of 1.91× instead of 3.31×,
which moves levetiracetam from 6th to 5th of six. It does not reach the top three under any sourced assumption.

Known gap: because the source reports a single block value rather than a dose–response curve, the model cannot
represent what a higher levetiracetam dose would do (effect gain 1.0× in the therapeutic-ceiling table).

### Independent cross-check of modelled concentrations (Wright et al. Front Neurol 2013;4:192)

Otoul et al., reported in that review, measured steady-state concentrations in 187 children on these drugs
while adding levetiracetam. Those are real patient levels, measured for a different purpose, so they are an
independent check on the mid-therapeutic-range concentrations ARIADNE assumes:

| Drug | Real patients (geometric mean, µg/mL) | ARIADNE total conc. (mg/L) | Agreement |
|---|---|---|---|
| Carbamazepine | 8.4 (and 8.1 in the control arm) | 8.0 (mid of 4–12) | within 5% |
| Lamotrigine | 8.2 (and 7.7) | 8.75 (mid of 2.5–15) | within 7% |
| Valproate | 83.8 (and 82.5) | 75 (mid of 50–100) | ARIADNE ~11% low |

The same source reports **no pharmacokinetic interaction** between levetiracetam and carbamazepine, valproate,
topiramate or lamotrigine, which matters because ARIADNE models each drug as monotherapy.

Also confirmed there: levetiracetam protein binding <10% (so fu ≥ 0.90; ARIADNE uses 1.0, overstating free
levetiracetam by at most 10%), half-life 6–8 h in adults, licensed dose 1,000–3,000 mg/day.

## The GABA weight is no longer assumed (derived 31 Oct 2026)

Previously `gaba_weight = 0.5`, flagged as an invented constant and the project's single biggest
vulnerability. It is now **derived**.

A region's excitability tracks net synaptic drive. A sodium-channel blocker scales excitatory drive by
(1 − f_Na); a GABA-A potentiator scales inhibitory drive by (1 + g_GABA):

```
net drive  ∝  g_E (1 − f_Na) − g_I (1 + g_GABA)
```

Normalising by resting excitatory conductance, the fractional reduction in net drive is

```
Δ = f_Na + ρ · g_GABA        with  ρ = g_I / g_E
```

This is the functional form ARIADNE already used — but the weight is **ρ, the cortical
inhibitory-to-excitatory conductance ratio**, a measured quantity rather than a free parameter.

| Source | Preparation | g_I/g_E |
|---|---|---|
| Haider et al., J Neurosci 2006 | in vivo, cortical up states | ≈ 1 (balanced) |
| Rudolph et al. 2007 | awake cat cortex | inhibition "several times" excitatory |
| Haider et al. 2013 | awake mouse | balance "profoundly skewed towards inhibition" |

**ρ = 1.0 adopted** as the most conservative sourced value; 1–4 as the sensitivity range. The previous
0.5 was *below* the lowest measured ratio. The ranking flip point is 0.05 — **20× below the lowest
plausible value** — so the cenobamate-first result is far more robust than the old "we assumed 0.5"
framing suggested.

Effect: cenobamate's excitability effect rises from 1.789 to 2.54, threshold 0.475 → 0.335. Order unchanged.

## What is still arbitrary, and exactly how much it matters

After this, the scaling constants `NAV_REF`, `X0_UNIT` and `COUPLING_UNIT` reduce to **one** number that
can affect the ranking: the ratio of the coupling scaling to the excitability scaling
(`COUPLING_UNIT / X0_UNIT`, currently 15). Everything else rescales all drugs equally and cannot reorder them.

Sweeping that ratio over a 250-fold range:

| Ratio (× current) | Ranking |
|---|---|
| 0.1× – 0.5× | cenobamate > lamotrigine > carbamazepine > valproate … |
| **1× (current)** | **cenobamate > lamotrigine > valproate > carbamazepine > pregabalin ≈ phenytoin > levetiracetam** |
| 2× – 4× | cenobamate > lamotrigine > pregabalin … |
| 10× – 25× | pregabalin takes first |

**Cenobamate ranks first across a 40-fold range of the only remaining arbitrary parameter, and
lamotrigine is second across the entire 250-fold range tested.**
