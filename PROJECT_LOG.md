# ARIADNE: project context for Claude Code

Owner: Anvi (B.Pharm student, NMIMS Mumbai). Entry for the Khoj healthcare innovation fair, **3 Nov 2026**.
She writes and presents the project; be a candid mentor: stress-test claims, never overstate results.

## What this is
Proof-of-concept research software (NOT a clinical or diagnostic tool) that ranks antiseizure drugs by how easily each one
stops a simulated focal seizure from spreading through a brain's structural connectome.
Pipeline: diffusion MRI -> connectome (DIPY) -> network Epileptor-2D seizure model (numba, verified against TVB) ->
calibrate to the brain's own spread threshold -> for each drug, bisection on drug strength until spread stops -> rank.
Optional inputs: scalp EEG (CHB-MIT), intracranial EEG (Penn 3T-7T, BIDS), genetics (gene -> drug safety rules).

## Run
```
pip install -r requirements.txt
python run_pipeline.py                      # Stanford scan connectome
python run_pipeline.py --sensitivity        # + mechanism re-weighting
python robustness.py                        # ranking vs model settings (~25 min on 2 cores)
python verify_model.py                      # 11 internal-consistency checks
python derive_parameters.py                 # drug_sources.yaml -> drugs.yaml (never edit drugs.yaml by hand)
streamlit run app.py                        # demo UI
python tests/test_against_tvb.py ; python tests/test_eeg.py
```
Penn iEEG case and genetics demos: see README.md (needs the Penn starter pack unzipped; path is an argument).

## Key files
- `ariadne/model.py` Epileptor-2D network (single-node threshold x0 between -2.07 and -2.06, verified numerically)
- `ariadne/engine.py` calibration, drug-threshold bisection, ranking, units (X0_UNIT, COUPLING_UNIT, K_COUPLING=-5)
- `drug_sources.yaml` every drug input with status (sourced / assumed) and citation; `derive_parameters.py` turns it into effects
- `genetics_rules.yaml` HLA-B*15:02 -> CBZ, Dravet/SCN1A LoF -> CBZ+LTG, POLG -> VPA
- `ariadne/eeg.py` (scalp, CHB-MIT), `ariadne/ieeg.py` (BIDS iEEG, label crosswalk to TVB template regions)
- `VALUE_VERIFICATION.md` source check of every number; `results/` precomputed outputs + reports

## Current results (30 Oct)
Stanford scan: cenobamate 0.47 > lamotrigine 0.72 > valproate 1.53 ~ carbamazepine 1.59 > pregabalin 2.16 (strength to stop spread; 1.0 = nominal).
Cenobamate top in 5/6 robustness variants; with GABA weight 0, lamotrigine is top. Thresholds resolve to ~0.06-0.13, so near-equal values are ties.

## Open decisions / next tasks
1. Lamotrigine binding RESOLVED: keep 55% (FDA label + Patsalos 2013); 66% (Patsalos 2018) as sensitivity. Phenytoin: range 10-20 mcg/mL verified (label); need binding % source before adding.
2. Still ASSUMED: GABA weight (0.5), NAV_REF (0.30), release reduction (LTG 0.3, VPA 0.2, PGB 0.6), valproate other excitability (0.4). Find quantitative sources.
3. Not openable from the cloud session: Vernillet 2020 steady-state AUCs (cenobamate free Cave 30.1 uM) and the pregabalin monograph PK table: verify from the papers.
4. Run real scalp EEG: download CHB-MIT chb01_03/15/16/26 .edf (PhysioNet) and run with `--eeg-edf` and the `.edf.seizures` files.
5. Real patient connectome: a patient with diffusion MRI + known onset (check full Penn dataset / OpenNeuro). Atlas registration needs a T1.
6. Multi-subject variability (5-10 diffusion scans).

## Honesty rules for this project
- Never call it diagnostic or clinically validated. "Proof-of-concept, hypothesis-generating, focal seizures."
- Every number needs a source or an "assumed" label. Do not fill gaps from memory.
- Healthy-volunteer scans, geometric parcels, template connectome for the Penn case: state these limitations.
- Data credits: Stanford HARDI (Rokem 2015), Sherbrooke 3-shell (Rokem 2017), Penn 3T-7T (Lucas 2026; Imtiaz 2026 Sci Data), CHB-MIT (PhysioNet).

## Update 31 Oct (IDEAS II real-patient data)
- `ariadne/ideas.py`: loads IDEAS II (OpenNeuro ds007401) 82-region Lausanne-36 Count connectomes; region order inferred from table_resected.csv.
- `data/ideas/`: metadata, resection table, Count CSVs for sub-1,2,3,7,13,15,108,143,145,210 (+ control 4097). Need the full batch (PowerShell zip).
- `ariadne/surgery.py` + `virtual_surgery_study.py`: virtual surgery vs ILAE outcome. Pilot (n=6): AUC 0.33, inconclusive.
  Known flaws to fix next: ceiling effect (baseline spread 93-100%; tune per-patient excitability blind to outcome) and partial resections.
- `CLINICAL_BENCHMARKS.md`: Cochrane 2022, SANAD II, WHO 2015, EpiPick 2022, Xie 2025, Zhang 2024, Löscher & Klein 2026 vs ARIADNE.
  Cenobamate-first matches placebo-subtracted seizure freedom (20%); pregabalin-vs-lamotrigine order depends on metric; ARIADNE ignores tolerability.

## Update 31 Oct (clinical gate: ILAE 2025 seizure type)
- `clinical_rules.yaml` + `ariadne/clinical.py`: ILAE 2025 classes (focal / generalized / unknown / unclassified) and types
  (FPC, FIC, FBTC; GTC, absence, myoclonic, tonic, atonic). Each modelled drug gets a guideline status per type
  (first-line / later-line / restricted / not listed / not licensed / avoid) citing NICE NG217 sec 5, NICE TA753, FDA labels.
- Only focal types are simulated; generalized/unknown/unclassified return guideline status only.
- Valproate restriction (NG217 Box 2): age < 55 or able to have children.
- `run_pipeline.py --seizure-type <id> [--age N] [--childbearing]`; app has a Seizure type tab; Explorer step 1.
- `check: true` rules (pregabalin focal add-on, myoclonic first-line wording) came from partly rendered NICE pages: re-verify.
- Rules are UK (NICE). Indian guidance and availability (e.g. cenobamate: CDSCO SEC recommendation Jun 2026 only) may differ.
- `python tests/test_clinical.py`

## Update 31 Oct, overnight (rigour pass: response surface, uncertainty, pre-registration)
Context: fair is **Tue 3 Nov**. Anvi asked to raise the project to PhD-level rigour. The ceiling is NOT effort,
it is that no public dataset links a brain scan to which drug worked. So the work targets *scientific craft*:
claim matched to evidence, limitations quantified not listed, a falsifiable prediction on record.

- `response_surface.py` -> `results/surface_<case>.json`. Maps the stop/no-stop boundary over
  (excitability reduction x coupling scale), 81x21 grid = 1701 sims (~71 min on 2 cores), boundary cells
  re-tested with extra noise seeds. Every drug is a straight ray in this plane, so thresholds become a lookup.
- `ariadne/surface.py` lookup; `ariadne/pharm.py` concentration->effect in one place (Hill n=1, reconstructs an
  implied IC50 for override-entered blocks, extrapolates ASSUMED effects with the same Hill form).
  Verified to reproduce drugs.yaml exactly at standard exposure (0 mismatches).
- `uncertainty.py` -> therapeutic-ceiling analysis, GABA-weight curve, levetiracetam inverse question,
  and TWO Monte Carlos: measurement-uncertainty-only vs full (also sampling unsourced parameters).
  The gap between them is the honest measure of how much rests on assumption.
- `make_report.py` -> `results/figures/*.png` + `results/uncertainty_report.md`. Drug palette is fixed and
  CVD-validated on all pairs (worst dE 9.9 deutan, 18.4 normal): CNB #6D3FD6, LTG #009E73, VPA #D55E00,
  CBZ #56B4E9, PGB #B5176B. Never recolour per chart.
- `PREREGISTRATION.md` - written BEFORE the surgery re-run. Hypothesis, analysis, decision rule, falsifier.
  Key finding: the n=6 pilot had power 0.00 (with 3 vs 3 the smallest attainable one-sided p IS 0.05).
- `virtual_surgery2.py` + `surgery.resect(graded=True)` - per-patient calibration to capacity ~0.50 (fixes the
  ceiling effect; uses pre-surgery connectome only, blind to outcome) and graded partial resections.
  Writes results/virtual_surgery2.csv; old pilot results kept for comparison.
- `PRESENTATION.md` - reframed claim, 60-second pitch, booth order, six Q&A drills.
- Drug sources gained sourced therapeutic ceilings: CBZ 12, LTG 15, VPA 100 mg/L (tops of reference ranges);
  PGB 600 mg/day and CNB 400 mg/day (FDA label maxima, both verified from the labels on 31 Oct).

### VERIFIED FINDING to state openly (a judge will otherwise find it)
For drugs with no coupling effect, threshold = boundary / (x0_effect x unit), so their RANKING is fixed by
pharmacology and is identical in every brain. Confirmed: CBZ/CNB threshold ratio is 3.36-3.40 observed across
Stanford, Sherbrooke, IDEAS sub-1 and Penn, against 3.350 predicted from x0_effects alone.
**Only the absolute dose is patient-specific** (CNB 0.47x Stanford vs 1.28x Sherbrooke). Order can only become
patient-specific via drugs mixing both mechanisms. This is a real limitation of the current formulation.

### BLOCKER for the surgery study (needs Anvi, cannot be done from the cloud session)
427 IDEAS II subjects have BOTH a resection map and a year-1 outcome (292 seizure-free, 135 not).
Only 6 are usable locally. 4 already-downloaded files (sub-13, 15, 108, 210) have NO outcome: do not re-download.
`results/ideas_download_list.md` lists 34 specific subject IDs to fetch from OpenNeuro ds007401
(file pattern `sub-<ID>_ses-1_Deterministic-Tractography_36-Lausanne_Count.csv`, drop into `data/ideas/`).
Network egress from this container is allowlisted and OpenNeuro is blocked, so she must download them.

### Overnight results (31 Oct, complete)
Response surface (1701 sims) + uncertainty + pre-registered surgery all finished. Headline numbers:
- Surface-derived thresholds match the independent bisection runs to within 0.07 (bisection resolution). VALIDATED.
- Boundary: needs x0 reduction 0.00338 at coupling 1.0, falling to 0 at coupling <= 0.92. 5/8 boundary cells
  agreed across extra noise seeds -> ~11% jitter on thresholds. Report it.
- **GABA weight flip is at 0.05, NOT 0.5.** The earlier robustness run only tested 0.0, which made the result look
  far more fragile than it is. Cenobamate wins for any GABA weight above ~0.05. This STRENGTHENS the headline.
- Monte Carlo (2000 draws x2): P(cenobamate first) = 98% measurement-only, 87% full uncertainty.
- Therapeutic ceiling: CNB and LTG reach containment within range; VPA, CBZ, PGB cannot. Dose doubling buys only
  1.19-1.43x effect (occupancy saturation).
- **Levetiracetam inverse: a pure release-reducer needs 175% (vs LTG) / 281% (vs CNB) release reduction = impossible.**
  Levetiracetam works clinically, so the coupling term is too weak or the mechanism is unmodelled. Best self-criticism
  in the project; it is on the poster and in the artifact.
- Virtual surgery, pre-registered: ceiling effect FIXED (pre-surgery capacity 0.54-0.69, was 0.93-1.00) but
  **AUC 0.44, p=0.65 -> H1 NOT supported.** Reported as a negative result. Binary-resection sensitivity AUC 0.11.
  Confound: patients who did worse had LARGER resections, so the measure partly tracks resection volume.
  Exploratory normalisation (reduction / resected fraction) gives AUC 0.78 but is driven by ONE patient (sub-3) = noise.
  Pre-registered as a co-primary for the next powered run (PREREGISTRATION.md Amendment 1, dated, written after seeing data).
- Artifact now has 8 steps (added "How sure?"), published: https://claude.ai/artifact/FWcKwmHHhNgUFEwgM3kYi6
- `ABSTRACT.md` and `PRESENTATION.md` are filled in with real numbers, no placeholders left.
- response_surface.py now checkpoints every 25 sims (a container restart killed a 97%-complete run once).

## Update 31 Oct midday: LEVETIRACETAM ADDED (derived, not assumed)
Anvi supplied the DrugBank entry; its target list includes CACNA1B (Cav2.2, N-type Ca channel) alongside SV2A.
That is quantifiable, so levetiracetam is now modelled - through calcium only, as an explicit LOWER BOUND.
- `drug_sources.yaml`: new `ca_current_block` (0.18 total HVA / 0.37 isolated N-type, Lukyanetz Epilepsia 2002;43:9-18,
  rat CA1, 1-200 uM, NO IC50 published) and cooperativity 2.6 (Reid/Bekkers/Clements J Neurosci 1998;18:2849-2855).
- `pharm.py` + `derive_parameters.py`: release_reduction = 1 - (1 - ca_block)^2.6 when ca_current_block is present,
  else the old assumed value. LEV coupling = 0.403 and is the ONLY non-assumed coupling value in the project.
- LEV PK: 12-46 mg/L (Patsalos & Gulcebi 2018; 3-34 and 20-40 also reported), not protein bound -> free 170.4 uM.
- `clinical_rules.yaml`: LEV is now a modelled drug, first_line for focal (NG217 5.2.1), GTC (5.1.1) and myoclonic
  (5.4.3, check:true), later_line absence (5.3.2), not_listed tonic/atonic. Removed from guideline_only lists.
- Stanford re-run with 6 drugs: CNB 0.47 > LTG 0.72 > VPA 1.53 > CBZ 1.59 > PGB 2.16 > **LEV 3.22**.
  LEV is last in ALL four sensitivity variants. Surface lookup gives 3.31 (diff +0.09, within resolution).
- Monte Carlo: LEV loses to LTG in 95% of draws, to CBZ in 74%, to PGB in 71%. Upper-bound (37% block) -> 1.91x, 5th.
- **THE HEADLINE SELF-CRITICISM**: ARIADNE ranks a NICE first-line focal drug last. Either the coupling term is far
  too weak or SV2A does nearly all the work. Stated on the poster, in the abstract and in Explorer step 6.
- Palette extended to 6: LEV = #00687A (light/report) and #a8d84a (artifact dark). 6-series CVD is WARN-band
  (worst pair dE 7.7 deutan), permitted only because every series is directly labelled.
- Known gap: source gives one block value, not a dose-response curve, so LEV's therapeutic-ceiling "effect gain" is 1.0x.
- Wright et al. Front Neurol 2013;4:192 (supplied 31 Oct midday): does NOT change any result, but adds
  (a) an independent real-patient cross-check of modelled concentrations via Otoul (CBZ 8.4 vs modelled 8.0;
  LTG 8.2 vs 8.75; VPA 83.8 vs 75 - ARIADNE ~11% low on VPA), (b) LEV protein binding <10% confirming fu,
  (c) no PK interaction LEV vs CBZ/VPA/TPM/LTG (matters: ARIADNE models monotherapy), (d) LEV half-life 6-8 h,
  (e) case reports of LEV AGGRAVATING myoclonus in JME -> added as a note on the myoclonic rule.

## PRIOR ART — must be cited (added 31 Oct, from PDFs Anvi supplied)
- **Sinha N. PhD thesis, Newcastle 2021** (supervisor Peter N Taylor): "Predicting seizure spread and
  neurosurgical outcomes in epilepsy by combining neuroimaging, machine learning, and computer modelling."
  Same question, same cohort lineage: Taylor's group released IDEAS; co-authors de Tisi/Vos/Duncan/Winston are
  the UCL Queen Square surgical cohort. Chapters published in Neurology, Epilepsia, NeuroImage:Clinical.
  KEY METHOD WE LACK: he standardises each patient's network against HEALTHY CONTROLS to get regional
  abnormality scores. We use raw spread capacity -> likely why our measure tracks resection volume.
  He FOUND an association between structural abnormality and post-surgical recurrence; our n=6 null is
  probably our implementation, not the premise. Recorded in PREREGISTRATION.md Amendment 2.
  POSITIONING: virtual surgery is a BORROWED validation test, not ARIADNE's contribution. The drug layer is.
- Gonen OM. PhD thesis, Melbourne 2020: "Multimodal Ultra-High-Field MRI Biomarkers of Epilepsy" (7T imaging
  methodology). Peripheral to ARIADNE; relevant only to lesion detection / better input imaging.
- A TDA dissertation abstract (persistent homology + zigzag persistence on iEEG) was also supplied: relevant as
  framing and as a possible source of a better onset estimate than our power-ratio heuristic. Do NOT bolt
  topology onto ARIADNE before the fair.
- **PHENYTOIN ADDED 31 Oct** (Anvi supplied Jun H et al. Drugs R D 2020;20:343-358): binding ~90% -> fu 0.10,
  range 10-20 mg/L (loss of effect <10, neurotoxicity >20) -> free 5.9 uM. Nav VDB IC50 47.78 uM from Kahlig 2022
  Table 1. **That fetch also re-verified CBZ 44.37 / LTG 39.09 / CNB 66.71 against the same table - all match.**
  x0_effect 0.369, threshold 2.29 (ties with pregabalin 2.225, inside the ~11% jitter).
  Caveat recorded: phenytoin has Michaelis-Menten elimination + concentration-dependent binding, so a single
  mid-range free concentration is cruder for phenytoin than for the others.
  Clinical gate: later_line focal (NG217 5.2.6 third-line add-on), avoid absence (5.3.3) and myoclonic (5.4.5).
- SEVEN drugs now: CNB 0.475 > LTG 0.76 > VPA 1.50 > CBZ 1.585 > PGB 2.225 ~ PHT 2.29 > LEV 3.31 (surface lookup).
- **CORRECTION 31 Oct: IDEAS is hosted on FIGSHARE, not OpenNeuro.** Index: https://sites.google.com/view/cnnp-lab/ideas-data
  Connectomes are in the "fully processed diffusion MRI ... and connectomes" release: https://figshare.com/s/2dbae1fbfe72f7e66e1a
  Metadata: .../bab70268afeb1071202b ; resection percentages: .../097ba0e254e36f0eee52 ;
  **Z-scored FreeSurfer vs normative controls: https://figshare.com/s/8c086fc295a75f85e628** (relevant to the
  Sinha control-normalisation fix). Earlier notes citing "OpenNeuro ds007401" sent Anvi to the wrong place.
  Session numbers vary: 6 of her 10 files are ses-1, 4 are ses-2. Do not filter on ses-1.

## Update 31 Oct ~13:05: PRE-REGISTERED NEGATIVE RESULT at n=405 (abnormality_study.py)
Anvi downloaded the IDEAS z-scored FreeSurfer release (figshare .../8c086fc295a75f85e628): thickness/area/volume
z-scores vs 100 normative controls for 428 patients -> `data/ideas/zscores/`. With the resection table and
outcomes, **415 patients have all three; 405 analysable**. FIRST ADEQUATELY POWERED SAMPLE IN THE PROJECT.
Pre-registered as Amendment 3 BEFORE running; result recorded immediately after.
- **H1 (removing more abnormal tissue -> seizure freedom): NOT SUPPORTED.** AUC 0.476 / 0.499 / 0.484
  (thickness / area / volume), all at or below chance. Adjusted logistic beta -0.277, p=0.15. Temporal-only
  AUC 0.512. **Amendment 1's proposed co-primary is formally WITHDRAWN.**
- H2 (whole-brain abnormality -> worse outcome, Sinha's direction): AUC 0.559, p=0.027 thickness only -
  does NOT clear the pre-specified Bonferroni 0.0167. Suggestive only.
- H3: **resected volume alone does NOT separate groups at n=405 (AUC 0.465)** - so the volume confound that
  dominated the n=6 pilot was itself a small-sample artefact. The n=6 exploratory AUC 0.78 was noise, confirmed.
- SCOPE: tests MORPHOMETRIC abnormality only. Does NOT test NETWORK abnormality (what Sinha used, what ARIADNE
  computes). The network version remains untested and still needs connectome files.
- `abnormality_study.py`, `results/abnormality_study.csv`, `results/abnormality_report.md`. Needs statsmodels.

## Update 31 Oct ~13:20: FULL COHORT ARRIVED (288 connectomes) + second pre-registered negative
Anvi uploaded the whole Lausanne-36 deterministic release as one zip: **288 Count matrices, 191 patients +
97 controls** -> all copied flat into `data/ideas/`. **183 patients analysable** (connectome + resection +
outcome; 128 SF / 55 NSF). Controls 4001-4099.
- **Analysis A (`network_abnormality.py`), pre-registered Amendment 4, NEGATIVE.** Nodal strength z-scored vs
  97 controls. H4 (network abnormality of resected tissue -> seizure freedom): n=182, AUC **0.535**, p=0.23;
  temporal-only AUC 0.499 (chance); adjusted logistic beta +0.512 p=0.30. H5 (global) AUC 0.523, p=0.31.
  Threshold was p<0.025. **Both morphometric (n=405) and network (n=182) abnormality hypotheses now fail.**
  Does NOT refute Sinha - our measure is a deliberately simple nodal-strength z-score, not his metric.
  Confound note: resected volume alone AUC 0.437 (worse outcomes had LARGER resections).
  -> `results/network_abnormality{,_tests}.csv`, `results/network_abnormality_report.md`
- **Analysis B (`virtual_surgery3.py`) RUNNING**: virtual surgery at n=183, ~298 s/patient (~15 h total).
  Patients processed INTERLEAVED (NSF, SF, NSF, SF...) so stopping early still leaves a balanced sample;
  resumable via results/virtual_surgery3.csv. Compute trade fixed in Amendment 4: 15-onset paired subsample,
  4x10 calibration, binary sensitivity dropped.
- Sherbrooke/IDEAS response-surface chain was STOPPED to free CPU (checkpointed, resumable later).
- **GABA WEIGHT NOW DERIVED (31 Oct ~13:45):** w = rho = g_I/g_E, the cortical inhibitory/excitatory CONDUCTANCE
  ratio, from net-drive algebra: reduction in net drive = f_Na + rho*g_GABA. Sourced: Haider J Neurosci 2006
  (in vivo up states, gI/gE ~1), Rudolph 2007 (awake cat, inhibition several-fold higher), Haider 2013.
  **rho = 1.0 adopted** (most conservative sourced value); old 0.5 was BELOW the lowest measurement.
  CNB x0_effect 1.789 -> 2.54, threshold 0.475 -> 0.335. ORDER UNCHANGED. Flip point 0.05 is 20x below the
  lowest plausible rho, so the result is far more robust than the old framing implied.
- **ONLY ONE ARBITRARY NUMBER LEFT that can reorder drugs**: COUPLING_UNIT/X0_UNIT (= 15). Swept 0.1x-25x in
  `results/unit_ratio_sensitivity.json`: **cenobamate first across 0.1x-4x (40-fold range)**; lamotrigine second
  across the whole 250-fold range. Pregabalin only takes first at >=10x.

## Update 31 Oct ~14:30: UI upgrade, 8th drug, and the clinical rules corrected against verbatim NG217

### Explorer (artifact) rewritten around mechanism class
- `scratchpad/story/build.py` is now the ONLY way to build `ariadne_story.html`: it injects `data.json` AND
  regenerates the embedded RULES payload from `clinical_rules.yaml`. A hand-pasted RULES copy went stale once
  (levetiracetam showed as "not a guideline option" after it had been made first-line). Never hand-edit the output.
- Drugs are grouped by **mechanism class** in a `<select>` with `<optgroup>` headers, replacing the flat chip row:
  Sodium-channel blockers (CBZ, LCM, LTG, PHT) / Sodium + GABA-A (CNB) / Calcium-transmitter release (PGB, LEV) /
  Broad spectrum (VPA). **Colour now encodes CLASS, with lightness steps inside a class.** This is a deliberate
  break from the report palette (which encodes drug): 8 categorical colours cannot pass CVD separation, and every
  element in the Explorer is directly labelled. Report/poster figures keep the per-drug palette.
- `drugsFor()` filters the selector by which drugs that brain actually has, so brains with different drug sets
  cannot throw.
- **Named-anatomy view** for the IDEAS patient: `data.json` case 3 has `kind:"named"`, 82 Lausanne-36 regions with
  name/lobe/hemisphere and recruitment `masks` per drug, instead of 3D coordinates. `drawNamed()` renders them as
  a lobe-grouped grid that lights up; `prep/draw/seizingCount/hud/strip/startPlay` are all guarded by `isNamed()`.
  There is no time course for this case (masks only), so it never animates — `startPlay` jumps straight to the
  end state rather than faking a spread.
- Onset chip is labelled "ONSET" in text because valproate's orange sits near the amber onset marker.
- Verify with `node /home/claude/check_story.js` (walks every step on every case and every drug option, reports
  page/console errors). The only expected console error is a blocked Google Fonts fetch offline.

### LACOSAMIDE ADDED (8th drug, fully sourced)
- `nav_vdb_ic50_uM: 269.3` — Kahlig 2022 Table 1, hNaV1.6, 269,300 nM slope 1.2. The SAME re-fetch of that table
  re-confirmed CNB 66,710 / PHT 47,780 / CBZ 44,370 / LTG 39,090 / OXC 42,000 nM, and revealed that
  **valproate appears only as "18% inhibition @ 1 mmol/L"** — which is where ARIADNE's 1.5% VPA Na block comes from.
- MW 250.3, protein binding <15% -> fu 0.85 (conservative end), range 10-20 mg/L -> free 50.9 uM.
- **x0_effect 0.53, identical to carbamazepine**, so the two tie (Stanford: both 1.594). Real and explainable.
- Caveats recorded: lacosamide acts on SLOW inactivation, which the Kahlig assay only partly captures; and the
  10-20 mg/L range is contested (Clin Neuropharmacol Dec 2022 proposes 10-12 mg/L, ~25% lower free concentration).
- **OXCARBAZEPINE DELIBERATELY NOT ADDED** despite having its IC50: Kahlig measured the parent, but the active
  circulating species is licarbazepine (MHD), which is what the reference range refers to. Pairing the parent's
  IC50 with MHD's concentration would be wrong. See DRUG_PANEL.md.

### CLINICAL RULES CORRECTED — three were wrong, both `check: true` items resolved
NG217 sections 5.1-5.5 read verbatim from the NCBI Bookshelf full text (NBK581165). `grep -c 'check: true'` is now 0.
- **valproate / focal: was `not_listed`. WRONG** — it is NG217 5.2.6, third-line add-on (Box 2 applies).
  Now `later_line`, so the childbearing rule correctly escalates it to `restricted`.
- **pregabalin / focal: was a guess from the Lyrica label.** It is NG217 5.2.5, second-line add-on.
- **myoclonic: both the valproate and levetiracetam rules invented sex-specific wording.** 5.4.2 simply offers
  "a choice of levetiracetam or sodium valproate" as first-line, with the Box 2 reference and no sex split.
  Levetiracetam was also mis-cited to 5.4.3; it is 5.4.2.
- Full focal ladder now recorded: 5.2.1 LTG/LEV -> 5.2.2 CBZ/OXC/ZNS -> 5.2.3 LCM mono -> 5.2.4 add-on
  CBZ/LCM/LTG/LEV/OXC/ZNS/TPM -> 5.2.5 BRV/CNB/ESL/perampanel/PGB -> 5.2.6 PB/PHT/VPA/tiagabine/vigabatrin.
  Cenobamate's 5.2.5 placement independently confirms TA753.
- Lacosamide placed: focal 5.2.3+5.2.4 later_line; GTC 5.1.3 later_line; not listed for absence/myoclonic/tonic/atonic.
- `python3 tests/test_clinical.py` passes.

### All four brains brought onto ONE pharmacology
`recompute_brains.py` + `refresh_sensitivity.py` + `add_lacosamide.py`. The stored results predated the derived
GABA weight and the last three drugs, so cenobamate's threshold disagreed between files. Cenobamate now:
Stanford 0.344, IDEAS sub-1 0.469, Penn 0.656, Sherbrooke 0.906 — all consistent with the exact prediction that
threshold scales as 1/x0_effect (ratio 0.704 from the weight change; observed 0.680-0.724, within bisection
resolution). The cross-brain table in the Explorer is rebuilt from these files, 7-8 columns.
NOTE: the CBZ/CNB ratio quoted earlier in this file (3.36-3.40) was computed at GABA weight 0.5; at the derived
weight 1.0 it is ~4.6-4.9. The CLAIM is unchanged — the ratio is constant across brains — only the number moved.

### DRUG_PANEL.md — the answer to "what new data do you need" for a 20-drug panel
Honest ceiling with today's model is 12-13 fully sourced drugs, not 20. The file lists, per drug, exactly which of
MW / range / binding / mechanism-number / ceiling is missing, which drugs need a NEW MODEL TERM (perampanel=AMPA,
vigabatrin=GABA-T, tiagabine=GAT-1, felbamate=NMDA), and two mapping decisions that must be made first
(T-type calcium is somatodendritic, NOT presynaptic release, so it must not go through the coupling term;
and prodrugs/active metabolites need the metabolite's own IC50 at the metabolite's own concentration).
Four papers would unlock nearly all of it; Patsalos 2018 TDM guidelines alone covers range+binding for every drug.
Highest-value additions: BRIVARACETAM (turns the levetiracetam-ranked-last embarrassment into a designed
experiment) and PHENOBARBITAL (a pure GABA-A drug, so it tests the derived rho independently of cenobamate).

### Compute notes
- virtual_surgery3.py STOPPED at 10/183 rows to free CPU (resumable, interleaved so the partial sample stays balanced).
- The Sherbrooke response surface reached 1700/1701 grid cells and was killed during boundary confirmation by a
  pkill whose pattern matched more than intended. Checkpoint intact; restarting costs ~1 cell plus confirmation.
- **pgrep/pkill WARNING, hit twice now**: a pattern matching the script's own command line kills the shell or waits
  forever. Use a lock file or a distinctive literal, never a substring of the running script.
