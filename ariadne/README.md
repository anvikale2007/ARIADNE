# ARIADNE: connectome-based antiseizure drug ranking (proof of concept)

Research prototype. **Not a diagnostic or clinical tool.** It simulates how antiseizure drugs, represented by
their mechanisms, change seizure spread across a brain's wiring, and ranks which drug stops the spread most easily
in *that* connectome and *that* onset site.

## What is in here
```
dwi_to_connectome.py   diffusion MRI (.nii.gz + .bval + .bvec) -> structural connectome (DIPY tractography)
ariadne/model.py       network Epileptor-2D seizure model (numba), verified against The Virtual Brain
ariadne/engine.py      calibration + "drug strength needed to stop spread" search + ranking
drugs.yaml             drug panel: mechanism -> model effect (ASSUMPTIONS, add citations!)
run_pipeline.py        runs everything, writes results/*.json + a ranking chart
app.py                 Streamlit app (ranking, brain view, sensitivity table)
led.py + arduino/      drives the 12-LED mini-brain (green / amber / red)
connectomes/           connectomes already built from the two real scans you supplied (+ QC images)
results/               precomputed results, so the app opens instantly
tests/                 equation check against TVB
```

## Setup
```
pip install -r requirements.txt
```

## Run
```
streamlit run app.py                         # demo UI on the precomputed results
python run_pipeline.py                       # recompute (about 5-6 min on 2 CPU cores)
python run_pipeline.py --onset 12            # different seizure onset region (index or label like P12)
python run_pipeline.py --sensitivity         # also test ranking robustness (adds ~20 min)
python run_pipeline.py --connectome connectomes/sherbrooke_3shell
```

### Use a new scan
```
python dwi_to_connectome.py --dwi scan.nii.gz --bval scan.bval --bvec scan.bvec --out connectomes/my_scan
python run_pipeline.py --connectome connectomes/my_scan
```
Any tractography pipeline works as long as you provide `weights.csv` (n x n streamline counts) and optionally
`labels.csv` in a folder. Passing `--atlas labels.nii.gz` (integer atlas already registered to the DWI) gives a
connectome with named anatomical regions.

### LED demo
```
python led.py --result results/<file>.json --state baseline --dry-run     # prints the frame
python led.py --result results/<file>.json --state cenobamate --port COM3 # real Arduino
```
Flash `arduino/ariadne_leds/ariadne_leds.ino` (Adafruit NeoPixel library, strip on pin 6).
Each LED = a group of nearby brain regions; red = seizing, amber = partly, green = quiet.

## How the model works
1. DWI -> tensor fit -> deterministic tractography -> regions -> streamline counts between regions = connectome.
2. Each region is an Epileptor-2D neural mass (fast activity x1, slow permittivity z). Regions are coupled through the
   connectome, so a seizing region pushes its neighbours towards seizure.
3. **Calibration:** background excitability is set to the lowest level where a seizure from the chosen onset recruits most
   of the network while a healthy network stays quiet, i.e. the brain's own seizure-spread threshold.
4. **Drugs** act through two knobs: lower excitability (Nav block, GABA-A) and weaker coupling (less transmitter release).
   For each drug the pipeline dials up its strength (bisection) until the spread stops. Lower strength needed = better
   for this brain. Differences between drugs therefore come from *how the brain's threshold responds to each mechanism*.
5. **Sensitivity analysis** re-weights the two mechanisms (x0.5 / x2) and recomputes the ranking, because the
   drug effect sizes are assumptions.

## Honest limitations (put these on your poster)
- Scans are healthy volunteers (Stanford HARDI, Sherbrooke 3-shell), downsampled to 4 mm: no patient, no epilepsy, no outcomes.
- Regions are geometric parcels (k-means), not named anatomy, because no T1/atlas was available. "Seizure onset P39" is a
  placeholder location; with an atlas you could place the onset in e.g. hippocampus.
- Drug potencies in `drugs.yaml` are mechanism-based assumptions, not measured. Add a citation for each.
- Spread is all-or-none near the threshold in this reduced model; the ranking reflects relative thresholds, not seizure frequency.
- Nothing is validated against clinical response. Next steps: (1) reproduce published drug-induced connectivity changes,
  (2) retrospective test on a hospital cohort with ethics approval, (3) prospective study. Not for clinical use.

## Data credits
- Stanford HARDI: Rokem et al. 2015, PLOS ONE 10(4):e0123272; https://purl.stanford.edu/yx282xq2090 (ODC-PDDL).
- Sherbrooke 3-shell: Rokem 2017, Dipy Data, https://hdl.handle.net/1773/38475 (CC0).
- Methods: DIPY (Garyfallidis et al. 2014); Epileptor (Jirsa et al. 2014, Brain); reduced Epileptor / VEP (Proix et al. 2014; Jirsa et al. 2017, NeuroImage);
  The Virtual Brain (Sanz-Leon et al. 2013).

## Update: sourced drug parameters + verification
- `drug_sources.yaml` holds the inputs (IC50s, therapeutic concentration ranges, free fractions) with a status on each number
  (`sourced` / `reference` / `assumed`). `python derive_parameters.py` turns them into `drugs.yaml`. Edit the sources file, not `drugs.yaml`.
- Anchored to a paper: persistent-Na IC50 for carbamazepine (77.5 uM) and lamotrigine (78.5 uM) and the negligible Nav effect of valproate
  (Kahlig et al. 2022, Epilepsia 63:697-708); cenobamate persistent-Na IC50 53.1 uM (Nakamura et al. 2019).
- Still assumptions: all concentrations (standard TDM ranges, confirm in a textbook), cenobamate's exposure, every GABA / T-type / alpha2-delta / release effect.
  Pregabalin and cenobamate rankings in particular are driven by assumed numbers.
- `python verify_model.py` runs 11 checks (healthy quiet, untreated spreads, dose-response monotone, noise and time-step robustness, onset dependence).
  Report: `results/verification_report.md`. This is verification of the simulator, not clinical validation.

## Multimodal inputs (diffusion MRI + EEG + genetics)
```
python run_pipeline.py --genetics example_patient.yaml                       # gene -> drug safety rules (edit genetics_rules.yaml)
python run_pipeline.py --eeg-edf chb01_03.edf --seizures chb01_03.edf.seizures   # EEG suggests the onset region
python run_pipeline.py --eeg-edf chb01_03.edf --seizures chb01_03.edf.seizures --genetics example_patient.yaml
```
- **EEG** (`ariadne/eeg.py`): finds the scalp channels whose line length rises most at seizure onset and maps their electrode positions to the nearest connectome region. This is only a *prior* for the onset zone: no head model and no registration between the EEG subject and the MRI subject.
- **Genetics** (`ariadne/genetics.py`, `genetics_rules.yaml`): published gene->drug rules. HLA-B*15:02 -> exclude carbamazepine (CPIC 2018) is sourced; the SCN1A and POLG rules still need their citations confirmed.
- The EEG, MRI and genetics in any demo come from different people unless you have a matched patient. Say so.
- `python tests/test_eeg.py` checks the annotation parser (on the CHB-MIT files) and the EEG onset finder (on synthetic data).

## Real epilepsy case: Penn 3T-7T dataset (sub-PID02, CC-BY-4.0)
```
python run_pipeline.py --template --ieeg-edf <ictal .edf> --ieeg-channels <channels.tsv> --ieeg-events <events.tsv> --ieeg-electrodes <electrodes.tsv>
```
(Paths are inside the Penn starter pack: `sample_sub-PID02/sub-PID02/ses-clinical01/ieeg/`.) `ariadne/ieeg.py` ranks intracranial contacts by early-ictal / baseline power,
reads each contact's anatomical label from the electrodes table, and maps labels to the TVB template's named regions with an approximate crosswalk.
Findings on the sample epoch: the power rise is widespread across the right hemisphere (grid contacts postcentral / precentral / frontal / temporal / hippocampal all high, ~3.6x median),
consistent with the recorded lateralisation (right only) but NOT a sharp onset zone. The 15 s epoch is short and probably already spread. The suggested onset (rS1) is therefore a weak prior.
Limitations: the connectome is the TVB template (the dataset sample has no diffusion MRI); the label crosswalk is approximate; one patient, one epoch; no genetics in this dataset.
Cite: Lucas et al. 2026 (Penn 3T-7T Paired Epilepsy Imaging, doi:10.26275/AOSZ-HJ9U); Imtiaz et al. 2026, Sci Data 13:1320.

## Update: user-supplied cited drug numbers (29 Oct)
`drug_sources.yaml` now uses:
- **Cenobamate** free Cave,ss 30.1 uM at 200 mg/day (54.4 uM at 400 mg/day, sensitivity) from Vernillet et al. 2020 / FDA review, fu 0.40 (FDA label);
  persistent-Na IC50 71.7 uM (Kahlig 2022, same assay as CBZ/LTG); tonic GABA-A EC50 36.6 uM (Sharma et al. 2020).
- **Pregabalin** Cave,ss 19.8 uM at 300 mg/day (32.7 uM at 600 mg/day), unbound; alpha2-delta-1 KD 6.0 nM (Li et al. 2011): binding is saturating at clinical levels.
- **Carbamazepine / lamotrigine / valproate** TDM ranges and protein binding from Patsalos et al. 2008/2018 and Reimers et al. 2018; valproate uses the concentration-dependent free fraction from the FDA label anchors.
- **Gene rules**: Dravet / SCN1A loss-of-function -> avoid sodium-channel blockers (Wirrell et al. 2022 international consensus); POLG -> valproate contraindicated (FDA label; Stewart et al. 2010). The SCN1A rule no longer fires for gain-of-function variants.

Still ASSUMED (the next numbers to find): the size of pregabalin's and lamotrigine's reduction in transmitter release, valproate's GABA / T-type effect,
the GABA weight (how GABA-A potentiation translates to excitability), and NAV_REF. `robustness.py` tests the ranking against these choices.

## Update: use-dependent (inactivated-state) sodium block + robustness (29-30 Oct)
- Sodium block now uses the inactivated-state (voltage-dependent) IC50s from Kahlig 2022 (CBZ 44.4, LTG 39.1, cenobamate 66.7 uM): the state
  channels are in during seizure depolarisation. 10 Hz use-dependent IC50s (CBZ 1418, LTG 516, cenobamate 749 uM) give <4% block at
  therapeutic free levels and are recorded but have little effect. `drug_sources.yaml: nav_mode` switches back to persistent current.
- Main ranking (Stanford scan): cenobamate > lamotrigine > valproate > carbamazepine > pregabalin.
- `results/robustness_report.md`: cenobamate is top in 5/6 variants (coupling K -4/-6, GABA weight 0/1, 400 mg/day exposure, old sodium model).
  With GABA weight = 0, lamotrigine is top: cenobamate's lead depends on how its GABA-A effect is weighted, which is still an assumption.
- Thresholds are resolved to ~0.06-0.13 (bisection), so near-equal values (e.g. carbamazepine vs valproate 1.53 vs 1.59) are ties.
- Pregabalin's effect does not change with exposure because its release reduction is a fixed assumption (binding is saturating, KD 6 nM).

## Clinical gate: seizure type first (ILAE 2025)
Doctors classify the seizure before choosing a drug, because some drugs worsen certain seizure types.
```
python run_pipeline.py --seizure-type focal_impaired_consciousness --age 24 --childbearing   # focal: simulated, status shown per drug
python run_pipeline.py --seizure-type absence                                              # generalized: guideline status only, no simulation
```
Rules and sources: `clinical_rules.yaml` (NICE NG217 section 5, NICE TA753, FDA labels; classification Beniczky et al., Epilepsia 2025).
Example: for absence, myoclonic, tonic and atonic seizures carbamazepine and pregabalin are flagged "avoid" (NG217 5.3.3 / 5.4.5 / 5.5.6);
cenobamate is "not licensed" outside focal seizures; valproate is "restricted" under 55 or if able to have children (NG217 Box 2).

## Rigour pass: how much of the answer is measurement, and how much is assumption?

Every drug enters the model through just two numbers (an excitability reduction and a coupling
reduction), so the boundary between "seizure spreads" and "seizure contained" can be mapped once and
reused. That turns each drug threshold into a lookup and makes uncertainty propagation affordable.

```
python response_surface.py --result results/stanford_hardi_onsetP39.json   # ~71 min on 2 cores
python uncertainty.py --surface results/surface_stanford_hardi_onsetP39.json
python make_report.py          # figures + results/uncertainty_report.md
```

This produces four things:

1. **Therapeutic ceiling.** Not "needs 1.59x strength" but "at the top of its therapeutic range, does
   this drug contain the seizure?" Occupancy saturates, so a doubled dose delivers less than a doubled
   effect; the analysis uses the actual concentration, not a linear multiplier.
2. **GABA-weight curve.** One unsourced constant decides which drug ranks first. The curve shows
   exactly where the ranking flips instead of hiding the dependence.
3. **Uncertainty propagation.** Two Monte Carlo runs - one sampling only ranges that have published
   sources, one also sampling the parameters with no source. The ranking becomes a probability, and
   the gap between the two runs measures how much rests on assumption.
4. **The levetiracetam question, inverted.** Levetiracetam is not modelled because no published number
   converts SV2A binding into an effect size. Instead the analysis reports how large a
   release-reducing effect would have to be to match each modelled drug - a falsifiable number.

### A limitation this analysis exposed
For drugs acting only on excitability, threshold = boundary / (x0_effect x unit). Their **ranking is
therefore fixed by pharmacology and identical in every brain**; only the absolute dose is
patient-specific. Verified: the carbamazepine/cenobamate threshold ratio is 3.36-3.40 across four
connectomes, against 3.35 predicted from the drug parameters alone.

## Pre-registered virtual surgery

`PREREGISTRATION.md` fixes the hypothesis, analysis and decision rule before the analysis is run.

```
python virtual_surgery2.py --workers 4
```

Per-patient calibration (blind to outcome) fixes the pilot's ceiling effect, and partial resections are
graded rather than counted all-or-nothing. See `results/ideas_download_list.md` for the specific
IDEAS II subjects to download: the test is written and runs, and is limited only by missing connectome
files.
