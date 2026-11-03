# IDEAS II: which connectome files to download next

The pre-registered virtual-surgery test is written and runs. It is limited by one thing only: missing
connectome files. The outcome data is already here.

- Subjects in IDEAS II with BOTH a resection map and a year-1 outcome: **427** (292 seizure-free, 135 not seizure-free)
- Connectome files present locally: **10**, of which **6** are usable
- Already downloaded but UNUSABLE (no year-1 outcome recorded, do not re-download): [13, 15, 108, 210]

## Power

80% power needs about **15 per group at AUC 0.80**, or **25 per group at AUC 0.70** (simulated, see
PREREGISTRATION.md). The current 3 vs 3 has power 0.00: with that sample the smallest attainable
one-sided p is exactly 0.05, so significance is arithmetically impossible.

## Download these 34 files

### WHERE THE FILES ACTUALLY ARE (corrected 31 Oct)

**Not OpenNeuro — figshare.** The index is the CNNP Lab IDEAS page:
<https://sites.google.com/view/cnnp-lab/ideas-data>

The connectomes live in the release described as *"Raw T1w and FLAIR, plus fully processed diffusion MRI
(tensor maps in native and MNI-152 space, and connectomes)"*:

**<https://figshare.com/s/2dbae1fbfe72f7e66e1a>**

Inside it, look for the deterministic-tractography connectomes. You already have the other two pieces:
`Metadata_Release_Anon.csv` (clinical metadata release) and `table_resected.csv` (percentage of each
Desikan-Killiany region resected).

Also on that index page and worth grabbing for the control-normalisation work:
**Z-scored FreeSurfer measures against normative controls** — <https://figshare.com/s/8c086fc295a75f85e628>

Each subject needs ONE file, named:

```
sub-<ID>_ses-<N>_Deterministic-Tractography_36-Lausanne_Count.csv
```

**Important: the session number is not always 1.** Of the ten files already downloaded, six are `ses-1`
and four are `ses-2` (sub-13, sub-15, sub-143 and sub-3 are all `ses-2`). Take whichever session the
subject has; the loader prefers `ses-1` when both exist. Do not skip a subject because there is no `ses-1`.

Put them in `data/ideas/`. Any subset helps; a balanced subset helps most.

**Seizure-free (ILAE 1-2), 17 subjects**

```
5, 8, 9, 11, 12, 16, 18, 19, 21, 22, 25, 26, 27, 30, 33, 34, 35
```

**Not seizure-free (ILAE 3-5), 17 subjects**

```
4, 6, 10, 14, 17, 23, 24, 28, 29, 31, 32, 38, 40, 42, 49, 50, 53
```

## Then run

```
python virtual_surgery2.py --workers 4
```

It resumes automatically and only processes subjects it has not already done.

## Why these were chosen

Lowest subject IDs with a usable outcome, split to give a balanced 20 vs 20. Balance matters because
the cohort is 292 seizure-free vs 135 not; an unbalanced sample wastes power.
