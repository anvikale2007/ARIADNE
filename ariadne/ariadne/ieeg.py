"""Intracranial EEG (iEEG) input for ARIADNE: find where a recorded seizure is strongest and name the anatomical region.

Built for BIDS iEEG such as the Penn 3T-7T epilepsy dataset (CC-BY-4.0): EDF + channels.tsv + events.tsv + electrodes.tsv
(electrodes carry MNI coordinates and an anatomical label per contact).

Method (simple on purpose): 60 Hz notch + 1-120 Hz band-pass; per channel, power in an early-ictal window divided by power
in a pre-ictal baseline; the strongest channels are mapped to their anatomical labels and then to connectome regions
through a hand-made label crosswalk (CROSSWALK below).

!! LIMITATIONS !!
 - Signals are only as good as the baseline: short baselines (30 s) make the ratios noisy.
 - The crosswalk from Desikan/AAL-style labels to the TVB template regions is approximate.
 - The connectome (template) is NOT the patient's own; there is no diffusion MRI in this dataset sample.
"""
import re
import numpy as np
import pandas as pd

ICTAL_WINDOW_S = 10.0
BASELINE_END_MARGIN_S = 5.0     # baseline ends this long before onset

# anatomical label (lower case, hemisphere prefix removed) -> TVB template region suffix.  Approximate.
CROSSWALK = {
    "hippocampus": "HC", "parahippocampal": "PHC", "amygdala": "AMYG", "entorhinal": "PHC", "fusiform": "TCV",
    "superior temporal": "TCS", "middle temporal": "TCC", "inferior temporal": "TCI", "temporal pole": "TCPOL",
    "transverse temporal": "A1", "insula": "IA", "precentral": "M1", "postcentral": "S1", "supramarginal": "S2",
    "superior frontal": "PFCDM", "rostral middle frontal": "PFCDL", "caudal middle frontal": "PFCCL",
    "pars opercularis": "PFCVL", "pars triangularis": "PFCVL", "pars orbitalis": "PFCORB",
    "lateral orbitofrontal": "PFCORB", "medial orbitofrontal": "PFCM", "frontal pole": "PFCPOL",
    "inferior parietal": "IP", "superior parietal": "PCS", "precuneus": "PCM", "cuneus": "V2", "lateral occipital": "V2",
    "lingual": "V1", "pericalcarine": "V1", "posterior cingulate": "CCP", "isthmus cingulate": "CCP",
    "caudal anterior cingulate": "CCA", "rostral anterior cingulate": "CCR",
}


def electrode_to_edf_name(name, edf_names):
    if name in edf_names:
        return name
    m = re.fullmatch(r"G(\d+)", name)                 # grid contacts: G12 in electrodes.tsv = RG12 in the EDF
    if m and f"RG{int(m.group(1)):02d}" in edf_names:
        return f"RG{int(m.group(1)):02d}"
    return None


def ictal_channel_scores(raw, onset_s, picks):
    """Early-ictal / baseline broadband power ratio per channel (1-120 Hz, 60 Hz notched).
    Line length was tried first and showed no ictal change on the sample epoch, while power rose ~3.6x on most channels."""
    r = raw.copy().pick(picks).load_data(verbose="ERROR")
    r.notch_filter(60.0, verbose="ERROR")
    r.filter(1.0, 120.0, verbose="ERROR")
    fs = r.info["sfreq"]; X = r.get_data()
    a0, b1 = int(2.0 * fs), int(max(onset_s - BASELINE_END_MARGIN_S, 3) * fs)   # skip 2 s of filter edge artefact
    i0, i1 = int(onset_s * fs), int((onset_s + ICTAL_WINDOW_S) * fs)
    pw = lambda a, b: (X[:, a:b] ** 2).mean(1)
    return dict(zip(r.ch_names, pw(i0, i1) / np.maximum(pw(a0, b1), 1e-30)))


def analyse(edf, channels_tsv, events_tsv, electrodes_tsv, top_n=8):
    import mne
    raw = mne.io.read_raw_edf(edf, preload=False, verbose="ERROR")
    ch = pd.read_csv(channels_tsv, sep="\t")
    keep = ch[ch["type"].isin(["ECOG", "SEEG"]) & (ch["status"] == "good")]["name"].tolist()
    ev = pd.read_csv(events_tsv, sep="\t"); onset = float(ev[ev["trial_type"] == "ictal"]["onset"].iloc[0])
    el = pd.read_csv(electrodes_tsv, sep="\t")
    edf_names = set(raw.ch_names)
    el["edf"] = [electrode_to_edf_name(n, edf_names) for n in el["name"]]
    el = el[el["edf"].notna() & el["edf"].isin(keep)]
    scores = ictal_channel_scores(raw, onset, el["edf"].tolist())
    el["score"] = el["edf"].map(scores)
    el = el.sort_values("score", ascending=False)
    top = el.head(top_n)[["name", "anat_label", "region_coarse", "tissue_atropos", "x", "y", "z", "score"]]
    # region score = sum of the 3 strongest channels carrying that anatomical label (named labels only)
    lab = el[el["anat_label"].str.lower() != "emptylabel"].copy()
    agg = lab.groupby("anat_label")["score"].apply(lambda s: s.nlargest(3).sum()).sort_values(ascending=False)
    return {"onset_s": onset, "n_channels": len(el), "top_channels": top.round(3).to_dict("records"),
            "label_scores": agg.round(3).to_dict()}


def to_template_region(anat_label, labels):
    """'right middle temporal' -> 'rTCC' if that region exists in the connectome labels."""
    s = anat_label.lower().strip()
    hemi = "r" if s.startswith("right") else ("l" if s.startswith("left") else None)
    core = re.sub(r"^(right|left)\s+", "", s)
    suf = CROSSWALK.get(core)
    if hemi and suf and hemi + suf in labels:
        return hemi + suf
    return None
