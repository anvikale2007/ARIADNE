"""IDEAS II (OpenNeuro ds007401; Taylor et al., Epilepsia 2026) loader: real epilepsy-patient connectomes + surgery data.

Connectomes: <sub>_ses-1_Deterministic-Tractography_36-Lausanne_Count.csv (82 x 82, no labels).
Region order is taken from table_resected.csv (same 82-region Lausanne-36 / Desikan-Killiany atlas). This order is
INFERRED (not stated in the CSV); it passes anatomical sanity checks (hemispheric block structure; hippocampus ->
amygdala / parahippocampal / thalamus) but confirm it against the dataset documentation.
"""
import os, re, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data", "ideas")


def labels():
    r = pd.read_csv(os.path.join(DATA, "table_resected.csv"))
    return [s.strip("'") for s in r.iloc[:, 0]]


def resected_fraction(sub_id):
    r = pd.read_csv(os.path.join(DATA, "table_resected.csv"))
    return r[str(sub_id)].to_numpy(dtype=float)


def metadata():
    return pd.read_csv(os.path.join(DATA, "Metadata_Release_Anon.csv"))


def count_file(sub_id, folder=None):
    pat = os.path.join(folder or DATA, f"*sub-{sub_id}_ses-*_Deterministic-Tractography_36-Lausanne_Count.csv")
    f = sorted(glob.glob(pat))      # prefer ses-1 if both sessions exist
    return f[0] if f else None


def to_connectome_folder(sub_id, out_root=None, folder=None):
    """Write connectomes/ideas_sub-<id>/weights.csv + labels.csv so run_pipeline.py can use it."""
    f = count_file(sub_id, folder)
    if f is None:
        raise FileNotFoundError(f"No Count CSV for sub-{sub_id}")
    C = np.loadtxt(f, delimiter=",")
    out = os.path.join(out_root or os.path.join(HERE, "connectomes"), f"ideas_sub-{sub_id}")
    os.makedirs(out, exist_ok=True)
    np.savetxt(os.path.join(out, "weights.csv"), C, delimiter=",", fmt="%.0f")
    np.savetxt(os.path.join(out, "labels.csv"), labels(), fmt="%s")
    return out


def available_subjects(folder=None):
    fs = glob.glob(os.path.join(folder or DATA, "*_36-Lausanne_Count.csv"))
    return sorted({int(re.search(r"sub-(\d+)_", os.path.basename(f)).group(1)) for f in fs})
