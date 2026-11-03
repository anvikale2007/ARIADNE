"""Connectome loading and preparation.

Inputs are plain CSV files, so a connectome from ANY tractography pipeline can be used:
  weights.csv  n x n streamline counts (or any non-negative connection strength)
  labels.csv   n region names (optional)
Build them from a diffusion MRI scan with dwi_to_connectome.py.
"""
import os
import numpy as np


def prepare(C):
    """Counts -> model-ready weights: log-compress, symmetrise, zero diagonal, scale so the strongest
    region's total input is 1 (keeps hubs from being over-driven)."""
    C = np.asarray(C, dtype=float)
    W = np.log1p(np.clip(C, 0, None))
    W = (W + W.T) / 2.0
    np.fill_diagonal(W, 0.0)
    return W / W.sum(1).max()


def load(folder):
    """Load weights.csv (+ labels.csv) from a folder produced by dwi_to_connectome.py."""
    C = np.loadtxt(os.path.join(folder, "weights.csv"), delimiter=",")
    lp = os.path.join(folder, "labels.csv")
    labels = [str(x) for x in np.loadtxt(lp, dtype=str, ndmin=1)] if os.path.exists(lp) else [f"R{i}" for i in range(len(C))]
    return prepare(C), labels


def load_tvb_template():
    """Optional: 76-region template connectome shipped with The Virtual Brain (needs tvb-library + tvb-data)."""
    import tvb_data, warnings, logging
    warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
    from tvb.simulator.lab import connectivity
    p = os.path.join(os.path.dirname(tvb_data.__file__), "connectivity", "connectivity_76.zip")
    c = connectivity.Connectivity.from_file(p)
    return prepare(c.weights), [str(x) for x in c.region_labels]
