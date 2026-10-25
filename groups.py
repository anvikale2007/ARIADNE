"""Merge the 76 template regions into 10 LED groups."""
import numpy as np
GROUP_NAMES = ["Hippocampus", "Amygdala", "L temporal", "R temporal", "L frontal",
               "R frontal", "L parietal", "R parietal", "Occipital", "Insula/cingulate"]
def group_of(lab):
    if lab in ("rHC", "rPHC"): return 0
    if lab == "rAMYG": return 1
    side, rest = lab[0], lab[1:]
    if rest in ("V1", "V2"): return 8
    if lab in ("lHC", "lAMYG", "lPHC") or rest.startswith(("TC", "A1", "A2")):
        return 2 if side == "l" else 3
    if rest.startswith(("PFC", "PMC", "M1", "FEF")): return 4 if side == "l" else 5
    if rest.startswith(("IP", "PC", "S1", "S2")): return 6 if side == "l" else 7
    return 9   # IA, CC*, G, etc.
def coarse(C, labels):
    g = np.array([group_of(l) for l in labels])
    M = np.zeros((10, 10))
    for a in range(10):
        for b in range(10):
            if a != b:
                M[a, b] = C[np.ix_(g == a, g == b)].sum()
    M = (M + M.T) / 2
    return M / M.max(), g
