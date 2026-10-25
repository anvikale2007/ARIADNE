"""Panel with real occupancies where sourced. f = state-dependence enhancement of Na+ blocker affinity (unknown; swept)."""
import numpy as np, warnings, sys
warnings.filterwarnings("ignore")
from epileptor_fast import simulate
from connectome import load_template
from groups import coarse
C, labels = load_template(); M, g = coarse(C, labels); M = np.minimum(M/M[0].max(), 1.0)
x0 = -2.4 + 0.3*M[0]; x0[0] = -1.6
def occ(C_free, Ki, f=1.0): return C_free / (C_free + Ki / f)
def drug_table(f):
    # name: (occ_on_x0, occ_on_K, status)
    return {
        "No drug":       (0.0, 0.0, "-"),
        "Carbamazepine": (occ(4.2, 140, f), 0.0, "real C & tonic IC50 x f"),
        "Lamotrigine":   (occ(5.3, 91, f), 0.0, "real C & tonic IC50 x f"),
        "Valproate":     (0.4, 0.3, "PLACEHOLDER"),
        "Cenobamate":    (0.7, 0.0, "PLACEHOLDER"),
        "Pregabalin":    (0.0, occ(12.6, 0.019), "real"),
        "Levetiracetam": (0.0, occ(53, 0.79), "real"),
        "Ethosuximide":  (0.0, 0.0, "assumed zero"),
    }
def run(dt_, SX, SK, K0, trials=4, seed=5):
    out = {}
    for n, (ox, oK, _) in dt_.items():
        xd = x0 - SX*ox; xd[0] = x0[0]
        o = simulate(M, xd, K=K0*(1 - SK*oK), n_trials=trials, seed=seed)
        out[n] = (int(((o > 0).any(axis=0).mean(0)[1:] >= 0.5).sum()), float((o[:, :, 1:] > 0).mean()))
    return out
if __name__ == "__main__":
    for f in (1, 10, 30, 100):
        t = drug_table(f); r = run(t, 0.10, 0.60, 0.7)
        print(f"\n--- f = {f}  (SX 0.10, SK 0.60, K0 0.7) ---")
        print(f"{'drug':15s}{'occ_x0':>8s}{'occ_K':>8s}{'groups':>8s}{'burden':>8s}  basis")
        for n in t:
            print(f"{n:15s}{t[n][0]:8.3f}{t[n][1]:8.3f}{r[n][0]:>6d}/9{r[n][1]:8.3f}  {t[n][2]}")
