"""Untreated seizure on the 10-group template network, at several coupling strengths K.
Research prototype. Settings are modelling choices, not clinical facts."""
import numpy as np, warnings
warnings.filterwarnings("ignore")
from connectome import load_template
from groups import coarse, GROUP_NAMES
from epileptor import simulate

C, labels = load_template()
M, g = coarse(C, labels)
M = np.minimum(M / M[0].max(), 1.0)           # strongest tie from the origin = 1
x0 = -2.4 + 0.3 * M[0]                        # excitability rises with tie to the origin
x0[0] = -1.6                                  # hippocampus = seizure origin

print("x0 per group:", dict(zip(GROUP_NAMES, x0.round(2))))
for K in [0.3, 0.5, 0.7, 1.0]:
    o = simulate(M, x0, K=K, T=6000, n_trials=4, seed=5)
    frac = (o > 0).any(axis=0).mean(0)
    rec = [GROUP_NAMES[i] for i in range(10) if frac[i] >= 0.5]
    print(f"K={K}: {len(rec)} of 10 groups seize -> {rec}", flush=True)
