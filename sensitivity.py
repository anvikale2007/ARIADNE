"""Parameter sensitivity: do the tier orderings hold when scale factors, baseline K and occupancies vary?
Occupancies are PLACEHOLDERS (jittered +/-30%). Replace with real values from occupancy_table.csv."""
import numpy as np, warnings, sys
warnings.filterwarnings("ignore")
from epileptor_fast import simulate
from connectome import load_template
from groups import coarse
from drugs import DRUGS
C, labels = load_template(); M, g = coarse(C, labels); M = np.minimum(M/M[0].max(), 1.0)
x0 = -2.4 + 0.3*M[0]; x0[0] = -1.6
rng = np.random.default_rng(0)
NA = ["Carbamazepine", "Lamotrigine", "Valproate", "Cenobamate"]
REL = ["Pregabalin", "Levetiracetam"]
N_SAMPLES = int(sys.argv[1]) if len(sys.argv) > 1 else 40
tests = {"every drug <= no drug": 0, "Na-class mean burden < release-class mean": 0,
         "every Na drug < every release drug": 0, "ethosuximide ~ no drug (built in)": 0}
kept = 0
for k in range(N_SAMPLES):
    SX, SK = rng.uniform(0.05, 0.20), rng.uniform(0.3, 0.9)
    K0 = rng.choice([0.5, 0.7, 1.0])
    b = {}
    for name, (ox, oK, _) in DRUGS.items():
        jx, jK = ox * rng.uniform(0.7, 1.3), oK * rng.uniform(0.7, 1.3)
        xd = x0 - SX * jx; xd[0] = x0[0]
        o = simulate(M, xd, K=K0 * (1 - SK * jK), n_trials=3, seed=1)
        b[name] = float((o[:, :, 1:] > 0).mean())
    base = b["No drug"]
    if base < 0.01:        # baseline barely spreads: sample uninformative
        continue
    kept += 1
    tests["every drug <= no drug"] += all(v <= base + 1e-9 for v in b.values())
    tests["Na-class mean burden < release-class mean"] += np.mean([b[n] for n in NA]) < np.mean([b[n] for n in REL])
    tests["every Na drug < every release drug"] += max(b[n] for n in NA) < min(b[n] for n in REL)
    tests["ethosuximide ~ no drug (built in)"] += abs(b["Ethosuximide"] - base) < 1e-9
print(f"informative samples: {kept} of {N_SAMPLES}")
for t, c in tests.items():
    print(f"{t:45s}{100*c/max(kept,1):6.0f}%")
