"""Run the whole drug panel, print the ranking, and write drugs.h for the LED model.
PLACEHOLDER occupancies (see drugs.py): ranking is NOT yet meaningful."""
import numpy as np, warnings
warnings.filterwarnings("ignore")
from epileptor_fast import simulate
from connectome import load_template
from groups import coarse, GROUP_NAMES
from drugs import DRUGS, apply

C, labels = load_template()
M, g = coarse(C, labels)
M = np.minimum(M / M[0].max(), 1.0)
x0 = -2.4 + 0.3 * M[0]; x0[0] = -1.6
K0, NEVER = 0.7, 255
results = {}
burdens = {}
for name in DRUGS:
    x0d, Kd = apply(name, x0, K0)
    o = simulate(M, x0d, K=Kd, n_trials=6, seed=5)           # (time, trials, 10)
    seized = (o > 0)
    hit = seized.any(axis=0).mean(0) >= 0.5                  # recruited in at least half the trials
    first = np.where(seized.any(axis=0), seized.argmax(axis=0), np.inf)
    t = np.where(hit, np.median(first, axis=0), np.inf)      # time each group is reached
    results[name] = (hit, t)
    burdens[name] = float((o[:, :, 1:] > 0).mean())
base = int(results["No drug"][0][1:].sum())
print(f"{'Drug':15s}{'groups reached':>15s}{'burden':>9s}{'suppression(burden)':>21s}  verdict")
for name, (hit, t) in sorted(results.items(), key=lambda kv: kv[1][0][1:].sum()):
    n = int(hit[1:].sum()); s = 1 - n / base
    v = "green" if n == 0 else ("amber" if n <= base / 2 else "red")
    sb = 1 - burdens[name] / burdens["No drug"]
    print(f"{name:15s}{n:>12d}/9{burdens[name]:>9.3f}{sb:>21.2f}  {v}")

tmax = max(t[np.isfinite(t)].max() for _, t in results.values())
with open("drugs.h", "w") as f:
    f.write(f"#define NUM_DRUGS {len(results)}\n#define MAX_STEP 12\n")
    f.write("const char* drugNames[] = {" + ",".join(f'"{n}"' for n in results) + "};\n")
    f.write("const uint8_t recruit[NUM_DRUGS][10] = {\n")
    for _, (hit, t) in results.items():
        f.write("  {" + ",".join(str(NEVER if np.isinf(x) else int(round(10 * x / tmax))) for x in t) + "},\n")
    f.write("};\n")
print("wrote drugs.h")
