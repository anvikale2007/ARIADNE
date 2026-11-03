"""Does the ranking survive changes to the model's own settings and to the assumed mapping choices?

  python robustness.py            # ~25 min on 2 cores -> results/robustness_report.md
Varies: coupling strength K (re-calibrating each time) and the assumed GABA weight used to turn GABA-A potentiation into
excitability reduction. A ranking that only holds for one setting should not be presented as a finding.
"""
import os, json, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
import numpy as np
from ariadne import connectome, engine
from derive_parameters import derive

HERE = os.path.dirname(os.path.abspath(__file__))


def run(W, onset, K, drugs):
    cal = engine.calibrate(W, onset, K=K)
    if not cal["ok"]:
        return None, None
    th = engine.drug_thresholds(W, onset, cal, drugs, iters=5)
    return engine.rank(th), {k: th[k]["threshold"] for k in th}


def main():
    W, labels = connectome.load(os.path.join(HERE, "connectomes", "stanford_hardi"))
    onset = int(np.argsort(W.sum(1))[len(labels) // 2])
    rows = []
    base_drugs = engine.load_drugs()
    for K in (-4.0, -6.0):   # K = -5 is the main analysis
        print("K", K, flush=True); r, t = run(W, onset, K, base_drugs); rows.append((f"coupling K = {K}", r, t))
    for gw in (0.0, 1.0):
        print("gaba weight", gw, flush=True)
        d, _, _ = derive(gw); r, t = run(W, onset, engine.K_COUPLING, d); rows.append((f"GABA weight = {gw}", r, t))
    print("high exposure", flush=True)
    d, _, _ = derive(high_exposure=True); r, t = run(W, onset, engine.K_COUPLING, d); rows.append(("400 mg/day exposure (cenobamate, pregabalin)", r, t))
    print("persistent mode", flush=True)
    d, _, _ = derive(nav_mode="persistent"); r, t = run(W, onset, engine.K_COUPLING, d); rows.append(("Na block = persistent current (old model)", r, t))
    tops = [r[0] for _, r, _ in rows if r]
    lines = ["# ARIADNE robustness to model settings", f"Connectome stanford_hardi, onset {labels[onset]}.", "",
             "| Variant | Ranking (best to worst) | Thresholds |", "|---|---|---|"]
    for name, r, t in rows:
        lines.append(f"| {name} | {' > '.join(r) if r else 'calibration failed'} | " +
                     (", ".join(f"{k} {('>4' if v is None else round(v,2))}" for k, v in t.items()) if t else "") + " |")
    lines += ["", f"Top drug in {tops.count(tops[0])}/{len(tops)} variants: {tops[0]}." if tops else "No variant calibrated."]
    open(os.path.join(HERE, "results", "robustness_report.md"), "w").write("\n".join(lines))
    json.dump([{"variant": n, "ranking": r, "thresholds": t} for n, r, t in rows], open(os.path.join(HERE, "results", "robustness.json"), "w"), indent=1)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
