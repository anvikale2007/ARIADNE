"""Verification suite: does the simulator behave the way a correct model must? (This is VERIFICATION, not clinical validation.)

  python verify_model.py [--connectome connectomes/stanford_hardi]   # ~5 min, writes results/verification_report.md
Checks: healthy control quiet | untreated seizure spreads | no-drug = baseline | dose-response is monotone |
noise-seed robustness | time-step convergence | seizure behaviour depends on onset site.
"""
import argparse, os, json, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from ariadne import connectome, engine, model, metrics

HERE = os.path.dirname(os.path.abspath(__file__))


def job(a):
    W, onset, sh, cs, seed, dt, x0h, K = a
    model.DT = dt
    r = metrics.analyse(model.run_simulation(W, onset, sh, cs, engine.LENGTH_MS, seed, K=K, x0_healthy=x0h))
    return r["recruited_regions"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--connectome", default=os.path.join(HERE, "connectomes", "stanford_hardi"))
    a = ap.parse_args()
    W, labels = connectome.load(a.connectome); n = len(labels)
    onset = int(np.argsort(W.sum(1))[n // 2])
    drugs = engine.load_drugs()
    print("calibrating ..."); cal = engine.calibrate(W, onset)
    assert cal["ok"], "calibration failed"
    x0h, K = cal["x0_healthy"], cal["K"]
    base = lambda **k: (W, k.get("onset", onset), k.get("sh", 0.0), k.get("cs", 1.0), k.get("seed", 0), k.get("dt", 0.1), x0h, K)

    jobs = {"healthy": base(onset=None)}
    for s in range(1, 4): jobs[f"baseline_seed{s}"] = base(seed=s)
    jobs["baseline_dt0.05"] = base(dt=0.05)
    jobs["baseline_dt0.1"] = base()
    jobs["no_drug_s0"] = (W, onset, *engine.drug_params(list(drugs.values())[0], 0.0), 0, 0.1, x0h, K)
    grid = [0.0, 0.5, 1.0, 2.0, 4.0]
    for nm, d in drugs.items():
        for s in grid: jobs[f"dose_{nm}_{s}"] = (W, onset, *engine.drug_params(d, s), 0, 0.1, x0h, K)
    order = np.argsort(W.sum(1)); others = [int(order[2]), int(order[n // 4]), int(order[3 * n // 4]), int(order[-3])]
    for o in others: jobs[f"onset_{o}"] = base(onset=o)

    keys = list(jobs)
    print(f"running {len(keys)} simulations ...")
    with ProcessPoolExecutor(2) as ex:
        res = dict(zip(keys, ex.map(job, [jobs[k] for k in keys])))

    checks = []
    def add(name, ok, detail): checks.append((name, bool(ok), detail))
    add("Healthy network stays quiet", res["healthy"] == 0, f"{res['healthy']} regions recruited")
    add("Untreated seizure spreads", res["baseline_seed1"] > engine.SPREAD_FRACTION * n, f"{res['baseline_seed1']}/{n} recruited")
    add("Zero-strength drug equals baseline", res["no_drug_s0"] == res["baseline_seed1"] or abs(res["no_drug_s0"] - res["baseline_seed1"]) <= 2,
        f"{res['no_drug_s0']} vs {res['baseline_seed1']}")
    seeds = [res[f"baseline_seed{s}"] for s in range(1, 4)]
    add("Result robust to noise seed", max(seeds) - min(seeds) <= 0.1 * n, f"recruited across 3 seeds: {seeds}")
    add("Result converged in time step (dt 0.1 vs 0.05)", abs(res["baseline_dt0.05"] - res["baseline_dt0.1"]) <= 0.1 * n,
        f"dt0.1 -> {res['baseline_dt0.1']}, dt0.05 -> {res['baseline_dt0.05']}")
    for nm in drugs:
        seq = [res[f"dose_{nm}_{s}"] for s in grid]
        mono = all(seq[i] >= seq[i + 1] - 2 for i in range(len(seq) - 1))
        add(f"Dose-response monotone: {nm}", mono, f"recruited at strengths {grid}: {seq}")
    ons = [res[f"onset_{o}"] for o in others]
    add("Spread depends on onset site (not identical everywhere)", len(set(ons)) > 1, f"recruited from 4 other onsets: {ons}")

    lines = ["# ARIADNE verification report", f"Connectome: {os.path.basename(a.connectome)} ({n} regions), onset {labels[onset]}, "
             f"calibrated x0={x0h}", "", "| Check | Result | Detail |", "|---|---|---|"]
    lines += [f"| {c} | {'PASS' if ok else '**FAIL**'} | {d} |" for c, ok, d in checks]
    npass = sum(ok for _, ok, _ in checks)
    lines += ["", f"**{npass}/{len(checks)} checks passed.**", "",
              "This verifies the simulator is internally consistent. It does NOT show the model predicts patient drug response."]
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    open(os.path.join(HERE, "results", "verification_report.md"), "w").write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
