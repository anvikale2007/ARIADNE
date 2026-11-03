"""Virtual surgery vs real outcome, per the plan fixed in PREREGISTRATION.md (31 Oct 2026).

    python virtual_surgery2.py [--workers 2] [--limit N]

Differences from the pilot (virtual_surgery_study.py), both decided from pre-surgery diagnostics with
no reference to outcome labels:

  1. Background excitability is calibrated PER PATIENT so that pre-surgery spread capacity sits near
     0.50 instead of saturating at 0.93-1.00. Calibration uses only the pre-surgery connectome.
  2. Partial resections are graded: a region's excitability moves toward the silenced value in
     proportion to the fraction removed. The original binary model is kept as a sensitivity analysis.

Primary outcome: reduction in spread capacity, compared between ILAE 1-2 and ILAE 3-5 by one-sided
Mann-Whitney U, reported with AUC. Results append to results/virtual_surgery2.csv so the run resumes.
"""
import argparse, json, os, time, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
import numpy as np, pandas as pd
from ariadne import ideas, connectome, surgery, engine

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "virtual_surgery2.csv")
TARGET_CAPACITY = 0.50      # pre-registered calibration target
CAL_ONSETS = 12             # pre-registered subsample size for calibration
CAL_ITERS = 6
X0_LO, X0_HI = -2.40, -1.95


def calibrate(W, K, hemi, workers, rng):
    """Lowest-error background excitability bringing PRE-surgery capacity near TARGET_CAPACITY.
    Uses the pre-surgery connectome only: no resection, no outcome label."""
    sub = sorted(rng.choice(hemi, size=min(CAL_ONSETS, len(hemi)), replace=False).tolist())
    lo, hi, trace = X0_LO, X0_HI, []
    best = None
    for _ in range(CAL_ITERS):
        mid = (lo + hi) / 2
        cap, _ = surgery.spread_capacity(W, mid, K, sub, workers=workers)
        trace.append((round(mid, 4), round(cap, 3)))
        if best is None or abs(cap - TARGET_CAPACITY) < abs(best[1] - TARGET_CAPACITY):
            best = (mid, cap)
        if cap > TARGET_CAPACITY:
            hi = mid          # too excitable
        else:
            lo = mid
    return round(best[0], 4), round(best[1], 3), trace


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=ideas.DATA)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    meta = ideas.metadata().drop_duplicates("ID").set_index("ID")
    lab = ideas.labels()
    res_table = pd.read_csv(os.path.join(ideas.DATA, "table_resected.csv"))
    subs = [s for s in ideas.available_subjects(a.counts) if s in meta.index and str(s) in res_table.columns
            and not np.isnan(meta.loc[s, "ILAE_Year1"])]
    done = set(pd.read_csv(OUT)["ID"]) if os.path.exists(OUT) else set()
    todo = [s for s in subs if s not in done][: a.limit]
    print(f"{len(subs)} eligible, {len(done)} done, {len(todo)} to run", flush=True)

    excluded = []
    for s in todo:
        t = time.time()
        frac = ideas.resected_fraction(s)
        if np.all(np.isnan(frac)):
            excluded.append({"ID": int(s), "reason": "no resection data"}); continue
        C = np.loadtxt(ideas.count_file(s, a.counts), delimiter=",")
        W = connectome.prepare(C)
        side = meta.loc[s, "Op_Side"]
        hemi = [i for i, l in enumerate(lab) if (l.startswith("Left") or "ctx-lh" in l) == (side == "L")]
        if len(hemi) < 4:
            excluded.append({"ID": int(s), "reason": f"operated hemisphere has {len(hemi)} regions"}); continue

        rng = np.random.default_rng(0)
        x0h, cal_cap, trace = calibrate(W, engine.K_COUPLING, hemi, a.workers, rng)
        pre, _ = surgery.spread_capacity(W, x0h, engine.K_COUPLING, hemi, workers=a.workers)
        post_g, _ = surgery.spread_capacity(W, x0h, engine.K_COUPLING, hemi, frac, workers=a.workers, graded=True)
        post_b, _ = surgery.spread_capacity(W, x0h, engine.K_COUPLING, hemi, frac, workers=a.workers, graded=False)
        row = {"ID": s, "ILAE_Year1": meta.loc[s, "ILAE_Year1"], "seizure_free": meta.loc[s, "ILAE_Year1"] <= 2,
               "Op_Side": side, "Op_Type": meta.loc[s, "Op_Type"], "Pathology": meta.loc[s, "Pathology"],
               "x0_calibrated": x0h, "cal_capacity": cal_cap,
               "capacity_pre": round(pre, 4), "capacity_post_graded": round(post_g, 4),
               "capacity_post_binary": round(post_b, 4),
               "reduction": round(pre - post_g, 4), "reduction_binary": round(pre - post_b, 4),
               "resected_regions_ge50pct": int(np.sum(np.nan_to_num(frac) >= 0.5)),
               "resected_fraction_total": round(float(np.nansum(np.nan_to_num(frac))), 3),
               "hemi_regions": len(hemi)}
        pd.DataFrame([row]).to_csv(OUT, mode="a", header=not os.path.exists(OUT), index=False)
        print(f"sub-{s}: ILAE {row['ILAE_Year1']:.0f}  x0 {x0h}  pre {pre:.3f} -> post {post_g:.3f} "
              f"(reduction {row['reduction']:.3f})  [{time.time()-t:.0f}s]", flush=True)
    if excluded:
        json.dump(excluded, open(os.path.join(HERE, "results", "virtual_surgery2_excluded.json"), "w"), indent=1)
    summarise()


def summarise():
    if not os.path.exists(OUT):
        print("no results yet"); return
    d = pd.read_csv(OUT).drop_duplicates("ID")
    from scipy.stats import mannwhitneyu
    L = ["# Virtual surgery vs real outcome (IDEAS II) - pre-registered analysis",
         "Plan fixed in PREREGISTRATION.md before this analysis was run.", "",
         f"Patients: {len(d)} (seizure-free {int(d.seizure_free.sum())}, not seizure-free {int((~d.seizure_free).sum())})",
         f"Per-patient calibration brought pre-surgery spread capacity to "
         f"{d.capacity_pre.min():.2f}-{d.capacity_pre.max():.2f} (pilot was 0.93-1.00, saturated).", ""]

    def test(col, label):
        sf, nsf = d.loc[d.seizure_free, col], d.loc[~d.seizure_free, col]
        if len(sf) < 2 or len(nsf) < 2:
            return [f"**{label}**: need both groups."]
        u, p = mannwhitneyu(sf, nsf, alternative="greater")
        auc = u / (len(sf) * len(nsf))
        verdict = ("H1 supported at this sample size" if auc > 0.5 and p < 0.05 else
                   "direction consistent with H1, not statistically supported" if auc > 0.5 else
                   "H1 not supported")
        return [f"**{label}**", f"- mean: seizure-free {sf.mean():.3f}, not seizure-free {nsf.mean():.3f}",
                f"- Mann-Whitney U={u:.0f}, one-sided p={p:.3g}, **AUC={auc:.2f}**",
                f"- Pre-registered decision rule: {verdict}.", ""]

    L += test("reduction", "PRIMARY: reduction in spread capacity (graded resection)")
    L += test("reduction_binary", "SENSITIVITY: same, original binary resection model")
    nsf_n = int((~d.seizure_free).sum()); sf_n = int(d.seizure_free.sum())
    if min(sf_n, nsf_n) >= 1:
        from math import comb
        mn = 1 / comb(sf_n + nsf_n, min(sf_n, nsf_n)) if sf_n + nsf_n < 30 else 0
        L += [f"**Power.** Smallest attainable one-sided p at n={sf_n} vs {nsf_n} is {mn:.3f}. "
              f"80% power needs about 15 per group at AUC 0.80, 25 per group at AUC 0.70 "
              f"(simulated, see PREREGISTRATION.md). This sample is a pilot and cannot confirm the hypothesis.", ""]
    L += ["**Confounds (reported regardless of significance)**", "",
          "| Variable | Seizure-free | Not seizure-free |", "|---|---|---|"]
    for col, nm in [("resected_regions_ge50pct", "Regions >=50% resected"),
                    ("resected_fraction_total", "Total resected fraction"),
                    ("capacity_pre", "Pre-surgery capacity"), ("x0_calibrated", "Calibrated excitability")]:
        L.append(f"| {nm} | {d.loc[d.seizure_free, col].mean():.3f} | {d.loc[~d.seizure_free, col].mean():.3f} |")
    L += ["", "| Patient | ILAE | Surgery | Pre | Post | Reduction |", "|---|---|---|---|---|---|"]
    for _, r in d.sort_values("ID").iterrows():
        L.append(f"| {r.ID:.0f} | {r.ILAE_Year1:.0f} | {r.Op_Type} | {r.capacity_pre:.3f} | "
                 f"{r.capacity_post_graded:.3f} | {r.reduction:.3f} |")
    L += ["", "Exploratory beyond the primary test. Model parameters were never tuned to outcomes.",
          "The earlier pilot (fixed excitability, binary resection, n=6) gave AUC 0.33 and is reported alongside."]
    open(os.path.join(HERE, "results", "virtual_surgery2_report.md"), "w").write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
