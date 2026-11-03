"""Analysis B of PREREGISTRATION.md Amendment 4: virtual surgery vs outcome at n~183.

Same hypothesis and test as the n=6 pre-registered run; the only change is the compute trade fixed in
Amendment 4 (15-onset paired subsample, 4x10 calibration, binary sensitivity dropped).

Patients are processed in an INTERLEAVED order (seizure-free, not, seizure-free, ...) so that if the run is
stopped early the completed sample is still balanced. Resumable: results append per patient.

    python virtual_surgery3.py [--workers 2] [--sweep 15]
"""
import argparse, os, time, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
import numpy as np, pandas as pd
from ariadne import ideas, connectome, surgery, engine

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "virtual_surgery3.csv")
TARGET, CAL_ITERS, CAL_ONSETS = 0.50, 4, 10


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--sweep", type=int, default=15)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    meta = ideas.metadata().drop_duplicates("ID").set_index("ID")
    lab = ideas.labels()
    rt = pd.read_csv(os.path.join(ideas.DATA, "table_resected.csv"))
    subs = [s for s in ideas.available_subjects()
            if s < 4000 and s in meta.index and str(s) in rt.columns and not np.isnan(meta.loc[s, "ILAE_Year1"])]
    sf = [s for s in subs if meta.loc[s, "ILAE_Year1"] <= 2]
    nsf = [s for s in subs if meta.loc[s, "ILAE_Year1"] > 2]
    order = [x for pair in zip(nsf, sf) for x in pair] + sf[len(nsf):] + nsf[len(sf):]   # interleaved
    done = set(pd.read_csv(OUT)["ID"]) if os.path.exists(OUT) else set()
    todo = [s for s in order if s not in done][: a.limit]
    print(f"{len(subs)} eligible ({len(sf)} SF, {len(nsf)} NSF); {len(done)} done; {len(todo)} to run", flush=True)

    for i, s in enumerate(todo, 1):
        t0 = time.time()
        try:
            frac = ideas.resected_fraction(s)
            if np.all(np.isnan(frac)):
                continue
            W = connectome.prepare(np.loadtxt(ideas.count_file(s), delimiter=","))
            side = meta.loc[s, "Op_Side"]
            hemi = [k for k, l in enumerate(lab) if (l.startswith("Left") or "ctx-lh" in l) == (side == "L")]
            if len(hemi) < 4:
                continue
            rng = np.random.default_rng(0)
            sweep = sorted(rng.choice(hemi, size=min(a.sweep, len(hemi)), replace=False).tolist())
            cal = sweep[: min(CAL_ONSETS, len(sweep))]
            lo, hi, best = -2.40, -1.95, None
            for _ in range(CAL_ITERS):
                mid = (lo + hi) / 2
                cap, _ = surgery.spread_capacity(W, mid, engine.K_COUPLING, cal, workers=a.workers)
                if best is None or abs(cap - TARGET) < abs(best[1] - TARGET):
                    best = (mid, cap)
                lo, hi = (lo, mid) if cap > TARGET else (mid, hi)
            x0h = round(best[0], 4)
            pre, _ = surgery.spread_capacity(W, x0h, engine.K_COUPLING, sweep, workers=a.workers)
            post, _ = surgery.spread_capacity(W, x0h, engine.K_COUPLING, sweep, frac, workers=a.workers, graded=True)
            row = {"ID": s, "ILAE_Year1": meta.loc[s, "ILAE_Year1"], "seizure_free": meta.loc[s, "ILAE_Year1"] <= 2,
                   "Op_Side": side, "Op_Type": meta.loc[s, "Op_Type"], "Pathology": meta.loc[s, "Pathology"],
                   "x0_calibrated": x0h, "cal_capacity": round(best[1], 3),
                   "capacity_pre": round(pre, 4), "capacity_post": round(post, 4),
                   "reduction": round(pre - post, 4),
                   "resected_fraction_total": round(float(np.nansum(np.nan_to_num(frac))), 3),
                   "n_onsets": len(sweep)}
            pd.DataFrame([row]).to_csv(OUT, mode="a", header=not os.path.exists(OUT), index=False)
            print(f"[{i}/{len(todo)}] sub-{s}: ILAE {row['ILAE_Year1']:.0f} "
                  f"pre {pre:.3f} -> post {post:.3f}  red {row['reduction']:.3f}  [{time.time()-t0:.0f}s]", flush=True)
        except Exception as e:
            print(f"[{i}/{len(todo)}] sub-{s} FAILED: {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    main()
