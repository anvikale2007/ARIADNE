"""Cohort test on IDEAS II: does ARIADNE's virtual surgery predict who became seizure-free?

  python virtual_surgery_study.py --counts <folder with *_36-Lausanne_Count.csv> [--workers 8]

For each patient with a connectome, a resection table entry and an ILAE year-1 outcome:
  spread capacity before vs after virtual resection (seizure seeded in every region of the operated hemisphere in turn).
Then compares the reduction between seizure-free (ILAE 1-2) and not seizure-free (ILAE 3-5) patients
(Mann-Whitney U, AUC). Results are appended per patient to results/virtual_surgery.csv, so the run can be resumed.
Pre-registered expectation: larger reduction in seizure-free patients. A null result is a valid, reportable outcome.
~3 min per patient on 2 cores; scales with --workers.
"""
import argparse, os, time, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
import numpy as np, pandas as pd
from ariadne import ideas, connectome, surgery, engine

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "virtual_surgery.csv")
X0H = -2.10


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=ideas.DATA); ap.add_argument("--workers", type=int, default=os.cpu_count() or 2)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    meta = ideas.metadata().drop_duplicates("ID").set_index("ID")  # metadata contains 1 duplicated ID
    lab = ideas.labels()
    res_table = pd.read_csv(os.path.join(ideas.DATA, "table_resected.csv"))
    subs = [s for s in ideas.available_subjects(a.counts) if s in meta.index and str(s) in res_table.columns
            and not np.isnan(meta.loc[s, "ILAE_Year1"])]
    done = set(pd.read_csv(OUT)["ID"]) if os.path.exists(OUT) else set()
    todo = [s for s in subs if s not in done][: a.limit]
    print(f"{len(subs)} eligible patients, {len(done)} done, {len(todo)} to run")
    for s in todo:
        t = time.time()
        frac = ideas.resected_fraction(s)
        if np.all(np.isnan(frac)):
            continue
        C = np.loadtxt(ideas.count_file(s, a.counts), delimiter=",")
        W = connectome.prepare(C)
        side = meta.loc[s, "Op_Side"]
        hemi = [i for i, l in enumerate(lab) if (l.startswith("Left") or "ctx-lh" in l) == (side == "L")]
        pre, _ = surgery.spread_capacity(W, X0H, engine.K_COUPLING, hemi, workers=a.workers)
        post, _ = surgery.spread_capacity(W, X0H, engine.K_COUPLING, hemi, frac, workers=a.workers)
        row = {"ID": s, "ILAE_Year1": meta.loc[s, "ILAE_Year1"], "seizure_free": meta.loc[s, "ILAE_Year1"] <= 2,
               "Op_Side": side, "Op_Type": meta.loc[s, "Op_Type"], "Pathology": meta.loc[s, "Pathology"],
               "Number_ASMs": meta.loc[s, "Number_ASMs"], "capacity_pre": pre, "capacity_post": post,
               "reduction": pre - post, "resected_regions_ge50pct": int(np.sum(np.nan_to_num(frac) >= 0.5))}
        pd.DataFrame([row]).to_csv(OUT, mode="a", header=not os.path.exists(OUT), index=False)
        print(f"sub-{s}: ILAE {row['ILAE_Year1']:.0f}  reduction {row['reduction']:.3f}  ({time.time()-t:.0f}s)", flush=True)
    summarise()


def summarise():
    if not os.path.exists(OUT):
        return
    d = pd.read_csv(OUT)
    sf, nsf = d[d.seizure_free].reduction, d[~d.seizure_free].reduction
    lines = ["# Virtual surgery vs real outcome (IDEAS II)", f"Patients: {len(d)} (seizure-free {len(sf)}, not {len(nsf)})",
             f"Mean reduction in spread capacity: seizure-free {sf.mean():.3f}, not seizure-free {nsf.mean():.3f}" if len(sf) and len(nsf) else "Need both outcome groups."]
    if len(sf) >= 3 and len(nsf) >= 3:
        from scipy.stats import mannwhitneyu
        u, p = mannwhitneyu(sf, nsf, alternative="greater")
        auc = u / (len(sf) * len(nsf))
        lines += [f"Mann-Whitney U (seizure-free > not): U={u:.0f}, one-sided p={p:.3g}, AUC={auc:.2f}",
                  "Confound check: number of resected regions also differs by group? -> compare 'resected_regions_ge50pct'.",
                  f"Resected regions >=50%: seizure-free {d[d.seizure_free].resected_regions_ge50pct.mean():.1f}, not {d[~d.seizure_free].resected_regions_ge50pct.mean():.1f}"]
    lines.append("Exploratory analysis; model parameters not tuned to outcomes. AUC 0.5 = no better than chance.")
    open(os.path.join(HERE, "results", "virtual_surgery_report.md"), "w").write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
