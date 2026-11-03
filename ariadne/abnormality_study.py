"""Does removing MORE ABNORMAL tissue predict seizure freedom? (IDEAS, n~415)

Pre-registered in PREREGISTRATION.md Amendment 3, written before this was run.

NOTE: this analysis does NOT use the ARIADNE model - no connectome, no Epileptor, no drugs. It tests the
principle behind Amendment 1's proposed fix (what matters is how abnormal the removed tissue was, not how
much) on the first adequately powered sample available. Methodologically it is close to Sinha (2021).

    python abnormality_study.py
"""
import os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.api as sm

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data", "ideas")
Z = os.path.join(D, "zscores")
METRICS = {"thickness": "thick", "area": "area", "volume": "vol"}


def load(metric):
    zp = pd.read_csv(os.path.join(Z, f"{METRICS[metric]}_zscores_patients.csv"))
    rt = pd.read_csv(os.path.join(D, "table_resected.csv"))
    rt.index = [str(x).strip().strip("'") for x in rt.iloc[:, 0]]
    rt = rt.drop(columns=rt.columns[0])
    # map ctx-lh-bankssts <-> lh_bankssts ; keep only cortical regions present in both
    zcols = {c: "ctx-" + c.split("_")[0] + "-" + c.split("_")[1] for c in zp.columns if c != "ID"}
    shared = [(zc, rc) for zc, rc in zcols.items() if rc in rt.index]
    meta = pd.read_csv(os.path.join(D, "Metadata_Release_Anon.csv")).drop_duplicates("ID").set_index("ID")
    rows = []
    for _, r in zp.iterrows():
        pid = int(r["ID"])
        if str(pid) not in rt.columns or pid not in meta.index:
            continue
        ilae = meta.loc[pid, "ILAE_Year1"]
        if pd.isna(ilae):
            continue
        f = rt[str(pid)].reindex([rc for _, rc in shared]).astype(float).fillna(0).clip(0, 1).values
        z = np.abs(pd.to_numeric(r[[zc for zc, _ in shared]], errors="coerce").astype(float).fillna(0).values)
        tot = f.sum()
        rows.append({"ID": pid, "ILAE": ilae, "seizure_free": ilae <= 2,
                     "vol_resected": tot,
                     "abnormality_removed": float((f * z).sum()),
                     "abnormality_density": float((f * z).sum() / tot) if tot > 0 else np.nan,
                     "global_abnormality": float(z.mean()),
                     "side": meta.loc[pid, "Op_Side"], "op": meta.loc[pid, "Op_Type"]})
    return pd.DataFrame(rows).dropna(subset=["abnormality_density"]), len(shared)


def test(d, col, greater_in_sf, label):
    sf, nsf = d.loc[d.seizure_free, col], d.loc[~d.seizure_free, col]
    a, b = (sf, nsf) if greater_in_sf else (nsf, sf)
    u, p = mannwhitneyu(a, b, alternative="greater")
    auc = u / (len(a) * len(b))
    return {"measure": label, "col": col, "n_sf": len(sf), "n_nsf": len(nsf),
            "mean_sf": sf.mean(), "mean_nsf": nsf.mean(), "AUC": auc, "p": p,
            "sig_bonferroni_0.0167": p < 0.0167}


def main():
    out = []
    for metric in METRICS:
        d, nreg = load(metric)
        print(f"\n=== {metric.upper()}  (n={len(d)}, {int(d.seizure_free.sum())} seizure-free, "
              f"{int((~d.seizure_free).sum())} not; {nreg} cortical regions matched)")
        res = [test(d, "abnormality_density", True, "H1 abnormality density of resected tissue"),
               test(d, "global_abnormality", False, "H2 whole-brain abnormality (higher = worse outcome)"),
               test(d, "vol_resected", True, "H3 resected volume alone (confound check)"),
               test(d, "abnormality_removed", True, "total abnormality removed (secondary)")]
        for r in res:
            print(f"  {r['measure']:52s} SF {r['mean_sf']:7.3f} vs NSF {r['mean_nsf']:7.3f}   "
                  f"AUC {r['AUC']:.3f}  p={r['p']:.4g}{'  *' if r['sig_bonferroni_0.0167'] else ''}")
            r["metric"] = metric
            out.append(r)
        # adjusted model
        X = pd.get_dummies(d[["abnormality_density", "vol_resected", "side"]], columns=["side"], drop_first=True)
        X = sm.add_constant(X.astype(float))
        m = sm.Logit(d.seizure_free.astype(int), X).fit(disp=0)
        print(f"  adjusted logistic (density | volume, side): beta={m.params['abnormality_density']:+.3f}  "
              f"p={m.pvalues['abnormality_density']:.4g}")
        temporal = d[d.op.astype(str).str.startswith("T")]
        if len(temporal) > 30:
            t = test(temporal, "abnormality_density", True, "H1 restricted to temporal surgery")
            print(f"  temporal only (n={len(temporal)}): AUC {t['AUC']:.3f}  p={t['p']:.4g}")
            t["metric"] = metric; out.append(t)
    pd.DataFrame(out).to_csv(os.path.join(HERE, "results", "abnormality_study.csv"), index=False)
    print("\nSaved results/abnormality_study.csv")


if __name__ == "__main__":
    main()
