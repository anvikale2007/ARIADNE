"""Analysis A of PREREGISTRATION.md Amendment 4: does removing more NETWORK-abnormal tissue predict
seizure freedom? Pre-registered before running.

Unlike abnormality_study.py (morphometry, refuted at n=405), this uses the structural CONNECTOME:
each patient's regional nodal strength is z-scored against 97 control connectomes.

    python network_abnormality.py
"""
import os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
from ariadne import ideas, connectome

HERE = os.path.dirname(os.path.abspath(__file__))


def strengths(sub):
    C = np.loadtxt(ideas.count_file(sub), delimiter=",")
    return connectome.prepare(C).sum(axis=1)          # nodal strength after ARIADNE's own normalisation


def main():
    subs = sorted(ideas.available_subjects())
    controls = [s for s in subs if s >= 4000]
    meta = ideas.metadata().drop_duplicates("ID").set_index("ID")
    rt = pd.read_csv(os.path.join(ideas.DATA, "table_resected.csv"))
    rt.index = [str(x).strip().strip("'") for x in rt.iloc[:, 0]]
    rt = rt.drop(columns=rt.columns[0])

    print(f"building normative distribution from {len(controls)} controls ...")
    S = np.vstack([strengths(s) for s in controls])
    mu, sd = S.mean(0), S.std(0, ddof=1)
    sd[sd == 0] = np.nan

    rows = []
    for s in subs:
        if s >= 4000 or s not in meta.index or str(s) not in rt.columns:
            continue
        ilae = meta.loc[s, "ILAE_Year1"]
        if pd.isna(ilae):
            continue
        z = np.abs((strengths(s) - mu) / sd)
        f = rt[str(s)].astype(float).fillna(0).clip(0, 1).values
        if len(f) != len(z):
            continue
        tot = f.sum()
        if tot <= 0:
            continue
        rows.append({"ID": s, "ILAE": ilae, "seizure_free": ilae <= 2,
                     "vol_resected": tot,
                     "net_abnormality_density": float(np.nansum(f * z) / tot),
                     "net_abnormality_removed": float(np.nansum(f * z)),
                     "net_global_abnormality": float(np.nanmean(z)),
                     "side": meta.loc[s, "Op_Side"], "op": str(meta.loc[s, "Op_Type"])})
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(HERE, "results", "network_abnormality.csv"), index=False)
    n_sf = int(d.seizure_free.sum())
    print(f"\nn = {len(d)}  ({n_sf} seizure-free, {len(d)-n_sf} not)\n")

    def test(col, greater_in_sf, label, dd=None):
        dd = d if dd is None else dd
        sf, nsf = dd.loc[dd.seizure_free, col], dd.loc[~dd.seizure_free, col]
        a, b = (sf, nsf) if greater_in_sf else (nsf, sf)
        u, p = mannwhitneyu(a, b, alternative="greater")
        auc = u / (len(a) * len(b))
        star = "  *" if p < 0.025 else ""
        print(f"  {label:54s} SF {sf.mean():7.3f} vs NSF {nsf.mean():7.3f}   AUC {auc:.3f}  p={p:.4g}{star}")
        return {"measure": label, "n": len(dd), "AUC": auc, "p": p, "mean_sf": sf.mean(), "mean_nsf": nsf.mean()}

    out = [test("net_abnormality_density", True,  "H4 network abnormality of resected tissue (PRIMARY)"),
           test("net_global_abnormality",  False, "H5 whole-brain network abnormality (higher = worse)"),
           test("vol_resected",            True,  "confound: resected volume alone"),
           test("net_abnormality_removed", True,  "secondary: total network abnormality removed")]
    X = pd.get_dummies(d[["net_abnormality_density", "vol_resected", "side"]], columns=["side"], drop_first=True)
    m = sm.Logit(d.seizure_free.astype(int), sm.add_constant(X.astype(float))).fit(disp=0)
    print(f"\n  adjusted logistic (density | volume, side): beta={m.params['net_abnormality_density']:+.3f}  "
          f"p={m.pvalues['net_abnormality_density']:.4g}")
    t = d[d.op.str.startswith("T")]
    if len(t) > 30:
        out.append(test("net_abnormality_density", True, f"H4 restricted to temporal surgery (n={len(t)})", t))
    pd.DataFrame(out).to_csv(os.path.join(HERE, "results", "network_abnormality_tests.csv"), index=False)
    print("\nSaved results/network_abnormality*.csv")


if __name__ == "__main__":
    main()
