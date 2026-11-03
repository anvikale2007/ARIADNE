"""Compute lacosamide's containment threshold on every brain already in results/."""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ariadne import connectome, engine, ideas  # noqa: E402

FILES = ["results/stanford_hardi_onsetP39.json", "results/sherbrooke_3shell_onsetP58.json",
         "results/penn_pid02_ieeg_template.json", "results/ideas_sub-1_onset_rh-entorhinal.json"]
DRUG = "lacosamide"


def load_W(R):
    c = R["connectome"]
    if c == "tvb_template":
        return connectome.load_tvb_template()[0]
    if c.startswith("ideas"):
        return connectome.prepare(np.loadtxt(ideas.count_file(1), delimiter=","))
    return connectome.load(os.path.join(HERE, "connectomes", c))[0]


def main():
    drugs = engine.load_drugs()
    sub = {DRUG: drugs[DRUG]}
    for rel in FILES:
        p = os.path.join(HERE, rel)
        R = json.load(open(p))
        W = load_W(R)
        th = engine.drug_thresholds(W, R["onset_index"], R["calibration"], sub, seeds=(0,), iters=6, workers=2)
        sn = engine.snapshot(W, R["onset_index"], R["calibration"], sub, th)
        R["drugs"][DRUG] = {**th[DRUG], "targets": sub[DRUG]["targets"], "rationale": sub[DRUG]["rationale"]}
        R.setdefault("snapshot", {})[DRUG] = sn[DRUG]
        R["ranking"] = engine.rank(R["drugs"])
        json.dump(R, open(p, "w"), indent=1)
        print(f"{os.path.basename(rel):42s} {DRUG}={th[DRUG]['threshold']}   {' > '.join(R['ranking'])}", flush=True)


if __name__ == "__main__":
    main()
