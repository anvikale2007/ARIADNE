"""Re-run the mechanism re-weighting sensitivity for every brain, on the current pharmacology.

The stored `sensitivity` blocks predate both the derived GABA weight and the addition of
levetiracetam and phenytoin, so they disagreed with the thresholds sitting next to them in the
same file. This recomputes all four variants for all drugs.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ariadne import connectome, engine  # noqa: E402

VARIANTS = {"excitability effects x0.5": (0.5, 1.0), "excitability effects x2": (2.0, 1.0),
            "coupling effects x0.5": (1.0, 0.5), "coupling effects x2": (1.0, 2.0)}
TARGETS = ["stanford_hardi_onsetP39", "sherbrooke_3shell_onsetP58", "penn_pid02_ieeg_template"]


def load_W(R):
    if R["connectome"] == "tvb_template":
        return connectome.load_tvb_template()
    return connectome.load(os.path.join(HERE, "connectomes", R["connectome"]))


def main():
    drugs = engine.load_drugs()
    for name in TARGETS:
        path = os.path.join(HERE, "results", name + ".json")
        R = json.load(open(path))
        W, _ = load_W(R)
        cal, onset = R["calibration"], R["onset_index"]
        R["sensitivity"] = {}
        for label, (xw, cw) in VARIANTS.items():
            t = engine.drug_thresholds(W, onset, cal, drugs, seeds=(0,), iters=6,
                                       x0_weight=xw, coupling_weight=cw, workers=2)
            rank = engine.rank(t)
            R["sensitivity"][label] = {"ranking": rank, "thresholds": {k: t[k]["threshold"] for k in t}}
            print(f"[{name}] {label:26s} -> {rank[0]} first", flush=True)
        json.dump(R, open(path, "w"), indent=1)
        top = {v["ranking"][0] for v in R["sensitivity"].values()}
        print(f"[{name}] top drug across variants: {', '.join(sorted(top))}", flush=True)


if __name__ == "__main__":
    main()
