"""Bring every brain's result file onto the SAME pharmacology.

Why this exists: the GABA weight used to be an assumed 0.5 and is now derived from the
excitatory/inhibitory conductance ratio (gaba_weight = g_I/g_E = 1.0, see drug_sources.yaml).
Cenobamate is the only drug with a GABA term, so its threshold changed and had to be re-run
everywhere; levetiracetam and phenytoin were added after these files were first written.

Running this leaves all four brains comparable, which is the only way the cross-brain table
in the paper and the Explorer means anything.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ariadne import connectome, engine  # noqa: E402

TARGETS = ["stanford_hardi_onsetP39", "sherbrooke_3shell_onsetP58", "penn_pid02_ieeg_template"]
# cenobamate is re-run because its x0 effect changed; the other two were simply absent before.
REDO = ["cenobamate"]
ADD = ["levetiracetam", "phenytoin"]


def load_W(R):
    if R["connectome"] == "tvb_template":
        return connectome.load_tvb_template()
    return connectome.load(os.path.join(HERE, "connectomes", R["connectome"]))


def main():
    drugs = engine.load_drugs()
    for name in TARGETS:
        path = os.path.join(HERE, "results", name + ".json")
        R = json.load(open(path))
        W, labels = load_W(R)
        cal, onset = R["calibration"], R["onset_index"]
        want = [d for d in REDO + ADD if d in drugs]
        sub = {k: drugs[k] for k in want}
        print(f"[{name}] n={len(labels)} onset={R['onset']} recomputing {want}", flush=True)
        th = engine.drug_thresholds(W, onset, cal, sub, seeds=(0,), iters=6, workers=1)
        snap = engine.snapshot(W, onset, cal, sub, th)
        R.setdefault("snapshot", {})
        for k in want:
            R["drugs"][k] = {**th[k], "targets": sub[k]["targets"], "rationale": sub[k]["rationale"]}
            R["snapshot"][k] = snap[k]
            rec = sum(snap[k]["recruited_mask"])
            print(f"   {k:16s} threshold={th[k]['threshold']}  recruited at standard dose {rec}/{len(labels)}", flush=True)
        R["ranking"] = engine.rank(R["drugs"])
        R["pharmacology_note"] = ("thresholds use the derived GABA weight g_I/g_E = 1.0 "
                                  "(drug_sources.yaml); cenobamate re-run, levetiracetam and phenytoin added")
        json.dump(R, open(path, "w"), indent=1)
        print(f"[{name}] ranking: " + " > ".join(R["ranking"]), flush=True)


if __name__ == "__main__":
    main()
