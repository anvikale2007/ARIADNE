"""Containment boundary for ONE brain at a few chosen coupling scales.

combination.py needs required_x0(coupling_scale) for pairs that include a release-reducing drug
(pregabalin, levetiracetam). The full 81x21 surface for the IDEAS case was never finished, and the
whole grid is not needed: only a handful of coupling scales are, because the coupling reduction a
realistic pair reaches is small. This bisects the boundary at each requested scale directly.

    python combo_boundary.py --result results/ideas_sub-1_onset_rh-entorhinal.json
"""
import argparse, json, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                                          # noqa: E402
from ariadne import connectome, engine, model, metrics       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def spreads(W, onset, x0_red, c_scale, x0h, K, seed=0):
    r = metrics.analyse(model.run_simulation(W, onset, -x0_red, c_scale, engine.LENGTH_MS,
                                             seed, K=K, x0_healthy=x0h))
    return r["recruited_regions"] > engine.SPREAD_FRACTION * W.shape[0]


def boundary_at(W, onset, c_scale, x0h, K, hi=0.012, iters=7):
    """Smallest excitability reduction that contains the seizure at this coupling scale."""
    if not spreads(W, onset, hi, c_scale, x0h, K):
        lo, h = 0.0, hi
        for _ in range(iters):
            mid = (lo + h) / 2
            if spreads(W, onset, mid, c_scale, x0h, K):
                lo = mid
            else:
                h = mid
        return round((lo + h) / 2, 6)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", default=os.path.join(HERE, "results",
                                                     "ideas_sub-1_onset_rh-entorhinal.json"))
    ap.add_argument("--scales", default="1.0,0.99,0.98,0.97,0.96,0.94,0.92")
    a = ap.parse_args()

    R = json.load(open(a.result))
    name = os.path.basename(a.result).replace(".json", "")
    cal = R["calibration"]
    x0h, K, onset = cal["x0_healthy"], cal["K"], R["onset_index"]
    if R["connectome"] == "tvb_template":
        W, labels = connectome.load_tvb_template()
    else:
        W, labels = connectome.load(os.path.join(HERE, "connectomes", R["connectome"]))

    out = {}
    for s in [float(x) for x in a.scales.split(",")]:
        t0 = time.time()
        b = boundary_at(W, onset, s, x0h, K)
        out[f"{s:.3f}"] = b
        print(f"coupling {s:.3f}  boundary x0 = {b}   ({time.time() - t0:.0f}s)", flush=True)

    p = os.path.join(HERE, "results", f"combo_boundary_{name}.json")
    json.dump({"brain": name, "onset_index": onset, "x0_healthy": x0h, "K": K,
               "boundary_by_coupling_scale": out}, open(p, "w"), indent=1)
    print("wrote", p)


if __name__ == "__main__":
    main()
