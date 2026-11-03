"""Map the stop/no-stop boundary in (excitability reduction, coupling reduction) space for one brain.

Every drug in ARIADNE acts through exactly two model quantities:
    x0_shift       = -strength * x0_effect * X0_UNIT          (how much excitability is lowered)
    coupling_scale = 1 - strength * coupling_effect * COUPLING_UNIT

So the whole drug model lives on a 2D plane. Mapping the boundary once turns every later
threshold calculation into a lookup, which is what makes Monte Carlo uncertainty propagation
affordable (a single simulation takes ~5 s).

    python response_surface.py --result results/stanford_hardi_onsetP39.json

Writes results/surface_<name>.json: the grid, the recruited counts, and the monotone boundary
(minimum excitability reduction needed to stop spread, as a function of coupling scale).
"""
import argparse, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ariadne import connectome, engine, model, metrics

HERE = os.path.dirname(os.path.abspath(__file__))
X0_MAX = 0.030      # largest excitability reduction mapped (cenobamate at strength 4 is ~0.029)
C_MIN = 0.80        # smallest coupling scale mapped (pregabalin at strength 4 is ~0.86)


def _job(a):
    W, onset, x0s, cs, seed, x0h, K = a
    r = metrics.analyse(model.run_simulation(W, onset, x0s, cs, engine.LENGTH_MS, seed, K=K, x0_healthy=x0h))
    return x0s, cs, seed, r["recruited_regions"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", default=os.path.join(HERE, "results", "stanford_hardi_onsetP39.json"))
    ap.add_argument("--nx", type=int, default=81, help="excitability-reduction grid points")
    ap.add_argument("--nc", type=int, default=21, help="coupling-scale grid points")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--confirm-seeds", type=int, default=3, help="extra noise seeds used to re-test the boundary cells")
    a = ap.parse_args()

    R = json.load(open(a.result))
    name = os.path.basename(a.result).replace(".json", "")
    cal = R["calibration"]
    x0h, K, onset = cal["x0_healthy"], cal["K"], R["onset_index"]
    if R["connectome"] == "tvb_template":
        W, labels = connectome.load_tvb_template()
    else:
        W, labels = connectome.load(os.path.join(HERE, "connectomes", R["connectome"]))
    n = len(labels)

    x0_grid = np.round(np.linspace(0.0, X0_MAX, a.nx), 5)       # excitability REDUCTION (positive)
    c_grid = np.round(np.linspace(1.0, C_MIN, a.nc), 4)         # coupling scale (1 = untouched)

    # Checkpoint so a container restart cannot cost the whole run.
    part_path = os.path.join(HERE, "results", f"surface_{name}.partial.json")
    rec = {}
    if os.path.exists(part_path):
        for k_, v in json.load(open(part_path)).items():
            xs, cs_ = k_.split("|")
            rec[(float(xs), float(cs_))] = v
        print(f"resuming: {len(rec)} cells already done", flush=True)

    jobs = [(W, onset, -float(x), float(c), 0, x0h, K)
            for c in c_grid for x in x0_grid if (float(x), float(c)) not in rec]
    total = len(x0_grid) * len(c_grid)
    print(f"{name}: {len(jobs)} of {total} simulations on {a.workers} workers "
          f"(~{len(jobs)*2.5/a.workers/60:.0f} min)", flush=True)

    def save_partial():
        os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
        json.dump({f"{x}|{c}": v for (x, c), v in rec.items()}, open(part_path + ".tmp", "w"))
        os.replace(part_path + ".tmp", part_path)

    t0 = time.time()
    if jobs:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for i, (x0s, cs, seed, k) in enumerate(ex.map(_job, jobs), 1):
                rec[(round(-x0s, 5), round(cs, 4))] = k
                if i % 25 == 0:
                    save_partial()
                    print(f"  {i}/{len(jobs)}  ({time.time()-t0:.0f}s)", flush=True)
        save_partial()

    grid = np.array([[rec[(float(x), float(c))] for x in x0_grid] for c in c_grid])   # [coupling, x0]
    stops = grid <= engine.SPREAD_FRACTION * n

    # Boundary: smallest excitability reduction that stops spread, per coupling scale.
    # Enforced monotone in x0 (more reduction can only help) to absorb noise jitter.
    mono = np.maximum.accumulate(stops, axis=1)
    boundary = []
    for ci, c in enumerate(c_grid):
        idx = np.argmax(mono[ci]) if mono[ci].any() else None
        boundary.append(None if idx is None else float(x0_grid[idx]))

    # Re-test the boundary cells with extra noise seeds (the grid itself used one seed).
    conf = []
    if a.confirm_seeds > 1:
        cells = [(ci, b) for ci, b in enumerate(boundary) if b is not None and b > 0]
        jobs2 = [(W, onset, -float(b), float(c_grid[ci]), s, x0h, K)
                 for ci, b in cells for s in range(1, a.confirm_seeds)]
        print(f"confirming {len(cells)} boundary cells with {a.confirm_seeds-1} extra seeds ...", flush=True)
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            got = {}
            for x0s, cs, seed, k in ex.map(_job, jobs2):
                got.setdefault((round(-x0s, 5), round(cs, 4)), []).append(k)
        for ci, b in cells:
            ks = got.get((round(b, 5), round(float(c_grid[ci]), 4)), [])
            agree = all(k <= engine.SPREAD_FRACTION * n for k in ks)
            conf.append({"coupling_scale": float(c_grid[ci]), "x0_reduction": b, "recruited_other_seeds": ks, "agrees": bool(agree)})

    out = {"name": name, "connectome": R["connectome"], "onset": R["onset"], "n_regions": n,
           "x0_healthy": x0h, "K": K, "spread_fraction": engine.SPREAD_FRACTION,
           "x0_grid": [float(x) for x in x0_grid], "coupling_grid": [float(c) for c in c_grid],
           "recruited": grid.tolist(), "stops": stops.tolist(), "boundary_x0_reduction": boundary,
           "boundary_confirmation": conf,
           "units": {"X0_UNIT": engine.X0_UNIT, "COUPLING_UNIT": engine.COUPLING_UNIT, "S_MAX": engine.S_MAX},
           "note": "x0_reduction is a POSITIVE excitability reduction in model units; coupling_scale 1.0 = untouched. "
                   "Grid used noise seed 0; boundary cells re-tested with further seeds (see boundary_confirmation)."}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    p = os.path.join(HERE, "results", f"surface_{name}.json")
    json.dump(out, open(p, "w"))
    if os.path.exists(part_path):
        os.remove(part_path)
    ok = sum(c["agrees"] for c in conf)
    print(f"\nSaved {p}  ({time.time()-t0:.0f}s)")
    print(f"Boundary confirmed on {ok}/{len(conf)} cells across extra seeds")
    for c, b in zip(c_grid, boundary):
        print(f"  coupling {c:.3f} -> needs x0 reduction {'none in range' if b is None else f'{b:.4f}'}")


if __name__ == "__main__":
    main()
