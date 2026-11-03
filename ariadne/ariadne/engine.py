"""ARIADNE engine: calibrate a connectome, then find how much of each drug is needed to stop seizure spread."""
import os, json
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import yaml
from . import model, metrics

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
K_COUPLING = -5.0        # coupling strength (negative = seizing regions pull neighbours' slow variable towards seizure)
X0_UNIT = 0.004          # model units of excitability reduction at x0_effect = 1, strength 1
COUPLING_UNIT = 0.06     # fractional coupling reduction at coupling_effect = 1, strength 1
S_MAX = 4.0              # strongest drug strength tested (in units of "nominal dose")
SPREAD_FRACTION = 0.5    # seizure "spreads" if more than this fraction of regions are recruited
LENGTH_MS = 60000


def load_drugs(path=None):
    return yaml.safe_load(open(path or os.path.join(HERE, "drugs.yaml")))


def drug_params(d, s, x0_weight=1.0, coupling_weight=1.0):
    """Model parameters for drug d at strength s."""
    x0_shift = -s * d["x0_effect"] * X0_UNIT * x0_weight
    coupling_scale = max(0.0, 1.0 - s * d["coupling_effect"] * COUPLING_UNIT * coupling_weight)
    return x0_shift, coupling_scale


def _sim(W, onset, x0_shift, cscale, seed, x0_healthy, K):
    return metrics.analyse(model.run_simulation(W, onset, x0_shift, cscale, LENGTH_MS, seed, K=K, x0_healthy=x0_healthy))


def spreads(res, n):
    return res["recruited_regions"] > SPREAD_FRACTION * n


# ---------- calibration ----------
def _cal_job(a):
    W, onset, x0h, K = a
    return x0h, _sim(W, onset, 0.0, 1.0, 0, x0h, K)


def calibrate(W, onset, K=K_COUPLING, grid=None, workers=2):
    """Find the brain's own seizure-spread threshold: the lowest background excitability at which a seizure
    starting at `onset` recruits most of the network while a healthy (no-onset) network stays quiet."""
    n = W.shape[0]
    grid = grid if grid is not None else np.round(np.arange(-2.20, -2.059, 0.01), 3)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        out = list(ex.map(_cal_job, [(W, onset, float(x), K) for x in grid]))
    chosen = None
    for x0h, res in out:
        if spreads(res, n):
            chosen = x0h; break
    table = [(float(x), r["recruited_regions"]) for x, r in out]
    if chosen is None:
        return {"ok": False, "table": table, "K": K}
    healthy = _sim(W, None, 0.0, 1.0, 0, chosen, K)
    return {"ok": True, "x0_healthy": chosen, "K": K, "healthy_recruited": healthy["recruited_regions"], "table": table}


# ---------- drug strength needed to stop spread ----------
def _drug_job(a):
    W, onset, name, d, cal, seeds, iters, xw, cw = a
    n = W.shape[0]
    x0h, K = cal["x0_healthy"], cal["K"]

    def spread_at(s):
        votes = []
        for seed in seeds:
            sh, cs = drug_params(d, s, xw, cw)
            votes.append(spreads(_sim(W, onset, sh, cs, seed, x0h, K), n))
        return np.mean(votes) > 0.5

    history = []
    if spread_at(S_MAX):
        s_star = None; history.append((S_MAX, True))
    else:
        lo, hi = 0.0, S_MAX
        for _ in range(iters):
            mid = (lo + hi) / 2
            sp = spread_at(mid); history.append((round(mid, 3), bool(sp)))
            lo, hi = (mid, hi) if sp else (lo, mid)
        s_star = round((lo + hi) / 2, 3)
    return name, s_star, history


def drug_thresholds(W, onset, cal, drugs, seeds=(0,), iters=5, x0_weight=1.0, coupling_weight=1.0, workers=2):
    jobs = [(W, onset, nm, d, cal, list(seeds), iters, x0_weight, coupling_weight) for nm, d in drugs.items()]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return {nm: {"threshold": s, "history": h} for nm, s, h in ex.map(_drug_job, jobs)}


def rank(th):
    """Lower threshold = less drug needed to stop spread = better. Drugs that never stop it rank last."""
    return sorted(th, key=lambda k: (th[k]["threshold"] is None, th[k]["threshold"] or 1e9))


def snapshot(W, onset, cal, drugs, th, seed=0):
    """Per-region state: untreated vs each drug at nominal dose (strength 1), for visualisation / LEDs."""
    n = W.shape[0]; x0h, K = cal["x0_healthy"], cal["K"]
    base = _sim(W, onset, 0.0, 1.0, seed, x0h, K)
    out = {"baseline": {"recruited_mask": base["recruited_mask"], "recruited_regions": base["recruited_regions"]}}
    for nm, d in drugs.items():
        sh, cs = drug_params(d, 1.0)
        r = _sim(W, onset, sh, cs, seed, x0h, K)
        out[nm] = {"recruited_mask": r["recruited_mask"], "recruited_regions": r["recruited_regions"]}
    return out
