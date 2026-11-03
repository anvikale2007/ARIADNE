"""Virtual surgery: does removing the resected tissue make this brain network less able to generate seizures?

Non-circular design (no onset is chosen from the resection):
  1. Network ictogenicity threshold x_c: give EVERY region the same excitability x0 and find (bisection) the lowest x0 at which
     the network produces a seizure on its own (noise-driven, >= MIN_RECRUITED regions recruited). Well-connected,
     tightly coupled networks ignite at lower x0.
  2. Virtual resection: scale each region's connections by (1 - fraction resected) and silence fully removed tissue.
  3. Delta x_c = x_c(after) - x_c(before). Larger Delta = surgery removed more of the network's seizure-generating capacity.
  4. Across patients, compare Delta between seizure-free (ILAE 1-2) and not seizure-free (ILAE 3-5) groups.
This is the logic of 'brain network ictogenicity' studies (e.g. the Newcastle CNNP group that released IDEAS); the exact
implementation here is ARIADNE's own and has not been validated: treat results as exploratory.
"""
import numpy as np
from . import model, metrics

MIN_RECRUITED = 1
LENGTH_MS = 60000
REMOVED_X0 = -3.5          # fully resected tissue cannot seize


def _simulate(W, x0, K, seed):
    rec_every = int(model.RECORD_EVERY_MS / model.DT)
    return model._simulate(np.ascontiguousarray(W, dtype=np.float64), np.asarray(x0, dtype=np.float64), K, model.R,
                           model.IEXT, float(LENGTH_MS), model.DT, rec_every, model.NOISE, int(seed))


def resect(W, frac, x0_base, graded=False):
    """Virtual resection. Connections are scaled by the fraction of tissue left in each region.

    graded=False (original): only regions at least 95% removed are silenced, so a region 60% resected
    keeps full seizure-generating capacity - which under-counts partial resections.
    graded=True: excitability moves toward the silenced value in proportion to the fraction removed."""
    f = np.clip(np.nan_to_num(np.asarray(frac, float)), 0, 1)
    keep = 1.0 - f
    W2 = W * keep[:, None] * keep[None, :]
    x0 = np.full(W.shape[0], x0_base, dtype=float)
    if graded:
        x0 = x0_base + f * (REMOVED_X0 - x0_base)
    else:
        x0[f >= 0.95] = REMOVED_X0
    return W2, x0


def ignites(W, x0, K, seeds=(0,)):
    votes = [metrics.analyse(_simulate(W, x0, K, s))["recruited_regions"] >= MIN_RECRUITED for s in seeds]
    return np.mean(votes) > 0.5


def threshold(W, K, frac=None, lo=-2.30, hi=-2.00, iters=6, seeds=(0,)):
    """Lowest homogeneous x0 at which the (optionally resected) network ignites."""
    def ign(x):
        if frac is None:
            return ignites(W, np.full(W.shape[0], x), K, seeds)
        W2, x0 = resect(W, frac, x)
        return ignites(W2, x0, K, seeds)
    if ign(lo):
        return lo
    if not ign(hi):
        return None
    for _ in range(iters):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if ign(mid) else (mid, hi)
    return round((lo + hi) / 2, 4)


def delta_threshold(W, frac, K, **kw):
    pre = threshold(W, K, None, **kw)
    post = threshold(W, K, frac, **kw)
    return pre, post, (None if pre is None or post is None else round(post - pre, 4))


# ---------------------------------------------------------------------------------------------------------------
# Design 2 (replaces the homogeneous-threshold test, which proved insensitive: with every region equally excitable
# the network sits at a common rest state, coupling carries no signal, and ignition happens at the single-node
# threshold regardless of surgery; delta = 0 for sub-1).
# Spread capacity: seed a seizure in EACH candidate region in turn (others at the calibrated background excitability)
# and measure the mean fraction of the network recruited. Do it before and after virtual resection; onsets that were
# removed by surgery count as 'no seizure' afterwards. Reduction in capacity = how much the surgery cut the network's
# ability to generate AND spread seizures. No onset is chosen from the resection, so the test is not circular.
# ---------------------------------------------------------------------------------------------------------------
from concurrent.futures import ProcessPoolExecutor

SWEEP_LENGTH_MS = 40000


def _sweep_job(args):
    W, x0, K, seed = args
    rec_every = int(model.RECORD_EVERY_MS / model.DT)
    x1 = model._simulate(np.ascontiguousarray(W, dtype=np.float64), np.asarray(x0, dtype=np.float64), K, model.R,
                         model.IEXT, float(SWEEP_LENGTH_MS), model.DT, rec_every, model.NOISE, int(seed))
    return metrics.analyse(x1)["recruited_regions"]


def spread_capacity(W, x0h, K, onsets, frac=None, x0_onset=-1.6, workers=2, seed=0, graded=False):
    n = W.shape[0]
    jobs, removed = [], []
    for o in onsets:
        if frac is None:
            W2, x0 = W, np.full(n, x0h)
        else:
            W2, x0 = resect(W, frac, x0h, graded=graded)
        if frac is not None and np.nan_to_num(frac[o]) >= 0.95:
            removed.append(o); continue
        x0 = x0.copy(); x0[o] = x0_onset
        jobs.append((W2, x0, K, seed))
    with ProcessPoolExecutor(workers) as ex:
        rec = list(ex.map(_sweep_job, jobs))
    per_onset = np.array(rec + [0] * len(removed), float) / n
    return float(per_onset.mean()), per_onset
