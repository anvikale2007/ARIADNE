"""Network Epileptor-2D simulation (Proix et al. 2014; the reduced Epileptor used in the Virtual Epileptic Patient).

Each connectome region is one neural-mass node with a fast variable x1 (electrical activity) and a slow
variable z (permittivity). Regions are coupled through the structural connectome: a seizing region pulls the
slow variable of its neighbours towards seizure, so activity spreads along white-matter connections.

The equations are re-implemented in numba for speed and verified against The Virtual Brain's reference
implementation in tests/test_against_tvb.py.
"""
import numpy as np
from numba import njit

# --- parameters --------------------------------------------------------------------------------------
X0_HEALTHY = -2.3        # excitability of healthy regions (single-node seizure threshold verified numerically between -2.07 and -2.06)
X0_EPILEPTOGENIC = -1.6  # excitability of the seizure onset zone
K_COUPLING = 6.0         # global coupling strength (re-calibrated per connectome by run_pipeline.py)
R = 0.00035              # slow time scale
IEXT = 3.1
DT = 0.1                 # ms
NOISE = 1e-4
RECORD_EVERY_MS = 10.0


@njit(cache=True)
def _dfun(x1, z, x0, Iext, r, K, W, rsum, dx1, dz):
    n = x1.shape[0]
    for i in range(n):
        c = -rsum[i] * x1[i]
        for j in range(n):
            c += W[i, j] * x1[j]
        # fast variable
        if x1[i] < 0.0:
            f = 1.0 * x1[i] ** 2 + (5.0 - 3.0) * x1[i]
        else:
            f = -0.0 - 0.6 * (z[i] - 4.0) ** 2 + 5.0 * x1[i]
        dx1[i] = 1.0 - z[i] + Iext - f * x1[i]
        # slow variable ("energy" term only active for z < 0)
        e = -0.1 * z[i] ** 7 if z[i] < 0.0 else 0.0
        h = 4.0 * (x1[i] - x0[i]) + e
        dz[i] = r * (h - z[i] + K * c)


@njit(cache=True)
def _simulate(W, x0, K, r, Iext, T, dt, rec_every, noise, seed):
    n = W.shape[0]
    rsum = np.zeros(n)
    for i in range(n):
        rsum[i] = W[i].sum()
    np.random.seed(seed)
    x1 = np.full(n, -1.8)
    z = np.full(n, 3.0)
    nsteps = int(T / dt)
    nrec = nsteps // rec_every
    out = np.empty((nrec, n))
    a1 = np.empty(n); a2 = np.empty(n); b1 = np.empty(n); b2 = np.empty(n)
    sq = np.sqrt(dt)
    for k in range(nsteps):
        _dfun(x1, z, x0, Iext, r, K, W, rsum, a1, a2)
        x1p = x1 + dt * a1
        zp = z + dt * a2
        _dfun(x1p, zp, x0, Iext, r, K, W, rsum, b1, b2)
        x1 = x1 + 0.5 * dt * (a1 + b1) + noise * sq * np.random.randn(n)
        z = z + 0.5 * dt * (a2 + b2)
        if k % rec_every == 0 and k // rec_every < nrec:
            out[k // rec_every] = x1
    return out


def run_simulation(W, onset_idx=None, x0_shift=0.0, coupling_scale=1.0, length_ms=80000, seed=0,
                   K=None, x0_healthy=None):
    """Simulate one seizure scenario. Returns x1[t, region] sampled every 10 ms.

    W              connectivity matrix (normalised, max = 1)
    onset_idx      index of the seizure onset zone (None = healthy control)
    x0_shift       drug effect on excitability, applied to ALL regions (negative = less excitable)
    coupling_scale drug effect on inter-regional coupling (<1 = weaker propagation)
    """
    n = W.shape[0]
    K = K_COUPLING if K is None else K
    x0_healthy = X0_HEALTHY if x0_healthy is None else x0_healthy
    x0 = np.full(n, x0_healthy)
    if onset_idx is not None:
        x0[onset_idx] = X0_EPILEPTOGENIC
    x0 = x0 + x0_shift
    rec_every = int(RECORD_EVERY_MS / DT)
    return _simulate(np.ascontiguousarray(W, dtype=np.float64), x0, K * coupling_scale, R, IEXT,
                     float(length_ms), DT, rec_every, NOISE, int(seed))
