"""Check our numba Epileptor-2D equations match The Virtual Brain's reference implementation."""
import sys, os, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from tvb.simulator.lab import models
from ariadne import model

rng = np.random.default_rng(1)
n = 12
x1 = rng.uniform(-2.5, 1.0, n); z = rng.uniform(-0.5, 5.0, n)
x0 = rng.uniform(-2.5, -1.5, n)
W = rng.random((n, n)); W = (W + W.T) / 2; np.fill_diagonal(W, 0)
K = 0.7
rsum = W.sum(1)
dx1 = np.empty(n); dz = np.empty(n)
model._dfun(x1, z, x0, model.IEXT, model.R, K, W, rsum, dx1, dz)

m = models.Epileptor2D(x0=x0, Ks=np.array([K]), r=np.array([model.R]), Iext=np.array([model.IEXT]))
state = np.stack([x1, z])[:, :, None]
c = ((W @ x1) - rsum * x1)[None, :, None]
ref = m.dfun(state, c, 0.0)[:, :, 0]
err = max(np.abs(ref[0] - dx1).max(), np.abs(ref[1] - dz).max())
print("max abs difference vs TVB:", err)
assert err < 1e-10, "Equations do not match TVB!"
print("PASS")
