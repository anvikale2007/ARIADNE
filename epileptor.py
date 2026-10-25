"""Epileptor network model (Jirsa et al. 2014; VEP, Jirsa et al. 2017).
Verify constants against the papers before citing."""
import numpy as np

P = dict(I1=3.1, I2=0.45, c=1.0, d=5.0, r=0.00035, tau2=10.0, gamma=0.01, a=6.0)

def simulate(C, x0, K=0.0, T=6000, dt=0.05, n_trials=1, noise=0.0005, seed=0, p=P, every=20):
    """C: (N,N) symmetric connectome. x0: scalar or (N,). Returns x1 series (steps, trials, N)."""
    rng = np.random.default_rng(seed)
    N = C.shape[0]
    x0 = np.broadcast_to(np.asarray(x0, float), (N,))
    shape = (n_trials, N)
    x1 = np.full(shape, -1.5); y1 = np.full(shape, -5.0); z = np.full(shape, 3.0)
    x2 = np.zeros(shape); y2 = np.zeros(shape); u = np.zeros(shape)
    rs = C.sum(1)
    steps = int(T / dt)
    out = np.empty((steps // every, n_trials, N), dtype=np.float32)
    k = 0
    for s in range(steps):
        coup = x1 @ C.T - x1 * rs
        f1 = np.where(x1 < 0, x1**3 - 3*x1**2, (x2 - 0.6*(z - 4)**2) * x1)
        f2 = np.where(x2 < -0.25, 0.0, p['a'] * (x2 + 0.25))
        dx1 = y1 - f1 - z + p['I1']
        dy1 = p['c'] - p['d'] * x1**2 - y1
        dz = p['r'] * (4*(x1 - x0) - z - K*coup)
        dx2 = -y2 + x2 - x2**3 + p['I2'] + 2*u - 0.3*(z - 3.5)
        dy2 = (-y2 + f2) / p['tau2']
        du = -p['gamma'] * (u - 0.1*x1)
        sq = np.sqrt(dt)
        x1 = x1 + dt*dx1 + noise*sq*rng.standard_normal(shape)
        y1 = y1 + dt*dy1
        z = z + dt*dz
        x2 = x2 + dt*dx2 + noise*sq*rng.standard_normal(shape)
        y2 = y2 + dt*dy2
        u = u + dt*du
        if s % every == 0 and k < out.shape[0]:
            out[k] = x1; k += 1
    return out
