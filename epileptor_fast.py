"""Fast Epileptor network (numba). Same equations as epileptor.py; different random numbers.
Returns x1 series (steps, trials, regions)."""
import numpy as np
from numba import njit

@njit(cache=True)
def _run(C, x0, K, T, dt, n_trials, noise, seed, every):
    np.random.seed(seed)
    N = C.shape[0]
    steps = int(T / dt)
    out = np.empty((steps // every, n_trials, N), dtype=np.float32)
    rs = np.zeros(N)
    for i in range(N):
        for j in range(N):
            rs[i] += C[i, j]
    I1, I2, c, d, r, tau2, gam, a = 3.1, 0.45, 1.0, 5.0, 0.00035, 10.0, 0.01, 6.0
    sq = np.sqrt(dt)
    for tr in range(n_trials):
        x1 = np.full(N, -1.5); y1 = np.full(N, -5.0); z = np.full(N, 3.0)
        x2 = np.zeros(N); y2 = np.zeros(N); u = np.zeros(N)
        k = 0
        for s in range(steps):
            nx1 = np.empty(N); nx2 = np.empty(N); ny1 = np.empty(N)
            ny2 = np.empty(N); nz = np.empty(N); nu = np.empty(N)
            for i in range(N):
                coup = 0.0
                for j in range(N):
                    coup += C[i, j] * (x1[j] - x1[i])
                if x1[i] < 0:
                    f1 = x1[i]**3 - 3*x1[i]**2
                else:
                    f1 = (x2[i] - 0.6*(z[i]-4)**2) * x1[i]
                f2 = 0.0 if x2[i] < -0.25 else a*(x2[i] + 0.25)
                dx1 = y1[i] - f1 - z[i] + I1
                dy1 = c - d*x1[i]**2 - y1[i]
                dz = r*(4*(x1[i]-x0[i]) - z[i] - K*coup)
                dx2 = -y2[i] + x2[i] - x2[i]**3 + I2 + 2*u[i] - 0.3*(z[i]-3.5)
                dy2 = (-y2[i] + f2)/tau2
                du = -gam*(u[i] - 0.1*x1[i])
                nx1[i] = x1[i] + dt*dx1 + noise*sq*np.random.standard_normal()
                ny1[i] = y1[i] + dt*dy1
                nz[i] = z[i] + dt*dz
                nx2[i] = x2[i] + dt*dx2 + noise*sq*np.random.standard_normal()
                ny2[i] = y2[i] + dt*dy2
                nu[i] = u[i] + dt*du
            x1, y1, z, x2, y2, u = nx1, ny1, nz, nx2, ny2, nu
            if s % every == 0 and k < out.shape[0]:
                for i in range(N):
                    out[k, tr, i] = x1[i]
                k += 1
    return out

def simulate(C, x0, K=0.0, T=6000, dt=0.05, n_trials=1, noise=0.0005, seed=0, every=20):
    C = np.ascontiguousarray(C, dtype=np.float64)
    x0 = np.ascontiguousarray(np.broadcast_to(np.asarray(x0, float), (C.shape[0],)))
    return _run(C, x0, float(K), float(T), float(dt), int(n_trials), float(noise), int(seed), int(every))
