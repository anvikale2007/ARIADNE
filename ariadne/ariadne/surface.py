"""Lookup against a precomputed stop/no-stop response surface (see response_surface.py).

Every ARIADNE drug acts through two model quantities only: an excitability reduction and a coupling
reduction. Once the boundary between "seizure spreads" and "seizure contained" has been mapped in that
2D plane, a drug's threshold is a geometry problem rather than a simulation, which is what makes
uncertainty propagation over thousands of parameter draws affordable.
"""
import json
import numpy as np


class Surface:
    def __init__(self, d):
        self.d = d
        self.x0_grid = np.array(d["x0_grid"])
        self.c_grid = np.array(d["coupling_grid"])
        self.X0_UNIT = d["units"]["X0_UNIT"]
        self.COUPLING_UNIT = d["units"]["COUPLING_UNIT"]
        self.S_MAX = d["units"]["S_MAX"]
        b = [np.inf if v is None else float(v) for v in d["boundary_x0_reduction"]]
        order = np.argsort(self.c_grid)                      # ascending coupling scale
        self._c = self.c_grid[order]
        self._b = np.array(b)[order]
        self.c_min, self.c_max = float(self._c.min()), float(self._c.max())
        self.clamped = 0                                     # how many lookups fell outside the mapped range

    # ---- core geometry -------------------------------------------------
    def required_x0(self, coupling_scale):
        """Minimum excitability reduction needed to contain the seizure at this coupling scale."""
        c = float(coupling_scale)
        if c < self.c_min:
            self.clamped += 1
            c = self.c_min
        if c > self.c_max:
            c = self.c_max
        finite = np.isfinite(self._b)
        if not finite.any():
            return np.inf
        return float(np.interp(c, self._c[finite], self._b[finite]))

    def stops(self, x0_reduction, coupling_scale):
        return float(x0_reduction) >= self.required_x0(coupling_scale)

    def required_x0_vec(self, coupling_scale):
        """Vectorised required_x0. Values below the mapped range are clamped to its lower edge."""
        c = np.clip(np.asarray(coupling_scale, dtype=float), self.c_min, self.c_max)
        finite = np.isfinite(self._b)
        if not finite.any():
            return np.full(c.shape, np.inf)
        return np.interp(c, self._c[finite], self._b[finite])

    def threshold(self, x0_effect, coupling_effect, s_max=None, step=0.005):
        """Smallest drug strength that contains the seizure, or None if it never does up to s_max.

        Strength scales both mechanisms linearly, exactly as engine.drug_params does."""
        s_max = self.S_MAX if s_max is None else s_max
        if x0_effect <= 0 and coupling_effect <= 0:
            return None
        s = np.arange(step, s_max + step, step)
        x0r = s * float(x0_effect) * self.X0_UNIT
        cs = np.maximum(0.0, 1.0 - s * float(coupling_effect) * self.COUPLING_UNIT)
        hit = x0r >= self.required_x0_vec(cs)
        if not hit.any():
            return None
        return round(float(s[int(np.argmax(hit))]), 4)

    def contained_at(self, x0_effect, coupling_effect):
        """Does a drug with these effects contain the seizure at exposure as given (strength 1.0)?"""
        return self.stops(x0_effect * self.X0_UNIT, max(0.0, 1.0 - coupling_effect * self.COUPLING_UNIT))


def load(path):
    return Surface(json.load(open(path)))


def rank(thresholds):
    """Lower threshold is better; drugs that never contain the seizure rank last (ties keep input order)."""
    return sorted(thresholds, key=lambda k: (thresholds[k] is None, thresholds[k] if thresholds[k] is not None else 0))
