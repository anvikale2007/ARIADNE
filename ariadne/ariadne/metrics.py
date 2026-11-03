"""Quantifying seizure spread from simulated activity x1 (regions x time), after discarding a burn-in."""
import numpy as np

BURN_IN_MS = 15000
RECORD_MS = 10.0
RECRUIT_AMPLITUDE = 1.0   # a region is "recruited" if its x1 swings by more than this (ictal swings are ~2.5)
ICTAL_LEVEL = -0.5        # x1 above this = actively seizing


def analyse(x1):
    y = x1[int(BURN_IN_MS / RECORD_MS):]
    amp = y.max(axis=0) - y.min(axis=0)
    recruited = amp > RECRUIT_AMPLITUDE
    burden = (y > ICTAL_LEVEL).mean(axis=0)          # fraction of time each region is ictal
    return {"recruited_regions": int(recruited.sum()),
            "mean_burden": float(burden.mean()),
            "recruited_mask": recruited.tolist(),
            "burden_by_region": burden.tolist()}
