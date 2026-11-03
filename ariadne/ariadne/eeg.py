"""EEG input: estimate where a seizure starts, as a PRIOR for the onset zone in ARIADNE.

Works with CHB-MIT style data: an EDF recording (bipolar scalp channels such as 'F7-T7') plus its
'.edf.seizures' annotation file (seizure start/end sample indices at 256 Hz).

Method (deliberately simple, so it can be explained and checked):
  1. band-pass 1-40 Hz
  2. for each channel, line length (sum of |sample-to-sample change|) in the first ICTAL_S seconds of the seizure
     divided by the line length in a pre-seizure baseline window
  3. the channels with the biggest increase are the earliest/strongest seizure channels
  4. each channel's electrode-pair midpoint (10-20 positions) is mapped to the nearest connectome region

!! LIMITATION !!  Scalp electrodes -> brain regions is crude: no head model, no registration between this EEG
subject and the MRI. Treat the output as a suggestion for the onset zone, not a localisation.
"""
import re, struct
import numpy as np

ICTAL_S = 10.0
BASE_START_S, BASE_END_S = 70.0, 10.0   # baseline = [onset-70 s, onset-10 s]
FS_DEFAULT = 256


def read_seizure_annotation(path, fs=FS_DEFAULT):
    """Parse a CHB-MIT '.seizures' WFDB annotation file -> list of (start_s, end_s)."""
    b = open(path, "rb").read(); i = 0; t = 0; ev = []
    while i + 2 <= len(b):
        w = struct.unpack_from("<H", b, i)[0]; i += 2
        if w == 0:
            break
        typ, d = w >> 10, w & 0x3FF
        if typ == 59:                                  # SKIP: 32-bit time increment follows
            hi, lo = struct.unpack_from("<HH", b, i); i += 4
            t += (hi << 16) | lo; continue
        if typ in (60, 61, 62, 63):
            if typ == 63:
                i += d + (d & 1)
            continue
        t += d; ev.append((typ, t))
    t0 = 1 << 32                                       # the first SKIP in these files carries a 2^32 offset
    marks = [(ty, (tt - t0) / fs) for ty, tt in ev if tt >= t0 and ty in (32, 33)]
    starts = [s for ty, s in marks if ty == 32]; ends = [s for ty, s in marks if ty == 33]
    return list(zip(starts, ends))


def _line_length(x):
    return np.abs(np.diff(x, axis=-1)).sum(axis=-1)


def channel_scores(raw, onset_s):
    """Ratio of ictal to baseline line length for every channel (higher = stronger seizure activity)."""
    r = raw.copy().filter(1.0, 40.0, verbose="ERROR")
    fs = r.info["sfreq"]
    data = r.get_data()
    a0, a1 = int(max(0, onset_s - BASE_START_S) * fs), int(max(1, onset_s - BASE_END_S) * fs)
    i0, i1 = int(onset_s * fs), int((onset_s + ICTAL_S) * fs)
    base = _line_length(data[:, a0:a1]) / max(a1 - a0 - 1, 1)
    ict = _line_length(data[:, i0:i1]) / max(i1 - i0 - 1, 1)
    scores = ict / np.maximum(base, 1e-12)
    return dict(zip(r.ch_names, scores))


def _electrode_names(ch):
    """'T8-P8-0' -> ['T8','P8'];  'FP1-F7' -> ['FP1','F7']."""
    parts = [p for p in re.split(r"[-]", ch) if p]
    parts = [p for p in parts if not p.isdigit()]
    return parts[:2]


def _positions():
    import mne
    mont = mne.channels.make_standard_montage("standard_1020")
    pos = mont.get_positions()["ch_pos"]
    return {k.upper(): np.asarray(v) for k, v in pos.items()}


def channel_position(ch, P):
    names = [n for n in _electrode_names(ch) if n.upper() in P]
    if not names:
        return None
    return np.mean([P[n.upper()] for n in names], axis=0)


def map_to_regions(positions, centres):
    """Scale head coordinates onto the connectome's extent and return the nearest region index per position."""
    C = np.asarray(centres, float)
    c0 = C.mean(0)
    P = np.array(positions, float)
    P0 = P.mean(0)
    scale = (C.max(0) - C.min(0)) / np.maximum(P.max(0) - P.min(0), 1e-9)
    scale = scale / 1.0
    Ps = (P - P0) * scale.mean() + c0                  # isotropic scaling keeps head shape
    d = np.linalg.norm(Ps[:, None, :] - C[None], axis=-1)
    return d.argmin(1)


def suggest_onset(raw, onset_s, centres, top_channels=3):
    """Return (ranked channel scores, region index votes, suggested region index)."""
    sc = channel_scores(raw, onset_s)
    P = _positions()
    ranked = sorted(sc.items(), key=lambda kv: -kv[1])
    chosen, pos = [], []
    for ch, s in ranked:
        p = channel_position(ch, P)
        if p is None:
            continue
        chosen.append((ch, float(s))); pos.append(p)
        if len(chosen) == top_channels:
            break
    reg = map_to_regions(pos, centres) if pos else []
    votes = {}
    for (ch, s), r in zip(chosen, reg):
        votes[int(r)] = votes.get(int(r), 0.0) + s
    best = max(votes, key=votes.get) if votes else None
    return {"top_channels": chosen, "region_votes": votes, "suggested_region": best}
