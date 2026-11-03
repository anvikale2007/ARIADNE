"""Synthetic check of the EEG onset module (no real EEG needed). Also verifies the .seizures parser on real CHB-MIT annotations."""
import sys, os, warnings, glob
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, mne
from ariadne import eeg

# 1. annotation parser against the 4 files supplied by the user (known CHB-MIT chb01 seizure times)
expected = {"chb01_03": (2996, 3036), "chb01_15": (1732, 1772), "chb01_16": (1015, 1066), "chb01_26": (1862, 1963)}
d = "/mnt/user-data/attachments/4dd36763-85b0-510b-8425-34461c4e8d03"
for rec, (s, e) in expected.items():
    fs_ = glob.glob(f"{d}/*{rec}.edf.seizures")
    if fs_:
        got = eeg.read_seizure_annotation(fs_[0])
        print(rec, got); assert got and abs(got[0][0] - s) < 1 and abs(got[0][1] - e) < 1, "annotation parse mismatch"

# 2. synthetic recording: 400 s of noise, seizure burst at 200 s concentrated in T7-P7
rng = np.random.default_rng(0)
chs = ["FP1-F7", "F7-T7", "T7-P7", "P7-O1", "FP1-F3", "F3-C3", "C3-P3", "P3-O1", "FZ-CZ", "CZ-PZ", "T8-P8-0", "F8-T8"]
fs = 256; n = 400 * fs
X = rng.standard_normal((len(chs), n)) * 10e-6
t = np.arange(int(40 * fs)) / fs
burst = 120e-6 * np.sin(2 * np.pi * 8 * t)
X[chs.index("T7-P7"), 200 * fs:200 * fs + len(t)] += burst
X[chs.index("F7-T7"), 200 * fs:200 * fs + len(t)] += 0.4 * burst
raw = mne.io.RawArray(X, mne.create_info(chs, fs, "eeg"), verbose="ERROR")
sc = eeg.channel_scores(raw, 200.0)
top = max(sc, key=sc.get)
print("top channel:", top); assert top == "T7-P7"
centres = np.random.default_rng(1).uniform(-60, 60, (70, 3))
res = eeg.suggest_onset(raw, 200.0, centres)
print(res); assert res["suggested_region"] is not None
print("PASS")
