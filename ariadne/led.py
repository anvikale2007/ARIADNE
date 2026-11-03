"""Drive the mini-brain LED demo from an ARIADNE result.

Regions are grouped into N_LEDS groups by position (front-to-back / left-to-right / up-down via k-means on region
centres). Each LED shows its group's state:   0 = green (quiet)   1 = amber (partly recruited)   2 = red (seizing)

  python led.py --result results/stanford_hardi_onsetP39.json --state baseline --port /dev/ttyACM0
  python led.py --result results/stanford_hardi_onsetP39.json --state cenobamate --dry-run   # print only
Windows ports look like COM3. See arduino/ariadne_leds/ariadne_leds.ino for the matching sketch.
"""
import argparse, json, os, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
N_LEDS = 12


def led_states(result, state, n_leds=N_LEDS):
    mask = np.array(result["snapshot"][state]["recruited_mask"], dtype=bool)
    cpath = os.path.join(HERE, "connectomes", result["connectome"], "centres_mm.csv")
    if os.path.exists(cpath):
        from sklearn.cluster import KMeans
        c = np.loadtxt(cpath, delimiter=",")
        g = KMeans(n_clusters=n_leds, n_init=4, random_state=0).fit_predict(c)
    else:
        g = np.arange(len(mask)) % n_leds
    states = []
    for k in range(n_leds):
        frac = mask[g == k].mean() if (g == k).any() else 0.0
        states.append(0 if frac == 0 else (1 if frac < 0.5 else 2))
    return states


def send(states, port=None, baud=115200):
    frame = "S:" + "".join(str(s) for s in states) + "\n"
    if port is None:
        print("[dry run]", frame.strip()); return
    import serial
    with serial.Serial(port, baud, timeout=1) as ser:
        time.sleep(2)  # Arduino resets when the port opens
        ser.write(frame.encode()); print("sent", frame.strip())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", required=True); ap.add_argument("--state", default="baseline",
                    help="'baseline' or a drug name")
    ap.add_argument("--port", default=None); ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    r = json.load(open(a.result))
    send(led_states(r, a.state), None if a.dry_run else a.port)
