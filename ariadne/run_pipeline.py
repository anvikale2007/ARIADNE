"""ARIADNE command-line pipeline.

  python run_pipeline.py                                   # Stanford scan connectome, default onset
  python run_pipeline.py --connectome connectomes/sherbrooke_3shell --onset 12
  python run_pipeline.py --sensitivity --seeds 0 1         # also test how robust the ranking is
  python run_pipeline.py --template                        # TVB 76-region template instead of your scans

Output: results/<name>.json  (read by app.py) and results/<name>.png
"""
import argparse, json, os, sys, time, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
import numpy as np
from ariadne import connectome, engine, genetics, eeg, ieeg, clinical

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--connectome", default=os.path.join(HERE, "connectomes", "stanford_hardi"))
    ap.add_argument("--template", action="store_true", help="use the TVB template connectome")
    ap.add_argument("--onset", default="auto", help="region index/label of the seizure onset zone ('auto' = median-strength region)")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0], help="noise seeds (majority vote)")
    ap.add_argument("--iters", type=int, default=6, help="bisection steps (precision = %.2f / 2^iters)" % engine.S_MAX)
    ap.add_argument("--genetics", default=None, help="patient YAML with a 'variants:' list (applies gene-drug safety rules)")
    ap.add_argument("--eeg-edf", default=None, help="EDF recording; with --seizures, suggests the onset region from EEG")
    ap.add_argument("--seizures", default=None, help="CHB-MIT .edf.seizures annotation file for --eeg-edf")
    ap.add_argument("--seizure-index", type=int, default=0, help="which annotated seizure to use")
    ap.add_argument("--ieeg-edf", default=None, help="BIDS iEEG EDF (ictal epoch); needs --template and the three files below")
    ap.add_argument("--ieeg-channels", default=None); ap.add_argument("--ieeg-events", default=None); ap.add_argument("--ieeg-electrodes", default=None)
    ap.add_argument("--seizure-type", default=None, choices=list(clinical.seizure_types()),
                    help="ILAE 2025 seizure type: applies guideline rules; only focal types are simulated")
    ap.add_argument("--age", type=float, default=None, help="patient age (valproate restriction under 55, NICE NG217 Box 2)")
    ap.add_argument("--childbearing", action="store_true", help="patient is able to have children (valproate restriction)")
    ap.add_argument("--sensitivity", action="store_true")
    ap.add_argument("--name", default=None)
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    a = ap.parse_args()

    clin = None
    if a.seizure_type:
        clin = clinical.assess(a.seizure_type, list(engine.load_drugs()), age=a.age, childbearing=a.childbearing)
        print(f"Seizure type: {clin['type_label']} ({clin['class_label']})")
        for d_, r_ in clin["drugs"].items():
            print(f"  {d_:14s} {r_['label']:30s} {r_['source'] or ''} {r_['rec']}")
        for g_ in clin["guideline_only"]:
            print(f"  {g_['drug']:14s} {g_['status']:30s} (guideline option, not modelled) {g_['rec']}")
        if not clin["simulate"]:
            os.makedirs(os.path.join(HERE, "results", "clinical"), exist_ok=True)
            out = os.path.join(HERE, "results", "clinical", (a.name or a.seizure_type) + ".json")
            json.dump(clin, open(out, "w"), indent=1)
            sys.exit(f"\n{clin['why']}\nNo simulation run. Guideline assessment saved to {out}")

    if a.template:
        W, labels = connectome.load_tvb_template(); src = "tvb_template"
    else:
        W, labels = connectome.load(a.connectome); src = os.path.basename(a.connectome.rstrip("/"))
    n = len(labels)
    eeg_info = None
    ieeg_info = None
    if a.ieeg_edf:
        if not a.template:
            sys.exit("--ieeg-edf maps anatomical labels to the named TVB template regions: add --template")
        ieeg_info = ieeg.analyse(a.ieeg_edf, a.ieeg_channels, a.ieeg_events, a.ieeg_electrodes)
        mapped = [(lab_, ieeg.to_template_region(lab_, labels), sc_) for lab_, sc_ in ieeg_info["label_scores"].items()]
        mapped = [m for m in mapped if m[1]]
        ieeg_info["mapped_regions"] = [{"label": l_, "region": r_, "score": s_} for l_, r_, s_ in mapped]
        print("iEEG strongest labels -> template regions:", [(l_, r_, round(s_, 1)) for l_, r_, s_ in mapped[:5]])
        if a.onset == "auto" and mapped:
            a.onset = mapped[0][1]
    if a.eeg_edf:
        import mne
        cpath = os.path.join(a.connectome, "centres_mm.csv")
        if a.template or not os.path.exists(cpath):
            sys.exit("--eeg-edf needs a connectome folder with centres_mm.csv")
        ann = eeg.read_seizure_annotation(a.seizures)
        raw = mne.io.read_raw_edf(a.eeg_edf, preload=True, verbose="ERROR")
        eeg_info = eeg.suggest_onset(raw, ann[a.seizure_index][0], np.loadtxt(cpath, delimiter=","))
        eeg_info["seizure_window_s"] = ann[a.seizure_index]
        print(f"EEG: top channels {[(c, round(s,1)) for c, s in eeg_info['top_channels']]} -> suggested onset region index {eeg_info['suggested_region']}")
        if a.onset == "auto":
            a.onset = str(eeg_info["suggested_region"])
    if a.onset == "auto":
        onset = int(np.argsort(W.sum(1))[n // 2])
    elif a.onset.isdigit():
        onset = int(a.onset)
    elif a.onset in labels:
        onset = labels.index(a.onset)
    else:
        sys.exit(f"Onset not found. Regions: {labels}")
    name = a.name or f"{src}_onset{labels[onset]}"
    drugs = engine.load_drugs()
    excluded, flagged = {}, {}
    if a.genetics:
        variants = genetics.load_patient(a.genetics)
        drugs, excluded, flagged = genetics.apply(drugs, variants)
        print(f"Genetics: variants {variants} -> excluded {list(excluded)} flagged {list(flagged)}")
    if clin:
        dropped = [d_ for d_ in list(drugs) if d_ not in clin["eligible_for_simulation"]]
        drugs = {k: v for k, v in drugs.items() if k not in dropped}
        if dropped:
            print(f"Clinical gate removed from simulation: {dropped}")
    t0 = time.time()
    print(f"Connectome: {src} ({n} regions) | seizure onset: {labels[onset]} | drugs: {', '.join(drugs)}")

    print("1/3 Calibrating the network's seizure-spread threshold ...")
    cal = engine.calibrate(W, onset, workers=a.workers)
    if not cal["ok"]:
        sys.exit("Calibration failed: a seizure from this onset never spreads in the tested range. "
                 "Try another --onset (a better-connected region).")
    print(f"    background excitability x0 = {cal['x0_healthy']} | healthy control recruited: {cal['healthy_recruited']} regions")

    print("2/3 Finding the drug strength needed to stop spread (lower = more potent here) ...")
    th = engine.drug_thresholds(W, onset, cal, drugs, seeds=a.seeds, iters=a.iters, workers=a.workers)
    ranking = engine.rank(th)
    res = {"connectome": src, "labels": labels, "onset": labels[onset], "onset_index": onset, "calibration": cal,
           "drugs": {k: {**th[k], "targets": drugs[k]["targets"], "rationale": drugs[k]["rationale"]} for k in th},
           "ranking": ranking, "sensitivity": {}, "excluded": excluded, "flagged": flagged, "eeg": eeg_info, "ieeg": ieeg_info,
           "inputs": {"diffusion_mri": True, "eeg": bool(eeg_info), "ieeg": bool(ieeg_info), "genetics": bool(a.genetics), "seizure_type": bool(clin)},
           "clinical": clin}

    if a.sensitivity:
        print("    sensitivity: re-weighting the two mechanisms (x0 and coupling effects) x0.5 / x2 ...")
        for label, (xw, cw) in {"excitability effects x0.5": (0.5, 1.0), "excitability effects x2": (2.0, 1.0),
                                "coupling effects x0.5": (1.0, 0.5), "coupling effects x2": (1.0, 2.0)}.items():
            t = engine.drug_thresholds(W, onset, cal, drugs, seeds=a.seeds, iters=a.iters, x0_weight=xw, coupling_weight=cw, workers=a.workers)
            res["sensitivity"][label] = {"ranking": engine.rank(t), "thresholds": {k: t[k]["threshold"] for k in t}}

    print("3/3 Saving per-region snapshot ...")
    res["snapshot"] = engine.snapshot(W, onset, cal, drugs, th)

    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    out = os.path.join(HERE, "results", name + ".json")
    json.dump(res, open(out, "w"), indent=1)

    for d_, rs in excluded.items():
        print(f"  EXCLUDED {d_}: {rs[0]['reason']} ({rs[0]['source']})")
    for d_, rs in flagged.items():
        print(f"  FLAG {d_}: {rs[0]['reason']}")
    print(f"\nRanking for seizure onset {labels[onset]} in {src}  (drug strength needed to stop spread; 1.0 = nominal dose)")
    for i, k in enumerate(ranking, 1):
        t = th[k]["threshold"]
        tag = f"   [{clin['drugs'][k]['label']}]" if clin else ""
        print(f"  {i}. {k:14s} {'no suppression up to %.1f' % engine.S_MAX if t is None else '%.2f' % t}{tag}")
    for lab, v in res["sensitivity"].items():
        print(f"  [{lab}] {' > '.join(v['ranking'])}")
    stable = all(v["ranking"][0] == ranking[0] for v in res["sensitivity"].values()) if res["sensitivity"] else None
    if stable is not None:
        print("  Top drug is", "STABLE" if stable else "NOT stable", "across the sensitivity variants.")
    print(f"\nSaved {out}  ({time.time()-t0:.0f}s)")
    print("Proof-of-concept simulation with ASSUMED drug parameters. Not clinical advice.")

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.6))
    vals = [(th[k]["threshold"] if th[k]["threshold"] is not None else engine.S_MAX) for k in ranking]
    ax.barh(ranking[::-1], vals[::-1], color="#7a4fd0")
    ax.set_xlabel("drug strength needed to stop seizure spread (lower = better)")
    ax.set_title(f"ARIADNE ranking: onset {labels[onset]}, {src}")
    plt.tight_layout(); plt.savefig(os.path.join(HERE, "results", name + ".png"), dpi=120)


if __name__ == "__main__":
    main()
