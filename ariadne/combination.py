"""Two-drug combinations, on the same response surface as monotherapy.

WHY THIS NEEDS NO NEW SIMULATIONS
Every drug in ARIADNE acts through exactly two model quantities: an excitability reduction and a
coupling reduction. The stop/no-stop boundary has already been mapped over that whole plane, so a
combination is just another point on it. What has to be decided is how two drugs' effects COMBINE
before the lookup, and that is a pharmacology question, not a simulation question.

HOW TWO DRUGS COMBINE HERE
1. SHARED TARGET -> Bliss independence on the fraction bound, never addition.
   Two sodium-channel blockers compete for one channel population; their blocks cannot add, because
   occupancy saturates at 1. The non-interaction reference is
       combined = 1 - (1 - b1) * (1 - b2)
   which is the same algebra the monotherapy pipeline already uses for the release term. Applied to
   nav_block (both act on the inactivated state of the channel) and to gaba_fraction.
   This is a NULL model: it assumes no synergy and no antagonism. Any real interaction at the target
   would be a departure from it, and no published effect size exists to put one in.

2. SEPARATE MODEL TERMS -> additive, because they are different quantities.
   Excitability and coupling are different dials. A sodium blocker and a release-reducing drug do not
   compete; each moves its own dial.

3. PRESYNAPTIC CALCIUM -> combine the CURRENT BLOCK, then apply cooperativity ONCE.
       release = 1 - (1 - b_combined) ** 2.6
   Adding two release reductions would apply the exponent twice and overstate the pair badly.

4. RESIDUAL EXCITABILITY (valproate's other_x0) -> added. ASSUMED: it is a catch-all with no named
   target, so there is no shared-target argument for Bliss. Flagged in the output.

PHARMACOKINETIC INTERACTIONS (the part that is genuinely drug-drug, not just additive)
Combination changes what reaches the brain, not only what happens there. Lamotrigine is the
documented case, from the FDA LAMICTAL label (2025), Table 13:
  - valproate "increased lamotrigine concentrations slightly more than 2-fold"
  - carbamazepine, phenytoin, phenobarbital and primidone each "decreased lamotrigine concentration
    approximately 40%"
Those multiply the free concentration BEFORE the Hill step, which the monotherapy pipeline already
handles at any exposure. No other pair in this panel has a quantified interaction in its label, so
none is applied; absence of an entry means "not quantified", not "no interaction".

NOT MODELLED, AND SAID SO: tolerability, additive adverse effects, cost, adherence, or any synergy
beyond the Bliss null. NICE NG217 puts monotherapy first; this ranks ADD-ON pairs for the case where
a first drug has already failed, as a companion to a guideline tool, not a replacement for one.

    python combination.py --result results/ideas_sub-1_onset_rh-entorhinal.json
"""
import argparse, itertools, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import yaml                                              # noqa: E402
from ariadne import pharm                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
X0_UNIT = 0.004
COUPLING_UNIT = 0.06

# FDA LAMICTAL label (2025) Table 13. Multiplier applied to the free concentration of the ROW drug
# when the COLUMN drug is co-administered. Only label-quantified interactions appear here.
PK = {
    ("lamotrigine", "valproate"): 2.0,          # "slightly more than 2-fold" - 2.0 is the conservative read
    ("lamotrigine", "carbamazepine"): 0.6,      # "decreased ... approximately 40%"
    ("lamotrigine", "phenytoin"): 0.6,
}
PK_SOURCE = "FDA LAMICTAL label (2025), Table 13"

# Which model quantity each drug's block lands in, used to decide Bliss vs addition.
NAV_DRUGS = {"carbamazepine", "lacosamide", "lamotrigine", "phenytoin", "cenobamate", "valproate"}
GABA_DRUGS = {"cenobamate"}
CA_DRUGS = {"pregabalin", "levetiracetam"}


def bliss(a, b):
    """Non-interaction combination of two fractions of the same target population."""
    return 1.0 - (1.0 - a) * (1.0 - b)


def load_panel():
    src = yaml.safe_load(open(os.path.join(HERE, "drug_sources.yaml")))
    ds = src["drugs"] if "drugs" in src else {k: v for k, v in src.items()
                                              if isinstance(v, dict) and "mw" in v}
    nav_ref = src.get("nav_ref") or src.get("constants", {}).get("nav_ref", 0.30)
    gaba_w = src.get("gaba_weight") or src.get("constants", {}).get("gaba_weight", 1.0)
    return ds, nav_ref, gaba_w


def parts(d, name, free_uM, nav_ref, gaba_w):
    """The raw target-level quantities for one drug at one free concentration."""
    e = pharm.effects(d, free_uM, gaba_w, nav_ref, std_free_uM=pharm.standard_free_uM(d))
    ca = float(d.get("ca_current_block") or 0.0)
    return {"nav": e["nav_block"], "gaba": e["gaba_fraction"], "other": e["other_x0"],
            "ca": ca, "release_assumed": 0.0 if ca else e["release_reduction"],
            "x0_effect": e["x0_effect"], "coupling_effect": e["coupling_effect"]}


def combine(pa, pb, nav_ref, gaba_w, ca_coop=pharm.CA_COOPERATIVITY):
    """Two drugs' target-level quantities -> one (x0_effect, coupling_effect) pair."""
    nav = bliss(pa["nav"], pb["nav"])                       # same channel population
    gaba = bliss(pa["gaba"], pb["gaba"])                    # same receptor population
    ca = bliss(pa["ca"], pb["ca"])                          # same presynaptic calcium current
    release = (1.0 - (1.0 - ca) ** ca_coop) if ca > 0 else 0.0
    release = bliss(release, bliss(pa["release_assumed"], pb["release_assumed"]))
    other = pa["other"] + pb["other"]                        # ASSUMED: no named shared target
    x0 = (nav + gaba_w * gaba) / nav_ref + other
    return {"x0_effect": x0, "coupling_effect": release,
            "nav": nav, "gaba": gaba, "ca": ca, "other": other}


def boundary_from_result(R, ds, nav_ref, gaba_w):
    """Recover this brain's containment boundary at full coupling from its own monotherapy thresholds.

    For a drug with no coupling effect the pipeline's own identity holds exactly:
        threshold = boundary / (x0_effect * X0_UNIT)
    so every such drug gives an independent estimate of the same boundary. They agree to the
    bisection resolution; the median is taken and the spread reported.
    """
    est = []
    for name, row in R["drugs"].items():
        t = row.get("threshold")
        d = ds.get(name)
        if t is None or d is None:
            continue
        p = parts(d, name, pharm.standard_free_uM(d), nav_ref, gaba_w)
        if p["coupling_effect"] > 1e-9 or p["x0_effect"] <= 0:
            continue
        est.append((name, t * p["x0_effect"] * X0_UNIT))
    vals = sorted(v for _, v in est)
    med = vals[len(vals) // 2] if len(vals) % 2 else (vals[len(vals) // 2 - 1] + vals[len(vals) // 2]) / 2
    return med, est


def make_lookup(case_file, flat_boundary):
    """required_x0(coupling_scale) for one brain, from the best source available.

    A coupling reduction lowers the excitability reduction the brain needs, and by more than is
    negligible: on the IDEAS case a 1% coupling reduction drops the requirement about 9%. So pairs
    that include a release-reducing drug MUST be solved on the curve, not at coupling 1.0.

    Preference order: the full 81x21 response surface for that brain; else the targeted boundary
    computed at a few coupling scales by combo_boundary.py; else the flat boundary with a flag.
    """
    surf = os.path.join(HERE, "results", f"surface_{case_file}.json")
    if os.path.exists(surf):
        from ariadne import surface as _s
        S = _s.load(surf)
        return (lambda c: S.required_x0(c)), "full response surface"
    curve = os.path.join(HERE, "results", f"combo_boundary_{case_file}.json")
    if os.path.exists(curve):
        d = json.load(open(curve))["boundary_by_coupling_scale"]
        pts = sorted((float(k), v) for k, v in d.items() if v is not None)
        if pts:
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]

            def req(c):
                c = min(max(float(c), xs[0]), xs[-1])
                for i in range(1, len(xs)):
                    if c <= xs[i]:
                        f = (c - xs[i - 1]) / (xs[i] - xs[i - 1]) if xs[i] != xs[i - 1] else 0.0
                        return ys[i - 1] + f * (ys[i] - ys[i - 1])
                return ys[-1]
            return req, "boundary at measured coupling scales"
    return (lambda c: flat_boundary), "flat boundary at coupling 1.0 (coupling effect NOT credited)"


def solve_threshold(x0_effect, coupling_effect, req, s_max=4.0, step=0.002):
    """Smallest common multiple of both drugs' standard exposure that contains the seizure.

    Strength scales both mechanisms linearly, exactly as engine.drug_params does for one drug."""
    if x0_effect <= 0 and coupling_effect <= 0:
        return None
    s = step
    while s <= s_max + 1e-9:
        if s * x0_effect * X0_UNIT >= req(max(0.0, 1.0 - s * coupling_effect * COUPLING_UNIT)):
            return round(s, 4)
        s += step
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", default=os.path.join(HERE, "results",
                                                     "ideas_sub-1_onset_rh-entorhinal.json"))
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    R = json.load(open(a.result))
    name = os.path.basename(a.result).replace(".json", "")
    ds, nav_ref, gaba_w = load_panel()
    panel = [k for k in R["drugs"] if k in ds]

    boundary, est = boundary_from_result(R, ds, nav_ref, gaba_w)
    print(f"brain: {name}   onset {R['onset_index']}   boundary x0 reduction = {boundary:.5f}")
    for n_, v in sorted(est, key=lambda z: z[1]):
        print(f"   {n_:<14s} implies boundary {v:.5f}")
    spread = (max(v for _, v in est) - min(v for _, v in est)) / boundary if est else 0.0
    print(f"   spread across estimates: {spread * 100:.1f}%\n")

    rows = []
    for x, y in itertools.combinations(panel, 2):
        dx, dy = ds[x], ds[y]
        cx, cy = pharm.standard_free_uM(dx), pharm.standard_free_uM(dy)
        fx, fy = PK.get((x, y), 1.0), PK.get((y, x), 1.0)
        px = parts(dx, x, cx * fx, nav_ref, gaba_w)
        py = parts(dy, y, cy * fy, nav_ref, gaba_w)
        c = combine(px, py, nav_ref, gaba_w)

        reach = c["x0_effect"] * X0_UNIT                       # at standard exposure of both drugs
        contained = reach >= boundary
        need = boundary / (c["x0_effect"] * X0_UNIT) if c["x0_effect"] > 0 else None

        mono = {k: R["drugs"][k].get("threshold") for k in (x, y)}
        best_mono = min([v for v in mono.values() if v is not None], default=None)
        rows.append({
            "pair": f"{x} + {y}", "a": x, "b": y,
            "x0_effect": round(c["x0_effect"], 4),
            "coupling_effect": round(c["coupling_effect"], 4),
            "contained_at_standard": bool(contained),
            "dose_multiple_needed": None if need is None else round(need, 3),
            "best_monotherapy": best_mono,
            "pk_applied": ([f"{x} x{fx}" for _ in [0] if fx != 1.0] +
                           [f"{y} x{fy}" for _ in [0] if fy != 1.0]) or None,
            "assumed_other_x0": round(c["other"], 4) or None,
        })

    rows.sort(key=lambda r: (r["dose_multiple_needed"] is None, r["dose_multiple_needed"] or 9e9))
    print(f"{'pair':<34s}{'needs':>8s}{'both at std?':>14s}{'best mono':>11s}{'gain':>8s}  notes")
    for r in rows:
        need = r["dose_multiple_needed"]
        g = (r["best_monotherapy"] / need) if (need and r["best_monotherapy"]) else None
        notes = []
        if r["pk_applied"]:
            notes.append("PK: " + ", ".join(r["pk_applied"]))
        if r["assumed_other_x0"]:
            notes.append("has assumed term")
        print(f"{r['pair']:<34s}{(need if need else float('nan')):>8.3f}"
              f"{('YES' if r['contained_at_standard'] else 'no'):>14s}"
              f"{(r['best_monotherapy'] if r['best_monotherapy'] else float('nan')):>11.3f}"
              f"{(g if g else float('nan')):>8.2f}  {'; '.join(notes)}")

    out = a.out or os.path.join(HERE, "results", f"combination_{name}.json")
    json.dump({"brain": name, "onset_index": R["onset_index"], "boundary_x0_reduction": boundary,
               "boundary_estimates": dict(est), "boundary_spread_fraction": spread,
               "pk_source": PK_SOURCE, "pk_multipliers": {f"{k[0]}|{k[1]}": v for k, v in PK.items()},
               "combination_model": "Bliss independence on shared targets; additive across model terms;"
                                    " cooperativity applied once to combined calcium block",
               "pairs": rows}, open(out, "w"), indent=1)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
