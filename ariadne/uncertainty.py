"""Three analyses that turn ARIADNE's single ranking into a quantified one.

    python uncertainty.py --surface results/surface_stanford_hardi_onsetP39.json

1. THERAPEUTIC CEILING. "Strength 1.59x" is not a clinical instruction. This asks the clinical
   question instead: at the top of its therapeutic range (or its maximum licensed dose), does the drug
   contain the seizure? Because receptor occupancy saturates, doubling the dose does not double the
   effect, so the equivalent strength multiplier is reported too.

2. GABA-WEIGHT CURVE. One unsourced constant (how much a fractional GABA-A potentiation counts
   relative to sodium-channel block) decides whether cenobamate or lamotrigine ranks first. This maps
   the ranking across the whole 0-1 range and finds the exact value where it flips.

3. UNCERTAINTY PROPAGATION. Two Monte Carlo runs:
     (A) measurement uncertainty only - sampling the ranges that have sources (therapeutic ranges,
         protein-binding disagreements, IC50 reproducibility)
     (B) full uncertainty - also sampling the parameters that have no source at all
   The ranking becomes a probability rather than a claim, and the gap between A and B measures how
   much of the answer rests on assumptions rather than on measurements.

All outputs: results/uncertainty_<surface name>.json + a markdown report.
"""
import argparse, json, os, sys
import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ariadne import surface as surf, pharm

HERE = os.path.dirname(os.path.abspath(__file__))
IC50_FOLD = 2.0      # ASSUMED: log-uniform +/- 2-fold spread for a published IC50 reproduced elsewhere
DRUGS_ORDER = ["cenobamate", "lamotrigine", "valproate", "carbamazepine", "pregabalin", "levetiracetam"]


def load_sources():
    S = yaml.safe_load(open(os.path.join(HERE, "drug_sources.yaml")))
    mode = S.pop("nav_mode"); nav_ref = S.pop("nav_ref"); gw = S.pop("gaba_weight")
    return S, mode, nav_ref, gw


def rank_of(S, Sf, gw, nav_ref, mode, free=None, coupling=None):
    """Thresholds + ranking for one parameter set.

    S may be a perturbed copy of the drug dictionary (IC50s moved, assumed effects replaced).
    `free` overrides each drug's exposure; `coupling` sets the transmitter-release reduction outright."""
    th = {}
    for n, d in S.items():
        f = (free or {}).get(n)
        if f is None:
            f = pharm.standard_free_uM(d)
        e = pharm.effects(d, f, gw, nav_ref, mode, std_free_uM=pharm.standard_free_uM(d),
                          ca_cooperativity=d.get("_coop", pharm.CA_COOPERATIVITY))
        cp = e["coupling_effect"] if not coupling or n not in coupling else coupling[n]
        th[n] = Sf.threshold(e["x0_effect"], cp)
    return th, surf.rank(th)


# ----------------------------------------------------------------- 1. therapeutic ceiling
def ceiling_analysis(S, Sf, gw, nav_ref, mode):
    rows = []
    for n, d in S.items():
        std, ceil = pharm.standard_free_uM(d), pharm.ceiling_free_uM(d)
        e_s = pharm.effects(d, std, gw, nav_ref, mode)
        e_c = pharm.effects(d, ceil, gw, nav_ref, mode)
        t_s = Sf.threshold(e_s["x0_effect"], e_s["coupling_effect"])
        # the dose increase expressed as an equivalent strength multiplier (sub-linear, because occupancy saturates)
        equiv = e_c["x0_effect"] / e_s["x0_effect"] if e_s["x0_effect"] > 0 else (
            e_c["coupling_effect"] / e_s["coupling_effect"] if e_s["coupling_effect"] > 0 else 1.0)
        rows.append({"drug": n, "basis": d.get("ceiling_basis", ""),
                     "free_standard_uM": round(std, 1), "free_ceiling_uM": round(ceil, 1),
                     "conc_ratio": round(ceil / std, 2), "effect_ratio": round(equiv, 2),
                     "x0_standard": round(e_s["x0_effect"], 3), "x0_ceiling": round(e_c["x0_effect"], 3),
                     "coupling_standard": round(e_s["coupling_effect"], 3), "coupling_ceiling": round(e_c["coupling_effect"], 3),
                     "threshold_strength": t_s,
                     "contained_at_standard": bool(Sf.contained_at(e_s["x0_effect"], e_s["coupling_effect"])),
                     "contained_at_ceiling": bool(Sf.contained_at(e_c["x0_effect"], e_c["coupling_effect"]))})
    return rows


# ----------------------------------------------------------------- 2. GABA weight curve
def gaba_curve(S, Sf, nav_ref, mode, n=51):
    out = []
    for gw in np.round(np.linspace(0, 1, n), 3):
        th, rk = rank_of(S, Sf, float(gw), nav_ref, mode)
        out.append({"gaba_weight": float(gw), "ranking": rk,
                    "thresholds": {k: th[k] for k in DRUGS_ORDER if k in th}})
    flip = None
    for a, b in zip(out, out[1:]):
        if a["ranking"][0] != b["ranking"][0]:
            flip = {"between": [a["gaba_weight"], b["gaba_weight"]], "from": a["ranking"][0], "to": b["ranking"][0]}
    return out, flip


# ----------------------------------------------------------------- 2b. the levetiracetam question, inverted
def release_needed(S, Sf, gw, nav_ref, mode, targets=("lamotrigine", "cenobamate")):
    """Levetiracetam is the most-prescribed focal drug and is NOT modelled, because no published number
    converts SV2A binding into a size of effect. Rather than invent one, invert the question: how large
    would a purely release-reducing drug's effect have to be, at a standard dose, to match each modelled
    drug? That is a falsifiable number an experiment could later check."""
    th, _ = rank_of(S, Sf, gw, nav_ref, mode)
    out = []
    # the coupling scale at which the seizure is contained with no excitability effect at all
    c_star = None
    for c, b in sorted(zip(Sf._c, Sf._b)):
        if np.isfinite(b) and b <= 0:
            c_star = max(c_star or 0, float(c))
    for tgt in targets:
        t = th.get(tgt)
        if t is None or c_star is None:
            out.append({"to_match": tgt, "target_threshold": t, "release_reduction_needed": None,
                        "feasible": None, "interpretable_as": "target never contains the seizure"})
            continue
        # need 1 - t * r * COUPLING_UNIT <= c_star  ->  r >= (1 - c_star) / (t * COUPLING_UNIT)
        need = (1.0 - c_star) / (t * Sf.COUPLING_UNIT)
        feasible = need <= 1.0
        out.append({"to_match": tgt, "target_threshold": t, "release_reduction_needed": round(float(need), 3),
                    "feasible": bool(feasible), "coupling_scale_for_containment": round(float(c_star), 3),
                    "interpretable_as": (f"a {need*100:.0f}% reduction in transmitter release at a standard dose"
                                         if feasible else
                                         f"a {need*100:.0f}% reduction in transmitter release - impossible, "
                                         f"since the maximum conceivable reduction is 100%")})
    return out


# ----------------------------------------------------------------- 3. Monte Carlo
IC50_KEYS = ("nav_vdb_ic50_uM", "nav_persistent_ic50_uM", "nav_udb10_ic50_uM", "gaba_ec50_uM")


def draw(S, rng, full):
    """One parameter draw. Returns (perturbed drug dicts, free concentrations, coupling overrides, gaba_weight, nav_ref)."""
    P, free, coupling = {}, {}, {}
    for n, d in S.items():
        p = dict(d)
        # IC50 reproducibility: each measured potency is moved log-uniformly within +/- IC50_FOLD.
        # Perturbing the IC50 itself (rather than the final effect) keeps the sodium and GABA
        # contributions separate and lets occupancy saturate correctly.
        for k in IC50_KEYS:
            if p.get(k):
                p[k] = float(p[k]) * float(np.exp(rng.uniform(-np.log(IC50_FOLD), np.log(IC50_FOLD))))
        for k in ("nav_block_override", "nav_vdb_block_override"):
            if p.get(k):
                p[k] = float(np.clip(p[k] * np.exp(rng.uniform(-np.log(IC50_FOLD), np.log(IC50_FOLD))), 0, 1))

        # calcium-current block and release cooperativity are both MEASURED, so their uncertainty belongs
        # in the measurement tier: block between the paper's two reported values (total HVA vs isolated
        # N-type), cooperativity across the Hill coefficients Reid 1998 reports (2.5-3.3).
        if p.get("ca_current_block") is not None:
            lo_b = float(p["ca_current_block"])
            hi_b = float(p.get("ca_current_block_high") or lo_b)
            p["ca_current_block"] = float(rng.uniform(min(lo_b, hi_b), max(lo_b, hi_b)))
            p["_coop"] = float(rng.uniform(2.5, 3.3))

        # exposure: sampled across the sourced therapeutic range / licensed dose range
        lo, hi = pharm.standard_free_uM(d), pharm.ceiling_free_uM(d)
        if d.get("therapeutic_mg_per_L"):
            a, _ = d["therapeutic_mg_per_L"]
            fu_lo = pharm.valproate_free_fraction(a) if d.get("concentration_dependent_binding") else d["free_fraction"]
            lo = pharm.MG_L_TO_UM(a, d["mw"]) * fu_lo
        f = float(rng.uniform(min(lo, hi), max(lo, hi)))
        # lamotrigine protein binding: two sourced values disagree (55% vs 66% bound)
        if n == "lamotrigine":
            f *= float(rng.uniform(0.34, 0.45)) / 0.45
        free[n] = f
        P[n] = p

    gw, nav_ref = 0.5, 0.30
    if full:
        gw = float(rng.uniform(0.0, 1.0))                                   # unsourced
        nav_ref = float(np.exp(rng.uniform(np.log(0.15), np.log(0.60))))    # unsourced scaling constant
        for n, d in S.items():
            if (d.get("release_reduction") or 0) > 0:                       # unsourced transmitter-release reduction
                P[n]["release_reduction"] = float(rng.uniform(0.0, 1.0))
            if (d.get("other_x0") or 0) > 0:                                # unsourced lumped excitability term
                P[n]["other_x0"] = float(rng.uniform(0.0, 1.0))
    return P, free, coupling, gw, nav_ref


def monte_carlo(S, Sf, mode, n_draws, full, seed=0):
    rng = np.random.default_rng(seed)
    firsts, orders, th_all = [], [], {k: [] for k in S}
    for _ in range(n_draws):
        P, free, coupling, gw, nav_ref = draw(S, rng, full)
        th, rk = rank_of(P, Sf, gw, nav_ref, mode, free=free, coupling=coupling)
        firsts.append(rk[0]); orders.append(tuple(rk))
        for k, v in th.items():
            th_all[k].append(v)
    p_first = {k: firsts.count(k) / n_draws for k in sorted(set(firsts), key=lambda x: -firsts.count(x))}
    top_orders = sorted({o: orders.count(o) / n_draws for o in set(orders)}.items(), key=lambda x: -x[1])[:5]
    pct = {}
    for k, v in th_all.items():
        got = [x for x in v if x is not None]
        pct[k] = {"p_contains_within_4x": len(got) / n_draws,
                  "median": None if not got else round(float(np.median(got)), 3),
                  "p10": None if not got else round(float(np.percentile(got, 10)), 3),
                  "p90": None if not got else round(float(np.percentile(got, 90)), 3)}
    # pairwise: how often does cenobamate beat lamotrigine, etc.
    pair = {}
    for i, a in enumerate(DRUGS_ORDER):
        for b in DRUGS_ORDER[i + 1:]:
            if a in th_all and b in th_all:
                wins = sum(1 for x, y in zip(th_all[a], th_all[b])
                           if (x is not None and y is None) or (x is not None and y is not None and x < y))
                pair[f"{a} < {b}"] = round(wins / n_draws, 3)
    return {"n_draws": n_draws, "full_uncertainty": full, "p_ranks_first": p_first,
            "most_likely_orders": [{"order": list(o), "p": round(p, 3)} for o, p in top_orders],
            "thresholds": pct, "pairwise_better": pair, "lookups_outside_mapped_range": Sf.clamped}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--surface", default=os.path.join(HERE, "results", "surface_stanford_hardi_onsetP39.json"))
    ap.add_argument("--draws", type=int, default=2000)
    a = ap.parse_args()

    Sf = surf.load(a.surface)
    S, mode, nav_ref, gw = load_sources()
    S = {k: S[k] for k in DRUGS_ORDER if k in S}
    name = os.path.basename(a.surface).replace("surface_", "").replace(".json", "")
    print(f"Surface: {Sf.d['name']} ({Sf.d['n_regions']} regions, onset {Sf.d['onset']})")

    base_th, base_rk = rank_of(S, Sf, gw, nav_ref, mode)
    print("Baseline from surface:", {k: base_th[k] for k in base_rk})

    print("\n1/3 therapeutic ceiling ...")
    ceil = ceiling_analysis(S, Sf, gw, nav_ref, mode)
    print("\n2/3 GABA-weight curve ...")
    curve, flip = gaba_curve(S, Sf, nav_ref, mode)
    print("   flip:", flip)
    inv = release_needed(S, Sf, gw, nav_ref, mode)
    for r in inv:
        print(f"   to match {r['to_match']}: {r['interpretable_as'] if r['release_reduction_needed'] is not None else 'unreachable'}")
    print(f"\n3/3 Monte Carlo ({a.draws} draws x2) ...")
    Sf.clamped = 0
    mc_sourced = monte_carlo(S, Sf, mode, a.draws, full=False, seed=1)
    Sf.clamped = 0
    mc_full = monte_carlo(S, Sf, mode, a.draws, full=True, seed=2)

    out = {"surface": os.path.basename(a.surface), "brain": Sf.d["name"], "onset": Sf.d["onset"],
           "baseline": {"gaba_weight": gw, "nav_ref": nav_ref, "nav_mode": mode,
                        "thresholds": base_th, "ranking": base_rk},
           "therapeutic_ceiling": ceil, "gaba_curve": curve, "gaba_flip": flip,
           "unmodelled_drug_inverse": inv,
           "monte_carlo_sourced": mc_sourced, "monte_carlo_full": mc_full,
           "ic50_fold_assumed": IC50_FOLD}
    p = os.path.join(HERE, "results", f"uncertainty_{name}.json")
    json.dump(out, open(p, "w"), indent=1)
    print("\nSaved", p)
    print("\nP(ranks first), measurement uncertainty only:", mc_sourced["p_ranks_first"])
    print("P(ranks first), full uncertainty:            ", mc_full["p_ranks_first"])


if __name__ == "__main__":
    main()
