"""Concentration -> model effect, in one place.

Everything uses Hill n=1 occupancy at the free (unbound) concentration:
    occupancy(C) = C / (C + IC50)
    x0_effect       = (nav_block + GABA_WEIGHT * gaba_fraction) / NAV_REF + other_x0
    coupling_effect = release_reduction

Two things this module adds over the original derive_parameters.py:

1. It can evaluate a drug at ANY exposure, not just the one modelled. Because occupancy saturates,
   doubling the dose does not double the effect, so "2x the strength" and "2x the dose" are different
   questions. The therapeutic-ceiling analysis needs the second one.

2. Where a block fraction was entered as a fixed override measured at one concentration (valproate),
   the implied IC50 is reconstructed so the same drug can be evaluated elsewhere. At the original
   concentration the result is unchanged by construction.

Effects marked ASSUMED in drug_sources.yaml (other_x0, release_reduction) have no measured EC50. They
are defined as "the effect at the standard therapeutic exposure", and are extrapolated to other
exposures with the same Hill form by back-calculating an implied EC50. That is a modelling choice, not
a measurement, and is reported as such.
"""

MG_L_TO_UM = lambda mg_per_L, mw: mg_per_L / mw * 1000.0


def occupancy(free_uM, ic50_uM):
    if ic50_uM is None or free_uM is None:
        return 0.0
    return float(free_uM) / (float(free_uM) + float(ic50_uM))


def implied_ic50(free_uM, fraction):
    """IC50 that reproduces this fraction at this concentration under Hill n=1."""
    f = float(fraction)
    if f <= 0:
        return None
    if f >= 1:
        return 0.0
    return float(free_uM) * (1.0 - f) / f


def valproate_free_fraction(total_mg_per_L):
    """FDA label anchors: 10% free at 40 mg/L, 18.5% at 130 mg/L -> fu = 0.0622 + 0.000944*C."""
    return 0.0622 + 0.000944 * float(total_mg_per_L)


def standard_free_uM(d):
    """Free concentration the drug is modelled at (mid therapeutic range, or stated maintenance dose)."""
    if d.get("free_conc_uM") is not None:
        return float(d["free_conc_uM"])
    lo, hi = d["therapeutic_mg_per_L"]
    mid = (lo + hi) / 2.0
    fu = valproate_free_fraction(mid) if d.get("concentration_dependent_binding") else d["free_fraction"]
    return MG_L_TO_UM(mid, d["mw"]) * fu


def ceiling_free_uM(d):
    """Free concentration at the therapeutic ceiling: top of the reference range, or max licensed dose."""
    if d.get("ceiling_free_uM") is not None:
        return float(d["ceiling_free_uM"])
    if d.get("ceiling_mg_per_L") is None:
        return None
    c = float(d["ceiling_mg_per_L"])
    fu = valproate_free_fraction(c) if d.get("concentration_dependent_binding") else d["free_fraction"]
    return MG_L_TO_UM(c, d["mw"]) * fu


CA_COOPERATIVITY = 2.6   # sourced: N-type Hill coefficient, Reid/Bekkers/Clements J Neurosci 1998;18:2849-2855


def effects(d, free_uM, gaba_weight, nav_ref, nav_mode="inactivated", std_free_uM=None,
            ca_cooperativity=CA_COOPERATIVITY, high_ca=False):
    """Model effects for drug dict `d` at free concentration `free_uM`.

    std_free_uM anchors the overrides and the assumed effects; defaults to the drug's standard exposure."""
    std = standard_free_uM(d) if std_free_uM is None else float(std_free_uM)

    # --- sodium-channel block ---
    if nav_mode == "inactivated":
        ic50 = d.get("nav_vdb_ic50_uM")
        override = d.get("nav_vdb_block_override")
    elif nav_mode == "udb10":
        ic50 = d.get("nav_udb10_ic50_uM")
        override = None if d.get("nav_udb10_ic50_uM") else d.get("nav_block_override")
    else:
        ic50 = d.get("nav_persistent_ic50_uM")
        override = d.get("nav_block_override")
    if ic50 is not None:
        nav = occupancy(free_uM, ic50)
    elif override is not None:
        nav = occupancy(free_uM, implied_ic50(std, override)) if override > 0 else 0.0
    else:
        nav = 0.0

    # --- GABA-A potentiation ---
    gaba = occupancy(free_uM, d["gaba_ec50_uM"]) if d.get("gaba_ec50_uM") else 0.0

    # --- assumed effects, extrapolated with the same Hill form ---
    def assumed(value):
        v = float(value or 0.0)
        if v <= 0:
            return 0.0
        return occupancy(free_uM, implied_ic50(std, v))

    other = assumed(d.get("other_x0"))

    # --- transmitter-release reduction ---
    # Preferred: DERIVE it from a measured calcium-current block, because release is a supralinear
    # function of presynaptic calcium entry. Reid, Bekkers & Clements (J Neurosci 1998) measured the
    # cooperativity at rat hippocampal autapses: Hill coefficient 2.6 for N-type, 2.5 for P/Q-type
    # (3.1-3.3 for the unblocked preparation, which the authors argue underestimates the true value).
    #     release_remaining = (1 - ca_block) ** n     ->     reduction = 1 - (1 - ca_block) ** n
    # Fallback: an assumed release_reduction, extrapolated with the Hill form like any other assumed term.
    if d.get("ca_current_block") is not None:
        blk = float(d["ca_current_block_high"] if (high_ca and d.get("ca_current_block_high") is not None)
                    else d["ca_current_block"])
        release = 1.0 - (1.0 - blk) ** float(ca_cooperativity)
    else:
        release = assumed(d.get("release_reduction"))

    x0 = (nav + gaba_weight * gaba) / nav_ref + other
    return {"free_uM": float(free_uM), "nav_block": nav, "gaba_fraction": gaba,
            "other_x0": other, "release_reduction": release,
            "x0_effect": x0, "coupling_effect": release}
