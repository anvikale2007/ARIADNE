"""Drug layer: each drug lowers excitability (x0 shift) and/or weakens coupling (K).
PLACEHOLDER: occupancy values below are NOT measured. Replace them with values
computed from ChEMBL / IUPHAR (occupancy = C / (C + Ki)) before citing anything.
effect on model:  x0 -> x0 - SX * occ_x0     K -> K * (1 - SK * occ_K)"""
SX, SK = 0.10, 0.60          # the single global scale factors (calibration knobs)

DRUGS = {                    # name: (occ_x0, occ_K, mechanism)
    "No drug":        (0.0, 0.0, "baseline"),
    "Carbamazepine":  (0.6, 0.0, "Na+ channel block"),
    "Lamotrigine":    (0.5, 0.0, "Na+ channel block"),
    "Valproate":      (0.4, 0.3, "Na+, T-type Ca2+, GABA"),
    "Cenobamate":     (0.7, 0.0, "persistent Na+ + GABA-A"),
    "Pregabalin":     (0.0, 0.6, "alpha2delta (less release)"),
    "Levetiracetam":  (0.0, 0.5, "SV2A (less release)"),
    "Ethosuximide":   (0.0, 0.0, "T-type Ca2+ (negative control)"),
}

def apply(drug, x0, K):
    ox, oK, _ = DRUGS[drug]
    x0d = x0.copy() - SX * ox
    x0d[0] = x0[0]            # origin kept as is: the drug must stop SPREAD, not the origin
    return x0d, K * (1 - SK * oK)
