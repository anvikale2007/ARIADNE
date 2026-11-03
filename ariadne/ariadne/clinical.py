"""Clinical gate: ILAE 2025 seizure type + simple patient factors -> guideline status of each drug, and whether
ARIADNE's focal-spread simulation applies. Rules and sources live in clinical_rules.yaml (no numbers invented here)."""
import os, yaml

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATUS_ORDER = ["first_line", "later_line", "restricted", "not_listed", "not_indicated", "avoid"]
LABEL = {"first_line": "First-line", "later_line": "Later-line / add-on", "restricted": "Restricted (valproate rules)",
         "not_listed": "Not a guideline option", "not_indicated": "Not licensed for this", "avoid": "Avoid: may worsen"}


def load_rules():
    return yaml.safe_load(open(os.path.join(HERE, "clinical_rules.yaml")))


def seizure_types():
    R = load_rules()
    return {k: {"label": v["label"], "class": v["class"]} for k, v in R["types"].items()}


def assess(seizure_type, drugs, age=None, childbearing=False):
    """Return a dict: class info, simulate flag, per-drug status (for the modelled drugs), guideline-only options,
    and which drugs may enter the simulation. Unknown drug names get status 'not_listed' with no source."""
    R = load_rules()
    if seizure_type not in R["types"]:
        raise ValueError(f"unknown seizure type '{seizure_type}'. Options: {', '.join(R['types'])}")
    T = R["types"][seizure_type]
    C = R["classes"][T["class"]]
    vr = R["valproate_restriction"]
    vpa_restricted = childbearing or (age is not None and age < 55)
    per = {}
    for d in drugs:
        r = dict(T["drugs"].get(d, {"status": "not_listed", "source": None, "rec": "no rule for this seizure type"}))
        r["source_text"] = R["sources"].get(r.get("source"), "")
        if d == "valproate" and vpa_restricted and r["status"] in ("first_line", "later_line"):
            r["guideline_status"] = r["status"]
            r["status"] = "restricted"
            r["rec"] += "; " + vr["rule"]
        r["label"] = LABEL[r["status"]]
        per[d] = r
    simulate = bool(C["simulate"])
    eligible = [d for d in drugs if simulate and per[d]["status"] not in ("avoid", "not_indicated")]
    return {"seizure_type": seizure_type, "type_label": T["label"], "class": T["class"], "class_label": C["label"],
            "simulate": simulate, "why": C["why"], "drugs": per, "eligible_for_simulation": eligible,
            "guideline_only": T.get("guideline_only", []),
            "patient": {"age": age, "childbearing": bool(childbearing), "valproate_restricted": vpa_restricted},
            "classification_source": R["sources"]["ILAE2025"]}
