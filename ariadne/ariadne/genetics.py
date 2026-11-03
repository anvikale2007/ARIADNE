"""Genetic-history input: apply published gene -> drug safety rules to the drug panel."""
import os, yaml

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_rules():
    return yaml.safe_load(open(os.path.join(HERE, "genetics_rules.yaml")))


def load_patient(path):
    p = yaml.safe_load(open(path)) or {}
    return [str(v).upper() for v in p.get("variants", [])]


def apply(drugs, variants):
    """Return (allowed_drugs, excluded, flagged). excluded/flagged: {drug: [ {reason, source, status, variant} ]}."""
    excluded, flagged = {}, {}
    for r in load_rules():
        hit = [v for v in variants if any(m.upper() in v for m in r["match"])]
        if not hit:
            continue
        for dname in ([r["drug"]] if isinstance(r["drug"], str) else r["drug"]):
            if dname not in drugs:
                continue
            rec = {"variant": hit[0], "reason": r["reason"], "source": r["source"], "status": r["status"]}
            (excluded if r["action"] == "exclude" else flagged).setdefault(dname, []).append(rec)
    allowed = {k: v for k, v in drugs.items() if k not in excluded}
    return allowed, excluded, flagged
