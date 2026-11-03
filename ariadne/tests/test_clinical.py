"""Clinical gate: every drug status must cite a source, and the key guideline rules must hold."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ariadne import clinical

DRUGS = ["cenobamate", "lamotrigine", "valproate", "carbamazepine", "pregabalin"]
R = clinical.load_rules()
for t, T in R["types"].items():
    for d, r in T["drugs"].items():
        assert r["source"] in R["sources"], (t, d)
        assert r["status"] in clinical.STATUS_ORDER, (t, d)
f = clinical.assess("focal_impaired_consciousness", DRUGS)
assert f["simulate"] and set(f["eligible_for_simulation"]) == set(DRUGS)
assert f["drugs"]["lamotrigine"]["status"] == "first_line"
for gen in ["absence", "myoclonic", "tonic", "atonic"]:
    g = clinical.assess(gen, DRUGS)
    assert not g["simulate"] and g["drugs"]["carbamazepine"]["status"] == "avoid" and g["drugs"]["pregabalin"]["status"] == "avoid", gen
assert clinical.assess("generalized_tonic_clonic", DRUGS)["drugs"]["cenobamate"]["status"] == "not_indicated"
assert clinical.assess("generalized_tonic_clonic", DRUGS, age=30)["drugs"]["valproate"]["status"] == "restricted"
assert clinical.assess("generalized_tonic_clonic", DRUGS, age=60)["drugs"]["valproate"]["status"] == "first_line"
assert clinical.assess("generalized_tonic_clonic", DRUGS, age=60, childbearing=True)["drugs"]["valproate"]["status"] == "restricted"
assert not clinical.assess("unknown", DRUGS)["simulate"]
print("clinical gate: all checks pass;", len(R["types"]), "seizure types")
