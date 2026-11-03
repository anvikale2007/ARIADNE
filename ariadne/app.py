"""ARIADNE web app.   streamlit run app.py"""
import glob, json, os
import numpy as np
import pandas as pd
import yaml
import streamlit as st
from ariadne import clinical
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
P = lambda *a: os.path.join(HERE, *a)
ACCENT, QUIET, SEIZE, ONSET = "#5b3fb8", "#2e9e5b", "#d63030", "#e0a800"

st.set_page_config(page_title="ARIADNE", page_icon="🧠", layout="wide")
st.markdown(f"""
<style>
  .block-container {{padding-top: 1.6rem; max-width: 1200px;}}
  .ariadne-title {{font-size: 2.3rem; font-weight: 800; letter-spacing: -0.02em; margin-bottom: 0;}}
  .ariadne-sub {{color: #6b6b80; margin-top: .2rem;}}
  .pill {{display:inline-block; padding:.15rem .6rem; border-radius:999px; font-size:.78rem; margin-right:.35rem;
          border:1px solid {ACCENT}33; background:{ACCENT}12; color:{ACCENT};}}
  .note {{border-left:3px solid {ONSET}; padding:.5rem .8rem; background:#e0a80012; border-radius:4px; font-size:.9rem;}}
</style>""", unsafe_allow_html=True)

st.markdown('<div class="ariadne-title">ARIADNE</div>', unsafe_allow_html=True)
st.markdown('<div class="ariadne-sub">Connectome-based simulation of how antiseizure drugs change seizure spread, '
            'ranked on real brain wiring with published drug potencies.</div>', unsafe_allow_html=True)
st.markdown('<span class="pill">Proof of concept</span><span class="pill">Focal seizures</span>'
            '<span class="pill">Not a clinical tool</span>', unsafe_allow_html=True)


def load(path, default=None):
    try:
        return open(path).read()
    except OSError:
        return default


def doc(path, default="Not found."):
    """Render a markdown report inside a tab: drop its H1 title, keep single line breaks."""
    t = load(path)
    if t is None:
        return default
    lines = [l for l in t.splitlines() if not l.startswith("# ")]
    return "\n".join(l if l.startswith("|") or not l.strip() else l + "  " for l in lines)


FRIENDLY = {
    "stanford_hardi": "Healthy volunteer: Stanford HARDI scan",
    "sherbrooke_3shell": "Healthy volunteer: Sherbrooke 3-shell scan",
    "ideas_sub-1": "Epilepsy patient: IDEAS II sub-1",
    "penn_pid02": "Epilepsy patient: Penn PID02 (iEEG onset, template brain)",
    "demo_genetics_hlab1502": "Genetics demo: HLA-B*15:02 carrier",
    "demo_genetics_dravet": "Genetics demo: Dravet / SCN1A loss of function",
}


def friendly(path):
    b = os.path.basename(path).replace(".json", "")
    for k, v in FRIENDLY.items():
        if b.startswith(k):
            return v
    return b


files = sorted(f for f in glob.glob(P("results", "*.json")) if not os.path.basename(f).startswith(("robustness",)))
if not files:
    st.warning("No results yet. Run `python run_pipeline.py` first."); st.stop()

with st.sidebar:
    st.header("Case")
    main = [i for i, f in enumerate(files) if os.path.basename(f).startswith("stanford_hardi")]
    pick = st.selectbox("Choose a brain / scenario", files, index=main[0] if main else 0, format_func=friendly)
    R = json.load(open(pick))
    labels, ranking = R["labels"], R["ranking"]
    st.markdown(f"**Connectome:** `{R['connectome']}`  \n**Regions:** {len(labels)}  \n**Seizure onset:** `{R['onset']}`")
    inp = R.get("inputs", {"diffusion_mri": True})
    st.markdown("**Inputs used:** " + ", ".join(k.replace("_", " ") for k, v in inp.items() if v))
    st.divider()
    st.header("Patient")
    TYPES = clinical.seizure_types()
    stype = st.selectbox("Seizure type (ILAE 2025)", list(TYPES), index=list(TYPES).index("focal_impaired_consciousness"),
                         format_func=lambda k: f"{TYPES[k]['class'].capitalize()} · {TYPES[k]['label']}")
    age = st.number_input("Age (years)", min_value=0, max_value=110, value=35)
    childbearing = st.checkbox("Able to have children", value=False)
    CLIN = clinical.assess(stype, list(R["drugs"]), age=age, childbearing=childbearing)
    st.divider()
    st.caption("Drug strengths are relative to a nominal therapeutic dose (1.0). Lower = less drug needed to stop the "
               "seizure spreading in this brain. Differences below ~0.1 are ties.")

tab_type, tab_rank, tab_patients, tab_evidence, tab_checks, tab_about = st.tabs(
    ["Seizure type", "Drug ranking", "Real patients (IDEAS II)", "Clinical evidence", "Verification", "About & limitations"])

STATUS_ICON = {"first_line": "🟢", "later_line": "🔵", "restricted": "🟠", "not_listed": "⚪", "not_indicated": "⚫", "avoid": "🔴"}

# --------------------------------------------------------------------------------------------- seizure type gate
with tab_type:
    st.subheader(f"{CLIN['type_label']}")
    st.markdown(f"**Class:** {CLIN['class_label']}. {CLIN['why']}")
    if CLIN["simulate"]:
        st.success("Focal seizure: ARIADNE's spread simulation applies. Drugs below are ranked in the next tab.")
    else:
        st.warning("Not a focal seizure: ARIADNE does not simulate this. Guideline status only.")
    rows = [{"Drug": d, "Status": f"{STATUS_ICON[r['status']]} {r['label']}", "Rule": r["rec"],
             "Source": r["source"] or "", "Modelled in ARIADNE": "yes"} for d, r in
            sorted(CLIN["drugs"].items(), key=lambda x: clinical.STATUS_ORDER.index(x[1]["status"]))]
    rows += [{"Drug": g["drug"], "Status": f"{STATUS_ICON[g['status']]} {clinical.LABEL[g['status']]}", "Rule": g["rec"],
              "Source": "NG217", "Modelled in ARIADNE": "no"} for g in CLIN["guideline_only"]]
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    else:
        st.info("No guideline drug list for this class. Classify the seizure first (EEG, MRI, history).")
    if CLIN["patient"]["valproate_restricted"]:
        st.caption("Valproate is restricted because the patient is under 55 or able to have children (NICE NG217 Box 2).")
    with st.expander("Sources"):
        R_ = clinical.load_rules()
        for k, v in R_["sources"].items():
            st.markdown(f"- **{k}**: {v}")
        st.caption("Rules marked as 'check' in clinical_rules.yaml were read from a partly rendered page and should be re-verified.")

# --------------------------------------------------------------------------------------------- drug ranking
with tab_rank:
    for d_, rs in R.get("excluded", {}).items():
        st.error(f"**{d_} excluded** ({rs[0]['variant']}): {rs[0]['reason']}  \n_{rs[0]['source']}_")
    for d_, rs in R.get("flagged", {}).items():
        st.warning(f"**{d_} flagged** ({rs[0]['variant']}): {rs[0]['reason']}")
    if R.get("ieeg"):
        e = R["ieeg"]
        st.info("iEEG-suggested onset (anatomical labels mapped to the template): "
                + ", ".join(f"{m['label']} → {m['region']}" for m in e.get("mapped_regions", [])[:3])
                + ". The seizure was widespread, so this onset is a weak prior.")

    if not CLIN["simulate"]:
        st.warning(f"{CLIN['type_label']} is {CLIN['class_label'].lower()}. ARIADNE only simulates focal seizures, so there is "
                   "no ranking for this seizure type. See the Seizure type tab for guideline options.")
    else:
        c1, c2 = st.columns([1, 1.15], gap="large")
        with c1:
            st.subheader("Which drug stops the spread most easily?")
            th = {k: R["drugs"][k]["threshold"] for k in ranking}
            smax = 4.0
            fig, ax = plt.subplots(figsize=(5.6, 0.55 * len(ranking) + 0.8))
            vals = [(th[k] if th[k] is not None else smax) for k in ranking[::-1]]
            bars = ax.barh(ranking[::-1], vals, color=[ACCENT if th[k] is not None else "#c9c4dd" for k in ranking[::-1]])
            for b, k in zip(bars, ranking[::-1]):
                ax.text(b.get_width() + 0.05, b.get_y() + b.get_height() / 2,
                        ">4 (no stop)" if th[k] is None else f"{th[k]:.2f}", va="center", fontsize=9)
            ax.axvline(1.0, color="gray", ls="--", lw=1)
            ax.set_xlim(0, smax + 0.8); ax.set_xlabel("drug strength needed to stop spread (dashed = nominal dose)")
            ax.spines[["top", "right"]].set_visible(False)
            st.pyplot(fig)
            for i, k in enumerate(ranking, 1):
                d = R["drugs"][k]
                cs = CLIN["drugs"].get(k, {"status": "not_listed", "label": ""})
                st.markdown(f"**{i}. {k}** &nbsp; <small>{STATUS_ICON[cs['status']]} {cs['label']}</small>  \n<small>{', '.join(d['targets'])}</small>", unsafe_allow_html=True)

        with c2:
            st.subheader("Seizure spread in the simulated brain")
            state = st.radio("Show", ["baseline"] + ranking, horizontal=True,
                             format_func=lambda s: "Untreated" if s == "baseline" else s)
            mask = np.array(R["snapshot"][state]["recruited_mask"], dtype=bool)
            cpath = P("connectomes", R["connectome"], "centres_mm.csv")
            if os.path.exists(cpath):
                C = np.loadtxt(cpath, delimiter=",")
                fig, axs = plt.subplots(1, 3, figsize=(9, 3))
                for ax, (i, j, t) in zip(axs, [(0, 1, "axial"), (0, 2, "coronal"), (1, 2, "sagittal")]):
                    ax.scatter(C[~mask, i], C[~mask, j], c=QUIET, s=55)
                    ax.scatter(C[mask, i], C[mask, j], c=SEIZE, s=55)
                    o = R["onset_index"]; ax.scatter(C[o, i], C[o, j], marker="*", s=260, c=ONSET, edgecolors="k")
                    ax.set_title(t); ax.set_aspect("equal"); ax.axis("off")
                st.pyplot(fig)
                st.caption("🟢 quiet   🔴 seizing   ⭐ onset")
            else:
                seiz = [labels[i] for i in np.where(mask)[0]]
                st.markdown('<div class="note">No 3D region positions for this connectome, so regions are listed instead.</div>',
                            unsafe_allow_html=True)
                st.write(", ".join(seiz[:25]) + (" …" if len(seiz) > 25 else "") if seiz else "No regions recruited.")
            m1, m2 = st.columns(2)
            base = int(np.sum(R["snapshot"]["baseline"]["recruited_mask"]))
            m1.metric("Regions seizing (this view)", f"{int(mask.sum())} / {len(mask)}",
                      None if state == "baseline" else f"{int(mask.sum()) - base} vs untreated", delta_color="inverse")
            m2.metric("Untreated", f"{base} / {len(mask)}")
            if state != "baseline":
                st.caption("Snapshot at the nominal therapeutic dose (strength 1.0). Drugs that need more than 1.0 still let the seizure spread here.")

        if R.get("sensitivity"):
            st.subheader("Does the ranking hold if the assumptions change?")
            rows = [{"Variant": "Main analysis", "Ranking (best → worst)": " › ".join(ranking)}]
            rows += [{"Variant": k, "Ranking (best → worst)": " › ".join(v["ranking"])} for k, v in R["sensitivity"].items()]
            rob = P("results", "robustness.json")
            if os.path.exists(rob) and R["connectome"] == "stanford_hardi":
                rows += [{"Variant": x["variant"], "Ranking (best → worst)": " › ".join(x["ranking"])}
                         for x in json.load(open(rob)) if x.get("ranking")]
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")

        with st.expander("Drug parameters used (and where each number comes from)"):
            try:
                D = yaml.safe_load(open(P("drugs.yaml")))
                st.dataframe(pd.DataFrame([{
                    "Drug": k, "Free conc. (µM)": v.get("free_conc_uM"),
                    "Na⁺ channel block": f"{100 * v.get('nav_block_fraction', 0):.0f}%",
                    "GABA-A boost": f"{100 * v.get('gaba_fraction', 0):.0f}%",
                    "Excitability effect": v["x0_effect"], "Coupling effect": v["coupling_effect"],
                    "Sources": v.get("citation", "")} for k, v in D.items()]), hide_index=True, width="stretch")
                st.caption("Release-reduction sizes, valproate's GABA/T-type effect and the GABA weight are assumptions "
                           "(see drug_sources.yaml and VALUE_VERIFICATION.md).")
            except Exception as ex:  # noqa: BLE001
                st.write(f"Could not read drugs.yaml: {ex}")

# --------------------------------------------------------------------------------------------- IDEAS patients
with tab_patients:
    st.subheader("Real epilepsy patients: IDEAS II (UCL / Queen Square, OpenNeuro ds007401)")
    st.markdown("Pre-surgical structural connectomes (82 named regions), surgical resections and ILAE outcomes after "
                "epilepsy surgery. Used here for **virtual surgery**: remove the tissue that was really resected and "
                "measure how much the simulated network's ability to generate and spread seizures drops.")
    vs = P("results", "virtual_surgery.csv")
    if os.path.exists(vs):
        d = pd.read_csv(vs)
        d["Outcome"] = np.where(d.seizure_free, "Seizure-free", "Not seizure-free")
        show = d[["ID", "Outcome", "ILAE_Year1", "Op_Type", "capacity_pre", "capacity_post", "reduction"]].copy()
        show.columns = ["Patient", "Outcome", "ILAE yr 1", "Surgery", "Spread before", "Spread after", "Reduction"]
        sf, nsf = d.loc[d.seizure_free, "reduction"], d.loc[~d.seizure_free, "reduction"]
        auc = float(np.mean([(a > b) + 0.5 * (a == b) for a in sf for b in nsf])) if len(sf) and len(nsf) else float("nan")
        k1, k2, k3, k4 = st.columns(4)
        k1.metric("Patients", f"{len(d)}", f"{len(sf)} seizure-free", delta_color="off")
        k2.metric("Reduction: seizure-free", f"{sf.mean():.3f}")
        k3.metric("Reduction: not seizure-free", f"{nsf.mean():.3f}")
        k4.metric("AUC (0.5 = chance)", f"{auc:.2f}")
        cA, cB = st.columns([1.7, 1], gap="medium")
        with cA:
            st.dataframe(show.round(3), hide_index=True, width="stretch")
        with cB:
            fig, ax = plt.subplots(figsize=(4.4, 3.2))
            for x, (lab, col) in enumerate([("Seizure-free", QUIET), ("Not seizure-free", SEIZE)]):
                y = d.loc[d.Outcome == lab, "reduction"]
                ax.scatter(np.full(len(y), x) + np.random.default_rng(0).uniform(-.08, .08, len(y)), y, s=70, c=col)
                if len(y):
                    ax.hlines(y.mean(), x - .25, x + .25, color="k", lw=2)
            ax.set_xticks([0, 1], ["Seizure-free", "Not seizure-free"]); ax.set_xlim(-.6, 1.6)
            ax.set_ylabel("Reduction in spread capacity"); ax.spines[["top", "right"]].set_visible(False)
            st.pyplot(fig)
        with st.expander("Full report"):
            st.markdown(doc(P("results", "virtual_surgery_report.md")))
        st.markdown('<div class="note"><b>Pilot only.</b> Too few patients to conclude anything. Known issues being fixed: '
                    'a ceiling effect (untreated seizures reach 93–100% of the brain) and partial resections being under-counted.</div>',
                    unsafe_allow_html=True)
    else:
        st.info("Run `python virtual_surgery_study.py` to populate this tab.")

# --------------------------------------------------------------------------------------------- clinical evidence
with tab_evidence:
    st.subheader("How does ARIADNE's ranking compare with clinical evidence?")
    st.markdown("ARIADNE simulates **efficacy only** (stopping spread). Clinical rankings also reflect side effects, and "
                "most trials studied newly diagnosed patients rather than drug-resistant focal epilepsy.")
    ev = pd.DataFrame([
        ["Cenobamate ranks first", "Placebo-subtracted seizure freedom 20% vs lamotrigine 1.6%, pregabalin 1.5% (Löscher & Klein 2026); "
         "highest seizure-freedom RR in 75-trial NMA (Zhang 2024)", "Consistent"],
        ["…but real-world data", "Cenobamate not clearly better for seizure freedom in emulated trials (Xie 2025)", "Mixed"],
        ["Lamotrigine effective in drug-resistant focal epilepsy", "Add-on vs placebo OR 2.51, 11 RCTs (WHO mhGAP 2015)", "Consistent"],
        ["Lamotrigine vs carbamazepine / valproate", "Lamotrigine better on treatment failure, largely tolerability; no difference on remission "
         "(Cochrane, Nevitt 2022)", "Partly"],
        ["Lamotrigine vs pregabalin", "Tie on placebo-subtracted seizure freedom; pregabalin higher on NMA risk ratio (Zhang 2024)",
         "Depends on metric"],
        ["Lamotrigine first-line for focal epilepsy", "Superior to levetiracetam and zonisamide for 12-month remission (SANAD II, 2021)",
         "Consistent (LEV/ZNS not modelled)"],
    ], columns=["ARIADNE result / claim", "Clinical evidence", "Verdict"])
    st.dataframe(ev, hide_index=True, width="stretch")
    st.markdown("**Bottom line:** the clinical literature has no head-to-head ranking of these drugs in drug-resistant "
                "epilepsy, so ARIADNE's ranking is a **mechanistic hypothesis**, not a validated prediction.")
    with st.expander("Full benchmark notes"):
        st.markdown(doc(P("CLINICAL_BENCHMARKS.md")))
    with st.expander("How ARIADNE differs from EpiPick"):
        st.markdown("EpiPick (validated in 425 patients, Hadady 2022) tailors drug choice to the **clinical profile**: "
                    "seizure type, age, sex, comorbidities. ARIADNE explores tailoring it to the **brain network** and drug "
                    "**mechanism**. They are complementary: EpiPick narrows the safe options, ARIADNE could rank them by "
                    "predicted effect on a patient's seizure network.")

# --------------------------------------------------------------------------------------------- verification
with tab_checks:
    st.subheader("Can the simulator be trusted to behave correctly?")
    v1, v2 = st.columns(2, gap="large")
    with v1:
        st.markdown("#### Model checks"); st.markdown(doc(P("results", "verification_report.md"), "Run `python verify_model.py`."))
    with v2:
        st.markdown("#### Ranking vs model settings"); st.markdown(doc(P("results", "robustness_report.md"), "Run `python robustness.py`."))
    with st.expander("Every number checked against its source (VALUE_VERIFICATION.md)"):
        st.markdown(doc(P("VALUE_VERIFICATION.md")))

# --------------------------------------------------------------------------------------------- about
with tab_about:
    st.subheader("How it works")
    st.markdown("""
1. **Brain wiring.** A diffusion-MRI connectome: either built from a scan with tractography, or a ready-made patient connectome (IDEAS II).
2. **Seizure model.** Every region runs a reduced **Epileptor** (verified against The Virtual Brain); regions are coupled through the
   connectome so a seizure spreads along white-matter tracts.
3. **Calibration.** Each brain is set just above its own seizure-spread threshold.
4. **Drugs as mechanisms.** Each drug lowers excitability (sodium-channel block, including inactivated-state block; GABA-A boost) and/or
   weakens coupling (reduced transmitter release), sized from published IC50/EC50 values at real free blood levels.
5. **Ranking.** Each drug is dialled up until the spread stops; less drug needed = better *for this brain and onset site*.
6. **Optional inputs.** Scalp EEG or intracranial EEG to suggest the onset; genetic variants to exclude unsafe drugs (HLA-B\\*15:02, Dravet/SCN1A, POLG).
""")
    st.subheader("Limitations")
    st.markdown("""
- **Not validated for drug choice.** No public dataset links brain scans to which drug worked; clinical trials also lack head-to-head rankings.
- **Efficacy only.** Side effects and tolerability (which drive many real-world choices) are not modelled.
- **Some effects are assumptions:** release-reduction sizes, valproate's GABA/T-type effects, and the GABA weight (which decides cenobamate vs lamotrigine).
- **All-or-none spread** near threshold; focal seizures only; five drugs (levetiracetam and topiramate not modelled because their mechanisms are not quantified).
- **Virtual surgery** on real patients is a pilot so far.
""")
    st.caption("Data: Stanford HARDI (Rokem 2015), Sherbrooke 3-shell (Rokem 2017), IDEAS II (Taylor et al., Epilepsia 2026; OpenNeuro ds007401), "
               "Penn 3T-7T (Lucas 2026; Imtiaz 2026 Sci Data), CHB-MIT (PhysioNet). Research prototype; not for clinical use.")
