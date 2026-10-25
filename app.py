"""ARIADNE prototype app. Research prototype - not for clinical use.
Run:  python -m streamlit run app.py --theme.base dark --theme.primaryColor "#22d3ee"
"""
import numpy as np, warnings
import streamlit as st
import streamlit.components.v1 as components
warnings.filterwarnings("ignore")
from epileptor_fast import simulate
from connectome import load_template
from groups import coarse, GROUP_NAMES
from run_real import drug_table

st.set_page_config(page_title="ARIADNE", layout="centered")
st.markdown("""
<style>
.stApp{background:radial-gradient(1200px 600px at 50% -10%,#1b2566 0%,#0a0f2a 45%,#060914 100%);}
.block-container{padding-top:4rem;max-width:980px}
h1.ari{font-size:3rem;margin:0;letter-spacing:.25em;background:linear-gradient(90deg,#22d3ee,#8b5cf6,#ff3b6b);
 -webkit-background-clip:text;-webkit-text-fill-color:transparent;font-weight:800}
p.sub{color:#9fb0ff;margin:.2rem 0 1rem 0}
div[role="radiogroup"]{gap:.4rem}
div[role="radiogroup"] label{background:#10183a;border:1px solid #2b3a8b;border-radius:999px;padding:.25rem .8rem}
div[role="radiogroup"] label:has(input:checked){background:linear-gradient(90deg,#0e7490,#6d28d9);border-color:#22d3ee}
</style>
<h1 class="ari">ARIADNE</h1>
<p class="sub">Watch how an antiseizure drug changes the way a seizure spreads through the brain.</p>
""", unsafe_allow_html=True)

@st.cache_data(show_spinner="Simulating the brain...")
def simulate_panel(f, SX, SK, K0):
    C, labels = load_template(); M, g = coarse(C, labels); M = np.minimum(M / M[0].max(), 1.0)
    x0 = -2.4 + 0.3 * M[0]; x0[0] = -1.6
    res = {}
    for n, (ox, oK, basis) in drug_table(f).items():
        xd = x0 - SX * ox; xd[0] = x0[0]
        o = simulate(M, xd, K=K0 * (1 - SK * oK), n_trials=4, seed=5)
        res[n] = dict(hit=(o > 0).any(axis=0).mean(0) >= 0.5, burden=float((o[:, :, 1:] > 0).mean()), basis=basis)
    return res

# node position, label position, label side  (viewBox 520 x 420, brain centred on x = 260)
NODES = {"Hippocampus": ((310, 255), (440, 255), "s"), "Amygdala": ((322, 215), (440, 215), "s"),
         "L temporal": ((135, 200), (95, 190), "e"), "R temporal": ((385, 200), (440, 175), "s"),
         "L frontal": ((180, 100), (95, 100), "e"), "R frontal": ((340, 100), (440, 100), "s"),
         "L parietal": ((180, 290), (95, 290), "e"), "R parietal": ((340, 300), (440, 300), "s"),
         "Occipital": ((260, 355), None, "m"), "Insula/cingulate": ((260, 175), None, "m")}
CORAL, MINT, SLATE = "#ff3b6b", "#2dffc4", "#3b4a7a"
SQUIGGLES = ["M200,70 C220,90 190,110 215,130 S200,160 225,175", "M150,130 C175,140 160,165 185,180",
             "M160,225 C185,215 175,250 203,245", "M155,275 C180,267 175,295 205,293",
             "M215,55 C235,75 225,100 245,115", "M225,145 C240,160 230,185 247,200",
             "M185,335 C205,325 210,350 230,345", "M125,190 C140,200 135,220 150,230"]

def brain_svg(hit, base_hit, tag):
    s = ['<svg viewBox="0 0 520 420" xmlns="http://www.w3.org/2000/svg" style="width:100%;display:block">',
         '<defs><radialGradient id="bg" cx="50%" cy="40%" r="65%"><stop offset="0" stop-color="#3b3f9a"/><stop offset="1" stop-color="#14183f"/></radialGradient>'
         '<filter id="glow" x="-80%" y="-80%" width="260%" height="260%"><feGaussianBlur stdDeviation="5" result="b"/><feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge></filter></defs>',
         '<path d="M260,28 C350,28 405,95 410,190 C415,290 370,380 260,392 C150,380 105,290 110,190 C115,95 170,28 260,28 Z" fill="url(#bg)" stroke="#22d3ee" stroke-opacity=".7" stroke-width="2.5" filter="url(#glow)"/>']
    for p in SQUIGGLES:
        s.append(f'<path d="{p}" fill="none" stroke="#9fb0ff" stroke-opacity=".22" stroke-width="3" stroke-linecap="round"/>')
        s.append(f'<path d="{p}" transform="translate(520,0) scale(-1,1)" fill="none" stroke="#9fb0ff" stroke-opacity=".22" stroke-width="3" stroke-linecap="round"/>')
    s.append('<path d="M260,32 C255,120 268,200 258,300 C254,340 262,370 260,388" fill="none" stroke="#060914" stroke-opacity=".8" stroke-width="3"/>')
    s.append('<text x="30" y="30" fill="#7c8bd6" font-size="13" font-family="Segoe UI,Arial">L</text><text x="480" y="30" fill="#7c8bd6" font-size="13" font-family="Segoe UI,Arial">R</text>')
    for i, nme in enumerate(GROUP_NAMES):
        (x, y), lab, side = NODES[nme]
        if i == 0 or hit[i]: c, state = CORAL, "seizing"
        elif base_hit[i]: c, state = MINT, "protected"
        else: c, state = SLATE, "quiet"
        if lab:
            s.append(f'<line x1="{x}" y1="{y}" x2="{lab[0] + (6 if side == "e" else -6)}" y2="{lab[1]}" stroke="#7c8bd6" stroke-opacity=".35" stroke-dasharray="3 3"/>')
        if state == "seizing":
            s.append(f'<circle cx="{x}" cy="{y}" r="14" fill="none" stroke="{c}" stroke-width="2"><animate attributeName="r" values="14;30;14" dur="1.8s" repeatCount="indefinite"/><animate attributeName="opacity" values=".9;0;.9" dur="1.8s" repeatCount="indefinite"/></circle>')
        ring = ' stroke="#fff" stroke-width="2.5"' if i == 0 else ""
        s.append(f'<circle cx="{x}" cy="{y}" r="14" fill="{c}"{ring} filter="url(#glow)"/>')
        if state == "protected":
            s.append(f'<path d="M{x-6},{y} l4,5 l8,-9" fill="none" stroke="#04281f" stroke-width="3" stroke-linecap="round" stroke-linejoin="round"/>')
        elif state == "seizing":
            s.append(f'<text x="{x}" y="{y+5}" text-anchor="middle" font-size="14" font-weight="bold" fill="#fff" font-family="Segoe UI,Arial">!</text>')
        col = "#ffffff" if state != "quiet" else "#8c9ad1"
        if lab:
            anchor = "end" if side == "e" else "start"
            s.append(f'<text x="{lab[0]}" y="{lab[1]+4}" text-anchor="{anchor}" font-size="12" fill="{col}" font-family="Segoe UI,Arial">{nme}</text>')
            s.append(f'<text x="{lab[0]}" y="{lab[1]+17}" text-anchor="{anchor}" font-size="9.5" fill="{c}" font-family="Segoe UI,Arial">{state}</text>')
        else:
            s.append(f'<text x="{x}" y="{y+32}" text-anchor="middle" font-size="12" fill="{col}" font-family="Segoe UI,Arial">{nme}</text>')
            s.append(f'<text x="{x}" y="{y+45}" text-anchor="middle" font-size="9.5" fill="{c}" font-family="Segoe UI,Arial">{state}</text>')
    s.append("</svg>")
    return "".join(s)

def card(title, tag, svg, accent):
    return (f'<div style="flex:1 1 400px;background:linear-gradient(160deg,#101846,#0a0f2a);border:1px solid {accent}55;'
            f'border-radius:18px;padding:14px 16px;box-shadow:0 0 30px {accent}22">'
            f'<div style="font:700 12px Segoe UI,Arial;letter-spacing:.2em;color:{accent}">{tag}</div>'
            f'<div style="font:700 20px Segoe UI,Arial;color:#fff;margin:2px 0 6px">{title}</div>{svg}</div>')

with st.expander("Advanced settings (leave alone for the demo)"):
    K0 = st.select_slider("Untreated spread strength", [0.5, 0.7, 1.0], value=0.7)
    f = st.select_slider("Sodium-blocker strength factor (unmeasured)", [1, 10, 30, 100], value=30)
    SX = st.slider("Excitability scale", 0.05, 0.20, 0.10, 0.01)
    SK = st.slider("Coupling scale", 0.3, 0.9, 0.6, 0.05)

res = simulate_panel(f, SX, SK, K0)
base = res["No drug"]
estimate = {n for n, r in res.items() if r["basis"] in ("PLACEHOLDER", "assumed zero")}
drug = st.radio("Choose a drug", [n for n in res if n != "No drug"], horizontal=True,
                format_func=lambda n: n + (" *" if n in estimate else ""), label_visibility="collapsed")
r = res[drug]
base_n = int(base["hit"][1:].sum()); n_hit = int(r["hit"][1:].sum()); protected = base_n - n_hit
if n_hit == 0: label, col = "SEIZURE SPREAD STOPPED", MINT
elif n_hit <= base_n / 2: label, col = "SEIZURE SPREAD REDUCED", "#ffb020"
else: label, col = "LITTLE OR NO EFFECT", CORAL
cut = max(0, 100 * (1 - r["burden"] / base["burden"]))
html = (f'<div style="display:flex;gap:18px;flex-wrap:wrap">'
        + card("Untreated seizure", "BEFORE", brain_svg(base["hit"], base["hit"], "before"), CORAL)
        + card(f"With {drug}{' *' if drug in estimate else ''}", "AFTER", brain_svg(r["hit"], base["hit"], "after"), "#22d3ee") + '</div>'
        f'<div style="margin-top:18px;background:{col}22;border:1.5px solid {col};border-radius:16px;padding:16px;text-align:center;'
        f'box-shadow:0 0 28px {col}33"><div style="font:800 24px Segoe UI,Arial;letter-spacing:.12em;color:{col}">{label}</div>'
        f'<div style="font:16px Segoe UI,Arial;color:#dfe5ff;margin-top:6px">{drug} protects <b>{protected} of {base_n}</b> regions '
        f'the untreated seizure reaches &nbsp;·&nbsp; seizure activity down <b>{cut:.0f}%</b></div></div>')
page = f'<body style="margin:0;background:transparent">{html}</body>'
if hasattr(st, 'iframe'):
    st.iframe(page, height=560)
else:
    components.html(page, height=560)

with st.expander("Compare all drugs"):
    order = sorted([n for n in res if n != "No drug"], key=lambda n: res[n]["burden"])
    st.bar_chart({n + (" *" if n in estimate else ""): 100 * (1 - res[n]["burden"] / base["burden"]) for n in order}, horizontal=True)
    st.caption("Bars show how much seizure activity each drug removes (%).")
st.caption("\\* estimate: no measured binding data yet. 'Stopped' means the model predicts suppression, not a clinical result. "
           "Template brain, not a patient. RESEARCH PROTOTYPE - NOT FOR CLINICAL USE.")
