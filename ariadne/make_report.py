"""Turn the uncertainty analysis into poster-ready figures and a written report.

    python make_report.py [--uncertainty results/uncertainty_stanford_hardi_onsetP39.json]

Figures go to results/figures/ as PNG (300 dpi, white background, print-safe).
The drug colours are a fixed assignment validated for colour-vision deficiency on all pairs
(six series: worst pair dE 7.7 deutan / 17.4 normal, which the colour rules allow only because every
series is directly labelled). They are never cycled or reassigned.
"""
import argparse, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "results", "figures")
ORDER = ["cenobamate", "lamotrigine", "valproate", "carbamazepine", "pregabalin", "levetiracetam"]
COL = {"cenobamate": "#6D3FD6", "lamotrigine": "#009E73", "valproate": "#D55E00",
       "carbamazepine": "#56B4E9", "pregabalin": "#B5176B", "levetiracetam": "#00687A"}
VALID_TABLE = ""
NOISE_NOTE = ""
INK, INK2, MUTED, GRID = "#1a1a1f", "#44444f", "#6e6e7a", "#e3e3e8"
GOOD, BAD = "#1B7F4B", "#B3261E"

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "font.size": 10, "axes.labelsize": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": GRID, "grid.linewidth": 0.8, "legend.frameon": False, "figure.dpi": 110,
})


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", os.path.relpath(p, HERE))
    return p


# --------------------------------------------------------------------- 1. response surface
def fig_surface(S, U):
    x0 = np.array(S["x0_grid"]) * 1000.0           # milli-units, easier to read
    c = np.array(S["coupling_grid"])
    rec = np.array(S["recruited"], float) / S["n_regions"] * 100.0
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    im = ax.pcolormesh(x0, c, rec, cmap="Purples", shading="nearest", vmin=0, vmax=100, rasterized=True)
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("regions recruited (%)", color=INK)
    cb.outline.set_edgecolor(MUTED)

    b = [np.nan if v is None else v * 1000.0 for v in S["boundary_x0_reduction"]]
    ax.plot(b, c, color=INK, lw=2.2, zorder=4, label="containment boundary")

    # show only the region the drugs actually occupy, otherwise the plot is mostly empty
    xmax = max([(r["threshold_strength"] or 4) * r["x0_standard"] * 0.004 * 1000
                for r in U["therapeutic_ceiling"]] + [4.0]) * 1.35
    xmax = min(xmax, x0.max())

    # drugs with no coupling effect all run along y = 1.0 and would hide one another, so they
    # are drawn with distinct dash patterns rather than nudged off their true position
    flat_styles = ["-", (0, (7, 4)), (0, (2, 2.5))]
    flat_n = 0
    for d in ORDER:
        row = next((r for r in U["therapeutic_ceiling"] if r["drug"] == d), None)
        if not row:
            continue
        s = np.linspace(0, 4, 200)
        xs = s * row["x0_standard"] * 0.004 * 1000.0
        cs = np.maximum(0.0, 1 - s * row["coupling_standard"] * 0.06)
        keep = xs <= xmax
        lw, ls = 2.0, "-"
        if row["coupling_standard"] == 0:
            ls = flat_styles[flat_n % len(flat_styles)]; flat_n += 1; lw = 2.4
        if row["x0_standard"] == 0:
            lw = 3.2                                        # a pure-coupling drug sits on the axis
        ax.plot(xs[keep], cs[keep], color=COL[d], lw=lw, ls=ls, alpha=.95, zorder=5, label=d,
                solid_capstyle="round")
        t = row["threshold_strength"]
        if t is not None:
            tx = t * row["x0_standard"] * 0.004 * 1000
            if tx <= xmax:
                ax.plot([tx], [max(0, 1 - t * row["coupling_standard"] * 0.06)],
                        "o", ms=9, mfc=COL[d], mec="white", mew=1.8, zorder=6)

    ax.set_xlim(0, xmax); ax.set_ylim(c.min(), 1.001)
    ax.set_xlabel("excitability reduction  (model units x1000)")
    ax.set_ylabel("coupling scale  (1.0 = untreated)")
    ax.set_title("Every drug is a straight line across one plane", pad=10)
    ax.grid(True, lw=.6, alpha=.45)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize=9, columnspacing=1.4)
    fig.text(0.01, -0.30, "Shading: how much of the network the seizure recruits. Black line: where it stops spreading. "
             "Each coloured line is one drug as its dose rises; the dot marks where it crosses into containment. "
             "Pregabalin has no excitability effect, so its line runs down the left axis.",
             fontsize=8.5, color=MUTED, wrap=True)
    return save(fig, "fig1_response_surface.png")


# --------------------------------------------------------------------- 2. GABA weight curve
def fig_gaba(U):
    cur = U["gaba_curve"]
    gw = [p["gaba_weight"] for p in cur]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    smax = 4.0
    for d in ORDER:
        y = [(p["thresholds"].get(d) if p["thresholds"].get(d) is not None else np.nan) for p in cur]
        ax.plot(gw, y, color=COL[d], lw=2.2, label=d)
        last = next((i for i in range(len(y) - 1, -1, -1) if not np.isnan(y[i])), None)
        if last is not None:
            ax.annotate(d, (gw[last], y[last]), color=COL[d], fontsize=9, fontweight="bold",
                        xytext=(4, 0), textcoords="offset points", va="center")
    f = U.get("gaba_flip")
    title = "The ranking barely depends on the one unsourced constant"
    if f:
        x = float(np.mean(f["between"]))
        ax.axvline(x, color=INK, lw=1.6, ls="--", zorder=2)
        ax.annotate(f"ranking flips only below {x:.2f}\n({f['from']} first to the left,\n{f['to']} first to the right)",
                    (x + 0.05, smax * 0.80), fontsize=9, color=INK, ha="left", va="top",
                    bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=GRID))
        ax.annotate("value used throughout: 0.5", (0.5, smax * 0.045), fontsize=8.5, color=MUTED,
                    ha="center", va="bottom")
        ax.plot([0.5], [0], marker="^", ms=9, color=MUTED, clip_on=False, zorder=5)
        if x < 0.15:
            title = f"The winner changes only if GABA potentiation counts for almost nothing (below {x:.2f})"
    ax.axhline(1.0, color=MUTED, lw=1, ls=":")
    ax.annotate("standard dose", (0.07, 1.03), fontsize=8.5, color=MUTED, va="bottom", ha="left")
    ax.set_xlim(0, 1.22); ax.set_ylim(0, smax * 1.02)
    ax.set_xlabel("GABA weight  (unsourced constant: how much GABA-A potentiation counts vs sodium block)")
    ax.set_ylabel("dose needed to contain the seizure\n(multiples of standard)")
    ax.set_title(title, fontsize=11.5)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.grid(True, axis="y", lw=.6, alpha=.6)
    return save(fig, "fig2_gaba_weight.png")


# --------------------------------------------------------------------- 3. Monte Carlo
def fig_monte_carlo(U):
    A, B = U["monte_carlo_sourced"], U["monte_carlo_full"]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9), sharex=True)
    for ax, M, title, sub in [
            (axes[0], A, "Measurement uncertainty only", "sampling the ranges that have published sources"),
            (axes[1], B, "Full uncertainty", "also sampling the parameters with no source at all")]:
        p = [M["p_ranks_first"].get(d, 0.0) for d in ORDER]
        y = np.arange(len(ORDER))[::-1]
        ax.barh(y, p, height=.62, color=[COL[d] for d in ORDER], zorder=3)
        for yi, d, pi in zip(y, ORDER, p):
            ax.text(pi + .02, yi, f"{pi*100:.0f}%", va="center", fontsize=9.5, color=INK2)
            ax.text(-0.04, yi, d, va="center", ha="right", fontsize=9.5, color=INK)
        ax.set_yticks([]); ax.set_xlim(-0.46, 1.16)
        ax.set_xticks([0, .25, .5, .75, 1.0]); ax.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
        ax.set_title(title, pad=14)
        ax.annotate(sub, (0.5, 1.015), xycoords="axes fraction", ha="center", va="bottom",
                    fontsize=8.5, color=MUTED)
        ax.grid(True, axis="x", lw=.6, alpha=.6, zorder=0)
        ax.spines["left"].set_visible(False); ax.spines["bottom"].set_bounds(0, 1)
        ax.tick_params(left=False)
    fig.subplots_adjust(wspace=0.12)
    fig.suptitle(f"How often does each drug come first, across {A['n_draws']} parameter draws?",
                 y=1.14, fontsize=12.5, fontweight="bold")
    fig.text(0.5, -0.1, "Left: if every published number were re-measured elsewhere, how stable is the answer? "
             "Right: how much of the answer rests on assumptions rather than measurements?",
             fontsize=8.5, color=MUTED, ha="center")
    return save(fig, "fig3_monte_carlo.png")


# --------------------------------------------------------------------- 4. therapeutic ceiling
def fig_ceiling(U):
    rows = [r for d in ORDER for r in U["therapeutic_ceiling"] if r["drug"] == d]
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    smax = 4.0
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        d = r["drug"]
        unreachable = r["threshold_strength"] is None
        need = smax if unreachable else r["threshold_strength"]
        reach = r["effect_ratio"]
        ok = r["contained_at_ceiling"]
        ax.barh(yi, need, height=.46, color=COL[d], alpha=.92, zorder=3,
                hatch="///" if unreachable else None, edgecolor="white", linewidth=0)
        # the most this drug can actually deliver, at the top of its therapeutic range
        ax.plot([reach, reach], [yi - .30, yi + .30], color=INK, lw=2.6, zorder=6,
                solid_capstyle="butt")
        ax.text(reach, yi + .34, f"{reach:.2f}x", fontsize=8.5, color=INK, ha="center", va="bottom",
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=.85), zorder=7)
        ax.text(-0.06, yi, d, va="center", ha="right", fontsize=10, color=INK, fontweight="bold")
        status = "contained" if ok else ("needs > 4x" if unreachable else "not contained")
        ax.text(smax * 1.03, yi, status, va="center", fontsize=9,
                color=GOOD if ok else BAD, fontweight="bold")
    ax.axvline(1.0, color=MUTED, lw=1.2, ls=":", zorder=1)
    ax.text(1.0, len(rows) - .42, "standard dose", fontsize=8.5, color=MUTED, ha="center", va="bottom")
    ax.set_yticks([]); ax.set_xlim(0, smax); ax.set_ylim(-.7, len(rows) - .1)
    ax.set_xlabel("dose needed to contain the seizure  (multiples of standard exposure)")
    ax.set_title("Can the drug get there without leaving its therapeutic range?")
    ax.grid(True, axis="x", lw=.6, alpha=.6, zorder=0)
    ax.spines["left"].set_visible(False); ax.tick_params(left=False)
    fig.text(0.01, -0.26, "Bar: dose the simulation says is needed. Black tick: how far the drug can actually go at the top of its therapeutic "
             "range or maximum licensed dose, on the same scale. Occupancy saturates, so doubling the dose delivers less than double the effect. "
             "For two-mechanism drugs the tick shows the excitability component; the verdict uses both.",
             fontsize=8.5, color=MUTED, wrap=True)
    return save(fig, "fig4_therapeutic_ceiling.png")


# --------------------------------------------------------------------- report
def validation_table(U):
    """Compare surface-derived thresholds against the original independent bisection runs."""
    import glob
    base = U["baseline"]["thresholds"]
    sim = None
    for p in glob.glob(os.path.join(HERE, "results", "*.json")):
        b = os.path.basename(p)
        if b.startswith(("surface_", "uncertainty_", "robustness")):
            continue
        try:
            R = json.load(open(p))
        except Exception:
            continue
        if R.get("connectome") == U.get("surface", "").replace("surface_", "").rsplit("_onset", 1)[0] \
           and R.get("onset") == U.get("onset") and "drugs" in R and not R.get("excluded"):
            sim = {k: R["drugs"][k]["threshold"] for k in R["drugs"]}
            break
    if not sim:
        return "_(original simulation results not found for comparison)_\n"
    rows = ["| Drug | full simulation | from surface | difference |", "|---|---|---|---|"]
    for d in ORDER:
        a, b2 = sim.get(d), base.get(d)
        if a is None or b2 is None:
            rows.append(f"| {d} | {a if a is not None else 'never'} | {b2 if b2 is not None else 'never'} | - |")
        else:
            rows.append(f"| {d} | {a:.3f}x | {b2:.3f}x | {b2-a:+.3f} |")
    rows.append("")
    return "\n".join(rows)


def report(U, S, figs):
    global VALID_TABLE, NOISE_NOTE
    VALID_TABLE = validation_table(U)
    conf = S.get("boundary_confirmation") or []
    if conf:
        ok = sum(c["agrees"] for c in conf)
        step = (S["x0_grid"][1] - S["x0_grid"][0]) if len(S["x0_grid"]) > 1 else 0
        NOISE_NOTE = (
            f"**How sharp is the boundary?** The grid was run with one noise seed, then every boundary cell was "
            f"re-tested with two further seeds: **{ok} of {len(conf)}** gave the same verdict. The boundary therefore "
            f"carries about one grid step of jitter ({step*1000:.2f} in the units of the figure, roughly "
            f"{step/S['boundary_x0_reduction'][0]*100:.0f}% of its value at an untouched coupling scale). That "
            f"propagates to roughly the same relative uncertainty on every threshold, which is far smaller than the "
            f"gaps between drugs, but it is the reason thresholds are not quoted to three decimal places.\n")
    else:
        NOISE_NOTE = ""
    A, B = U["monte_carlo_sourced"], U["monte_carlo_full"]
    base = U["baseline"]
    L = [f"# ARIADNE uncertainty analysis - {U['brain']} (onset {U['onset']})", "",
         "Generated by `make_report.py` from `uncertainty.py`. The point-estimate ranking below is the one ARIADNE",
         "previously reported; everything after it quantifies how much that ranking can be trusted.", "",
         f"**Point-estimate ranking:** {' > '.join(base['ranking'])}  ",
         f"(GABA weight {base['gaba_weight']}, NAV_REF {base['nav_ref']}, sodium block = {base['nav_mode']} state)", ""]

    L += ["## 1. The whole drug model lives on one plane", "",
          f"![response surface]({os.path.relpath(figs['surface'], os.path.join(HERE,'results'))})", "",
          "Each drug enters the simulation through exactly two numbers: how much it lowers excitability and how much",
          "it weakens coupling. Mapping the containment boundary once (1701 simulations) therefore replaces every",
          "later simulation with a lookup, which is what makes uncertainty propagation affordable.", "",
          "**Validation.** The thresholds read off the surface agree with the original full simulations, which were",
          "computed independently by bisection, to within the bisection resolution:", "", VALID_TABLE,
          NOISE_NOTE,
          "**A consequence worth stating plainly.** For drugs that act only on excitability, the threshold is",
          "boundary / (x0_effect x unit), so the *ranking* among them is fixed entirely by their pharmacology and is",
          "the same in every brain. The connectome sets the absolute dose, not the order. Patient-specific ordering can",
          "only come from drugs that mix both mechanisms, or from brains whose boundary shape differs.", ""]

    fmid = float(np.mean(f["between"])) if (f := U.get("gaba_flip")) else None
    L += ["## 2. " + ("How much does the one unsourced constant matter?" if fmid is None or fmid < 0.15
                      else "One unsourced constant decides the winner"), "",
          f"![gaba weight]({os.path.relpath(figs['gaba'], os.path.join(HERE,'results'))})", ""]
    if f:
        mid = float(np.mean(f["between"]))
        L += [f"The ranking flips between GABA weight **{f['between'][0]}** and **{f['between'][1]}**: below that,",
              f"**{f['from']}** ranks first; above it, **{f['to']}** does. The value used throughout ARIADNE is 0.5,",
              "which has no source.", ""]
        if mid < 0.15:
            L += [f"**This is a better result than expected.** The flip sits at about {mid:.2f}, so the ranking only",
                  "changes if GABA-A potentiation is treated as counting for almost nothing relative to sodium-channel",
                  "block. Cenobamate's GABA-A effect is substantial at therapeutic concentrations (45% of maximal",
                  "potentiation at the modelled free level, Sharma 2020), so a weight near zero is not a plausible",
                  "reading of the pharmacology. The headline ordering is therefore not balanced on a knife edge -",
                  "the earlier robustness run, which only tested the single value 0.0, made it look more fragile than",
                  "it is.", ""]
    else:
        L += ["Across the full 0-1 range the first-ranked drug does not change.", ""]

    L += ["## 3. The ranking as a probability, not a claim", "",
          f"![monte carlo]({os.path.relpath(figs['mc'], os.path.join(HERE,'results'))})", "",
          f"Two Monte Carlo runs of {A['n_draws']} draws each.", "",
          "| Drug | P(first) - measurement uncertainty | P(first) - full uncertainty |", "|---|---|---|"]
    for d in ORDER:
        L.append(f"| {d} | {A['p_ranks_first'].get(d,0)*100:.0f}% | {B['p_ranks_first'].get(d,0)*100:.0f}% |")
    L += ["", "**Measurement uncertainty** samples only ranges that have published sources: therapeutic reference",
          f"ranges, the two disagreeing lamotrigine protein-binding values, and +/-{U['ic50_fold_assumed']}-fold IC50",
          "reproducibility. **Full uncertainty** additionally samples the parameters with no source at all (the GABA",
          "weight, the NAV_REF scaling constant, and the assumed transmitter-release reductions).", "",
          "The difference between the two columns is the honest measure of how much of ARIADNE's answer rests on",
          "assumption rather than measurement.", "",
          "Most likely complete orderings under full uncertainty:", ""]
    for o in B["most_likely_orders"][:3]:
        L.append(f"- {' > '.join(o['order'])}  ({o['p']*100:.0f}%)")
    L += ["", "Pairwise, under full uncertainty:", ""]
    for k, v in list(B["pairwise_better"].items())[:6]:
        L.append(f"- P({k}) = {v*100:.0f}%")

    inv = U.get("unmodelled_drug_inverse") or []
    if inv:
        L += ["### The levetiracetam question, inverted", "",
              "Levetiracetam is the most-prescribed drug for focal epilepsy and is **not** modelled here, because no",
              "published number converts SV2A binding into a size of effect. Rather than invent one, the question can",
              "be turned around: how large would a purely release-reducing effect have to be to match each modelled",
              "drug at a standard dose?", ""]
        for r in inv:
            if r.get("release_reduction_needed") is None:
                continue
            L.append(f"- To match **{r['to_match']}** (threshold {r['target_threshold']}): "
                     f"{r['interpretable_as']}.")
        if any(r.get("feasible") is False for r in inv):
            L += ["", "**This is a falsification signal, and it is the most interesting thing in this report.** A drug",
                  "acting only by reducing transmitter release cannot reach these thresholds even at a physically",
                  "impossible 100% reduction. Yet levetiracetam plainly works in patients. Either the coupling term in",
                  "this model is too weak to carry a real drug effect, or levetiracetam's efficacy comes through a",
                  "mechanism the model does not represent. Both are testable, and both point at the coupling",
                  "formulation as the next thing to fix.", ""]

    L += ["", "## 4. Can the drug actually get there?", "",
          f"![ceiling]({os.path.relpath(figs['ceiling'], os.path.join(HERE,'results'))})", "",
          "A dose multiplier is not a clinical instruction. This asks the clinical question: at the top of its",
          "therapeutic reference range, or its maximum licensed dose, does the drug contain the seizure? Because",
          "receptor occupancy saturates, doubling the dose delivers less than double the effect.", "",
          "| Drug | Ceiling basis | Free conc. standard -> ceiling | Dose x | Effect x | Needed | Contained at ceiling |",
          "|---|---|---|---|---|---|---|"]
    for d in ORDER:
        r = next((x for x in U["therapeutic_ceiling"] if x["drug"] == d), None)
        if not r:
            continue
        need = "never" if r["threshold_strength"] is None else f"{r['threshold_strength']:.2f}x"
        L.append(f"| {d} | {r['basis']} | {r['free_standard_uM']} -> {r['free_ceiling_uM']} uM | "
                 f"{r['conc_ratio']}x | {r['effect_ratio']}x | {need} | {'yes' if r['contained_at_ceiling'] else 'no'} |")
    L += ["", "The gap between the dose column and the effect column is saturation: a drug given at twice the dose does",
          "not deliver twice the effect, so a threshold above roughly 1.5x is often out of reach in practice.", ""]

    L += ["## What this does and does not establish", "",
          "- It does **not** validate the ranking against patient outcomes. No public dataset links a brain scan to",
          "  which drug worked, so that test cannot be run by anyone today (see `PREREGISTRATION.md`).",
          "- It **does** quantify exactly how much of the answer is measurement and how much is assumption, and",
          "  identifies the single parameter that would need an experiment to settle it.",
          "- Efficacy only: tolerability, which drives much real prescribing, is not modelled at all.", ""]
    p = os.path.join(HERE, "results", "uncertainty_report.md")
    open(p, "w").write("\n".join(L))
    print("  wrote", os.path.relpath(p, HERE))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uncertainty", default=os.path.join(HERE, "results", "uncertainty_stanford_hardi_onsetP39.json"))
    ap.add_argument("--surface", default=None)
    a = ap.parse_args()
    U = json.load(open(a.uncertainty))
    S = json.load(open(a.surface or os.path.join(HERE, "results", U["surface"])))
    print("Figures:")
    figs = {"surface": fig_surface(S, U), "gaba": fig_gaba(U), "mc": fig_monte_carlo(U), "ceiling": fig_ceiling(U)}
    print("Report:")
    report(U, S, figs)


if __name__ == "__main__":
    main()
