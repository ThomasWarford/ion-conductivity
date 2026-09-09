"""Plots for the Onsager analysis (onsager.py output).

  python onsager_plots.py
  python onsager_plots.py --charge-balanced-only   # drop non-stoichiometric cells

Writes results/figs/onsager_curves_1500K.png   -- MSD_self, MSD_coll, C_dist vs lag
       results/figs/onsager_f_running.png       -- running f(t) = C_dist/MSD_self, per T
       results/figs/onsager_f_vs_T.png          -- f +/- block error vs T
and prints results/onsager_summary.md
"""
from __future__ import annotations

import gzip
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from charge_balance import classify

HERE = os.path.dirname(__file__)
RES = os.path.join(HERE, "results")
SPECIES = "Li"                      # set by --species
CB_ONLY = False                     # set by --charge-balanced-only
FIG = os.path.join(RES, "figs", "li")
TEMPS = [1000, 1500, 2000, 2500]
TCOL = dict(zip(TEMPS, plt.cm.viridis(np.linspace(0, .82, 4))))


def load():
    df = pd.read_csv(os.path.join(RES, "onsager.csv"))
    df = df[df.error.fillna("") == ""].copy()
    if "species" in df:
        df = df[df.species == SPECIES].copy()
    if CB_ONLY:
        # Most NCSD cells are packings at compositions that are not compounds
        # (Br-Li is Li3Br, not LiBr); their "ionic" transport is metallic
        # self-diffusion in a Li-rich melt.  See charge_balance.py.
        keep = df.formula.map(lambda f: classify(f).status == "balanced")
        dropped = sorted(df.loc[~keep, "system"].unique())
        df = df[keep].copy()
        if dropped:
            print(f"charge-balance filter dropped {len(dropped)} systems: "
                  + ", ".join(dropped))
    with gzip.open(os.path.join(RES, "onsager_curves.json.gz")) as fh:
        curves = json.load(fh)
    order = (df[df.temperature == 1500].sort_values("D_self_cm2s", ascending=False)
             .system.tolist())
    order += [s for s in df.system.unique() if s not in order]
    return df, curves, order


def grid(n):
    c = math.ceil(math.sqrt(n))
    r = math.ceil(n / c)
    return r, c


# one colour per curve; the in-panel number is drawn in a darker shade of the
# SAME hue, so each number is unambiguously tied to the curve it comes from.
CURVE = {"self": "#2a78d6", "coll": "#eb6834", "dist": "#1baf7a"}
LABEL = {"self": "#1a5fa8", "coll": "#c04a18", "dist": "#0f7d55"}


def _sci(v):
    """1.7e-4 style, compact enough for an in-panel label."""
    if v is None or not np.isfinite(v) or v == 0:
        return "n/a"
    e = int(np.floor(np.log10(abs(v))))
    return f"{v / 10 ** e:.2f}e{e}"


def fig_curves(df, curves, order, T=1500):
    n = len(order)
    r, c = grid(n)
    fig, axs = plt.subplots(r, c, figsize=(c * 3.0, r * 2.6), squeeze=False)
    for ax in axs.flat:
        ax.axis("off")
    for i, sys in enumerate(order):
        ax = axs.flat[i]; ax.axis("on")
        key = f"{sys}|{T}|{SPECIES}"
        row = df[(df.system == sys) & (df.temperature == T)]
        if key not in curves or row.empty:
            ax.set_title(sys, fontsize=8); continue
        cd = curves[key]; row = row.iloc[0]
        t = np.array(cd["t"])
        ax.plot(t, np.array(cd["MSD_self"]) / 1e3, lw=1.4, color=CURVE["self"],
                label=r"MSD$_{self}$ = $\langle\sum_a|\Delta r_a|^2\rangle$   $\rightarrow$  $D_{self}$")
        ax.plot(t, np.array(cd["MSD_coll"]) / 1e3, lw=1.4, color=CURVE["coll"],
                label=r"MSD$_{coll}$ = $\langle|\sum_a\Delta r_a|^2\rangle$   $\rightarrow$  $D_{coll}$")
        ax.plot(t, np.array(cd["C_dist"]) / 1e3, lw=1.4, color=CURVE["dist"],
                label=r"C$_{dist}$ = MSD$_{coll}$ $-$ MSD$_{self}$   $\rightarrow$  $f$")
        ax.axhline(0, color="#999", lw=.6)
        lo, hi = eval(row["window"]) if isinstance(row["window"], str) else row["window"]
        ax.axvspan(lo * t[-1], hi * t[-1], color="#000", alpha=.05)
        ax.set_title(f"{sys}", fontsize=8.5, pad=2)
        # --- annotate the panel with the transport numbers ---
        d, ds = row["D_self_cm2s"], row.get("D_self_cm2s_stat", np.nan)
        rel = 100 * ds / abs(d) if d else np.nan
        ax.text(.03, .97, f"$D_{{self}}$ = {_sci(d)} cm$^2$/s  (±{rel:.0f}%)",
                transform=ax.transAxes, va="top", fontsize=6.6, color=LABEL["self"])
        ax.text(.03, .885, f"$D_{{coll}}$ = {_sci(row['D_coll_cm2s'])} cm$^2$/s",
                transform=ax.transAxes, va="top", fontsize=6.6, color=LABEL["coll"])
        fe = row.get("f_err", np.nan)
        ax.text(.03, .80, f"$f$ = {row['f']:+.2f} ± {fe:.2f}",
                transform=ax.transAxes, va="top", fontsize=6.6, color=LABEL["dist"])
        ax.tick_params(labelsize=7)
        lo_y, hi_y = ax.get_ylim()
        ax.set_ylim(lo_y, hi_y + 0.42 * (hi_y - lo_y))
        if i % c == 0:
            ax.set_ylabel(r"$10^3\ \AA^2$", fontsize=7)
        if i // c == r - 1:
            ax.set_xlabel("lag time (ps)", fontsize=7)
    h, lab = axs.flat[0].get_legend_handles_labels()
    fig.suptitle(f"{SPECIES} displacement correlation functions at {T} K   —   "
                 r"$D$ = slope over the shaded window $/\,6N$", fontsize=12, y=0.998)
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, 0.978),
               ncol=3, fontsize=9, frameon=False, handlelength=1.8, columnspacing=3.0)
    fig.tight_layout(rect=[0, 0, 1, 0.945])
    fig.savefig(f"{FIG}/curves_{T}K.png", dpi=130)
    plt.close(fig)


def fig_f_running(df, curves, order):
    n = len(order)
    r, c = grid(n)
    fig, axs = plt.subplots(r, c, figsize=(c * 2.7, r * 2.3), squeeze=False)
    for ax in axs.flat:
        ax.axis("off")
    for i, sys in enumerate(order):
        ax = axs.flat[i]; ax.axis("on")
        ax.axhspan(-1, 1, color="#1baf7a", alpha=.06)
        ax.axhline(0, color="#999", lw=.6)
        for T in TEMPS:
            key = f"{sys}|{T}|{SPECIES}"
            if key not in curves:
                continue
            cd = curves[key]; t = np.array(cd["t"])
            fr = np.array(cd["f_running"], dtype=float)
            m = t > 0.15 * t[-1]
            ax.plot(t[m], np.clip(fr[m], -3, 3), lw=1.2, color=TCOL[T], label=f"{T} K")
        ax.set_ylim(-2.6, 2.6)
        ax.set_title(sys, fontsize=8.5, pad=2)
        ax.tick_params(labelsize=7)
        if i % c == 0:
            ax.set_ylabel("f(t) = C_dist/MSD_self", fontsize=7)
        if i // c == r - 1:
            ax.set_xlabel("lag time (ps)", fontsize=7)
    axs.flat[0].legend(fontsize=6.5, ncol=2, loc="upper right")
    fig.suptitle("Running correlation factor f(t).  A converged f is a plateau inside the green ±1 band.",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f"{FIG}/f_running.png", dpi=130)
    plt.close(fig)


def fig_f_vs_T(df, order):
    n = len(order)
    r, c = grid(n)
    fig, axs = plt.subplots(r, c, figsize=(c * 2.5, r * 2.1), squeeze=False, sharey=True)
    for ax in axs.flat:
        ax.axis("off")
    for i, sys in enumerate(order):
        ax = axs.flat[i]; ax.axis("on")
        g = df[df.system == sys].sort_values("temperature")
        ax.axhspan(-1, 1, color="#1baf7a", alpha=.06)
        ax.axhline(0, color="#999", lw=.6)
        ax.errorbar(g.temperature, g.f, yerr=g.f_err, fmt="o-", ms=4,
                    lw=1, color="#2a78d6", ecolor="#b0740f", capsize=2)
        ax.set_ylim(-3, 3)
        ax.set_title(sys, fontsize=8.5, pad=2)
        ax.tick_params(labelsize=7)
        if i % c == 0:
            ax.set_ylabel("f", fontsize=8)
        if i // c == r - 1:
            ax.set_xlabel("T (K)", fontsize=7)
    fig.suptitle("Correlation factor f vs temperature   (error bars = systematic + statistical, in quadrature)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{FIG}/f_vs_T.png", dpi=130)
    plt.close(fig)


def summary(df, order):
    L = [f"# {SPECIES} Onsager transport — electrolyte shortlist\n",
         f"Following Karan et al. (2026), Eqs 2–6. {SPECIES} sub-lattice, COM frame, "
         "slopes fitted over the 0.2–0.6·t_max lag window of the ~20–40 ps NCSD "
         "production trajectories. Errors on f are split into `f_sys` (spread over "
         "later fit windows — the curve has not straightened) and `f_stat` (spread "
         "over 3 blocks fitted at matched lags), added in quadrature.\n",
         "`MSD_self = Σ_a|Δr_a|²`, `MSD_coll = |Σ_a Δr_a|²`, `C_dist = MSD_coll − MSD_self`; "
         "`D_xx = slope/6N` (cm²/s); `L_xx = slope/(6 k_BT V)` "
         "(10⁻⁸ mol²·J⁻¹·m⁻¹·s⁻¹); `f = L_dist/L_self = D_coll/D_self − 1`.\n",
         "## Convergence (relative scatter over matched-lag blocks, median over all 100 system×T)\n",
         f"- **D_self** : {(df.D_self_cm2s_stat/df.D_self_cm2s.abs()).median()*100:.0f}%  → usable",
         f"- **D_coll / L_coll** : {(df.D_coll_cm2s_stat/df.D_coll_cm2s.abs()).median()*100:.0f}%  → not converged",
         f"- **f** : total error ≈ {df.f_err.median():.2f} (systematic {df.f_sys.median():.2f}, statistical {df.f_stat.median():.2f}); physical range is −1…+1  → magnitude not converged",
         f"- f with a sign stable across all fit windows: **{int(df.f_sign_stable.sum())} / {len(df)}**",
         f"- f resolved to better than ±0.3: **{(df.f_err<0.3).sum()} / {len(df)}**\n",
         "## Values at 1500 K\n",
         "| system | D_self | D_coll | D_dist | f ± block | L_self | L_coll | L_dist | σ_NE | σ_true |",
         "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
    g = df[df.temperature == 1500].set_index("system")
    for s in order:
        if s not in g.index:
            continue
        r = g.loc[s]
        L.append(
            f"| {s} | {r.D_self_cm2s:.2e} | {r.D_coll_cm2s:.2e} | {r.D_dist_cm2s:+.2e} "
            f"| {r.f:+.2f} ± {r.f_err:.2f} | {r.L_self*1e8:.2f} | {r.L_coll*1e8:.2f} "
            f"| {r.L_dist*1e8:+.2f} | {r.sigma_NE_mScm:.0f} | {r.sigma_true_mScm:.0f} |")
    L += ["\n(σ in mS/cm, z_Li = 1. Full table incl. all temperatures: results/onsager.csv)"]
    open(os.path.join(FIG, "summary.md"), "w").write("\n".join(L) + "\n")


def main():
    global SPECIES, FIG, CB_ONLY
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="Li")
    ap.add_argument("--charge-balanced-only", action="store_true",
                    help="keep only stoichiometric ionic compositions")
    a = ap.parse_args()
    SPECIES = a.species
    globals()["CB_ONLY"] = a.charge_balanced_only
    FIG = os.path.join(RES, "figs", SPECIES.lower() + ("_cb" if a.charge_balanced_only else ""))
    os.makedirs(FIG, exist_ok=True)
    df, curves, order = load()
    for T in TEMPS:
        fig_curves(df, curves, order, T)
    fig_f_running(df, curves, order)
    fig_f_vs_T(df, order)
    summary(df, order)
    print(f"wrote {FIG}/*.png and {FIG}/summary.md")


if __name__ == "__main__":
    main()
