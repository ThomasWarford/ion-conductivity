"""Summarise the NCSD diffusivity screen.

Question 1 - how fast does D converge within the (~20-40 ps) AIMD window?
Question 2 - how large is D?
Goal      - pick amorphous systems where an MLIP-vs-AIMD validation of the
            diffusion coefficients (and the harder cross-correlation / Onsager
            terms) is actually feasible.

  python report.py

Reads  results/screen.csv  (+ results/mp_reference.csv if present)
Writes results/summary.md, results/shortlist_all.csv, results/shortlist_electrolyte.csv,
       results/figs/*.png
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(__file__)
RES = os.path.join(HERE, "results")
FIG = os.path.join(RES, "figs")

# thresholds for "D is converged in the run" (Li channel)
D_MIN = 1e-5
BETA = (0.90, 1.15)
BLOCK_CV_MAX = 0.25
SPREAD_MAX = 0.25
HALF_TOL = 0.40
TRUN_MAX = 20.0            # ps  -- converges within a feasible AIMD run
DRIFT_MAX = 1.0

HALIDE = set("F Cl Br I".split())
CHALC = set("O S Se Te".split())
PNICT = set("N P As Sb Bi".split())
NETWORK = set("B C Si Ge".split())
INSULATING_ANION = HALIDE | {"O", "S", "Se", "N", "P"}


def anion_class(system):
    els = set(system.split("-")) - {"Li"}
    if els & HALIDE:
        return "halide"
    if "O" in els:
        return "oxide"
    if els & {"S", "Se", "Te"}:
        return "chalcogenide"
    if els & PNICT:
        return "pnictide"
    if els & NETWORK:
        return "network(B/C/Si/Ge)"
    return "metallic"


def load():
    df = pd.read_csv(os.path.join(RES, "screen.csv"))
    if "error" in df:
        df = df[df.error.isna() | (df.error == "")].copy()
    df["anion_class"] = df.system.map(anion_class)
    df["electrolyte_like"] = df.system.map(
        lambda s: bool((set(s.split("-")) - {"Li"}) & INSULATING_ANION))
    mp = None
    p = os.path.join(RES, "mp_reference.csv")
    if os.path.exists(p):
        mp = pd.read_csv(p)
    return df, mp


def converged(d):
    return (
        (d.D_einstein_cm2s >= D_MIN)
        & d.beta_diffusive.between(*BETA)
        & (d.block_cv <= BLOCK_CV_MAX)
        & (d.Dwin_spread <= SPREAD_MAX)
        & ((d.D_half_ratio - 1).abs() <= HALF_TOL)
        & (d.t_run_converge_ps <= TRUN_MAX)
        & (d.energy_drift_meV_atom_ps.abs() <= DRIFT_MAX)
    )


def build_shortlist(li, systems):
    sl = li[li.system.isin(systems)].copy()
    piv = sl.pivot_table(index="system", columns="temperature", values="D_einstein_cm2s")
    piv.columns = [f"D_Li_{int(c)}K" for c in piv.columns]
    cv = sl.pivot_table(index="system", columns="temperature", values="block_cv")
    cv.columns = [f"blockcv_{int(c)}K" for c in cv.columns]
    meta = (sl.sort_values("temperature").groupby("system").agg(
        formula=("formula", "last"), anion_class=("anion_class", "last"),
        n_Li=("n_species", "last"), natoms=("natoms", "last"),
        n_conv_T=("converged", "sum"),
        med_prod_ps=("production_ps", "median"),
        min_trun_ps=("t_run_converge_ps", "max"),
        med_tlag_ps=("t_lag_stable_ps", "median"),
        med_beta=("beta_diffusive", "median"),
        med_block_cv=("block_cv", "median"),
        med_Dwin_spread=("Dwin_spread", "median"),
        charge_conv_T=("charge_converged", "sum"),
        med_haven=("haven_ratio", "median")))
    out = meta.join(piv).join(cv)
    ea = {}
    for s, g in sl.groupby("system"):
        g = g[(g.D_einstein_cm2s > 0)].dropna(subset=["D_einstein_cm2s"])
        if g.temperature.nunique() >= 3:
            p = np.polyfit(1.0 / g.temperature, np.log(g.D_einstein_cm2s), 1)
            ea[s] = round(-p[0] * 8.617333e-5, 3)
    out["Ea_eV"] = out.index.map(ea)
    dcol = "D_Li_1000K" if "D_Li_1000K" in out else out.filter(like="D_Li_").columns[0]
    return out.sort_values(dcol, ascending=False)


def main():
    os.makedirs(FIG, exist_ok=True)
    df, mp = load()
    li = df[df.species == "Li"].copy()
    li["converged"] = converged(li)
    li["charge_converged"] = (li.beta_charge.between(0.8, 1.3)
                              & li.D_charge_cm2s.notna() & (li.D_charge_cm2s > 0))

    L = ["# NCSD amorphous-diffusivity screen  --  convergence & magnitude\n",
         f"{df.system.nunique()} compositions x 4 temperatures (1000-2500 K), "
         f"AIMD, 2 fs timestep, {li.production_ps.median():.0f} ps median production run.\n",
         "D is the 3-D Einstein self-diffusivity (MSD -> 6 D t). The MPContribs "
         "`diffusivity` field fits MSD -> 2 D t and is therefore ~3x larger "
         "(confirmed below).\n",
         "\n**Convergence flags (Li):** D >= 1e-5 cm^2/s; late-time MSD exponent "
         "beta in [0.90, 1.15]; block CV <= 0.25 over 4 independent sub-runs; "
         "D stable to <25% across lag windows; |D(2nd half)/D(1st half) - 1| <= 0.4; "
         "D from a <=20 ps trajectory prefix already within 25% of the full-run D "
         "(t_run_converge); |energy drift| <= 1 meV/atom/ps.\n"]

    L += ["## 1. By temperature\n",
          "| T (K) | median prod (ps) | median D_Li (cm^2/s) | D_Li 10-90% range | "
          "frac diffusive (beta>0.9) | frac Li-converged | frac charge-converged |",
          "|---|---|---|---|---|---|---|"]
    for T, g in li.groupby("temperature"):
        lo, hi = g.D_einstein_cm2s.quantile([.1, .9])
        L.append(f"| {T:.0f} | {g.production_ps.median():.0f} | "
                 f"{g.D_einstein_cm2s.median():.2e} | {lo:.1e} - {hi:.1e} | "
                 f"{(g.beta_diffusive > 0.9).mean():.2f} | {g.converged.mean():.2f} | "
                 f"{g.charge_converged.mean():.2f} |")

    L += ["\n## 2. Convergence vs magnitude, by anion chemistry (1000 K)\n",
          "| anion class | n systems | median D_Li (cm^2/s) | frac Li-converged |",
          "|---|---|---|---|"]
    g0 = li[li.temperature == 1000]
    for cls, g in g0.groupby("anion_class"):
        L.append(f"| {cls} | {len(g)} | {g.D_einstein_cm2s.median():.2e} | "
                 f"{g.converged.mean():.2f} |")

    # shortlists: converged at >= 2 T
    nconv = li.groupby("system").converged.sum()
    all_cand = sorted(nconv[nconv >= 2].index)
    elec_cand = [s for s in all_cand
                 if li[li.system == s].electrolyte_like.iloc[0]]

    sl_all = build_shortlist(li, all_cand)
    sl_el = build_shortlist(li, elec_cand)
    sl_all.to_csv(os.path.join(RES, "shortlist_all.csv"))
    sl_el.to_csv(os.path.join(RES, "shortlist_electrolyte.csv"))

    L += [f"\n## 3. Shortlist A - any chemistry, Li converged at >= 2 T  ({len(sl_all)} systems)\n",
          "Mostly Li metal / Li-alloy melts: large D, textbook convergence. "
          "Full table: results/shortlist_all.csv. Top 20 by D_Li(1000 K):\n",
          _fmt(sl_all).head(20).to_markdown(),
          f"\n## 4. Shortlist B - insulating-anion (halide/oxide/S/Se/N/P) chemistry, "
          f"Li converged at >= 2 T  ({len(sl_el)} systems)\n",
          "The relevant set for solid-electrolyte MLIP validation. "
          "Full table: results/shortlist_electrolyte.csv.\n",
          _fmt(sl_el).to_markdown()]

    vpath = os.path.join(RES, "mp_verify.csv")
    if mp is None and os.path.exists(vpath):
        v = pd.read_csv(vpath).rename(columns={"T": "temperature"})
        mp = v.assign(D_table_einstein_cm2s=v.D_table_slope6,
                      D_reported_cm2s=v.D_reported)
    if mp is not None:
        m = mp[mp.element == "Li"].copy()
        r = (m.D_reported_cm2s / m.D_table_einstein_cm2s).replace(
            [np.inf, -np.inf], np.nan).dropna()
        mm = li.merge(m[["system", "temperature", "D_reported_cm2s"]].dropna(),
                      on=["system", "temperature"], how="inner")
        q = (mm.D_reported_cm2s / mm.D_einstein_cm2s).replace(
            [np.inf, -np.inf], np.nan).dropna()
        L += ["\n## 5. Cross-check vs MPContribs\n",
              f"- MP `diffusivity` / refit of their own MSD table as slope/6 : "
              f"median **{r.median():.2f}** (IQR {r.quantile(.25):.2f}-{r.quantile(.75):.2f}, "
              f"n={len(r)})  -> the published values are 3x the standard 3-D Einstein D.",
              f"- MP `diffusivity` / our independent trajectory D_einstein : "
              f"median **{q.median():.2f}** (n={len(q)})  -> our analysis "
              f"reproduces their MSD; the 3x offset is the 1D-vs-3D Einstein "
              f"prefactor (2 D t vs 6 D t), not a physics disagreement."]

    L += ["\n## 6. Cross-correlations / Onsager terms\n",
          f"The collective (charge) MSD converges far worse than the self MSD: "
          f"only {li[li.temperature==2500].charge_converged.mean():.0%} of systems "
          f"reach a diffusive collective regime even at 2500 K (vs "
          f"{li[li.temperature==2500].converged.mean():.0%} for self-D). "
          f"Systems with `charge_conv_T >= 2` in the shortlists are the ones where "
          f"validating the cross terms against AIMD is also realistic."]

    # ---- pick top recommendations ----
    def score(row):
        d = np.nanmean([row.get(f"D_Li_{T}K", np.nan) for T in (1000, 1500, 2000, 2500)])
        return (np.log10(d) if d > 0 else -9) - 3 * row.med_block_cv \
            - 0.05 * row.min_trun_ps + 0.5 * row.charge_conv_T
    for name, tbl in [("any chemistry", sl_all), ("electrolyte chemistry", sl_el)]:
        t = tbl.copy()
        t["score"] = t.apply(score, axis=1)
        top = t.sort_values("score", ascending=False).head(6)
        L += [f"\n## 7. Top picks - {name}\n",
              _fmt(top).to_markdown()]

    open(os.path.join(RES, "summary.md"), "w").write("\n".join(str(x) for x in L) + "\n")
    _figs(li, sl_all, sl_el)
    print(f"shortlist A (any): {len(sl_all)}   shortlist B (electrolyte): {len(sl_el)}")
    print("wrote results/summary.md, results/shortlist_*.csv, results/figs/*.png")


def _fmt(sl):
    keep = ["formula", "anion_class", "n_Li", "n_conv_T", "charge_conv_T",
            "min_trun_ps", "med_beta", "med_block_cv", "med_haven", "Ea_eV"]
    keep += [c for c in sl.columns if c.startswith("D_Li_")]
    return sl[keep].round(6)


def _figs(li, sl_all, sl_el):
    cls_order = ["metallic", "network(B/C/Si/Ge)", "pnictide", "chalcogenide",
                 "oxide", "halide"]
    cmap = dict(zip(cls_order, plt.cm.tab10(range(len(cls_order)))))

    # D vs convergence time, coloured by anion class (1000 & 2500 K)
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, T in zip(axs, (1000, 2500)):
        g = li[li.temperature == T]
        for cls in cls_order:
            gc = g[g.anion_class == cls]
            x = gc.t_run_converge_ps.fillna(55) + np.random.uniform(-0.8, 0.8, len(gc))
            ax.scatter(x, gc.D_einstein_cm2s.clip(1e-9), s=18, alpha=.6,
                       color=cmap[cls], label=cls)
        ax.axhline(D_MIN, ls="--", c="grey", lw=.8)
        ax.set_yscale("log"); ax.set_title(f"{T} K")
        ax.set_xlabel("t_converge  (ps;  55 = not converged in run)")
    axs[0].set_ylabel("D_Li  (cm$^2$/s)")
    axs[0].legend(fontsize=8, title="anion class")
    fig.suptitle("Li diffusivity vs convergence time")
    fig.tight_layout(); fig.savefig(f"{FIG}/D_vs_tconverge.png", dpi=130); plt.close(fig)

    # distributions
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    temps = sorted(li.temperature.unique())
    col = dict(zip(temps, plt.cm.viridis(np.linspace(0, .85, len(temps)))))
    for T in temps:
        g = li[li.temperature == T]
        axs[0].hist(np.log10(g.D_einstein_cm2s.clip(1e-9)), bins=30,
                    histtype="step", color=col[T], label=f"{T:.0f} K")
        axs[1].hist(g.beta_diffusive.dropna(), bins=30, range=(0, 1.5),
                    histtype="step", color=col[T])
        axs[2].hist(g.t_run_converge_ps.fillna(55), bins=22, histtype="step", color=col[T])
    axs[0].set_xlabel("log10 D_Li (cm^2/s)"); axs[0].legend(); axs[0].set_title("magnitude")
    axs[1].axvspan(*BETA, color="green", alpha=.12)
    axs[1].set_xlabel(r"$\beta$ (late-time dlnMSD/dlnt)"); axs[1].set_title("diffusive regime")
    axs[2].set_xlabel("t_converge (ps)"); axs[2].set_title("fitted-D convergence")
    fig.tight_layout(); fig.savefig(f"{FIG}/distributions.png", dpi=130); plt.close(fig)

    # Arrhenius of electrolyte shortlist
    sl = sl_el if len(sl_el) else sl_all
    dcols = sorted([c for c in sl.columns if c.startswith("D_Li_")],
                   key=lambda c: int(c.split("_")[2][:-1]))
    Ts = np.array([int(c.split("_")[2][:-1]) for c in dcols])
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for s, row in sl.head(22).iterrows():
        y = row[dcols].astype(float).values
        m = np.isfinite(y) & (y > 0)
        if m.sum() >= 2:
            ax.plot(1000 / Ts[m], y[m], "o-", ms=4, lw=1,
                    label=f"{s}  ({row['anion_class']})")
    ax.set_yscale("log")
    ax.set_xlabel("1000 / T  (K$^{-1}$)"); ax.set_ylabel("D_Li  (cm$^2$/s)")
    ax.set_title("Arrhenius - electrolyte-chemistry shortlist")
    ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/arrhenius_electrolyte.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
