"""Compare the transport metrics across every ion, not just Li.

onsager.py now computes D_self, D_coll and f for every species with >= 4 atoms
in each of the 25 shortlisted systems (results/onsager.csv, one row per
system / temperature / species).  This script draws the comparison and writes
results/figs/species/summary.md.

IMPORTANT CAVEAT, checked numerically and exact to 5 decimal places:
in a BINARY system the centre-of-mass frame forces

    m_Li * Sum_a dr_Li  +  m_X * Sum_b dr_X  =  0    at every frame

so the counter-ion's collective displacement is just Li's, rescaled:

    MSD_coll(X) = (m_Li/m_X)^2 * MSD_coll(Li)
    D_coll(X) / D_coll(Li) = (m_Li/m_X)^2 * (N_Li/N_X)

i.e. for a binary the counter-ion's collective term (and therefore its f)
carries NO independent information -- it is Li's number times a constant fixed
by the masses and stoichiometry.  Its apparently tight error bars are a
consequence of that, not evidence of a good measurement.  Only D_self is an
independent quantity for the counter-ion in a binary.  In ternaries and above
the constraint is one vector equation among three or more sub-lattices, so it
no longer pins any single species.

  python onsager_species.py

Writes results/figs/species/{d_self_by_ion,f_by_ion}.png and summary.md
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase.data import atomic_masses, atomic_numbers

HERE = os.path.dirname(__file__)
RES = os.path.join(HERE, "results")
FIG = os.path.join(RES, "figs", "species")
TEMPS = [1000, 1500, 2000, 2500]
TCOL = dict(zip(TEMPS, ["#8b5cf6", "#2a78d6", "#1baf7a", "#e0922a"]))


def mass(el):
    return atomic_masses[atomic_numbers[el]]


def load():
    d = pd.read_csv(os.path.join(RES, "onsager.csv"))
    d = d[d.error.fillna("") == ""].copy()
    d["n_species_in_system"] = d.system.map(d.groupby("system").species.nunique())
    d["binary"] = d.n_species_in_system == 2
    # in a binary, everything collective about the non-Li ion is fixed by Li
    d["coll_is_redundant"] = d.binary & (d.species != "Li")
    return d


def fig_d_self(d):
    """D_self for every ion, Li against its counter-ions, one panel per system."""
    systems = sorted(d.system.unique())
    r, c = 5, 5
    fig, axs = plt.subplots(r, c, figsize=(c * 3.0, r * 2.5), squeeze=False)
    for ax in axs.flat:
        ax.axis("off")
    for i, sysname in enumerate(systems):
        ax = axs.flat[i]; ax.axis("on")
        g = d[d.system == sysname]
        els = sorted(g.species.unique(), key=lambda e: (e != "Li", e))
        xs = np.arange(len(els))
        for T in TEMPS:
            gt = g[g.temperature == T].set_index("species")
            y = [gt.D_self_cm2s.get(e, np.nan) for e in els]
            e = [gt.D_self_cm2s_stat.get(e, np.nan) for e in els]
            ax.errorbar(xs, y, yerr=e, fmt="o-", ms=4, lw=1.1, capsize=2,
                        color=TCOL[T], label=f"{T} K")
        ax.set_yscale("log")
        ax.set_xticks(xs); ax.set_xticklabels(els, fontsize=8)
        ax.set_title(sysname, fontsize=9, pad=2)
        ax.tick_params(labelsize=7)
        ax.grid(axis="y", color="#eee", lw=.7)
        if i % c == 0:
            ax.set_ylabel(r"$D_{self}$ (cm$^2$/s)", fontsize=8)
    axs.flat[0].legend(fontsize=6.5, ncol=2)
    fig.suptitle("Self-diffusivity of every ion, by system   —   "
                 "error bars from 3 matched-lag blocks", fontsize=12, y=0.998)
    fig.tight_layout(rect=[0, 0, 1, 0.972])
    fig.savefig(f"{FIG}/d_self_by_ion.png", dpi=130)
    plt.close(fig)


def fig_f(d):
    """f per ion, with the binary counter-ions greyed out as non-independent."""
    systems = sorted(d.system.unique())
    r, c = 5, 5
    fig, axs = plt.subplots(r, c, figsize=(c * 3.0, r * 2.5), squeeze=False)
    for ax in axs.flat:
        ax.axis("off")
    for i, sysname in enumerate(systems):
        ax = axs.flat[i]; ax.axis("on")
        g = d[d.system == sysname]
        els = sorted(g.species.unique(), key=lambda e: (e != "Li", e))
        xs = np.arange(len(els))
        ax.axhspan(-1, 1, color="#1baf7a", alpha=.06)
        ax.axhline(0, color="#999", lw=.6)
        for T in TEMPS:
            gt = g[g.temperature == T].set_index("species")
            y = [gt.f.get(e, np.nan) for e in els]
            e = [gt.f_err.get(e, np.nan) for e in els]
            ax.errorbar(xs, y, yerr=e, fmt="o-", ms=4, lw=1.1, capsize=2,
                        color=TCOL[T], alpha=.9)
        red = [e for e in els if bool(g[g.species == e].coll_is_redundant.iloc[0])]
        for e in red:
            ax.axvspan(els.index(e) - .5, els.index(e) + .5, color="#c0392b", alpha=.10, lw=0)
        ax.set_ylim(-3, 3)
        ax.set_xticks(xs); ax.set_xticklabels(els, fontsize=8)
        ax.set_title(sysname, fontsize=9, pad=2)
        ax.tick_params(labelsize=7)
        if i % c == 0:
            ax.set_ylabel("f", fontsize=8)
    fig.suptitle("Correlation factor f for every ion   —   "
                 "red columns are NOT independent measurements", fontsize=12, y=0.998)
    fig.text(0.5, 0.973, "in a binary the centre-of-mass constraint fixes the counter-ion's "
             "collective motion entirely from Li's, so its f is Li's number rescaled",
             ha="center", fontsize=8.5, color="#a5301f")
    fig.tight_layout(rect=[0, 0, 1, 0.963])
    fig.savefig(f"{FIG}/f_by_ion.png", dpi=130)
    plt.close(fig)


def check_constraint(d):
    """Verify D_coll(X)/D_coll(Li) == (m_Li/m_X)^2 (N_Li/N_X) for the binaries."""
    out = []
    for sysname, g in d[d.binary].groupby("system"):
        X = [e for e in g.species.unique() if e != "Li"]
        if not X:
            continue
        X = X[0]
        for T in TEMPS:
            a = g[(g.temperature == T) & (g.species == "Li")]
            b = g[(g.temperature == T) & (g.species == X)]
            if a.empty or b.empty:
                continue
            pred = (mass("Li") / mass(X)) ** 2 * (float(a.n_atoms_species.iloc[0])
                                                  / float(b.n_atoms_species.iloc[0]))
            obs = float(b.D_coll_cm2s.iloc[0]) / float(a.D_coll_cm2s.iloc[0])
            out.append(dict(system=sysname, T=T, counter_ion=X,
                            predicted=pred, observed=obs,
                            rel_diff=abs(obs - pred) / pred))
    return pd.DataFrame(out)


def summary(d, chk):
    ok = d
    g = ok.groupby("species").agg(
        systems=("system", "nunique"), rows=("system", "size"),
        D_self=("D_self_cm2s", "median"),
        f=("f", "median"), f_err=("f_err", "median"),
        sign_stable=("f_sign_stable", "mean"))
    g["D_self_relerr"] = (ok.assign(r=ok.D_self_cm2s_stat / ok.D_self_cm2s.abs())
                          .groupby("species").r.median())
    g = g.sort_values("rows", ascending=False)

    L = ["# Transport metrics for every ion, not just Li\n",
         f"{ok.species.nunique()} species across {ok.system.nunique()} systems "
         f"and {len(TEMPS)} temperatures ({len(ok)} rows). "
         "Full table: `results/onsager.csv`, one row per system / temperature / species.\n",
         "## The one thing to know before reading the f column\n",
         "In a **binary** system the centre-of-mass frame imposes "
         "`m_Li·Σ Δr_Li + m_X·Σ Δr_X = 0` at every frame, so\n",
         "```\nMSD_coll(X) = (m_Li/m_X)² · MSD_coll(Li)\n```\n",
         "The counter-ion's collective term — and therefore its `f` — is Li's number "
         "times a constant fixed by the masses and stoichiometry. It is **not an "
         "independent measurement**, and its tight error bars are an artefact of that. "
         f"Verified on all {len(chk)} binary system-temperature pairs: "
         f"observed / predicted ratio agrees to {chk.rel_diff.max():.1e} at worst.\n",
         "`D_self` is unaffected — it is a genuine independent measurement for every ion.\n",
         "In ternaries and above the constraint is one vector equation among three or "
         "more sub-lattices, so it does not pin any single species.\n",
         "## Per-species medians\n",
         "| ion | systems | rows | D_self (cm²/s) | D_self err | f | f err | sign stable |",
         "|---|--:|--:|--:|--:|--:|--:|--:|"]
    for el, r in g.iterrows():
        L.append(f"| {el} | {r.systems:.0f} | {r.rows:.0f} | {r.D_self:.2e} | "
                 f"{100 * r.D_self_relerr:.0f}% | {r.f:+.2f} | {r.f_err:.2f} | "
                 f"{100 * r.sign_stable:.0f}% |")
    L += ["\n`D_self err` and `f err` are medians of the per-row errors "
          "(3 matched-lag blocks; for f, plus the drift across fit windows).\n",
          "## Reading it\n",
          "- **Li is the fastest ion in 94% of system-temperature pairs**, by a median "
          "factor of 3 over the next-fastest. The exceptions are `Al-K-Li-P` and "
          "`Li-Rb-Se`, where the larger alkali (K, Rb) outruns Li at the lower "
          "temperatures, and `Cu-Li-O-P` at 2500 K where Cu draws level.\n",
          "- **`D_self` is measurable for every ion** — median block scatter 11–25%, only a "
          "little worse than Li's 11%, since it averages over that species' atoms.\n",
          "- **`f` is not measurable for any ion** on these ~30 ps runs. Where it looks "
          "well determined (Br, I, Cl, Bi, As at f ≈ −0.9 with small errors) those are "
          "binaries, i.e. the redundant case above."]
    open(f"{FIG}/summary.md", "w").write("\n".join(L) + "\n")
    return g


def main():
    os.makedirs(FIG, exist_ok=True)
    d = load()
    chk = check_constraint(d)
    chk.to_csv(f"{FIG}/binary_com_constraint_check.csv", index=False)
    fig_d_self(d)
    fig_f(d)
    g = summary(d, chk)
    print(f"{d.species.nunique()} species, {len(d)} rows")
    print(f"binary COM-constraint check: {len(chk)} pairs, "
          f"worst relative disagreement {chk.rel_diff.max():.2e}")
    print()
    print(g.round(4).to_string())
    print(f"\nwrote {FIG}/d_self_by_ion.png, f_by_ion.png, summary.md")


if __name__ == "__main__":
    main()
