"""Un-averaged version of the Onsager correlation functions.

The curves in results/figs/onsager_curves_*.png are averages over time origins:
at lag tau, MSD(tau) averages (n - tau) values, one per starting frame.  This
script plots the individual per-origin values behind the mean, so you can see
how noisy a single measurement is and how the average emerges -- and, crucially,
how few *independent* (non-overlapping) origins survive at long lag.

Per origin t0 and lag tau, for the Li sub-lattice:
    MSD_self : Sum_a |dr_a(t0+tau) - dr_a(t0)|^2
    MSD_coll : |Sum_a [dr_a(t0+tau) - dr_a(t0)]|^2
    C_dist   : MSD_coll - MSD_self
Averaging each over t0 reproduces exactly the curves in the main figure.

  python onsager_raw.py [--procs N] [--temps 1000,1500,2000,2500]

Writes results/figs/raw/raw_<quantity>_<T>K.png and results/figs/raw/origin_counts.csv
"""
from __future__ import annotations

import argparse
import glob
import math
import os
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ncsd_traj import load_production
from diffusivity import com_correct

HERE = os.path.dirname(__file__)
RES = os.path.join(HERE, "results")
FIG = os.path.join(RES, "figs", "li", "raw")
DATA = os.path.join(HERE, "..", "data_ncsd", "trajectories")
SPECIES = "Li"                      # set by --species
WINDOW = (0.2, 0.6)        # the shaded fit window, as a fraction of t_max
NLAG = 42                  # lags at which the raw cloud is drawn
NSHOW = 140                # origins drawn per lag (subsampled for legibility)

CURVE = {"MSD_self": "#2a78d6", "MSD_coll": "#eb6834", "C_dist": "#1baf7a"}
TITLE = {"MSD_self": r"MSD$_{self}$ = $\sum_a|\Delta r_a|^2$",
         "MSD_coll": r"MSD$_{coll}$ = $|\sum_a \Delta r_a|^2$",
         "C_dist": r"C$_{dist}$ = MSD$_{coll}$ $-$ MSD$_{self}$"}


def raw_clouds(args):
    """Per-origin values at a grid of lags, plus the mean and the origin counts."""
    system, T = args
    g = glob.glob(os.path.join(DATA, f"temperature={T}K",
                               f"chemical_system={system}", "*.jsonl.gz"))
    if not g:
        return None
    tr = load_production(g[0])
    dt = tr.ps_per_frame
    cart = com_correct(tr.cart_unwrapped(), tr.symbols)
    start = max(tr.n_equil_frames, int(round(2.0 / dt)))
    if (tr.nframes - start) * dt < 8:
        start = max(int(round(2.0 / dt)), tr.nframes - int(round(20.0 / dt)))
    sel = tr.symbols == SPECIES
    x = cart[start:, sel] - cart[start:, sel][0]        # (n, N, 3)
    n, N, _ = x.shape
    s = x.sum(axis=1)                                   # (n, 3) collective

    lags = np.unique(np.linspace(1, n - 2, NLAG).astype(int))
    out = {k: dict(t=[], mean=[], pts_t=[], pts_y=[]) for k in CURVE}
    counts = []
    for lag in lags:
        d_self = ((x[lag:] - x[:n - lag]) ** 2).sum(-1).sum(1)    # (n-lag,) per origin
        d_coll = ((s[lag:] - s[:n - lag]) ** 2).sum(-1)           # (n-lag,)
        vals = {"MSD_self": d_self, "MSD_coll": d_coll, "C_dist": d_coll - d_self}
        n_orig = n - lag
        step = max(1, n_orig // NSHOW)
        for k, v in vals.items():
            out[k]["t"].append(lag * dt)
            out[k]["mean"].append(float(v.mean()))
            sub = v[::step][:NSHOW]
            out[k]["pts_t"].append(np.full(len(sub), lag * dt))
            out[k]["pts_y"].append(sub)
        counts.append((lag * dt, n_orig, max(1, n_orig // lag)))

    for k in out:
        out[k] = dict(t=np.array(out[k]["t"]), mean=np.array(out[k]["mean"]),
                      pts_t=np.concatenate(out[k]["pts_t"]),
                      pts_y=np.concatenate(out[k]["pts_y"]))
    tmax = (n - 1) * dt
    cnt = np.array(counts)
    lo_ps, hi_ps = WINDOW[0] * tmax, WINDOW[1] * tmax
    def at(ps):
        i = int(np.argmin(np.abs(cnt[:, 0] - ps)))
        return int(cnt[i, 1]), int(cnt[i, 2])
    o_lo, i_lo = at(lo_ps)
    o_hi, i_hi = at(hi_ps)
    meta = dict(system=system, temperature=T, N_Li=N, nframes=n, total_ps=round(tmax, 1),
                fit_lo_ps=round(lo_ps, 1), fit_hi_ps=round(hi_ps, 1),
                origins_at_lo=o_lo, origins_at_hi=o_hi,
                independent_at_lo=i_lo, independent_at_hi=i_hi)
    return system, T, out, meta


def grid(n):
    c = math.ceil(math.sqrt(n))
    return math.ceil(n / c), c


def fig_raw(cache, order, quantity, T):
    keys = [s for s in order if (s, T) in cache]
    r, c = grid(len(keys))
    fig, axs = plt.subplots(r, c, figsize=(c * 3.0, r * 2.5), squeeze=False)
    for ax in axs.flat:
        ax.axis("off")
    col = CURVE[quantity]
    for i, sys in enumerate(keys):
        ax = axs.flat[i]; ax.axis("on")
        d, meta = cache[(sys, T)]
        q = d[quantity]
        ax.axvspan(meta["fit_lo_ps"], meta["fit_hi_ps"], color="#000", alpha=.05, lw=0)
        ax.scatter(q["pts_t"], q["pts_y"] / 1e3, s=1.2, color=col, alpha=.10,
                   linewidths=0, rasterized=True)
        ax.plot(q["t"], q["mean"] / 1e3, lw=1.8, color="#111", zorder=5)
        ax.plot(q["t"], q["mean"] / 1e3, lw=1.0, color=col, zorder=6)
        ax.axhline(0, color="#999", lw=.6)
        ax.set_title(sys, fontsize=8.5, pad=2)
        ax.text(.03, .965,
                f"{meta['independent_at_lo']} → {meta['independent_at_hi']} independent origins\n"
                f"across the shaded window",
                transform=ax.transAxes, va="top", fontsize=6.2, color="#666")
        ax.tick_params(labelsize=7)
        if i % c == 0:
            ax.set_ylabel(r"$10^3\ \AA^2$", fontsize=7)
        if i // c == r - 1:
            ax.set_xlabel("lag time (ps)", fontsize=7)
    fig.suptitle(f"{SPECIES}:  {TITLE[quantity]}  at {T} K  —  every time origin, un-averaged",
                 fontsize=12, y=1.002)
    fig.text(0.5, 0.9655,
             f"faint dots = one starting frame each (up to {NSHOW} drawn per lag);  "
             "black/coloured line = their mean, i.e. the curve in the main figure",
             ha="center", fontsize=8.5, color="#555")
    fig.tight_layout(rect=[0, 0, 1, 0.953])
    out = f"{FIG}/raw_{quantity}_{T}K.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=min(24, len(os.sched_getaffinity(0))))
    ap.add_argument("--temps", default="1000,1500,2000,2500")
    ap.add_argument("--species", default="Li")
    args = ap.parse_args()
    global SPECIES, FIG
    SPECIES = args.species
    FIG = os.path.join(RES, "figs", SPECIES.lower(), "raw")
    temps = [int(t) for t in args.temps.split(",")]
    os.makedirs(FIG, exist_ok=True)

    on = pd.read_csv(os.path.join(RES, "onsager.csv"))
    on = on[on.error.fillna("") == ""]
    if "species" in on:
        on = on[on.species == SPECIES]
    order = (on[on.temperature == 1500].sort_values("D_self_cm2s", ascending=False)
             .system.tolist())
    order += [s for s in on.system.unique() if s not in order]

    jobs = [(s, T) for s in order for T in temps]
    print(f"{len(jobs)} (system, T) on {args.procs} procs", flush=True)
    cache, metas = {}, []
    with Pool(args.procs) as pool:
        for k, res in enumerate(pool.imap_unordered(raw_clouds, jobs)):
            if res is None:
                continue
            sysname, T, out, meta = res
            cache[(sysname, T)] = (out, meta)
            metas.append(meta)
            if (k + 1) % 20 == 0:
                print(f"  {k + 1}/{len(jobs)}", flush=True)

    for T in temps:
        for q in CURVE:
            print("wrote", fig_raw(cache, order, q, T), flush=True)

    md = pd.DataFrame(metas).sort_values(["system", "temperature"])
    md.to_csv(os.path.join(FIG, "origin_counts.csv"), index=False)
    print("\nAveraging depth across the shaded fit window "
          f"(lags {WINDOW[0]}–{WINDOW[1]} × t_max), median over all runs:")
    print(f"  nominal time origins : {md.origins_at_lo.median():.0f}"
          f"  →  {md.origins_at_hi.median():.0f}")
    print(f"  INDEPENDENT windows  : {md.independent_at_lo.median():.0f}"
          f"  →  {md.independent_at_hi.median():.0f}")


if __name__ == "__main__":
    main()
