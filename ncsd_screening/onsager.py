"""Onsager transport analysis for the electrolyte shortlist, following
Karan et al., Nat. Mater. (2026), Eqs 2-6.

For the Li sub-lattice of each (system, temperature), all in Angstrom^2 vs lag:

  MSD_self(t) = < Sum_a |dr_a(t)|^2 >     mean-square displacement, summed over
                                          the N_Li ions (= N_Li x the per-ion MSD)
  MSD_coll(t) = < |Sum_a dr_a(t)|^2 >     mean-square displacement of the *net*
                                          Li displacement (the "charge" MSD)
  C_dist(t)   = MSD_coll(t) - MSD_self(t) the ion-ion correlation part

Naming: the first two really are mean-square displacements, so they say MSD.
C_dist is a *difference* of mean squares -- not positive-definite, routinely
negative -- so it keeps the correlation-function C.  In the general Onsager
formalism (Karan Eq 2) every one of these is a displacement correlation
function C_ij(t) = <Sum_a dr_i^a . Sum_b dr_j^b>; the i != j terms are not
mean squares at all.  Only the i == j ones collapse to an MSD, which is why
the diagonal gets the more legible name here.

  L_xx = slope(.) / (6 k_B T V)          Onsager coefficient  [mol^2 / (J m s)]
  D_xx = slope(.) / (6 N_Li)             effective diffusivity [cm^2/s]
  f    = L_dist / L_self = D_coll/D_self - 1        Karan correlation factor

Displacements are taken in the mass-weighted centre-of-mass frame (Eq 2 is
defined in the COM frame).  Slopes are OLS fits over a fractional lag window.

Uncertainties on f are reported as TWO separate numbers, because they have
different causes and different cures:

  f_stat  statistical.  Std over NBLOCK contiguous blocks, every block fitted
          over the SAME absolute lag window (MATCH_WIN x one block length).
          Shrinks with more independent trajectory.  Note: scaling each block's
          window to its own length instead -- the naive thing -- makes short
          blocks sample short lags, which mixes f_sys into f_stat and badly
          overstates the noise.
  f_sys   systematic.  Half the spread of f refitted over DRIFT_WINDOWS, i.e.
          progressively later lags.  Non-zero means the correlation function
          has not reached its asymptotic linear regime, so there is no single
          slope yet.  Shrinks only with a LONGER trajectory.
  f_err   the two added in quadrature.

MSD_self needs no such care: it is linear at every lag (R^2 >= 0.999), so the
fit window does not matter and its block scatter is a genuine statistical error.

  python onsager.py [--procs N] [--window 0.2 0.6]

Writes results/onsager.csv and results/onsager_curves.json.gz
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import os
from multiprocessing import Pool

import numpy as np

from ncsd_traj import load_production
from diffusivity import com_correct, msd_fft

KB = 1.380649e-23           # J/K
NA = 6.02214076e23          # 1/mol
QE = 1.602176634e-19        # C
ANG2PS_TO_M2S = 1e-8        # (1e-10 m)^2 / (1e-12 s)
ANG2PS_TO_CM2S = 1e-4

HERE = os.path.dirname(__file__)
RES = os.path.join(HERE, "results")
DATA = os.path.join(HERE, "..", "data_ncsd", "trajectories")
SPECIES = "all"           # "all", or a single element symbol
MIN_ATOMS = 4             # skip species with fewer atoms than this
TEMPS = (1000, 1500, 2000, 2500)
NCURVE = 220               # points kept for plotting
NBLOCK = 3                 # blocks for the statistical error
MATCH_WIN = (0.25, 0.60)   # fit window as a fraction of ONE BLOCK, applied
                           # as the same absolute lag range to every block
DRIFT_WINDOWS = ((0.10, 0.25), (0.18, 0.40), (0.25, 0.55),
                 (0.35, 0.70))   # fractions of the full trajectory, for f_sys


def _slope(t, y, lo, hi):
    tmax = t[-1]
    m = (t >= lo * tmax) & (t <= hi * tmax) & (t > 0)
    if m.sum() < 5:
        return np.nan
    A = np.vstack([t[m], np.ones(m.sum())]).T
    return np.linalg.lstsq(A, y[m], rcond=None)[0][0]     # Angstrom^2 / ps


def _abs_window_fit(cart_li, dt, lo_ps, hi_ps):
    """D_self, D_coll and f from a fit over an ABSOLUTE lag window [lo_ps, hi_ps].

    Used for the error estimates, where every block and every window must be
    fitted over the same lags for the comparison to mean anything.
    """
    n, N, _ = cart_li.shape
    tt = np.arange(n) * dt
    m = (tt >= lo_ps) & (tt <= hi_ps) & (tt > 0)
    if m.sum() < 5:
        return dict(D_self_cm2s=np.nan, D_coll_cm2s=np.nan, f=np.nan)
    msd_self = N * msd_fft(cart_li).sum(axis=2).mean(axis=1)
    msd_coll = msd_fft(cart_li.sum(axis=1)).sum(axis=1)
    A = np.vstack([tt[m], np.ones(m.sum())]).T
    s_self = np.linalg.lstsq(A, msd_self[m], rcond=None)[0][0]
    s_coll = np.linalg.lstsq(A, msd_coll[m], rcond=None)[0][0]
    d_self = s_self / (6 * N) * ANG2PS_TO_CM2S
    d_coll = s_coll / (6 * N) * ANG2PS_TO_CM2S
    return dict(D_self_cm2s=d_self, D_coll_cm2s=d_coll,
                f=(d_coll / d_self - 1.0) if s_self else np.nan)


def _coeffs(cart_li, dt, volume_A3, T, window):
    """cart_li: (nframes, N, 3) COM-corrected displacements of the Li sub-lattice."""
    n, N, _ = cart_li.shape
    tt = np.arange(n) * dt
    # per-ion MSD, averaged over ions and time origins, then scaled by N so that
    # it is on the same footing as the collective one and the difference is meaningful
    msd_per_ion = msd_fft(cart_li).sum(axis=2).mean(axis=1)
    msd_self = N * msd_per_ion
    # sum the displacement vectors over ions FIRST, then take one MSD of that resultant
    msd_coll = msd_fft(cart_li.sum(axis=1)).sum(axis=1)
    c_dist = msd_coll - msd_self

    lo, hi = window
    s_self = _slope(tt, msd_self, lo, hi)
    s_coll = _slope(tt, msd_coll, lo, hi)
    s_dist = s_coll - s_self

    kT = KB * T
    V = volume_A3 * 1e-30
    pref_L = 1.0 / (6.0 * kT * V) * ANG2PS_TO_M2S / NA**2         # -> mol^2/(J m s)
    n_dens = N / V                                                # m^-3

    D_self = s_self / (6 * N) * ANG2PS_TO_CM2S
    D_coll = s_coll / (6 * N) * ANG2PS_TO_CM2S
    D_dist = D_coll - D_self
    f = D_dist / D_self if D_self else np.nan

    # conductivity: sigma = z^2 e^2 n D / (kB T). D in m^2/s.
    sig_NE = QE**2 * n_dens * (D_self * 1e-4) / kT / 100.0         # S/cm  (assumes z = 1)
    sig_true = QE**2 * n_dens * (D_coll * 1e-4) / kT / 100.0

    curves = dict(t=tt, MSD_self=msd_self, MSD_coll=msd_coll, C_dist=c_dist)
    vals = dict(
        n_atoms_species=N, volume_A3=volume_A3,
        D_self_cm2s=D_self, D_coll_cm2s=D_coll, D_dist_cm2s=D_dist, f=f,
        L_self=s_self * pref_L, L_coll=s_coll * pref_L, L_dist=s_dist * pref_L,
        slope_MSD_self_A2ps=s_self, slope_MSD_coll_A2ps=s_coll, slope_C_dist_A2ps=s_dist,
        sigma_NE_mScm=sig_NE * 1e3, sigma_true_mScm=sig_true * 1e3,
        haven_ratio=(D_self / D_coll if D_coll else np.nan),
    )
    return vals, curves


def analyse(args):
    """One trajectory -> one row (and one curve set) per species present."""
    system, T = args
    g = glob.glob(os.path.join(DATA, f"temperature={T}K",
                               f"chemical_system={system}", "*.jsonl.gz"))
    if not g:
        return [dict(system=system, temperature=T, error="no file")], {}
    try:
        tr = load_production(g[0])
        dt = tr.ps_per_frame
        cart = com_correct(tr.cart_unwrapped(), tr.symbols)
        start = max(tr.n_equil_frames, int(round(2.0 / dt)))
        if (tr.nframes - start) * dt < 8:
            start = max(int(round(2.0 / dt)), tr.nframes - int(round(20.0 / dt)))
    except Exception as e:  # noqa: BLE001
        import traceback
        return [dict(system=system, temperature=T, error=str(e),
                     tb=traceback.format_exc())], {}

    present = sorted({s for s in tr.symbols.tolist()
                      if (tr.symbols == s).sum() >= MIN_ATOMS})
    if SPECIES != "all":
        present = [s for s in present if s == SPECIES]

    rows, curveset = [], {}
    for element in present:
        try:
            r, key, cd = _one_species(tr, cart, start, dt, element, system, T)
            rows.append(r)
            curveset[key] = cd
        except Exception as e:  # noqa: BLE001
            import traceback
            rows.append(dict(system=system, temperature=T, species=element,
                             error=str(e), tb=traceback.format_exc()))
    return rows, curveset


def _one_species(tr, cart, start, dt, element, system, T):
    if True:
        sel = tr.symbols == element
        use = cart[start:, sel] - cart[start:, sel][0]

        vals, curves = _coeffs(use, dt, tr.volume, float(T), tuple(WINDOW))
        tt = curves["t"]
        total_ps = tt[-1]

        # ---- statistical error: NB blocks, every one fitted over the SAME
        # absolute lag window, so we compare like with like.  (Scaling each
        # block's window to its own length instead would make short blocks
        # sample short lags and confuse lag-dependence with noise.)
        nb = use.shape[0] // NBLOCK
        blk_ps = (nb - 1) * dt
        abs_win = (MATCH_WIN[0] * blk_ps, MATCH_WIN[1] * blk_ps)
        bl = []
        for b in range(NBLOCK):
            seg = use[b * nb:(b + 1) * nb]
            seg = seg - seg[0]
            bl.append(_abs_window_fit(seg, dt, *abs_win))
        full_m = _abs_window_fit(use, dt, *abs_win)
        for k in ("D_self_cm2s", "D_coll_cm2s", "f"):
            arr = np.array([x[k] for x in bl], float)
            vals[k + "_stat"] = float(np.nanstd(arr))
        vals["match_window_ps"] = [round(abs_win[0], 2), round(abs_win[1], 2)]
        vals["f_matched"] = full_m["f"]

        # ---- systematic error: refit the full curve over progressively later
        # lag windows.  A settled quantity gives the same answer in each.
        fs = []
        for lo, hi in DRIFT_WINDOWS:
            m = _abs_window_fit(use, dt, lo * total_ps, hi * total_ps)
            if np.isfinite(m["f"]):
                fs.append(m["f"])
        if fs:
            vals["f_sys"] = float((max(fs) - min(fs)) / 2)
            vals["f_window_min"] = float(min(fs))
            vals["f_window_max"] = float(max(fs))
            vals["f_sign_stable"] = bool(min(fs) > 0 or max(fs) < 0)
        vals["f_err"] = float(np.hypot(vals.get("f_sys", np.nan),
                                       vals.get("f_stat", np.nan)))

        # downsample curves for plotting
        idx = np.unique(np.linspace(0, len(curves["t"]) - 1, NCURVE).astype(int))
        cd = {k: curves[k][idx].tolist() for k in curves}
        cd["f_running"] = (np.where(np.abs(curves["MSD_self"][idx]) > 0,
                                    curves["C_dist"][idx] / curves["MSD_self"][idx],
                                    np.nan)).tolist()

        row = dict(system=system, temperature=T, species=element,
                   formula=_formula(tr.symbols),
                   production_ps=round((use.shape[0] - 1) * dt, 1),
                   window=list(WINDOW), error="", **{k: _num(v) for k, v in vals.items()})
        return row, f"{system}|{T}|{element}", cd


def _formula(sym):
    from collections import Counter
    c = Counter(sym.tolist())
    return "".join(f"{k}{c[k]}" for k in sorted(c))


def _num(v):
    try:
        v = float(v)
        return v if np.isfinite(v) else None
    except (TypeError, ValueError):
        return v


WINDOW = [0.2, 0.6]


def main():
    global WINDOW, SPECIES
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=min(32, len(os.sched_getaffinity(0))))
    ap.add_argument("--window", type=float, nargs=2, default=[0.2, 0.6])
    ap.add_argument("--species", default="all",
                    help='"all" (default) or a single element symbol, e.g. Li')
    ap.add_argument("--systems", default="shortlist",
                    help='"shortlist" (results/shortlist_electrolyte.csv, default), '
                         '"shortlist_cb" (the charge-balanced one), "targets" '
                         '(md_targets.txt), or a comma-separated list of systems')
    ap.add_argument("--merge", action="store_true",
                    help="merge into the existing onsager.csv / onsager_curves.json.gz "
                         "instead of replacing them; rows for the systems being "
                         "recomputed are dropped first")
    args = ap.parse_args()
    WINDOW = args.window
    SPECIES = args.species

    import pandas as pd
    if args.systems == "shortlist":
        systems = list(pd.read_csv(os.path.join(RES, "shortlist_electrolyte.csv"))["system"])
    elif args.systems == "shortlist_cb":
        systems = list(pd.read_csv(os.path.join(RES, "shortlist_electrolyte_cb.csv"))["system"])
    elif args.systems == "targets":
        from fetch_ncsd import read_targets
        systems = read_targets()
    else:
        systems = [x.strip() for x in args.systems.split(",") if x.strip()]
    jobs = [(s, T) for s in systems for T in TEMPS]
    print(f"{len(jobs)} (system, T) pairs on {args.procs} procs, "
          f"window {WINDOW}, species={SPECIES}", flush=True)

    rows, curves = [], {}
    with Pool(args.procs) as p:
        for i, (rws, cvs) in enumerate(p.imap_unordered(analyse, jobs)):
            rows.extend(rws)
            curves.update(cvs)
            for r in rws:
                if r.get("error"):
                    print("  ERR", r["system"], r["temperature"],
                          r.get("species", ""), r["error"], flush=True)
            if (i + 1) % 20 == 0:
                print(f"  {i + 1}/{len(jobs)}", flush=True)

    df = pd.DataFrame(rows)
    if args.merge:
        # keep everything already computed for systems we did not just rerun
        old_csv = os.path.join(RES, "onsager.csv")
        old_gz = os.path.join(RES, "onsager_curves.json.gz")
        if os.path.exists(old_csv):
            old = pd.read_csv(old_csv)
            keep = old[~old.system.isin(systems)]
            df = pd.concat([keep, df], ignore_index=True)
            print(f"merged: kept {len(keep)} existing rows for "
                  f"{keep.system.nunique()} systems, added {len(rows)} new")
        if os.path.exists(old_gz):
            with gzip.open(old_gz, "rt") as fh:
                old_curves = json.load(fh)
            old_curves = {k: v for k, v in old_curves.items()
                          if k.split("|")[0] not in systems}
            old_curves.update(curves)
            curves = old_curves
    df = df.sort_values(["system", "temperature", "species"])
    cols = ["system", "temperature", "species", "formula", "n_atoms_species",
            "production_ps", "window",
            "D_self_cm2s", "D_self_cm2s_stat", "D_coll_cm2s", "D_coll_cm2s_stat",
            "D_dist_cm2s", "f", "f_err", "f_sys", "f_stat", "f_sign_stable",
            "f_window_min", "f_window_max", "f_matched", "match_window_ps",
            "L_self", "L_coll", "L_dist",
            "slope_MSD_self_A2ps", "slope_MSD_coll_A2ps", "slope_C_dist_A2ps",
            "sigma_NE_mScm", "sigma_true_mScm", "haven_ratio", "error"]
    df = df.reindex(columns=[c for c in cols if c in df.columns]
                    + [c for c in df.columns if c not in cols])
    df.to_csv(os.path.join(RES, "onsager.csv"), index=False)
    with gzip.open(os.path.join(RES, "onsager_curves.json.gz"), "wt") as fh:
        json.dump(curves, fh)
    ok = df[df.error.fillna("") == ""]
    print(f"wrote results/onsager.csv ({len(ok)} rows, "
          f"{ok.species.nunique()} species, {len(df) - len(ok)} errors)")
    print()
    print("median over each species:")
    g = ok.groupby("species").agg(
        n=("system", "size"),
        D_self=("D_self_cm2s", "median"),
        D_self_relerr=("D_self_cm2s_stat", lambda x: np.nan),
        f=("f", "median"), f_err=("f_err", "median"),
        sign_stable=("f_sign_stable", "mean"))
    g["D_self_relerr"] = (ok.assign(r=ok.D_self_cm2s_stat / ok.D_self_cm2s.abs())
                          .groupby("species").r.median())
    print(g.sort_values("n", ascending=False).round(4).to_string())


if __name__ == "__main__":
    main()
