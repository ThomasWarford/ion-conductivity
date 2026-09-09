"""Build the DATA blob embedded in results/aimd_targets.html.

Everything comes from the already-computed screen / Onsager outputs -- no
trajectory is re-read:

    results/screen.csv            per (system, T, species) diffusivity metrics
    results/onsager.csv           D_self / D_coll / D_dist / f per system x T
    results/onsager_curves.json.gz   the MSD_self / MSD_coll / C_dist curves

Only charge-balanced compositions are considered (charge_balance.py); the
previous version of the artifact recommended Br-Li = Li3Br and Li-S = Li3S,
which are not compounds.

    python aimd_targets_data.py --list                 # rank the candidates
    python aimd_targets_data.py --picks Cl-Li,Li-N,... # write results/aimd_targets_data.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import os

import numpy as np
import pandas as pd

from charge_balance import classify

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
TEMPS = [1000, 1500, 2000, 2500]
# absolute lag windows in ps -- fixed in ps (not as a fraction of the block) so
# that the spread measures lag-dependence, not block length.  See summary.md.
WINDOWS = [(2, 6), (3, 8), (5, 12), (7, 18), (10, 24)]


def _slope_r2(t, y, lo, hi):
    m = (t >= lo) & (t <= hi)
    if m.sum() < 4:
        return None, None
    p = np.polyfit(t[m], y[m], 1)
    resid = y[m] - np.polyval(p, t[m])
    ss = ((y[m] - y[m].mean()) ** 2).sum()
    return float(p[0]), (1.0 - float((resid ** 2).sum()) / ss if ss > 0 else None)


def drift_for(curves, system, T):
    c = curves.get(f"{system}|{T}|Li")
    if c is None:
        return None
    t = np.asarray(c["t"], float)
    self_, dist = np.asarray(c["MSD_self"], float), np.asarray(c["C_dist"], float)
    fs, r2d, r2s = [], None, None
    for lo, hi in WINDOWS:
        s_self, r2self = _slope_r2(t, self_, lo, hi)
        s_dist, r2dist = _slope_r2(t, dist, lo, hi)
        fs.append(None if not s_self else round(s_dist / s_self, 3))
        if (lo, hi) == (5, 12):
            r2d, r2s = r2dist, r2self
    ok = [v for v in fs if v is not None]
    return {
        "f_by_window": fs,
        "drift": round(max(ok) - min(ok), 3) if len(ok) > 1 else None,
        "sign_stable": bool(ok) and (min(ok) > 0 or max(ok) < 0),
        "r2_dist": None if r2d is None else round(r2d, 3),
        "r2_self": None if r2s is None else round(r2s, 4),
    }


def load():
    scr = pd.read_csv(os.path.join(RES, "screen.csv"))
    ons = pd.read_csv(os.path.join(RES, "onsager.csv"))
    ons = ons[(ons.get("error").fillna("") == "") & (ons.species == "Li")].copy()
    ons = ons[ons.formula.map(lambda f: classify(f).status == "balanced")].copy()
    with gzip.open(os.path.join(RES, "onsager_curves.json.gz")) as fh:
        curves = json.load(fh)
    return scr, ons, curves


def candidate_table(scr, ons, curves):
    li = scr[scr.species == "Li"]
    rows = []
    for s, g in ons.groupby("system"):
        d = {int(r.temperature): r for r in g.itertuples()}
        dr = {T: drift_for(curves, s, T) for T in TEMPS}
        stable = sum(1 for T in TEMPS if dr[T] and dr[T]["sign_stable"])
        drifts = [dr[T]["drift"] for T in TEMPS if dr[T] and dr[T]["drift"]]
        li_s = li[li.system == s]
        rows.append(dict(
            system=s, formula=g.formula.iloc[0], reduced=classify(g.formula.iloc[0]).reduced,
            D1500=d[1500].D_self_cm2s if 1500 in d else np.nan,
            f1500=d[1500].f if 1500 in d else np.nan,
            ferr1500=d[1500].f_err if 1500 in d else np.nan,
            sign_stable_T=stable, med_drift=round(float(np.median(drifts)), 2) if drifts else np.nan,
            med_block_cv=round(float(li_s.block_cv.median()), 3),
            n_Li=int(li_s.n_species.iloc[0]), natoms=int(li_s.natoms.iloc[0])))
    t = pd.DataFrame(rows).sort_values("D1500", ascending=False)
    return t


def build(scr, ons, curves, picks, meta):
    li = scr[scr.species == "Li"]
    out = {"picks": [], "background": [], "drift": {"windows": WINDOWS, "data": {}}}

    for s in picks:
        g = ons[ons.system == s]
        d = {int(r.temperature): r for r in g.itertuples()}
        li_s = li[li.system == s]
        by_T = lambda col: {str(T): (None if T not in d else
                                     round(float(getattr(d[T], col)), 10)) for T in TEMPS}
        ea = np.nan
        gg = li_s[li_s.D_einstein_cm2s > 0]
        if gg.temperature.nunique() >= 3:
            p = np.polyfit(1 / gg.temperature, np.log(gg.D_einstein_cm2s), 1)
            ea = round(float(-p[0] * 8.617333e-5), 3)
        pick = dict(meta[s])
        pick.update(
            system=s, formula=g.formula.iloc[0], reduced=classify(g.formula.iloc[0]).reduced,
            natoms=int(li_s.natoms.iloc[0]), n_Li=int(li_s.n_species.iloc[0]), Ea=ea,
            D_self=by_T("D_self_cm2s"),
            t_conv={str(T): (None if r.empty else float(r.t_run_converge_ps.iloc[0]))
                    for T in TEMPS for r in [li_s[li_s.temperature == T]]},
            block_cv={str(T): (None if r.empty else round(float(r.block_cv.iloc[0]), 3))
                      for T in TEMPS for r in [li_s[li_s.temperature == T]]},
            f=by_T("f"), f_err=by_T("f_err"),
            curves={str(T): _pick_curve(curves, s, T) for T in TEMPS})
        out["picks"].append(pick)
        out["drift"]["data"][s] = {str(T): drift_for(curves, s, T) for T in TEMPS}

    for s, g in ons[ons.temperature == 1500].groupby("system"):
        r = g.iloc[0]
        out["background"].append(dict(system=s, reduced=classify(r.formula).reduced,
                                      D=round(float(r.D_self_cm2s), 10),
                                      f=float(r.f), ferr=float(r.f_err),
                                      pick=s in picks))
    return out


def _pick_curve(curves, s, T, n=60):
    c = curves.get(f"{s}|{T}|Li")
    if c is None:
        return None
    idx = np.unique(np.linspace(0, len(c["t"]) - 1, n).astype(int))
    out = {k: [round(float(np.asarray(c[k], float)[i]), 4) for i in idx]
           for k in ("t", "MSD_self", "MSD_coll", "C_dist")}
    # running f(t), the key the artifact's chF() reads
    fr = np.asarray(c["f_running"], float)
    out["f"] = [None if not np.isfinite(fr[i]) else round(float(fr[i]), 4) for i in idx]
    return out


META = {
    "Cl-Li": dict(
        role="the calibration point", regime="f not resolved",
        analogue="molten LiCl -- decades of conductivity, tracer-D and Haven-ratio data",
        why="The only balanced alkali halide in the set and the fastest Li in it. Self-D is "
            "converged by 5 ps at every temperature, and molten LiCl has real experimental "
            "transport numbers to check the absolute value against rather than another "
            "simulation. Its correlation factor is NOT resolved (f moves by 2.5 across fit "
            "windows at 1500 K) -- this system is here for D_self, not for f."),
    "Li-N": dict(
        role="the one f we can nearly measure", regime="f < 0",
        analogue="Li3N -- a real, well characterised fast Li-ion conductor",
        why="The best-behaved correlation factor anywhere in the balanced set: negative at all "
            "four temperatures and in all five fit windows, drifting by only 0.24-0.55 against a "
            "set median of 0.96, with the tightest block CV (0.04-0.12) and C_dist fits at "
            "R2 = 0.96-0.998. If any of these systems yields an ab-initio f, it is this one."),
    "Bi-Li-S": dict(
        role="the sulfide, and the hardest case", regime="f > 0, but T-dependent",
        analogue="Li3BiS3 -- thio-bismuthate, sulfide-electrolyte chemistry",
        why="Second-fastest Li in the balanced set with clean self-D (block CV 0.07-0.14) and "
            "smooth correlation curves (R2 0.97-0.99). Its f holds one sign within each "
            "temperature but changes sign between them -- positive at 1000/1500 K, negative at "
            "2000 K. Either that is real physics or the curves have not straightened; only "
            "longer trajectories can tell, which makes it the sharpest test of the plan."),
    "Cu-Li-S": dict(
        role="two mobile cations", regime="f >> 0",
        analogue="Li-Cu thiolate / Li3CuS2",
        why="The only balanced system with a second genuinely mobile cation -- Cu diffuses at "
            "6.6e-5 cm2/s at 1500 K, half the Li rate -- so it is the one place we can look at "
            "cation-cation cross terms rather than only Li-Li. It also sits at the extreme "
            "positive-correlation end (f = +3.2 at 1500 K), so it stretches the potential "
            "furthest in that direction."),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--picks", default="Cl-Li,Li-N,Bi-Li-S,Cu-Li-S")
    ap.add_argument("--out", default=os.path.join(RES, "aimd_targets_data.json"))
    a = ap.parse_args()
    scr, ons, curves = load()
    if a.list:
        print(candidate_table(scr, ons, curves).to_string(index=False))
        return
    picks = a.picks.split(",")
    data = build(scr, ons, curves, picks, META)
    json.dump(data, open(a.out, "w"), separators=(",", ":"))
    print(f"wrote {a.out}  ({os.path.getsize(a.out)/1024:.0f} KB)  picks={picks}")


if __name__ == "__main__":
    main()
