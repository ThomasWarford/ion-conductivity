"""Per-file analysis worker: NCSD trajectory file -> list of result rows
(one row per diffusing species, plus an 'all-atoms' row).
"""
from __future__ import annotations

import traceback

import numpy as np

from ncsd_traj import load_production
from diffusivity import (com_correct, self_msd, collective_msd, fit_D,
                         diffusive_exponent, convergence_metrics,
                         prefix_convergence, nn_distance)

WINDOWS_PS = (1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 12.5, 15, 20, 25, 30, 40, 50, 60)
PREFIX_PS = (5, 7.5, 10, 15, 20, 25, 30, 40)
MIN_SPECIES_ATOMS = 4


def _energy_drift(energy, dt, natoms):
    """meV/atom/ps slope of the potential energy over the usable trajectory."""
    e = energy
    ok = np.isfinite(e)
    if ok.sum() < 10:
        return np.nan
    t = np.arange(len(e))[ok] * dt
    slope = np.polyfit(t, e[ok], 1)[0]        # eV / ps
    return slope / natoms * 1000.0            # meV / atom / ps


def analyse_file(path):
    try:
        tr = load_production(path)
    except Exception as e:  # noqa: BLE001
        return [dict(path=path, error=f"load: {e}")]

    dt = tr.ps_per_frame
    cart = com_correct(tr.cart_unwrapped(), tr.symbols)

    start = max(tr.n_equil_frames, int(round(2.0 / dt)))
    if (tr.nframes - start) * dt < 8.0:                    # keep >= 8 ps
        start = max(int(round(2.0 / dt)), tr.nframes - int(round(20.0 / dt)))
    use = cart[start:] - cart[start:][0]
    frac0 = tr.frac[start] if start < tr.nframes else tr.frac[-1]
    n = use.shape[0]
    t = np.arange(n) * dt
    prod_ps = (n - 1) * dt

    base = dict(
        system=tr.system, temperature=tr.temperature, contribs_id=tr.contribs_id,
        formula=_formula(tr.symbols), natoms=tr.natoms, volume_A3=tr.volume,
        total_ps=round(tr.total_ps, 2), equil_ps=round(tr.n_equil_frames * dt, 2),
        production_ps=round(prod_ps, 2),
        energy_drift_meV_atom_ps=round(_energy_drift(tr.energy[start:], dt, tr.natoms), 4),
        error="",
    )

    rows = []
    species = sorted(set(tr.symbols.tolist()))
    targets = [s for s in species if (tr.symbols == s).sum() >= MIN_SPECIES_ATOMS]
    for el in targets + ["*all*"]:
        sel = np.ones(tr.natoms, bool) if el == "*all*" else (tr.symbols == el)
        nsp = int(sel.sum())
        sub = use[:, sel]

        msd = self_msd(sub)
        D, slope = fit_D(t, msd)
        beta = diffusive_exponent(t, msd)
        cm = convergence_metrics(t, msd, sub, dt, WINDOWS_PS, tol=0.2)

        # collective / charge channel (self vs collective => cross-correlations).
        # Only meaningful for a genuine sub-lattice of a multi-component system.
        n_elements = len(set(tr.symbols.tolist()))
        if el == "*all*" or n_elements < 2 or nsp == tr.natoms:
            Dch, cbeta, haven = np.nan, np.nan, np.nan
        else:
            cmsd = collective_msd(sub)
            Dch, _ = fit_D(t, cmsd / nsp)
            cbeta = diffusive_exponent(t, cmsd)
            haven = D / Dch if (np.isfinite(Dch) and Dch > 0.05 * abs(D)) else np.nan

        dnn = nn_distance(frac0, tr.lattice, sel) if el != "*all*" else nn_distance(
            frac0, tr.lattice, np.ones(tr.natoms, bool))
        half_box = 0.5 * tr.volume ** (1 / 3)

        wins = cm["D_window"]
        big = [v for k, v in wins.items() if k >= 5 and np.isfinite(v)]
        spread = (max(big) - min(big)) / np.median(big) if big and np.median(big) else np.nan

        Dpref, t_run = prefix_convergence(sub, dt, PREFIX_PS, tol=0.25)

        row = dict(base)
        row.update(
            species=el, n_species=nsp,
            D_einstein_cm2s=_sig(D), D_mp_conv_cm2s=_sig(3 * D) if np.isfinite(D) else np.nan,
            msd_slope_A2ps=_sig(slope), beta_diffusive=_r(beta, 3),
            t_lag_stable_ps=cm["t_converge_ps"],
            t_run_converge_ps=t_run,
            Dwin_spread=_r(spread, 3), block_cv=_r(cm["block_cv"], 3),
            D_half_ratio=_r(cm["D_half_ratio"], 3),
            msd_tmax_A2=_r(float(msd[-1]), 1),
            msd_tmax_per_dNN2=_r(float(msd[-1]) / dnn ** 2, 1) if dnn else np.nan,
            disp_rms_tmax_A=_r(float(np.sqrt(msd[-1])), 2),
            disp_over_halfbox=_r(float(np.sqrt(msd[-1])) / half_box, 2),
            dNN_A=_r(dnn, 2),
            D_charge_cm2s=_sig(Dch), haven_ratio=_r(haven, 2), beta_charge=_r(cbeta, 3),
        )
        for W in WINDOWS_PS:
            row[f"Dlag_{W}"] = _sig(wins.get(float(W), np.nan))
        for P in PREFIX_PS:
            row[f"Dpre_{P}"] = _sig(Dpref.get(float(P), np.nan))
        rows.append(row)
    return rows


def _formula(symbols):
    from collections import Counter
    c = Counter(symbols.tolist())
    return "".join(f"{k}{c[k]}" for k in sorted(c))


def _sig(x, n=4):
    return float(f"{x:.{n}g}") if x is not None and np.isfinite(x) else np.nan


def _r(x, n):
    return round(float(x), n) if x is not None and np.isfinite(x) else np.nan


if __name__ == "__main__":
    import json
    import sys
    for p in sys.argv[1:]:
        for r in analyse_file(p):
            print(json.dumps(r))
