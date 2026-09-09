"""MSD, Einstein diffusivity, and convergence diagnostics for NCSD trajectories.

Answers two questions per (system, temperature, species):
  1. How large is the diffusivity?          -> D_einstein  (cm^2/s)
  2. How fast does it converge in the run?   -> t_converge_ps, beta, block_cv, ...

Conventions
-----------
* self-diffusion Einstein relation in 3D:   MSD_self(t) -> 6 D t
* D reported in cm^2/s.  1 Angstrom^2/ps = 1e-4 cm^2/s.
* The MP amorphous_diffusivity table fits slope/2 (a 1D relation applied to the
  full 3D MSD) so its values are 3x our D_einstein.  D_mp_conv = 3 * D_einstein
  is reported for direct comparison (see project memory).
"""
from __future__ import annotations

import numpy as np

ANG2_PS_TO_CM2_S = 1e-4


def masses(symbols):
    from ase.data import atomic_masses, atomic_numbers
    return np.array([atomic_masses[atomic_numbers[s]] for s in symbols])


def com_correct(cart, symbols):
    """Remove mass-weighted centre-of-mass drift at every frame."""
    m = masses(symbols)
    com = (cart * m[None, :, None]).sum(axis=1) / m.sum()
    return cart - com[:, None, :]


def msd_fft(x):
    """Vectorised FFT (FCA) MSD over all lag times.

    x : (nframes, ...) array.  Returns (nframes, ...) MSD averaged over time
    origins for each trailing-axis series.  MSD is summed over the trailing
    components? No -- it is computed per component; caller sums.
    """
    n = x.shape[0]
    rest = x.shape[1:]
    x2 = x.reshape(n, -1)
    # S2: autocorrelation per column via FFT
    f = np.fft.rfft(x2, 2 * n, axis=0)
    ac = np.fft.irfft(f * np.conjugate(f), axis=0)[:n].real
    ac /= np.arange(n, 0, -1)[:, None]
    # S1 recursion
    sq = x2 ** 2
    total = 2.0 * sq.sum(axis=0)
    s1 = np.empty_like(ac)
    q = total.copy()
    for m in range(n):
        if m > 0:
            q -= sq[m - 1] + sq[n - m]
        s1[m] = q / (n - m)
    msd = s1 - 2.0 * ac
    return msd.reshape((n,) + rest)


def self_msd(cart):
    """cart (nframes, natoms, 3) -> (nframes,) self-MSD averaged over atoms and
    time origins (sum over x,y,z)."""
    per = msd_fft(cart)              # (nframes, natoms, 3)
    return per.sum(axis=2).mean(axis=1)


def collective_msd(cart):
    """<|sum_i dr_i(t)|^2> summed over x,y,z.  /(6 N t) -> D_charge."""
    s = cart.sum(axis=1)            # (nframes, 3)
    return msd_fft(s).sum(axis=1)


def fit_D(t, msd, lo=0.2, hi=0.8):
    """OLS slope of msd vs t over fractional lag window [lo*tmax, hi*tmax].
    Returns (D_cm2_s, slope_ang2_ps)."""
    t = np.asarray(t)[:len(msd)]
    msd = np.asarray(msd)[:len(t)]
    if len(t) < 5:
        return np.nan, np.nan
    tmax = t[-1]
    m = (t >= lo * tmax) & (t <= hi * tmax) & (t > 0)
    if m.sum() < 5:
        return np.nan, np.nan
    A = np.vstack([t[m], np.ones(m.sum())]).T
    slope = np.linalg.lstsq(A, msd[m], rcond=None)[0][0]
    return slope / 6.0 * ANG2_PS_TO_CM2_S, slope


def diffusive_exponent(t, msd, lo=0.4, hi=0.9):
    """d ln MSD / d ln t over the late window; ~1 == diffusive, <1 == sub-diffusive."""
    tmax = t[-1]
    m = (t >= lo * tmax) & (t <= hi * tmax) & (msd > 0) & (t > 0)
    if m.sum() < 5:
        return np.nan
    return float(np.polyfit(np.log(t[m]), np.log(msd[m]), 1)[0])


def _D_fit_window(t, msd, tau_max, tau_min_frac=0.2):
    m = (t >= tau_min_frac * tau_max) & (t <= tau_max) & (t > 0)
    if m.sum() < 5:
        return np.nan
    A = np.vstack([t[m], np.ones(m.sum())]).T
    slope = np.linalg.lstsq(A, msd[m], rcond=None)[0][0]
    return slope / 6.0 * ANG2_PS_TO_CM2_S


def prefix_convergence(sub, dt, prefixes_ps, tol=0.25):
    """Recompute D from the first P ps of the trajectory for P in prefixes_ps.

    Directly answers "if the run were only P ps long, would D be right?".
    Returns (D_prefix dict, t_run_converge_ps) where t_run_converge_ps is the
    shortest prefix whose D -- and every longer prefix's D -- is within `tol`
    (relative) of the full-run D.
    """
    n = sub.shape[0]
    T = (n - 1) * dt
    D_full = fit_D(np.arange(n) * dt, self_msd(sub - sub[0]), 0.2, 0.8)[0]
    Dp = {}
    for P in prefixes_ps:
        if P > T * 1.001:
            continue
        k = min(max(int(round(P / dt)) + 1, 20), n)
        s = sub[:k] - sub[:k][0]
        Dp[float(P)] = fit_D(np.arange(s.shape[0]) * dt, self_msd(s), 0.2, 0.8)[0]
    t_run = np.nan
    if np.isfinite(D_full) and D_full > 0:
        Ps = sorted(Dp)
        for i, P in enumerate(Ps):
            if all(np.isfinite(Dp[p]) and abs(Dp[p] - D_full) / D_full <= tol
                   for p in Ps[i:]):
                t_run = P
                break
    return Dp, t_run


def convergence_metrics(t, msd, sub, dt, windows_ps, tol=0.2):
    """t, msd: full-trajectory lag axis (ps) and self-MSD (Angstrom^2), already
    computed once.  sub: (nframes, n_species, 3) com-corrected displacements, for
    the independent-block statistics.

    Returns:
      D_full          D fitted over lag window [0.2 T, 0.8 T]           (cm^2/s)
      D_window        {tau_max: D fitted over lags [0.2 tau_max, tau_max]}
      t_converge_ps   smallest tau_max s.t. every larger tau_max is within
                      `tol` (relative) of D_full  -- "the MSD has to reach this
                      lag before the fitted D stops moving"
      block_cv        coeff. of variation of D over 4 equal contiguous
                      sub-trajectories, each analysed on its own clock
      D_half_ratio    D(2nd half) / D(1st half)
    """
    T = t[-1]
    D_full = fit_D(t, msd, 0.2, 0.8)[0]

    D_window = {}
    for W in windows_ps:
        if W > T * 1.001:
            continue
        D_window[float(W)] = _D_fit_window(t, msd, min(W, T))

    t_conv = np.nan
    if np.isfinite(D_full) and D_full > 0:
        Ws = sorted(D_window)
        for i, W in enumerate(Ws):
            tail = [D_window[w] for w in Ws[i:]]
            if tail and all(np.isfinite(d) and abs(d - D_full) / D_full <= tol
                            for d in tail):
                t_conv = W
                break

    n = sub.shape[0]
    bl = n // 4
    Ds = []
    for b in range(4):
        seg = sub[b * bl:(b + 1) * bl]
        seg = seg - seg[0]
        tb = np.arange(seg.shape[0]) * dt
        Ds.append(fit_D(tb, self_msd(seg), 0.2, 0.8)[0])
    Ds = np.array(Ds, float)
    block_cv = (float(np.nanstd(Ds) / np.nanmean(Ds))
                if np.isfinite(np.nanmean(Ds)) and np.nanmean(Ds) else np.nan)

    h = n // 2
    t1 = np.arange(h) * dt
    d1 = fit_D(t1, self_msd(sub[:h] - sub[:h][0]), 0.2, 0.8)[0]
    d2 = fit_D(t1, self_msd(sub[h:2 * h] - sub[h:2 * h][0]), 0.2, 0.8)[0]
    half_ratio = float(d2 / d1) if d1 else np.nan

    return dict(D_full=D_full, D_window=D_window, t_converge_ps=t_conv,
                block_cv=block_cv, D_half_ratio=half_ratio, D_blocks=Ds.tolist())


def nn_distance(frac0, lattice, sel, other=None):
    """First-frame nearest-neighbour distance for atoms `sel` to atoms `other`
    (default: all other atoms).  Minimum-image, in Angstrom."""
    if other is None:
        other = np.ones(len(frac0), bool)
    fs = frac0[sel]
    fo = frac0[other]
    d = fs[:, None, :] - fo[None, :, :]
    d -= np.round(d)
    cart = d @ lattice
    r = np.sqrt((cart ** 2).sum(-1))
    r[r < 1e-6] = np.inf
    return float(np.median(r.min(axis=1)))
