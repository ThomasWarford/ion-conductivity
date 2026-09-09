"""Find chemical-trend series in the NCSD screen: keep every element but one
fixed, swap that one for a same-group congener.

A trend series is only useful if every member is a real compound, so this
reports charge balance (charge_balance.py) alongside the diffusivity and
convergence metrics.  Two grades of series:

  matched     every member reduces to the same stoichiometry with the swapped
              element renamed (Li3N / Li3P) -- an isovalent substitution and
              nothing else changes, which is what a clean trend needs
  unmatched   same elements-minus-one, but the packing used a different
              stoichiometry (Li3BiS3 / Li5SbS4) -- still a chemical trend, with
              a composition change confounded into it

It also reports BLOCKED series: a balanced member whose congener exists in the
dataset only at a composition that is not a compound.  Those are the cases where
a new ~100-atom cell at the balanced stoichiometry would restore the trend.

    python group_series.py
    python group_series.py --min-D 1e-5     # only series that are fast enough
"""
from __future__ import annotations

import argparse
import os
from collections import defaultdict
from functools import reduce
from math import gcd

import numpy as np
import pandas as pd

from charge_balance import classify, parse_formula

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

# main-group and the transition-metal triads that actually appear in the screen
GROUPS = {
    "1  alkali": "Li Na K Rb Cs", "2  alkaline earth": "Mg Ca Sr Ba",
    "13 triel": "B Al Ga In Tl", "14 tetrel": "C Si Ge Sn Pb",
    "15 pnictogen": "N P As Sb Bi", "16 chalcogen": "O S Se Te",
    "17 halogen": "F Cl Br I",
    "11 coinage": "Cu Ag Au", "12": "Zn Cd Hg",
    "3": "Sc Y La", "4": "Ti Zr Hf", "5": "V Nb Ta", "6": "Cr Mo W",
    "7": "Mn Tc Re", "8": "Fe Ru Os", "9": "Co Rh Ir", "10": "Ni Pd Pt",
}
GROUP_OF = {e: g for g, els in GROUPS.items() for e in els.split()}


def stoich(counts, rename_from=None):
    """Composition reduced to lowest terms, with one element anonymised as '*'."""
    c = defaultdict(int)
    for k, v in counts.items():
        c["*" if k == rename_from else k] += v
    g = reduce(gcd, c.values())
    return tuple(sorted((k, v // g) for k, v in c.items()))


def load():
    scr = pd.read_csv(os.path.join(RES, "screen.csv"))
    li = scr[scr.species == "Li"]
    form = scr[scr.species == "*all*"].groupby("system").formula.last()
    D = {int(T): li[li.temperature == T].set_index("system").D_einstein_cm2s
         for T in sorted(li.temperature.unique())}
    cv = li[li.temperature == 1500].set_index("system").block_cv
    nat = li.groupby("system").natoms.last()

    from report import converged            # the screen's own convergence flags
    nconv = li.assign(ok=converged(li)).groupby("system").ok.sum()

    info = {}
    for s, f in form.items():
        v = classify(f)
        info[s] = dict(system=s, formula=f, reduced=v.reduced, cb=v.status,
                       q=v.q_per_atom, counts=parse_formula(f), states=v.states,
                       natoms=int(nat[s]),
                       D1000=D[1000].get(s, np.nan), D1500=D[1500].get(s, np.nan),
                       cv=cv.get(s, np.nan), n_conv_T=int(nconv.get(s, 0)))
    return info


def families(info, isovalent=False):
    """Group systems by (fixed elements, swap axis).

    The swap axis is the periodic group by default.  With isovalent=True it is
    the oxidation state the balanced assignment gives the swapped element, which
    is what "isovalent substitution" actually means -- group is only a proxy for
    it, and misses Ga(3+) -> Sc/Y/La(3+) or Ti(4+) -> Sn(4+).  Only balanced
    systems have an assignment, so that mode is balanced-only by construction.
    """
    fam = defaultdict(list)
    for s, d in info.items():
        for e in d["counts"]:
            if isovalent:
                st = (d["states"] or {}).get(e)
                axis = None if st is None else f"{st:+g}"
            else:
                axis = GROUP_OF.get(e)
            if axis is None:
                continue
            fam[(frozenset(d["counts"]) - {e}, axis)].append((s, e))
    return {k: v for k, v in fam.items() if len(v) > 1}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-D", type=float, default=0.0,
                    help="require D_Li(1000 K) above this in every member")
    ap.add_argument("--isovalent", action="store_true",
                    help="swap axis is the assigned oxidation state, not the periodic "
                         "group (a superset: catches Ga3+ -> Sc/Y/La3+ etc.)")
    a = ap.parse_args()
    info = load()
    fam = families(info, isovalent=a.isovalent)

    usable, blocked = [], []
    for (fixed, g), members in sorted(fam.items(), key=lambda kv: -len(kv[1])):
        bal = [(s, e) for s, e in members if info[s]["cb"] == "balanced"]
        bad = [(s, e) for s, e in members if info[s]["cb"] != "balanced"]
        if len(bal) < 2:
            if bal and bad:
                blocked.append((fixed, g, bal, bad))
            continue
        if a.min_D and any(not (info[s]["D1000"] >= a.min_D) for s, _ in bal):
            continue
        keys = {stoich(info[s]["counts"], e) for s, e in bal}
        usable.append((fixed, g, bal, bad, len(keys) == 1))

    usable.sort(key=lambda r: (-len(r[2]), not r[4]))
    print(f"{len(usable)} usable trend series (>=2 charge-balanced members)\n")
    for fixed, g, bal, bad, matched in usable:
        head = " + ".join(sorted(fixed)) or "(binary)"
        axis = f"oxidation state {g}" if a.isovalent else f"group-{g}"
        print(f"[{'MATCHED ' if matched else 'unmatched'}] {head}  ·  swap the {axis} element"
              f"  ({len(bal)} balanced members)")
        for s, e in sorted(bal, key=lambda x: -info[x[0]]["D1000"]):
            d = info[s]
            print(f"     {e:2s}  {s:12s} {d['formula']:16s} {d['reduced']:11s} "
                  f"{d['natoms']:3d} at  D1000={d['D1000']:.2e}  D1500={d['D1500']:.2e}  "
                  f"cv={d['cv']:.2f}  conv@{d['n_conv_T']}T")
        for s, e in sorted(bad):
            d = info[s]
            q = f"q/atom {d['q']:+.2f}" if d["q"] is not None else "no anion"
            print(f"     {e:2s}  {s:12s} {d['formula']:16s} {d['reduced']:11s} "
                  f"-- EXCLUDED, {q}")
        print()

    print(f"\n{len(blocked)} BLOCKED series -- one balanced member, congeners present "
          f"only at non-compound compositions:\n")
    for fixed, g, bal, bad in blocked:
        head = " + ".join(sorted(fixed)) or "(binary)"
        s, e = bal[0]
        print(f"  {head} · {g}:  have {info[s]['reduced']} ({s}); "
              f"blocked " + ", ".join(
                  f"{info[b]['reduced']} ({b}, "
                  + (f"q/atom {info[b]['q']:+.2f}" if info[b]['q'] is not None else "alloy") + ")"
                  for b, _ in sorted(bad)))


if __name__ == "__main__":
    main()
