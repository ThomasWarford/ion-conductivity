"""Is an NCSD amorphous cell a charge-balanced ionic compound?

The MP amorphous cells are PACKMOL packings at a target composition, and many of
those compositions are not stoichiometric salts.  `Br-Li` is not LiBr -- it is
Br25Li75, i.e. Li3Br, which carries +2 per formula unit.  Physically that cell is
not "charged" (VASP runs a neutral cell); it is Li metal with Br dissolved in it,
and the Li mobility it reports is metallic self-diffusion, not ionic conduction.
For picking solid-electrolyte validation targets those systems are meaningless,
so this module flags them.

Method.  Each element gets a curated list of common oxidation states (cation
states, plus an anion state for the elements that form anions).  The single most
electronegative element present is FORCED into its anion state -- that is what
makes the assignment a salt rather than an arbitrary redox bookkeeping, and it is
what stops Ag(+3)/Si(-4) style "solutions" for what is really an alloy.  Every
other element ranges over all of its allowed states.  A composition is balanced
if some combination sums to zero.

Deliberately *permissive* on oxidation state (Sn is allowed +2 and +4, S is
allowed -2/+4/+6 so that sulfates pass, ...): the aim is to catch compositions no
oxidation-state assignment can rescue, not to pin down the true valence.

    from charge_balance import classify
    classify("Br25Li75")
    -> Verdict(status='unbalanced', q_excess=50.0, q_per_atom=0.5,
               reduced='Li3Br', states=None, anion='Br')

`status` is one of
    balanced    some assignment of common oxidation states sums to zero
    unbalanced  none does; `q_excess` is the smallest achievable |sum q|
    no_anion    no element present forms an anion -> an alloy/intermetallic,
                not an ionic compound at all (Ag25Li75 = Li3Ag)
"""
from __future__ import annotations

import itertools
from collections import namedtuple
from functools import lru_cache

# Pauling electronegativity, for picking the anion.
CHI = {
    "F": 3.98, "O": 3.44, "Cl": 3.16, "N": 3.04, "Br": 2.96, "I": 2.66,
    "S": 2.58, "C": 2.55, "Se": 2.55, "Au": 2.54, "W": 2.36, "Pb": 2.33,
    "Pt": 2.28, "Ir": 2.20, "Os": 2.20, "Pd": 2.20, "Ru": 2.20, "P": 2.19,
    "As": 2.18, "Mo": 2.16, "Te": 2.10, "Sb": 2.05, "B": 2.04, "Bi": 2.02,
    "Hg": 2.00, "Ge": 2.01, "Sn": 1.96, "Rh": 2.28, "Tc": 1.90, "Re": 1.90,
    "Ag": 1.93, "Si": 1.90, "Ni": 1.91, "Cu": 1.90, "Co": 1.88, "Fe": 1.83,
    "Ga": 1.81, "In": 1.78, "Cd": 1.69, "Cr": 1.66, "Zn": 1.65, "V": 1.63,
    "Tl": 1.62, "Al": 1.61, "Nb": 1.60, "Mn": 1.55, "Ti": 1.54, "Ta": 1.50,
    "Sc": 1.36, "Zr": 1.33, "Mg": 1.31, "Hf": 1.30, "Y": 1.22, "La": 1.10,
    "Ca": 1.00, "Li": 0.98, "Sr": 0.95, "Na": 0.93, "Ba": 0.89, "K": 0.82,
    "Rb": 0.82, "Cs": 0.79,
}

# Anion state of the elements that form one (the value forced on the most
# electronegative element present).
ANION = {
    "F": -1, "Cl": -1, "Br": -1, "I": -1,
    "O": -2, "S": -2, "Se": -2, "Te": -2,
    "N": -3, "P": -3, "As": -3, "Sb": -3, "Bi": -3,
    "C": -4, "Si": -4, "Ge": -4, "Sn": -4, "Pb": -4,   # Zintl anions
}

# Common positive oxidation states.
CATION = {
    "Li": (1,), "Na": (1,), "K": (1,), "Rb": (1,), "Cs": (1,),
    "Mg": (2,), "Ca": (2,), "Sr": (2,), "Ba": (2,), "Zn": (2,), "Cd": (2,),
    "Ag": (1,), "Cu": (1, 2), "Au": (1, 3), "Hg": (1, 2),
    "Al": (3,), "Ga": (3,), "In": (1, 3), "Tl": (1, 3), "B": (3,),
    "Sc": (3,), "Y": (3,), "La": (3,),
    "C": (4,), "Si": (4,), "Ge": (2, 4), "Sn": (2, 4), "Pb": (2, 4),
    "N": (3, 5), "P": (3, 5), "As": (3, 5), "Sb": (3, 5), "Bi": (3, 5),
    "S": (4, 6), "Se": (4, 6), "Te": (4, 6),
    "Ti": (2, 3, 4), "Zr": (4,), "Hf": (4,),
    "V": (2, 3, 4, 5), "Nb": (3, 4, 5), "Ta": (3, 4, 5),
    "Cr": (2, 3, 6), "Mo": (2, 3, 4, 5, 6), "W": (2, 3, 4, 5, 6),
    "Mn": (2, 3, 4, 6, 7), "Tc": (4, 6, 7), "Re": (4, 6, 7),
    "Fe": (2, 3), "Co": (2, 3), "Ni": (2, 3),
    "Ru": (2, 3, 4, 6, 8), "Os": (2, 3, 4, 6, 8),
    "Rh": (1, 3, 4), "Ir": (1, 3, 4), "Pd": (2, 4), "Pt": (2, 4),
}

Verdict = namedtuple("Verdict",
                     "status q_excess q_per_atom reduced states anion")


def parse_formula(formula):
    """'Ag17F51Li34' -> {'Ag': 17, 'F': 51, 'Li': 34}"""
    out, el, num = {}, "", ""
    for ch in formula + "$":
        if ch.isupper() or ch == "$":
            if el:
                out[el] = out.get(el, 0) + (int(num) if num else 1)
            el, num = ("" if ch == "$" else ch), ""
        elif ch.islower():
            el += ch
        elif ch.isdigit():
            num += ch
    return out


def _reduce(counts):
    from math import gcd
    from functools import reduce as _r
    g = _r(gcd, counts.values())
    return {el: n // g for el, n in counts.items()}, g


def _reduced_formula(counts):
    red, _ = _reduce(counts)
    order = sorted(red, key=lambda el: (CHI.get(el, 9.9), el))
    return "".join(f"{el}{red[el] if red[el] > 1 else ''}" for el in order)


@lru_cache(maxsize=4096)
def classify(formula):
    counts = parse_formula(formula)
    nat = sum(counts.values())
    els = list(counts)
    reduced = _reduced_formula(counts)

    anion = max(els, key=lambda e: CHI.get(e, 0.0))
    if anion not in ANION:
        return Verdict("no_anion", None, None, reduced, None, None)

    choices = []
    for el in els:
        if el == anion:
            choices.append((ANION[el],))
        else:
            st = list(CATION.get(el, ()))
            if el in ANION:
                st.append(ANION[el])
            if not st:
                return Verdict("no_anion", None, None, reduced, None, anion)
            choices.append(tuple(st))

    best, best_states = None, None
    for combo in itertools.product(*choices):
        q = sum(counts[el] * s for el, s in zip(els, combo))
        if best is None or abs(q) < abs(best):
            best, best_states = q, dict(zip(els, combo))
            if q == 0:
                break
    if best == 0:
        return Verdict("balanced", 0.0, 0.0, reduced, best_states, anion)
    return Verdict("unbalanced", float(best), best / nat, reduced,
                   None, anion)


if __name__ == "__main__":
    import sys
    for f in sys.argv[1:]:
        print(f, classify(f))
