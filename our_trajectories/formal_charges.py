"""Nominal (formal, integer) oxidation states per element, one dict per target
compound -- not derived algorithmically, since several elements here (P, Sb,
Nb, Ge, Sn) take different formal states in different compounds depending on
local bonding, and a single global per-element mapping wouldn't hold. Each
entry was checked to charge-balance to exactly zero against that compound's
real stoichiometry (read from the NCSD trajectory data itself, not guessed).

Li-Nb-S's Nb5+ charge-balances cleanly even though this is the one compound
needing ISPIN=2 -- the real DFT ground state carries unpaired spin somewhere
this nominal value doesn't capture (see REF_ddec6_spin_moments for that).

Ge-Li-Sn (Li6Ge2Sn) is a Zintl-type polyanion: charge balance alone can't
split the anionic charge between Ge and Sn (one equation, two unknowns). Ge
and Sn are treated identically (-2 each) as isovalent group-14 congeners
sharing one anionic sublattice equally -- an approximation, not a literature
value.
"""
from __future__ import annotations

FORMAL_CHARGES: dict[str, dict[str, int]] = {
    "Li-P": {"Li": 1, "P": -3},
    "Al-Li-Na-P": {"Al": 3, "Li": 1, "Na": 1, "P": -3},
    "Li-Rb-S": {"Li": 1, "Rb": 1, "S": -2},
    "Li-Rb-Se": {"Li": 1, "Rb": 1, "Se": -2},
    "F-Li-P": {"F": -1, "Li": 1, "P": 5},
    "Li-O-P-Ti": {"Li": 1, "O": -2, "P": 5, "Ti": 4},
    "Li-S-Sb": {"Li": 1, "S": -2, "Sb": 3},
    "Bi-Li-S": {"Bi": 3, "Li": 1, "S": -2},
    "Cl-Li": {"Cl": -1, "Li": 1},
    "Cu-Li-S": {"Cu": 1, "Li": 1, "S": -2},
    "Li-O-Zr": {"Li": 1, "O": -2, "Zr": 4},
    "Li-Nb-S": {"Li": 1, "Nb": 5, "S": -2},
    "Ge-Li-Sn": {"Ge": -2, "Li": 1, "Sn": -2},
}


def formal_charges_for(compound: str) -> dict[str, int]:
    """Indirection point so a future automatic oxidation-state guesser
    replaces one function, not every call site."""
    return FORMAL_CHARGES[compound]
