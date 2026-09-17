"""
Perf: use vasp_gam (not vasp_std) and module vasp/cpu_6.4.2-intel (needs
`sg vasp -c "sbatch ..."` for group access). Run on exactly 1 node -- for
this gamma-only, ~100-atom case, going to 2+ nodes is 10-25x slower, not
faster. Use 64 ranks/node with --cpus-per-task=2 (1 rank/physical core;
these nodes have 64 physical cores but expose 128 via SMT-2 -- packing
128 ranks/node is ~1.6x slower).
"""

from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp.sets import MatPESStaticSet
from pymatgen.io.vasp.inputs import Kpoints

incar_overrides = {
    # deviations from default MatPESStaticSet values
    "EFERMI": "Midgap",
    "ENAUG": None,
    "GGA_COMPAT": False,
    "ENCUT": 520, # old MPStaticSet value, so decent
    "ISPIN": 1, # doesn't matter for most systems, and if it does this affects accuracy rather than precision
    # NVT ensemble (fixed cell, velocity-scaling thermostat)
    "IBRION": 0,
    "ISIF": 0,
    "ISYM": 0,
    "SMASS": 0,
    "NSW": 2000,
    "POTIM": 2.0,
    "TEBEG": 1000,
    "TEEND": 1000,
    # performance
    "NCORE": 32,
}

if __name__ == "__main__":
    si = Structure(Lattice.cubic(5.43), ["Si"] * 2, [[0, 0, 0], [0.25, 0.25, 0.25]])

    mset = MatPESStaticSet(
        si,
        user_incar_settings=incar_overrides,
        user_kpoints_settings=Kpoints.gamma_automatic((1, 1, 1)),
    )
    mset.write_input("example_nvt")
