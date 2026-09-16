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
    "NCORE": 8, 
}

si = Structure(Lattice.cubic(5.43), ["Si"] * 2, [[0, 0, 0], [0.25, 0.25, 0.25]])

mset = MatPESStaticSet(
    si,
    user_incar_settings=incar_overrides,
    user_kpoints_settings=Kpoints.gamma_automatic((1, 1, 1)), # also make sure to use vasp_gam rather than vasp_std for speed
)
mset.write_input("example_nvt")
