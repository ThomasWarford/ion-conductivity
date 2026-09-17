"""Build VASP single-point directories that replicate the MP amorphous_diffusivity
(NCSD) AIMD settings, for a handful of frames pulled straight out of
../../data_ncsd.

Settings come from `MPMorphMDSetGenerator` -- the same generator
`our_trajectories/aimd_submit.py` uses -- turned into a single point:
IBRION=-1, NSW=0, nothing else touched.  For each frame we also write a
`tight/` copy whose only change is EDIFF 0.005 -> 1e-6, i.e. the fully
SCF-converged energy/forces under the same basis / FFT grid / projectors.

The reference e_0_energy / e_fr_energy / e_wo_entrp, forces and stress for the
frame are saved next to the inputs as ref.json.

    conda run -n dft python prep.py
"""
from __future__ import annotations

import gzip
import json
import os
import warnings
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure

warnings.simplefilter("ignore")
from atomate2.vasp.sets.mpmorph import MPMorphMDSetGenerator  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "data_ncsd" / "trajectories"
POTCAR_FUNCTIONAL = "PBE_54"

# (chemical_system, temperature K).  One production frame from each.
FRAMES = [
    ("Cl-Li", 1000),
    ("Cl-Li", 2500),
    ("Li-P", 1500),
    ("Ge-Li-Sn", 2000),
    ("Li-O-Zr", 1500),
    ("Al-K-Li-P", 1000),
]


def load_segments(path: Path) -> list[dict]:
    with gzip.open(path, "rb") as fh:
        recs = [json.loads(ln) for ln in fh if ln.strip()]
    recs.sort(key=lambda r: int(r["metadata"]["_id"][:8], 16))
    return recs


def pick_frame(path: Path) -> tuple[Structure, dict]:
    """Middle frame of the last (equilibrium-volume production) segment."""
    recs = load_segments(path)
    seg = recs[-1]["trajectory"]
    lat = np.asarray(seg["lattice"], dtype=float)
    lat = lat[0] if lat.ndim == 3 else lat
    syms = [s["element"] if isinstance(s, dict) else s for s in seg["species"]]
    fi = len(seg["coords"]) // 2
    frac = np.asarray(seg["coords"], dtype=float)[fi]
    fp = seg["frame_properties"][fi]
    struct = Structure(Lattice(lat), syms, frac, coords_are_cartesian=False)
    ref = dict(
        chemical_system=path.parent.name.split("chemical_system=")[-1],
        source_file=path.name,
        segment_index=len(recs) - 1,
        frame_index=fi,
        natoms=len(struct),
        formula=struct.composition.formula,
        e_0_energy=fp["e_0_energy"],
        e_fr_energy=fp["e_fr_energy"],
        e_wo_entrp=fp["e_wo_entrp"],
        kinetic=fp.get("kinetic"),
        forces=np.asarray(fp["forces"], dtype=float).tolist(),
        stress=np.asarray(fp["stress"], dtype=float).tolist(),
        lattice=lat.tolist(),
    )
    return struct, ref


def find_file(system: str, temp: int) -> Path:
    d = DATA / f"temperature={temp}K" / f"chemical_system={system}"
    files = sorted(d.glob("*.jsonl.gz"))
    if not files:
        raise FileNotFoundError(d)
    return files[0]


# Three tiers per frame:
#   exact  - the MPMorphMDSetGenerator INCAR verbatim, only IBRION=-1/NSW=0.
#            Loose EDIFF (~5e-3/cell), LREAL=.TRUE., PREC=Normal, ENCUT=ENMAX.
#            This is exactly what aimd_submit.py would run.
#   tight  - identical basis/grid/projectors, EDIFF=1e-7.  Isolates the SCF
#            convergence noise in the reference frames.
#   gold   - EDIFF=1e-7 AND LREAL=.FALSE., PREC=Accurate, ENCUT=700.  The
#            fully converged E/F for this geometry; the gap gold-tight is the
#            systematic cost of MP's own LREAL/PREC/ENCUT choices, which
#            "replicating MP" necessarily inherits.
TIERS = ("exact", "tight", "gold")


def write_case(struct: Structure, tag: str, tier: str) -> Path:
    gen = MPMorphMDSetGenerator(user_potcar_functional=POTCAR_FUNCTIONAL)
    vis = gen.get_input_set(struct)
    incar = vis.incar
    # MD -> single point.
    incar["IBRION"] = -1
    incar["NSW"] = 0
    incar["LWAVE"] = False
    incar["LCHARG"] = False
    incar["NCORE"] = 8
    if tier in ("tight", "gold"):
        incar["EDIFF"] = 1e-7
    if tier == "gold":
        incar["LREAL"] = False
        incar["PREC"] = "Accurate"
        incar["ENCUT"] = 700
        incar["ADDGRID"] = False  # not needed with the Accurate grid + exact projectors
    out = HERE / "runs" / tag / tier
    out.mkdir(parents=True, exist_ok=True)
    incar.write_file(out / "INCAR")
    vis.poscar.write_file(out / "POSCAR")
    vis.kpoints.write_file(out / "KPOINTS")
    vis.potcar.write_file(out / "POTCAR")
    return out


def main() -> None:
    (HERE / "runs").mkdir(exist_ok=True)
    index = []
    for system, temp in FRAMES:
        tag = f"{system}_{temp}K"
        path = find_file(system, temp)
        struct, ref = pick_frame(path)
        for tier in TIERS:
            out = write_case(struct, tag, tier)
            (out / "ref.json").write_text(json.dumps(ref, indent=1))
        index.append(tag)
        print(f"{tag:18s}  {ref['formula']:22s}  natoms={ref['natoms']:3d}  "
              f"seg{ref['segment_index']} frame{ref['frame_index']}  "
              f"E0={ref['e_0_energy']:.4f} eV  "
              f"max|F|={np.abs(np.array(ref['forces'])).max():.3f}")
    (HERE / "runs" / "cases.txt").write_text("\n".join(index) + "\n")
    print(f"\n{len(index)} cases -> {HERE/'runs'}")


if __name__ == "__main__":
    main()
