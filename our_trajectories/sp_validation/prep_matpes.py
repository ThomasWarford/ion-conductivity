"""Build short NVT-MD timing runs using the MatPESStaticSet-based recipe from
../nvt_md_input.py, for the same 6 frames used in prep.py, so we can compare
per-step wall time against the MPMorphMDSetGenerator-based exact/tight/gold
tiers and extrapolate a feasible trajectory length.

    conda run -n dft python prep_matpes.py
"""
from __future__ import annotations

from pathlib import Path

from pymatgen.io.vasp.sets import MatPESStaticSet
from pymatgen.io.vasp.inputs import Kpoints

from prep import FRAMES, find_file, pick_frame

HERE = Path(__file__).resolve().parent
NSW = 40  # enough MD steps for a stable per-step average, not a real trajectory


def write_case(tag: str, struct, temp: int) -> Path:
    incar_overrides = {
        "EFERMI": "Midgap",
        "ENAUG": None,
        "GGA_COMPAT": False,
        "ENCUT": 520,
        "IBRION": 0,
        "ISIF": 0,
        "ISYM": 0,
        "SMASS": 0,
        "NCORE": 8,
        "NSW": NSW,
        "POTIM": 2.0,
        "TEBEG": temp,
        "TEEND": temp,
    }
    mset = MatPESStaticSet(
        struct,
        user_incar_settings=incar_overrides,
        user_kpoints_settings=Kpoints.gamma_automatic((1, 1, 1)),
    )
    out = HERE / "runs" / tag / "matpes"
    out.mkdir(parents=True, exist_ok=True)
    mset.incar.write_file(out / "INCAR")
    mset.poscar.write_file(out / "POSCAR")
    mset.kpoints.write_file(out / "KPOINTS")
    mset.potcar.write_file(out / "POTCAR")
    return out


def main() -> None:
    for system, temp in FRAMES:
        tag = f"{system}_{temp}K"
        struct, ref = pick_frame(find_file(system, temp))
        out = write_case(tag, struct, temp)
        print(f"{tag:18s}  natoms={ref['natoms']:3d}  TEBEG=TEEND={temp}  -> {out}")


if __name__ == "__main__":
    main()
