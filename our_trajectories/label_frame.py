"""Run chargemol/DDEC6 on a finished single-point VASP directory and label it
with DFT + DDEC6 training targets, then clean up.

The DFT-mechanics part (running chargemol, merging vasprun.xml + DDEC6 into
one Atoms object, atom-order checks, volumetric compression) lives in
`volumetric_tools.vasp.label_vasp_dir` (see ../../volumetric-tools). This
module supplies the two things that are project-specific: a per-compound
formal-oxidation-state lookup (`formal_charges.py`) and a `config_type` tag,
threaded through `label_vasp_dir`'s `extra_arrays`/`extra_info`. All
training-target fields are prefixed `REF_`.

On top of `label_vasp_dir`'s own cleanup, this module also deletes
POTCAR/REPORT/WAVECAR/XDATCAR/IBZKPT/CONTCAR on success (nothing downstream
needs them once vasprun.xml and atoms.xyz exist), then marks the directory
OK and gzips everything remaining. On failure, everything is left in place
and marked failed.

    conda run -n dft python label_frame.py <case_dir> [<case_dir> ...]
"""
from __future__ import annotations

import gzip
import os
import shutil
import sys
import traceback
from pathlib import Path

import numpy as np
from ase.io import read as ase_read

# ~/bin (where Chargemol_09_26_2017_linux_parallel lives) isn't always on PATH in
# non-interactive/batch shells -- make sure it is before chargemol runs.
os.environ["PATH"] = f"{Path.home() / 'bin'}:{os.environ.get('PATH', '')}"

from pymatgen.io.vasp.inputs import Incar  # noqa: E402
from volumetric_tools.vasp import label_vasp_dir  # noqa: E402

from formal_charges import formal_charges_for  # noqa: E402

ATOMIC_DENSITIES = Path.home() / "chargemol" / "chargemol_09_26_2017" / "atomic_densities"

# Not needed downstream once vasprun.xml/atoms.xyz exist; label_vasp_dir
# handles AECCAR0/CHG/ddec-duplicates itself.
_EXTRA_CLEANUP_FILENAMES = ("POTCAR", "REPORT", "WAVECAR", "XDATCAR", "IBZKPT", "CONTCAR")


def gzip_all(case_dir: Path) -> None:
    for f in case_dir.rglob("*"):
        if f.is_file() and f.suffix not in (".gz", ".h5"):
            with open(f, "rb") as fin, gzip.open(f"{f}.gz", "wb") as fout:
                shutil.copyfileobj(fin, fout)
            f.unlink()


def process(case_dir: Path) -> bool:
    os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))
    compound = case_dir.resolve().parents[2].name
    try:
        # label_vasp_dir does its own vasprun.xml read internally; this extra
        # read (cheap, unlike the chargemol step) keeps it fully generic with
        # no project-specific concept baked in.
        preview = ase_read(case_dir / "vasprun.xml")
        symbols = preview.get_chemical_symbols()
        extra_arrays = {
            "REF_formal_charges": np.array(
                [formal_charges_for(compound)[sym] for sym in symbols], dtype=float
            )
        }
        extra_info = {"config_type": Incar.from_file(case_dir / "INCAR")["SYSTEM"]}

        label_vasp_dir(
            case_dir,
            atomic_densities_path=ATOMIC_DENSITIES,
            extra_info=extra_info,
            extra_arrays=extra_arrays,
            compress_log_path=case_dir / "voltools.log",
            delete_originals=True,
        )
    except Exception:
        (case_dir / "CHARGEMOL_FAILED").write_text(traceback.format_exc())
        gzip_all(case_dir)
        return False

    for name in _EXTRA_CLEANUP_FILENAMES:
        f = case_dir / name
        if f.exists():
            f.unlink()
    (case_dir / "CHARGEMOL_OK").touch()
    gzip_all(case_dir)
    return True


def main() -> None:
    ok = 0
    for arg in sys.argv[1:]:
        case_dir = Path(arg)
        success = process(case_dir)
        ok += success
        print(f"{'OK  ' if success else 'FAIL'} {case_dir}")
    print(f"{ok}/{len(sys.argv) - 1} succeeded")


if __name__ == "__main__":
    main()
