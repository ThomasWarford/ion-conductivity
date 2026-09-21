"""Run chargemol/DDEC6 on a finished single-point VASP directory (CHGCAR/AECCAR0/
AECCAR2/POTCAR present), label the frame with DFT + DDEC6 training targets as
an extended-xyz file, then clean up.

Labeling (`build_atoms`): merges three sources into one ASE `Atoms` object --
`vasprun.xml` (energy/forces/stress/Fermi level), the live `ChargemolAnalysis`
object (DDEC6 charges/spin moments/bond orders/dipoles -- not the r-moments,
unused by MACE-SCF, though the raw `ddec/` output files with those are still
kept on disk), and a static per-compound formal-oxidation-state lookup
(`formal_charges.py`). All fields are prefixed `REF_` (training targets); see
`formal_charges.py` for the oxidation-state rationale. Atom order across all
three sources is asserted to match before anything is attached -- see the
comment on `build_atoms` for why that's not something pymatgen itself checks.

Cleanup: on success, delete AECCAR0 and CHG (AECCAR0 is redundant once DDEC
charges exist; CHG is VASP's non-augmented pseudo density, redundant with
CHGCAR), plus POTCAR/REPORT/WAVECAR/XDATCAR/IBZKPT/CONTCAR (run bookkeeping/
restart files with nothing training needs, now that vasprun.xml and atoms.xyz
exist), lossily compress the remaining volumetric files (CHGCAR/AECCAR1/
AECCAR2) to zfp-compressed .h5 via `voltools compress --delete-originals` (see
../../volumetric-tools; verified error is far below DFT/AIMD noise, see
compression_test/RESULTS.md and volumetric-tools/benchmarks/vasp.md), and mark
the directory OK; on failure, leave everything in place (including the raw
volumetric files, uncompressed) and mark it failed, since we can't be sure
they won't be needed to retry. Chargemol's own output (DDEC6 charges, bond
orders, etc.) is written to a persistent `ddec/` subdirectory -- pymatgen's
default is to run it in a temporary directory that's deleted once parsed into
memory, which would silently discard the whole point of this pipeline. Either
way, gzip everything else remaining, recursively (voltools' own .h5 output is
already compressed, so it's skipped).

    conda run -n dft python label_frame.py <case_dir> [<case_dir> ...]
"""
from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read as ase_read
from ase.io import write as ase_write

# ~/bin (where Chargemol_09_26_2017_linux_parallel lives) isn't always on PATH in
# non-interactive/batch shells -- make sure it is before pymatgen looks for it.
os.environ["PATH"] = f"{Path.home() / 'bin'}:{os.environ.get('PATH', '')}"

from pymatgen.command_line.chargemol_caller import ChargemolAnalysis  # noqa: E402

from formal_charges import FORMAL_CHARGES  # noqa: E402

ATOMIC_DENSITIES = Path.home() / "chargemol" / "chargemol_09_26_2017" / "atomic_densities"
VOLTOOLS = Path.home() / "volumetric-tools" / ".venv" / "bin" / "voltools"

# VASP prints this as a single line, e.g.:
#   E-fermi :   0.6764     XC(G=0):  -7.1214     alpha+bet : -6.5542
# only in OUTCAR, not vasprun.xml.
FERMI_LINE_RE = re.compile(
    r"E-fermi\s*:\s*([-\d.]+)\s+XC\(G=0\):\s*([-\d.]+)\s+alpha\+bet\s*:\s*([-\d.]+)"
)


def compress_volumetric(case_dir: Path) -> None:
    """Compress CHGCAR/AECCAR1/AECCAR2 to zfp .h5 (default accuracy) and delete
    the raw files, logging the CLI's per-file ratio/error report."""
    with open(case_dir / "voltools.log", "wb") as log:
        subprocess.run(
            [str(VOLTOOLS), "compress", str(case_dir), "--delete-originals"],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


def gzip_all(case_dir: Path) -> None:
    for f in case_dir.rglob("*"):
        if f.is_file() and f.suffix not in (".gz", ".h5"):
            with open(f, "rb") as fin, gzip.open(f"{f}.gz", "wb") as fout:
                shutil.copyfileobj(fin, fout)
            f.unlink()


def parse_fermi_level(outcar_path: Path) -> tuple[float, float]:
    """(raw E-fermi, E-fermi + alpha+bet), both eV."""
    match = FERMI_LINE_RE.search(outcar_path.read_text())
    if match is None:
        raise ValueError(f"no E-fermi/XC(G=0)/alpha+bet line in {outcar_path}")
    efermi, _xc_g0, alpha_bet = (float(g) for g in match.groups())
    return efermi, efermi + alpha_bet


def build_atoms(case_dir: Path, chargemol: ChargemolAnalysis, compound: str) -> Atoms:
    """Merge vasprun.xml (energy/forces/stress/Fermi level), the DDEC6 result
    already held by `chargemol`, and `compound`'s formal-charge lookup into one
    ASE Atoms object. All three sources are keyed by plain atom order (no ids),
    so a silent mismatch would mislabel every array. Checked directly against a
    real case during planning: CHGCAR is written species-grouped, chargemol's
    DDEC output preserves that grouping with no reordering logic of its own,
    and VASP preserves input atom order in every output -- true in practice,
    but never verified by pymatgen or VASP itself, so it's asserted here
    instead of assumed.
    """
    atoms = ase_read(case_dir / "vasprun.xml")

    chargemol_symbols = [str(s) for s in chargemol.structure.species]
    if atoms.get_chemical_symbols() != chargemol_symbols:
        raise ValueError(
            f"atom order mismatch in {case_dir}: vasprun.xml="
            f"{atoms.get_chemical_symbols()} vs chargemol={chargemol_symbols}"
        )
    n = len(atoms)

    def as_array(values: list[float] | None) -> np.ndarray:
        return np.full(n, np.nan) if values is None else np.asarray(values, dtype=float)

    ddec_charges = as_array(chargemol.ddec_charges)
    ddec_dipoles = np.asarray(chargemol.dipoles, dtype=float)

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    efermi, efermi_plus_alpha_bet = parse_fermi_level(case_dir / "OUTCAR")
    config_type = atoms.calc.parameters["system"]  # INCAR SYSTEM tag, echoed into vasprun.xml

    atoms.calc = None  # drop vasprun.xml's own energy/forces/stress -- REF_* below replace them

    atoms.arrays["REF_forces"] = forces
    atoms.arrays["REF_ddec6_charges"] = ddec_charges
    atoms.arrays["REF_ddec6_spin_moments"] = as_array(chargemol.ddec_spin_moments)
    atoms.arrays["REF_ddec6_bond_order_sums"] = as_array(chargemol.bond_order_sums)
    atoms.arrays["REF_ddec6_dipoles"] = ddec_dipoles
    atoms.arrays["REF_formal_charges"] = np.array(
        [FORMAL_CHARGES[compound][sym] for sym in chargemol_symbols], dtype=float
    )

    multipoles = np.zeros((n, 4))
    multipoles[:, 0] = ddec_charges
    multipoles[:, 1:] = ddec_dipoles[:, [1, 2, 0]]  # Cartesian (x,y,z) -> e3nn/MACE 1o (y,z,x)
    atoms.arrays["REF_multipoles"] = multipoles

    atoms.info["REF_energy"] = energy
    atoms.info["REF_stress"] = stress
    atoms.info["REF_total_charge"] = 0.0
    atoms.info["REF_vasp_fermi_level"] = efermi
    atoms.info["REF_vasp_fermi_level_plus_alpha_bet"] = efermi_plus_alpha_bet

    # not a training target (no REF_ prefix) -- a grouping/weighting tag for
    # fitting frameworks, reusing the SYSTEM tag already written by
    # build_chargemol_frames.py (compound_temp_seg-label_frame).
    atoms.info["config_type"] = config_type

    return atoms


def process(case_dir: Path) -> bool:
    os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))
    compound = case_dir.resolve().parents[2].name
    try:
        # path must be absolute: pymatgen's _execute_chargemol chdir's into
        # run_dir (path/ddec) before re-resolving chargemol_output_path from
        # the *same* string it was given -- with a relative path that second
        # resolution is relative to the new cwd, doubling the path and making
        # every output file "not found" even though chargemol ran fine.
        chargemol = ChargemolAnalysis(
            path=str(case_dir.resolve()),
            atomic_densities_path=str(ATOMIC_DENSITIES),
            run_chargemol=True,
            run_dir="ddec",
        )
        atoms = build_atoms(case_dir, chargemol, compound)
    except Exception:
        (case_dir / "CHARGEMOL_FAILED").write_text(traceback.format_exc())
        gzip_all(case_dir)
        return False

    ase_write(case_dir / "atoms.xyz", atoms)

    for name in (
        "AECCAR0",
        "CHG",
        "POTCAR",
        "REPORT",
        "WAVECAR",
        "XDATCAR",
        "IBZKPT",
        "CONTCAR",
    ):
        f = case_dir / name
        if f.exists():
            f.unlink()
    # pymatgen copies uncompressed CHGCAR/AECCAR0/AECCAR2/POTCAR into ddec/ as
    # chargemol's working inputs -- pure duplicates of data we already have (or,
    # for AECCAR0, intentionally discard above) that would otherwise dwarf the
    # actual DDEC6 output (a few hundred KB) with ~150MB of redundant copies.
    for name in ("CHGCAR", "AECCAR0", "AECCAR2", "POTCAR"):
        f = case_dir / "ddec" / name
        if f.exists():
            f.unlink()
    compress_volumetric(case_dir)
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
