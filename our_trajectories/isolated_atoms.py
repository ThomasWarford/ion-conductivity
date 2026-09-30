"""Isolated-atom reference energies (E0s) for MLIP training.

E0s must come from the same Hamiltonian as the training labels, so every INCAR
here is the single-point recipe from build_chargemol_frames.py (MatPESStaticSet +
nvt_md_input.incar_overrides, Gamma-only, same POTCARs), with only
energy-neutral changes on top: ISYM=0 (let the atom break symmetry), a larger
NELM (atoms converge slowly), and no CHGCAR/AECCAR/WAVECAR output.

One atom sits in an 18 x 18.1 x 18.2 A orthorhombic box (unequal edges so open
p/d shells aren't pinned to cubic degeneracy). Box-size convergence is checked
on the diffuse/representative atoms (Li, Na, Rb, Cl, O) at 15 and 21 A: 15 A was
the first choice, but Rb_sv is 2.3 meV off the 21 A value there (18 A: 0.2 meV;
all other atoms < 0.3 meV at 15 A).

Both ISPIN=1 (same Hamiltonian as 12/13 systems' labels -- the default E0 set)
and ISPIN=2 (physical spin-polarized atom, Hund's-rule starting moment) are run.

Energies are read exactly like REF_energy in the training frames: ASE
get_potential_energy() on vasprun.xml, i.e. the sigma->0 energy (e_0_energy),
not the smeared free energy.

    conda run -n dft python isolated_atoms.py build
    sg vasp -c "sbatch our_trajectories/run_isolated_atoms.sbatch"   # from repo root
    conda run -n dft python isolated_atoms.py collect
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ase.io import read as ase_read
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MatPESStaticSet

from build_chargemol_frames import _MD_ONLY_KEYS
from build_runs import BY_SYSTEM
from nvt_md_input import incar_overrides as md_incar_overrides

HERE = Path(__file__).resolve().parent
OUT = HERE / "isolated_atoms"

ELEMENTS = sorted({el for system in BY_SYSTEM for el in system.split("-")})

BOX_OFFSETS = (0.0, 0.1, 0.2)
L_DEFAULT = 18
L_CONVERGENCE = (15, 21)
CONVERGENCE_ELEMENTS = ("Li", "Na", "Rb", "Cl", "O")

# Hund's-rule ground-state moments (mu_B), used as ISPIN=2 starting MAGMOM.
HUND_MOMENT = {
    "Al": 1, "Bi": 3, "Cl": 1, "Cu": 1, "F": 1, "Ge": 2, "Li": 1, "Na": 1, "Nb": 5,
    "O": 2, "P": 3, "Rb": 1, "S": 2, "Sb": 3, "Se": 2, "Sn": 2, "Ti": 2, "Zr": 2,
}

NELM = 400


def atom_incar(element: str, ispin: int) -> dict:
    incar = {k: v for k, v in md_incar_overrides.items() if k not in _MD_ONLY_KEYS}
    incar.update(
        IBRION=-1,
        NSW=0,
        ISPIN=ispin,
        ISYM=0,
        NELM=NELM,
        LAECHG=False,
        LCHARG=False,
        LWAVE=False,
        LORBIT=None,
        SYSTEM=f"{element}_atom_ispin{ispin}",
    )
    if ispin == 2:
        incar["MAGMOM"] = {element: HUND_MOMENT[element]}
    return incar


def case_dir(element: str, ispin: int, box: int) -> Path:
    return OUT / f"ispin{ispin}" / f"L{box}" / element


def write_case(element: str, ispin: int, box: int) -> Path:
    a, b, c = (box + d for d in BOX_OFFSETS)
    struct = Structure(Lattice.orthorhombic(a, b, c), [element], [[0, 0, 0]])
    mset = MatPESStaticSet(
        struct,
        user_incar_settings=atom_incar(element, ispin),
        user_kpoints_settings=Kpoints.gamma_automatic((1, 1, 1)),
    )
    out_dir = case_dir(element, ispin, box)
    mset.write_input(str(out_dir))
    return out_dir


def all_cases() -> list[tuple[str, int, int]]:
    cases = [(el, ispin, L_DEFAULT) for ispin in (1, 2) for el in ELEMENTS]
    cases += [(el, 1, box) for box in L_CONVERGENCE for el in CONVERGENCE_ELEMENTS]
    return cases


def build() -> None:
    OUT.mkdir(exist_ok=True)
    dirs = [write_case(*case) for case in all_cases()]
    (OUT / "cases.txt").write_text("\n".join(str(d.relative_to(HERE)) for d in dirs) + "\n")
    print(f"{len(dirs)} cases ({len(ELEMENTS)} elements: {' '.join(ELEMENTS)}) -> {OUT / 'cases.txt'}")


def read_case(d: Path) -> dict | None:
    vr_path = d / "vasprun.xml"
    if not vr_path.exists():
        return None
    try:
        vr = Vasprun(str(vr_path), parse_dos=False, parse_eigen=False, parse_potcar_file=False)
    except Exception:  # truncated vasprun.xml from a crashed run
        return None
    n_elec = len(vr.ionic_steps[-1]["electronic_steps"])
    result = {
        # identical to how REF_energy was produced for the training frames
        "energy": ase_read(vr_path).get_potential_energy(),
        "converged": vr.converged_electronic and n_elec < NELM,
        "n_elec": n_elec,
    }
    if vr.parameters.get("ISPIN") == 2:
        outcar_mag = [line for line in (d / "OSZICAR").read_text().splitlines() if "mag=" in line]
        result["magmom"] = float(outcar_mag[-1].split("mag=")[1]) if outcar_mag else None
    return result


def collect() -> None:
    results = {case: read_case(case_dir(*case)) for case in all_cases()}
    missing = [case for case, r in results.items() if r is None]
    unconverged = [case for case, r in results.items() if r and not r["converged"]]
    for case in missing:
        print(f"MISSING     {case}")
    for case in unconverged:
        print(f"UNCONVERGED {case} ({results[case]['n_elec']} electronic steps)")

    print("\nbox convergence (ISPIN=1), meV relative to largest box:")
    for el in CONVERGENCE_ELEMENTS:
        boxes = sorted((L_DEFAULT, *L_CONVERGENCE))
        energies = [results[(el, 1, box)] for box in boxes]
        if any(e is None for e in energies):
            continue
        ref = energies[-1]["energy"]
        row = "  ".join(f"L{box}: {1000 * (e['energy'] - ref):+7.2f}" for box, e in zip(boxes, energies))
        print(f"  {el:2s}  {row}")

    print(f"\nper-element E0 (eV), L{L_DEFAULT}:")
    print(f"  {'el':2s}  {'ISPIN=1':>12s}  {'ISPIN=2':>12s}  {'2-1':>8s}  {'mag':>5s}")
    for el in ELEMENTS:
        r1, r2 = results[(el, 1, L_DEFAULT)], results[(el, 2, L_DEFAULT)]
        e1 = f"{r1['energy']:12.6f}" if r1 else f"{'-':>12s}"
        e2 = f"{r2['energy']:12.6f}" if r2 else f"{'-':>12s}"
        de = f"{r2['energy'] - r1['energy']:+8.4f}" if r1 and r2 else f"{'-':>8s}"
        mag = f"{r2['magmom']:5.2f}" if r2 and r2.get("magmom") is not None else f"{'-':>5s}"
        print(f"  {el:2s}  {e1}  {e2}  {de}  {mag}")

    for ispin in (1, 2):
        rs = {el: results[(el, ispin, L_DEFAULT)] for el in ELEMENTS}
        if all(r and r["converged"] for r in rs.values()):
            path = OUT / f"e0s_ispin{ispin}.json"
            path.write_text(json.dumps({el: r["energy"] for el, r in rs.items()}, indent=2) + "\n")
            print(f"wrote {path}")
        else:
            print(f"not writing e0s_ispin{ispin}.json: missing/unconverged cases")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=("build", "collect"))
    args = ap.parse_args()
    {"build": build, "collect": collect}[args.command]()


if __name__ == "__main__":
    main()
