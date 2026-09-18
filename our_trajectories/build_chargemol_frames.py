"""Build charge-partitioning single-point VASP input directories from stratified
samples of the NCSD trajectories, for chargemol/DDEC6 training-data generation.

For each (compound, temperature): decode the *entire* raw trajectory (every
segment -- EOS volume-scan runs, iterative equilibration, and production; see
segment_labels.py), and take N evenly-spaced samples across that full frame pool.
Each sample becomes a single-point (NSW=0) VASP input directory at:

    our_trajectories/our_singlepoints/<compound>/<temp>K/seg<NN>-<label>/f<frame_idx>/

using the exact recipe in nvt_md_input.py (MatPESStaticSet + incar_overrides),
minus the MD-only tags, plus LAECHG/LCHARG for chargemol. The directory path plus
the INCAR SYSTEM tag fully encode provenance (compound, temperature, originating
MD segment/label, frame index) -- no separate manifest needed.

    conda run -n dft python build_chargemol_frames.py --compounds Li-P --n-per-temp 8
    conda run -n dft python build_chargemol_frames.py  # full 13-compound, 750/temp batch
"""
from __future__ import annotations

import argparse
import copy
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MatPESStaticSet

from build_runs import BY_SYSTEM
from nvt_md_input import incar_overrides as md_incar_overrides
from segment_labels import load_frame_pool

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data_ncsd" / "trajectories"
OUT = HERE / "our_singlepoints"

TEMPS = (1000, 1500, 2000, 2500)
N_PER_TEMP_DEFAULT = 750  # 3000/compound across the 4 temperatures
BATCH_SIZE_DEFAULT = 50  # cases per cases_batch_<i>.txt (tune from pilot timing)

# MD-only keys that don't apply to a single point (no ionic dynamics/thermostat).
_MD_ONLY_KEYS = ("IBRION", "ISIF", "ISYM", "SMASS", "NSW", "POTIM", "TEBEG", "TEEND")


def sp_incar_overrides(compound: str) -> dict:
    incar = {k: v for k, v in md_incar_overrides.items() if k not in _MD_ONLY_KEYS}
    incar.update(
        IBRION=-1,
        NSW=0,
        ISPIN=BY_SYSTEM[compound][0],
        LAECHG=True,
        LCHARG=True,
        LWAVE=False,
        LORBIT=None,  # we get charges from chargemol/DDEC6, not VASP's own PROCAR/DOS
    )
    return incar


def find_traj_file(compound: str, temp: int) -> Path:
    matches = sorted((DATA / f"temperature={temp}K" / f"chemical_system={compound}").glob("*.jsonl.gz"))
    if not matches:
        raise FileNotFoundError(f"no dt=*.jsonl.gz for {compound} @ {temp}K")
    return matches[0]


def stratified_sample_indices(n_pool: int, n_samples: int) -> np.ndarray:
    idx = np.unique(np.round(np.linspace(0, n_pool - 1, n_samples)).astype(int))
    if len(idx) != n_samples:
        raise ValueError(f"linspace collision: wanted {n_samples} distinct indices from pool of {n_pool}, got {len(idx)}")
    return idx


def write_case(compound: str, temp: int, frame, incar: dict) -> Path:
    struct = Structure(
        Lattice(frame.lattice),
        frame.species.tolist(),
        frame.frac,
        coords_are_cartesian=False,
    )
    case_incar = copy.deepcopy(incar)
    case_incar["SYSTEM"] = f"{compound}_{temp}K_seg{frame.segment_index:02d}-{frame.label}_f{frame.frame_index}"
    mset = MatPESStaticSet(
        struct,
        user_incar_settings=case_incar,
        user_kpoints_settings=Kpoints.gamma_automatic((1, 1, 1)),
    )
    out_dir = OUT / compound / f"{temp}K" / f"seg{frame.segment_index:02d}-{frame.label}" / f"f{frame.frame_index:05d}"
    mset.write_input(str(out_dir))
    return out_dir


def build_compound(compound: str, n_per_temp: int) -> list[Path]:
    incar = sp_incar_overrides(compound)
    case_dirs = []
    for temp in TEMPS:
        path = find_traj_file(compound, temp)
        pool = load_frame_pool(path)
        idx = stratified_sample_indices(len(pool), n_per_temp)
        for i in idx:
            case_dirs.append(write_case(compound, temp, pool[i], incar))
        print(f"{compound} @ {temp}K: {len(pool)} frames in pool -> {len(idx)} sampled")
    return case_dirs


def write_batches(case_dirs: list[Path], batch_size: int, tag: str) -> None:
    rel = [str(d.relative_to(HERE)) for d in case_dirs]
    n_batches = (len(rel) + batch_size - 1) // batch_size
    for i in range(n_batches):
        chunk = rel[i * batch_size : (i + 1) * batch_size]
        (OUT / f"cases_{tag}_batch_{i}.txt").write_text("\n".join(chunk) + "\n")
    print(f"{len(rel)} cases -> {n_batches} batch file(s) (cases_{tag}_batch_*.txt, batch_size={batch_size})")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--compounds", default=",".join(BY_SYSTEM), help="comma-separated compound list")
    ap.add_argument("--n-per-temp", type=int, default=N_PER_TEMP_DEFAULT)
    ap.add_argument("--batch-size", type=int, default=BATCH_SIZE_DEFAULT)
    ap.add_argument("--tag", default="full", help="prefix for the cases_<tag>_batch_*.txt files, e.g. 'pilot'")
    args = ap.parse_args()

    compounds = args.compounds.split(",")
    unknown = set(compounds) - set(BY_SYSTEM)
    if unknown:
        raise SystemExit(f"unknown compound(s): {sorted(unknown)}")

    OUT.mkdir(exist_ok=True)
    all_cases = []
    for compound in compounds:
        all_cases.extend(build_compound(compound, args.n_per_temp))
    write_batches(all_cases, args.batch_size, args.tag)


if __name__ == "__main__":
    main()
