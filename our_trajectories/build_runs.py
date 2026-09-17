"""Build production NVT-AIMD input directories for the NCSD-final-structure
campaign, at scale, using the exact recipe in `nvt_md_input.py`.

For each (system, temperature, replica) case: load the last frame of that
system's NCSD `*_last_quarter.traj` at that temperature as the starting
structure, override TEBEG/TEEND/ISPIN/NSW on top of
`nvt_md_input.incar_overrides`, and write VASP inputs into
`our_trajectories/runs/<system>_<temp>K_r<replica>/`.

Replicas are NOT differentiated via INCAR's RANDOM_SEED tag: setting it (even to
values inside VASP's own documented valid range) was found to hang vasp_gam
indefinitely on this cluster's vasp/cpu_6.4.2-intel build (confirmed by isolating
it as the only changed setting from an otherwise-working INCAR -- see git history
for the debugging session). Left unset, VASP seeds its own RNG per run, which is
enough to make independent replicas in practice.

Two manifests are written: `cases_base.txt` (13 systems x {1000,2500}K, replica 0
-- submit now) and `cases_extra.txt` (Cu-Li-S and Li-P at all 4 temperatures x 3
replicas, minus what's already in the base set -- build now, submit later).

    conda run -n dft python build_runs.py
"""
from __future__ import annotations

import copy
from pathlib import Path

from ase.io.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MatPESStaticSet

from nvt_md_input import incar_overrides

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data_ncsd" / "trajectories"
RUNS = HERE / "runs"

# (system, ISPIN, measured s/step at 2500K -- see sp_validation/logs/sp_perf2500_*.out)
TARGETS = [
    ("F-Li-P", 1, 9.38),
    ("Li-P", 1, 10.43),
    ("Li-O-Zr", 1, 10.47),
    ("Li-S-Sb", 1, 14.26),
    ("Cl-Li", 1, 15.13),
    ("Li-O-P-Ti", 1, 15.70),
    ("Al-Li-Na-P", 1, 18.65),
    ("Cu-Li-S", 1, 19.22),
    ("Bi-Li-S", 1, 19.36),
    ("Li-Nb-S", 2, 29.15),
    ("Li-Rb-S", 1, 29.72),
    ("Li-Rb-Se", 1, 31.96),
    ("Ge-Li-Sn", 1, 32.11),
]
BY_SYSTEM = {system: (ispin, s_per_step) for system, ispin, s_per_step in TARGETS}

BASE_TEMPS = (1000, 2500)
EXTRA_SYSTEMS = ["Cu-Li-S", "Li-P"]
EXTRA_TEMPS = (1000, 1500, 2000, 2500)
N_REPLICAS_EXTRA = 3

TARGET_WALL_SECONDS = 7 * 24 * 3600  # one sbatch link's --time
SAFETY = 0.85  # 15% margin so a link finishes comfortably inside its allocation


def nsw_for(s_per_step: float) -> int:
    return int(SAFETY * TARGET_WALL_SECONDS / s_per_step)


def build_manifest() -> tuple[set[tuple[str, int, int]], set[tuple[str, int, int]]]:
    base_cases = {(system, temp, 0) for system, _, _ in TARGETS for temp in BASE_TEMPS}
    all_extra = {
        (system, temp, rep)
        for system in EXTRA_SYSTEMS
        for temp in EXTRA_TEMPS
        for rep in range(N_REPLICAS_EXTRA)
    }
    extra_cases = all_extra - base_cases
    return base_cases, extra_cases


def case_name(system: str, temp: int, replica: int) -> str:
    return f"{system}_{temp}K_r{replica}"


def load_final_structure(system: str, temp: int):
    matches = sorted((DATA / f"temperature={temp}K" / f"chemical_system={system}").glob("*_last_quarter.traj"))
    if not matches:
        raise FileNotFoundError(f"no *_last_quarter.traj for {system} @ {temp}K")
    traj = Trajectory(str(matches[0]))
    return AseAtomsAdaptor.get_structure(traj[-1])


def incar_for(temp: int, ispin: int, nsw: int) -> dict:
    incar = copy.deepcopy(incar_overrides)
    incar.update(
        TEBEG=temp,
        TEEND=temp,
        ISPIN=ispin,
        NSW=nsw,
        # Disk hygiene only -- MatPESStaticSet defaults to writing Bader/DOS/PROCAR
        # data we don't need for a pure MD trajectory (see plan's Quotas section).
        LAECHG=False,
        LCHARG=False,
        LORBIT=None,
    )
    return incar


def write_case(system: str, temp: int, replica: int) -> tuple[Path, int]:
    ispin, s_per_step = BY_SYSTEM[system]
    nsw = nsw_for(s_per_step)
    struct = load_final_structure(system, temp)
    mset = MatPESStaticSet(
        struct,
        user_incar_settings=incar_for(temp, ispin, nsw),
        user_kpoints_settings=Kpoints.gamma_automatic((1, 1, 1)),
    )
    out_dir = RUNS / case_name(system, temp, replica)
    mset.write_input(str(out_dir))
    return out_dir, nsw


def main() -> None:
    RUNS.mkdir(exist_ok=True)
    base_cases, extra_cases = build_manifest()

    potim = incar_overrides["POTIM"]
    names = {}
    for system, temp, replica in sorted(base_cases | extra_cases):
        out_dir, nsw = write_case(system, temp, replica)
        ps = nsw * potim / 1000
        manifest = "base" if (system, temp, replica) in base_cases else "extra"
        names[(system, temp, replica)] = out_dir.name
        print(f"{out_dir.name:20s} ISPIN={BY_SYSTEM[system][0]} NSW={nsw:6d} ({ps:6.1f} ps) [{manifest}]")

    (RUNS / "cases_base.txt").write_text(
        "\n".join(names[c] for c in sorted(base_cases)) + "\n"
    )
    (RUNS / "cases_extra.txt").write_text(
        "\n".join(names[c] for c in sorted(extra_cases)) + "\n"
    )
    print(f"\n{len(base_cases)} base cases -> {RUNS/'cases_base.txt'}")
    print(f"{len(extra_cases)} extra cases -> {RUNS/'cases_extra.txt'}")


if __name__ == "__main__":
    main()
