"""For one NCSD .jsonl.gz trajectory file (selected by Slurm array index),
write the last quarter of its frames -- concatenating all segments in time
order -- as a sibling .traj and .xyz.gz file next to the source.

    conda run -n dft python last_quarter.py <RANK> [data_dir]
"""
from __future__ import annotations

import gzip
from pathlib import Path
from sys import argv

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write as ase_write
from ase.io.trajectory import Trajectory

try:
    import orjson

    def _loads(s):
        return orjson.loads(s)
except ImportError:  # pragma: no cover
    import json

    def _loads(s):
        return json.loads(s)


def load_segments(path: Path) -> list[dict]:
    with gzip.open(path, "rb") as fh:
        recs = [_loads(ln) for ln in fh if ln.strip()]
    recs.sort(key=lambda r: int(r["metadata"]["_id"][:8], 16))
    return recs


def build_atoms(tr: dict, fi: int, cell, symbols, fp: dict) -> Atoms:
    frac = np.asarray(tr["coords"][fi], dtype=float)
    atoms = Atoms(symbols=symbols, scaled_positions=frac, cell=cell, pbc=True)
    energy = fp.get("e_0_energy")
    forces = fp.get("forces")
    stress = fp.get("stress")
    if energy is not None or forces is not None or stress is not None:
        atoms.calc = SinglePointCalculator(
            atoms,
            energy=energy,
            forces=np.asarray(forces, dtype=float) if forces is not None else None,
            stress=np.asarray(stress, dtype=float) if stress is not None else None,
        )
    return atoms


def process(path: Path) -> tuple[int, int]:
    recs = load_segments(path)
    if not recs:
        raise ValueError(f"no segments in {path}")

    seg_lens = []
    for rec in recs:
        tr = rec["trajectory"]
        if tr.get("coords_are_displacement"):
            raise ValueError(f"unexpected displacement coords in {path}")
        seg_lens.append(len(tr["coords"]))
    total = sum(seg_lens)
    start_global = total - total // 4

    frames = []
    offset = 0
    for seg_idx, (rec, seg_len) in enumerate(zip(recs, seg_lens)):
        seg_end = offset + seg_len
        if seg_end <= start_global:
            offset = seg_end
            continue

        tr = rec["trajectory"]
        lat_arr = np.asarray(tr["lattice"], dtype=float)
        symbols = [s["element"] if isinstance(s, dict) else s for s in tr["species"]]
        frame_props = tr.get("frame_properties") or [{}] * seg_len
        local_start = max(0, start_global - offset)

        for fi in range(local_start, seg_len):
            cell = lat_arr[fi] if lat_arr.ndim == 3 else lat_arr
            fp = frame_props[fi] if fi < len(frame_props) else {}
            atoms = build_atoms(tr, fi, cell, symbols, fp)
            atoms.info["segment_index"] = seg_idx
            atoms.info["frame_index"] = fi
            if "temperature" in rec.get("metadata", {}):
                atoms.info["temperature"] = rec["metadata"]["temperature"]
            frames.append(atoms)

        offset = seg_end

    stem = path.name.split(".")[0]
    traj_dst = path.parent / f"{stem}_last_quarter.traj"
    xyz_dst = path.parent / f"{stem}_last_quarter.xyz.gz"

    with Trajectory(str(traj_dst), "w") as traj:
        for atoms in frames:
            traj.write(atoms)

    with gzip.open(xyz_dst, "wt") as fh:
        ase_write(fh, frames, format="extxyz")

    return total, len(frames)


def main() -> None:
    args = argv[1:]
    rank = int(args[0])
    data_dir = Path(args[1]) if len(args) > 1 else Path("data_ncsd/trajectories")

    paths = sorted(data_dir.glob("**/*.jsonl.gz"))
    path = paths[rank]

    total, n = process(path)
    stem = path.name.split(".")[0]
    print(f"[{rank}] {path}: {n}/{total} frames -> {stem}_last_quarter.{{traj,xyz.gz}}")


if __name__ == "__main__":
    main()
