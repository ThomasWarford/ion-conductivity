"""Decode every segment of an NCSD dt=*.jsonl.gz (not just the equilibrium-volume
production chain that ncsd_traj.load_production keeps) and label each one by its
role in the MPMorph pipeline, so every frame in the file can be attributed to a
specific MD run: an EOS volume-scan segment (labeled by its volume ratio to the
final, equilibrium-volume segment) or an equilibration/production segment at the
equilibrium volume (labeled by the existing <=2500-frame heuristic from
ncsd_traj.py, numbered in order of appearance).

Reuses ncsd_screening's own JSON decoding (_seg_records/_decode_segment) rather
than reimplementing it.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ncsd_screening"))
from ncsd_traj import _decode_segment, _seg_records  # noqa: E402

VOL_RTOL = 2e-3
EQUIL_MAX_FRAMES = 2500


@dataclass
class Frame:
    segment_index: int
    label: str  # e.g. "eos0.8", "equil0", "prod1"
    frame_index: int  # index within its segment
    species: np.ndarray
    frac: np.ndarray  # (natoms, 3)
    lattice: np.ndarray  # (3, 3)


def label_segments(segs: list[dict]) -> list[str]:
    """One label per segment, in the same order as `segs` (assumed time-sorted)."""
    eq_vol = float(abs(np.linalg.det(segs[-1]["lat"])))
    counts: dict[str, int] = defaultdict(int)
    labels = []
    for s in segs:
        vol = float(abs(np.linalg.det(s["lat"])))
        ratio = vol / eq_vol
        if abs(ratio - 1) > VOL_RTOL:
            base = f"eos{ratio:.1f}"
        elif s["frac"].shape[0] <= EQUIL_MAX_FRAMES:
            base = "equil"
        else:
            base = "prod"
        labels.append(f"{base}{counts[base]}")
        counts[base] += 1
    return labels


def load_frame_pool(path: Path) -> list[Frame]:
    """Every frame in the file, in chronological (segment, then in-segment) order,
    each tagged with its originating segment index/label/in-segment frame index."""
    segs = sorted((_decode_segment(r) for r in _seg_records(path)), key=lambda s: s["ts"])
    if not segs:
        raise ValueError(f"no segments in {path}")
    labels = label_segments(segs)

    pool = []
    for seg_idx, (seg, label) in enumerate(zip(segs, labels)):
        for fi in range(seg["frac"].shape[0]):
            pool.append(
                Frame(
                    segment_index=seg_idx,
                    label=label,
                    frame_index=fi,
                    species=seg["syms"],
                    frac=seg["frac"][fi],
                    lattice=seg["lat"],
                )
            )
    return pool


if __name__ == "__main__":
    # Spot-check: print each segment's index/label/nframes/volume for one file.
    import glob

    argv = sys.argv[1:]
    pattern = argv[0] if argv else "data_ncsd/trajectories/temperature=1000K/chemical_system=Li-P/*.jsonl.gz"
    path = glob.glob(pattern)[0]
    segs = sorted((_decode_segment(r) for r in _seg_records(path)), key=lambda s: s["ts"])
    labels = label_segments(segs)
    eq_vol = float(abs(np.linalg.det(segs[-1]["lat"])))
    print(f"{path}  (eq_vol={eq_vol:.2f})")
    for i, (s, label) in enumerate(zip(segs, labels)):
        vol = float(abs(np.linalg.det(s["lat"])))
        print(f"  seg{i:02d}  label={label:10s}  nframes={s['frac'].shape[0]:5d}  vol={vol:9.2f}  ratio={vol/eq_vol:.3f}")
    pool = load_frame_pool(path)
    print(f"total frames in pool: {len(pool)}")
