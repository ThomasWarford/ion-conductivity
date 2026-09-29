"""Concatenate the labeled single points in our_singlepoints/ into train/valid/test
extxyz files for MLIP training.

Frames live at <system>/<T>K/<segXX-label>/fNNNNN/atoms.xyz.gz, where each segment is
its own VASP MD run and NNNNN is the MD step within that run. The splits are:

- test:  the last `--test-frac` of sampled frames (by step) of every segment, i.e. a
         forward-in-time hold-out that covers every system, temperature and segment;
- valid: a random `--valid-frac` of each (system, T) group, drawn from what is left;
- train: everything else.

Frame text is copied verbatim (no ASE round-trip), so every REF_* field is preserved
exactly. Writes <out>/{train,valid,test}.xyz.gz, manifest.csv and split_info.json.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SPLITS = ("train", "valid", "test")


@dataclass(frozen=True)
class Frame:
    system: str
    temperature: int
    segment: str
    step: int
    path: Path  # relative to root


def discover(root: Path) -> list[Frame]:
    frames = []
    for p in sorted(root.glob("*/*K/*/f*/atoms.xyz.gz")):
        rel = p.relative_to(root)
        system, temp, segment, fdir = rel.parts[:4]
        if system.startswith("_"):
            continue
        frames.append(Frame(system, int(temp[:-1]), segment, int(fdir[1:]), rel))
    return frames


def assign(frames: list[Frame], test_frac: float, valid_frac: float,
           seed: int) -> dict[Frame, str]:
    split = {}
    by_seg = defaultdict(list)
    for f in frames:
        by_seg[(f.system, f.temperature, f.segment)].append(f)
    for seg_frames in by_seg.values():
        seg_frames.sort(key=lambda f: f.step)
        n_test = max(1, round(test_frac * len(seg_frames)))
        for f in seg_frames[-n_test:]:
            split[f] = "test"

    rng = np.random.default_rng(seed)
    by_group = defaultdict(list)
    for f in frames:
        by_group[(f.system, f.temperature)].append(f)
    for key in sorted(by_group):
        group = sorted(by_group[key], key=lambda f: (f.segment, f.step))
        rest = [f for f in group if f not in split]
        n_valid = round(valid_frac * len(group))
        order = rng.permutation(len(rest))
        for i, idx in enumerate(order):
            split[rest[idx]] = "valid" if i < n_valid else "train"
    return split


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=HERE / "our_singlepoints")
    ap.add_argument("--out", type=Path, default=HERE / "splits")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--test-frac", type=float, default=0.05,
                    help="fraction of each segment's last frames used as test")
    ap.add_argument("--valid-frac", type=float, default=0.05,
                    help="fraction of each (system, T) group drawn at random as valid")
    ap.add_argument("--expect", type=int, default=39000,
                    help="expected frame count; 0 disables the check")
    ap.add_argument("--force", action="store_true", help="overwrite an existing --out")
    args = ap.parse_args()

    frames = discover(args.root)
    if args.expect and len(frames) != args.expect:
        raise SystemExit(f"found {len(frames)} frames, expected {args.expect}")
    group_sizes = Counter((f.system, f.temperature) for f in frames)
    small = {k: n for k, n in group_sizes.items() if n < 20}
    if small:
        raise SystemExit(f"(system, T) groups with <20 frames: {small}")

    if args.out.exists():
        if not args.force:
            raise SystemExit(f"{args.out} exists; pass --force to overwrite")
        shutil.rmtree(args.out)
    args.out.mkdir(parents=True)

    split = assign(frames, args.test_frac, args.valid_frac, args.seed)
    ordered = sorted(frames, key=lambda f: (SPLITS.index(split[f]), f.system,
                                            f.temperature, f.segment, f.step))

    with open(args.out / "manifest.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["split", "system", "temperature", "segment", "frame", "path"])
        for f in ordered:
            w.writerow([split[f], f.system, f.temperature, f.segment, f.step, f.path])

    for name in SPLITS:
        tmp = args.out / f"{name}.xyz.gz.tmp"
        with gzip.open(tmp, "wt") as out:
            for f in ordered:
                if split[f] != name:
                    continue
                with gzip.open(args.root / f.path, "rt") as src:
                    text = src.read()
                out.write(text if text.endswith("\n") else text + "\n")
        tmp.rename(args.out / f"{name}.xyz.gz")

    counts = Counter(split.values())
    per_group = defaultdict(Counter)
    for f in frames:
        per_group[f"{f.system}_{f.temperature}K"][split[f]] += 1
    info = {
        "root": str(args.root), "seed": args.seed,
        "test_frac_per_segment": args.test_frac, "valid_frac_per_group": args.valid_frac,
        "counts": {s: counts[s] for s in SPLITS},
        "per_group": {k: {s: v[s] for s in SPLITS} for k, v in sorted(per_group.items())},
    }
    (args.out / "split_info.json").write_text(json.dumps(info, indent=2) + "\n")
    print(" ".join(f"{s}={counts[s]}" for s in SPLITS), f"-> {args.out}")


if __name__ == "__main__":
    main()
