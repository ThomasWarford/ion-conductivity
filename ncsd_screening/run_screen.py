"""Parallel driver: analyse every NCSD trajectory file and write results.

  python run_screen.py [--temps 1000,1500,2000,2500] [--procs N]
                       [--limit N] [--out results/screen]

Writes:
  <out>.jsonl   one JSON row per (system, temperature, species), appended live
  <out>.csv     the same, as a table (written at the end)
Re-running skips (system, temperature) pairs already present in <out>.jsonl.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from multiprocessing import Pool

DATA = os.path.join(os.path.dirname(__file__), "..", "data_ncsd", "trajectories")


def find_files(temps):
    out = []
    for T in temps:
        out += sorted(glob.glob(os.path.join(
            DATA, f"temperature={T}K", "chemical_system=*", "*.jsonl.gz")))
    return out


def _run(path):
    from worker import analyse_file
    t0 = time.time()
    try:
        rows = analyse_file(path)
    except Exception as e:  # noqa: BLE001
        import traceback
        rows = [dict(path=path, error=f"{e}", tb=traceback.format_exc())]
    for r in rows:
        r.setdefault("wall_s", round(time.time() - t0, 1))
    return path, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--temps", default="1000,1500,2000,2500")
    ap.add_argument("--procs", type=int, default=len(os.sched_getaffinity(0)))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "results", "screen"))
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    jsonl = args.out + ".jsonl"

    done = set()
    if os.path.exists(jsonl):
        for ln in open(jsonl):
            try:
                r = json.loads(ln)
                done.add((r.get("system"), r.get("temperature")))
            except json.JSONDecodeError:
                pass

    temps = [t.strip() for t in args.temps.split(",") if t.strip()]
    files = find_files(temps)

    def key(p):
        return (os.path.basename(os.path.dirname(p)).split("chemical_system=")[-1],
                float(os.path.basename(os.path.dirname(os.path.dirname(p))).split("temperature=")[-1].rstrip("K")))
    files = [p for p in files if key(p) not in done]
    if args.limit:
        files = files[:args.limit]

    print(f"{len(files)} files to process on {args.procs} procs "
          f"({len(done)} already done)", flush=True)
    t0 = time.time()
    n = 0
    with open(jsonl, "a", buffering=1) as fh, Pool(args.procs) as pool:
        for path, rows in pool.imap_unordered(_run, files, chunksize=1):
            for r in rows:
                fh.write(json.dumps(r) + "\n")
            n += 1
            if n % 25 == 0 or n == len(files):
                el = time.time() - t0
                print(f"  {n}/{len(files)}  {el:.0f}s  "
                      f"eta {el / n * (len(files) - n):.0f}s", flush=True)

    # consolidate to CSV
    import pandas as pd
    rows = [json.loads(ln) for ln in open(jsonl)]
    df = pd.DataFrame(rows).sort_values(["system", "temperature", "species"])
    df.to_csv(args.out + ".csv", index=False)
    nerr = df["error"].astype(bool).sum() if "error" in df else 0
    print(f"wrote {args.out}.csv  ({len(df)} rows, {nerr} errors)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
