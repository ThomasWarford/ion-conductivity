"""Download MP amorphous_diffusivity (NCSD) trajectories from S3.

The bucket is public and answers plain REST, so this uses only the standard
library -- no boto3, no AWS credentials, no `aws` CLI.  It mirrors the bucket
layout exactly:

    <dest>/trajectories/temperature=<T>K/chemical_system=<S>/dt=<...>.jsonl.gz

Downloads are resumable: a file whose size already matches the size reported by
the bucket listing is skipped, so re-running is a cheap integrity check.  Each
file lands as `<name>.part` and is renamed only once complete, so an interrupted
run never leaves something that looks finished.

  python fetch_ncsd.py --dry-run          # what would be fetched, and how much
  python fetch_ncsd.py                    # the 12 MD targets, 1000-2500 K
  python fetch_ncsd.py --systems Li-S,Li-Se --temps 1500
  python fetch_ncsd.py --systems all      # the whole dataset (74 GB w/o 5000 K)

Default systems come from md_targets.txt; default destination is
$SCRATCH/data/ncsd (symlinked in this repo as ../data_ncsd).

5000 K is supported but not in the default temperature list: those trajectories
are roughly 5x larger than the others (Li-S alone is 633 MB).
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed

BUCKET = "https://materialsproject-contribs.s3.amazonaws.com"
ROOT = "amorphous_diffusivity/trajectories"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"
HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS = os.path.join(HERE, "md_targets.txt")
DEFAULT_TEMPS = (1000, 1500, 2000, 2500)
RETRIES = 4


def read_targets(path=TARGETS):
    """The committed list of systems, ignoring blank lines and # comments."""
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                out.append(line)
    return out


def _get(url, tries=RETRIES):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 ** k)


def _list_prefix(prefix):
    """Every (key, size) under a bucket prefix, following continuation tokens."""
    out, token = [], None
    while True:
        q = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            q["continuation-token"] = token
        root = ET.fromstring(_get(f"{BUCKET}/?{urllib.parse.urlencode(q)}"))
        for c in root.findall(f"{NS}Contents"):
            out.append((c.findtext(f"{NS}Key"), int(c.findtext(f"{NS}Size"))))
        if root.findtext(f"{NS}IsTruncated") == "true":
            token = root.findtext(f"{NS}NextContinuationToken")
        else:
            return out


def list_systems(temp):
    """Every chemical_system present at one temperature."""
    out, token = [], None
    prefix = f"{ROOT}/temperature={temp}K/"
    while True:
        q = {"list-type": "2", "prefix": prefix, "delimiter": "/", "max-keys": "1000"}
        if token:
            q["continuation-token"] = token
        root = ET.fromstring(_get(f"{BUCKET}/?{urllib.parse.urlencode(q)}"))
        for p in root.findall(f"{NS}CommonPrefixes"):
            out.append(p.findtext(f"{NS}Prefix").rstrip("/").split("chemical_system=")[-1])
        if root.findtext(f"{NS}IsTruncated") == "true":
            token = root.findtext(f"{NS}NextContinuationToken")
        else:
            return sorted(out)


def list_keys(temps, systems, procs=8):
    """(key, size) for every requested (system, temperature).

    Raises if a requested system has no objects at some temperature -- that is
    how a typo, or a system that simply is not in the dataset, gets caught
    before anything is downloaded.
    """
    jobs = [(T, s) for T in temps for s in systems]
    found, missing = [], []
    with ThreadPoolExecutor(procs) as ex:
        futs = {ex.submit(_list_prefix,
                          f"{ROOT}/temperature={T}K/chemical_system={s}/"): (T, s)
                for T, s in jobs}
        for f in as_completed(futs):
            T, s = futs[f]
            got = f.result()
            if got:
                found.extend(got)
            else:
                missing.append((T, s))
    if missing:
        bad = sorted({s for _, s in missing})
        msg = ["not in the bucket: " + ", ".join(bad)]
        for s in bad:
            ts = sorted(T for T, x in missing if x == s)
            msg.append(f"  {s}: no objects at {', '.join(f'{T} K' for T in ts)}")
        msg.append(f"\nSystems available at {temps[0]} K:")
        msg.append("  " + " ".join(list_systems(temps[0])))
        raise SystemExit("\n".join(msg))
    return sorted(found)


def local_path(dest, key):
    return os.path.join(dest, key[len("amorphous_diffusivity/"):])


def download(dest, key, size):
    """Fetch one key; returns 'skip' if the local copy is already complete."""
    out = local_path(dest, key)
    if os.path.exists(out) and os.path.getsize(out) == size:
        return "skip"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    part = out + ".part"
    url = f"{BUCKET}/{urllib.parse.quote(key)}"
    for k in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(part, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            break
        except Exception:
            if k == RETRIES - 1:
                raise
            time.sleep(2 ** k)
    got = os.path.getsize(part)
    if got != size:
        os.remove(part)
        raise IOError(f"{key}: got {got} bytes, expected {size}")
    os.replace(part, out)
    return "ok"


def human(n):
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or u == "TB":
            return f"{n:.1f} {u}"
        n /= 1024


def default_dest():
    scratch = os.environ.get("SCRATCH")
    if scratch:
        return os.path.join(scratch, "data", "ncsd")
    return os.path.abspath(os.path.join(HERE, "..", "data_ncsd"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--systems", default=None,
                    help="comma-separated chemical_system tags, or 'all' "
                         f"(default: the {len(read_targets())} in md_targets.txt)")
    ap.add_argument("--temps", default=",".join(str(t) for t in DEFAULT_TEMPS),
                    help="comma-separated temperatures in K (5000 is supported "
                         "but not included by default)")
    ap.add_argument("--dest", default=default_dest())
    ap.add_argument("--procs", type=int, default=16, help="parallel downloads")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    temps = [int(t) for t in args.temps.split(",") if t.strip()]
    if args.systems is None:
        systems = read_targets()
        src = os.path.relpath(TARGETS, os.getcwd())
    elif args.systems.strip().lower() == "all":
        systems = list_systems(temps[0])
        src = f"every system at {temps[0]} K"
    else:
        systems = [s.strip() for s in args.systems.split(",") if s.strip()]
        src = "--systems"

    print(f"{len(systems)} systems ({src}) x {len(temps)} temperatures "
          f"{temps}\ndest: {args.dest}", flush=True)
    if 5000 in temps:
        print("NOTE: 5000 K trajectories are ~5x larger than the others.", flush=True)

    keys = list_keys(temps, systems)
    have = [(k, s) for k, s in keys
            if os.path.exists(local_path(args.dest, k))
            and os.path.getsize(local_path(args.dest, k)) == s]
    todo = [(k, s) for k, s in keys if (k, s) not in set(have)]
    print(f"{len(keys)} files, {human(sum(s for _, s in keys))} total; "
          f"{len(have)} already present, {len(todo)} to fetch "
          f"({human(sum(s for _, s in todo))})", flush=True)

    if args.dry_run:
        for k, s in todo:
            print(f"  {human(s):>9}  {k[len('amorphous_diffusivity/'):]}")
        return
    if not todo:
        print("nothing to do")
        return

    t0, done, failed = time.time(), 0, []
    with ThreadPoolExecutor(args.procs) as ex:
        futs = {ex.submit(download, args.dest, k, s): k for k, s in todo}
        for f in as_completed(futs):
            k = futs[f]
            try:
                f.result()
            except Exception as e:                       # keep going, report at end
                failed.append((k, repr(e)))
            done += 1
            if done % 5 == 0 or done == len(todo):
                print(f"  {done}/{len(todo)}  ({time.time() - t0:.0f}s)", flush=True)

    print(f"done in {time.time() - t0:.0f}s")
    if failed:
        print(f"{len(failed)} FAILED (re-run to retry):", file=sys.stderr)
        for k, e in failed:
            print(f"  {k}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
