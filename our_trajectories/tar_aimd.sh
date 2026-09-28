#!/bin/bash
# Pack the AIMD output into a single archive for transfer off-cluster:
#   runs/              production NVT-AIMD, 4.06 ns across 13 compounds (51 GB)
#   sp_validation/runs/ ISPIN=2 + parallelisation validation (26 GB)
# Originals are left untouched -- this only reads.
#
# Unlike the singlepoint tar (940k files, metadata-bound, needed a Slurm job)
# this is ~4,300 large mostly-gzipped files, so it streams and runs on the
# login node in minutes.
#
#   ~/ion-conductivity/our_trajectories/tar_aimd.sh

set -uo pipefail

ROOT="/home/twarford/ion-conductivity/our_trajectories"
OUT="${TAR_OUT:-$HOME/aimd_tars}"
NAME=aimd_runs
MEMBERS=(runs sp_validation/runs)

mkdir -p "$OUT"
cd "$ROOT" || exit 1

echo "=== $NAME: $(date) ==="
echo "members: ${MEMBERS[*]}"
du -shc "${MEMBERS[@]}" 2>/dev/null | tail -1

# One pass over the data yields archive, manifest and checksum together:
# --index-file gets the member listing for free, and the tee'd sha256 costs
# only idle CPU rather than a second read. -v is what generates the listing --
# --index-file merely redirects it, and without -v the file comes out empty.
tar -cvf - --index-file="$OUT/$NAME.manifest" -C "$ROOT" "${MEMBERS[@]}" \
  | tee >(sha256sum | awk '{print $1}' > "$OUT/$NAME.sha256") \
  > "$OUT/$NAME.tar"
rc=$?
[ $rc -eq 0 ] || { echo "tar pipeline failed (rc=$rc)"; exit 1; }

# A zero exit status would not catch a truncated archive, so compare what tar
# wrote against what is on disk. tar indexes every member, files and dirs,
# which is exactly what find counts.
n_tar=$(wc -l < "$OUT/$NAME.manifest")
n_disk=$(find "${MEMBERS[@]}" | wc -l)
echo "members archived=$n_tar on-disk=$n_disk"
if [ "$n_tar" -ne "$n_disk" ]; then
  echo "MISMATCH: archived $n_tar, expected $n_disk"
  exit 1
fi

# sha256sum -c needs "hash<space><space>filename"; the pipeline above writes a
# bare hash, so build the checkable form here rather than after the fact.
printf '%s  %s\n' "$(cat "$OUT/$NAME.sha256")" "$NAME.tar" > "$OUT/SHA256SUMS"

echo "size=$(du -h "$OUT/$NAME.tar" | cut -f1)  sha256=$(cat "$OUT/$NAME.sha256")"
echo "=== $NAME OK: $(date) ==="
