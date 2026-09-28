#!/bin/bash
#SBATCH --job-name=tar_singlepoints
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --qos=normal
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out

# Pack our_singlepoints/ into one tar per compound, staged for rsync off-cluster.
# Originals are left completely untouched -- this only reads.
#
# Tarring here is bound by per-file NFS metadata latency, not CPU or bandwidth
# (a timed Li-P run: 90 GiB / 72k files in 22m50s, of which only 6s was user
# CPU), so the win comes from running many tars at once. Partition `standard` is
# OverSubscribe=NO, so a task gets a whole node regardless -- rather than waste
# one node per archive, each of the 3 array tasks runs its share concurrently in
# the background, spreading the load over 3 NFS clients.
#
# Submit:
#   sbatch --array=0-2 our_trajectories/run_tar.sbatch

set -euo pipefail

ROOT="$SLURM_SUBMIT_DIR/our_trajectories"
OUT="${TAR_OUT:-$HOME/singlepoint_tars}"
NTASKS_TOTAL=3

mkdir -p "$OUT"
cd "$ROOT"

# Largest first (sizes measured by sampling, see plan) so the round-robin split
# below balances -- the biggest archive sets the critical path.
ARCHIVES=(
  Li-Rb-S Li-Rb-Se Al-Li-Na-P Bi-Li-S Cu-Li-S Li-Nb-S F-Li-P
  Ge-Li-Sn Cl-Li Li-S-Sb Li-P Li-O-P-Ti Li-O-Zr misc
)

members_for() {
  case "$1" in
    misc) printf '%s\n' our_singlepoints/cases_*.txt \
            our_singlepoints/_pilot_gzip_baseline pilot_archive ;;
    *)    printf '%s\n' "our_singlepoints/$1" ;;
  esac
}

make_archive() {
  local name="$1" log="$OUT/$1.log"
  local -a members
  mapfile -t members < <(members_for "$name")

  {
    echo "=== $name: $(date) ==="
    echo "members: ${members[*]}"

    # One pass over the data produces all three outputs: --index-file gives the
    # member listing for free, and the tee'd checksum costs only idle CPU rather
    # than a second read of ~100 GiB.
    # -v is what generates the listing; --index-file only redirects it (without
    # -v the index file comes out empty).
    tar -cvf - --index-file="$OUT/$name.manifest" -C "$ROOT" "${members[@]}" \
      | tee >(sha256sum | awk '{print $1}' > "$OUT/$name.sha256") \
      > "$OUT/$name.tar"

    # A zero exit status alone would not catch a truncated archive, so compare
    # what tar actually wrote against what is on disk. tar indexes every member
    # (files and dirs), which is exactly what find counts.
    local n_tar n_disk
    n_tar=$(wc -l < "$OUT/$name.manifest")
    n_disk=$(find "${members[@]}" | wc -l)
    echo "members archived=$n_tar on-disk=$n_disk"
    if [ "$n_tar" -ne "$n_disk" ]; then
      echo "MISMATCH $name: archived $n_tar, expected $n_disk"
      return 1
    fi

    echo "size=$(du -h "$OUT/$name.tar" | cut -f1)  sha256=$(cat "$OUT/$name.sha256")"
    echo "=== $name OK: $(date) ==="
  } >> "$log" 2>&1
}

declare -a NAMES=() PIDS=()
for i in "${!ARCHIVES[@]}"; do
  if [ $((i % NTASKS_TOTAL)) -eq "$SLURM_ARRAY_TASK_ID" ]; then
    name="${ARCHIVES[$i]}"
    make_archive "$name" &
    NAMES+=("$name")
    PIDS+=($!)
  fi
done

echo "task $SLURM_ARRAY_TASK_ID: ${#NAMES[@]} archives -> $OUT (${NAMES[*]}) $(date)"

# `set -e` does not fail the script on a background job, so check each one.
rc=0
for i in "${!PIDS[@]}"; do
  if wait "${PIDS[$i]}"; then
    echo "OK   ${NAMES[$i]}"
  else
    echo "FAIL ${NAMES[$i]} (see $OUT/${NAMES[$i]}.log)"
    rc=1
  fi
done

echo "done task $SLURM_ARRAY_TASK_ID ($(date))"
exit $rc
