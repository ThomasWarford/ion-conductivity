#!/bin/bash
#SBATCH --job-name=aimd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=long
#SBATCH --qos=normal
#SBATCH --time=7-00:00:00
#SBATCH --exclude=c11-[6-29,36-55]
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out

# Production NVT-AIMD, one system/temperature/replica per array task. Finalized
# recipe: vasp_gam (Gamma-only), Intel build, 64 ranks/node, NCORE=32 (INCAR itself
# is our_trajectories/nvt_md_input.py + build_runs.py; parallelization tuning came
# from sp_validation/run_perf_final.sbatch).
#
# One 7-day job can't cover the full 2-week window, so every case runs two links
# back-to-back (LINK=1, then LINK=2 continuing from link1's CONTCAR), chained with
# --dependency=aftercorr so link2's task i starts as soon as link1's task i ends.
# Each link gets its own subdirectory so link1's OUTCAR/XDATCAR history is never
# overwritten by link2.
#
# Submit (module is group `vasp`-restricted, must wrap sbatch itself):
#   N=$(wc -l < our_trajectories/runs/cases_base.txt)
#   JOB1=$(sg vasp -c "sbatch --array=0-$((N-1)) --parsable our_trajectories/run_aimd.sbatch")
#   sg vasp -c "sbatch --array=0-$((N-1)) --export=ALL,LINK=2 --dependency=aftercorr:$JOB1 our_trajectories/run_aimd.sbatch"
#
# For cases_extra.txt, add --export=ALL,CASES_FILE=our_trajectories/runs/cases_extra.txt
# (and LINK=2 on the second call), optionally --qos=high.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1

CASES_FILE="${CASES_FILE:-our_trajectories/runs/cases_base.txt}"
LINK="${LINK:-1}"

mapfile -t CASES < "$CASES_FILE"
CASE="${CASES[$SLURM_ARRAY_TASK_ID]}"
base="our_trajectories/runs/$CASE"
d="$base/link$LINK"

mkdir -p "$d"
cp "$base"/{INCAR,KPOINTS,POTCAR} "$d/"
ln -sf /dev/null "$d/vasprun.xml"  # ~600KB/step, not needed -- OUTCAR/XDATCAR have the same data

if [ "$LINK" -gt 1 ]; then
  prev="$base/link$((LINK-1))/CONTCAR"
  [ -s "$prev" ] || { echo "missing $prev, nothing to continue"; exit 1; }
  cp "$prev" "$d/POSCAR"
else
  cp "$base/POSCAR" "$d/POSCAR"
fi

echo "task $SLURM_ARRAY_TASK_ID -> $CASE link$LINK ($(date))"
( cd "$d"
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE link$LINK ($(date))"
