#!/bin/bash
#SBATCH --job-name=chargemol_sp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out

# Charge-partitioning single points: one array task processes a whole batch of
# case directories (from build_chargemol_frames.py) sequentially -- each case is
# a single-point VASP run (NSW=0, LAECHG/LCHARG on) followed by chargemol/DDEC6
# (run_chargemol.py), which lossily compresses CHGCAR/AECCAR1/AECCAR2 via voltools
# and gzips the rest of the case dir when done. Same node/MPI layout as
# run_aimd.sbatch (vasp_gam, Intel build, 64 ranks/node, NCORE=32) since it's the
# same ~100-atom Gamma-only cell size; --time/--partition are single-point-scale,
# not the 7-day AIMD scale, and will be retuned from the pilot's measured
# wall-time/frame (see plan).
#
# Submit (module is group `vasp`-restricted, must wrap sbatch itself):
#   N=$(ls our_trajectories/our_singlepoints/cases_pilot_batch_*.txt | wc -l)
#   sg vasp -c "sbatch --array=0-$((N-1)) our_trajectories/run_chargemol_sp.sbatch"
#
# For the full run, add --export=ALL,TAG=full (reads cases_full_batch_<i>.txt).

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1  # vasp_gam is MPI-only

TAG="${TAG:-pilot}"
CASES_FILE="our_trajectories/our_singlepoints/cases_${TAG}_batch_${SLURM_ARRAY_TASK_ID}.txt"
mapfile -t CASES < "$CASES_FILE"

echo "task $SLURM_ARRAY_TASK_ID: ${#CASES[@]} cases from $CASES_FILE ($(date))"
for CASE in "${CASES[@]}"; do
  CASE_PATH="our_trajectories/$CASE"
  echo "-- $CASE_PATH ($(date))"
  ( cd "$CASE_PATH"
    srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
    tail -n 3 OSZICAR
  )
  OMP_NUM_THREADS=64 conda run -p ~/.conda-envs/dft python our_trajectories/run_chargemol.py "$CASE_PATH"
done
echo "done task $SLURM_ARRAY_TASK_ID ($(date))"
