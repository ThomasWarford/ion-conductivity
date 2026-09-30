#!/bin/bash
#SBATCH --job-name=isolated_atoms
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out

# Isolated-atom E0 single points (isolated_atoms.py build), run sequentially in
# one job -- each is one atom, Gamma-only, minutes apiece. Same node/MPI layout
# as run_chargemol_sp.sbatch (vasp_gam, Intel build, 64 ranks/node). No
# labeling/cleanup step: outputs are left in place for isolated_atoms.py collect.
#
# Submit from the repo root (module is group `vasp`-restricted, must wrap sbatch):
#   sg vasp -c "sbatch our_trajectories/run_isolated_atoms.sbatch"

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1  # vasp_gam is MPI-only
ulimit -s unlimited  # the 21 A box segfaulted in MPI_Allreduce with the default stack

CASES_FILE="${CASES_FILE:-our_trajectories/isolated_atoms/cases.txt}"
mapfile -t CASES < "$CASES_FILE"

echo "${#CASES[@]} cases ($(date))"
for CASE in "${CASES[@]}"; do
  CASE_PATH="our_trajectories/$CASE"
  echo "-- $CASE_PATH ($(date))"
  ( cd "$CASE_PATH"
    srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1 || echo "vasp_gam failed"
    tail -n 2 OSZICAR
  )
done
echo "done ($(date))"
