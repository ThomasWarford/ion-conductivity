#!/bin/bash
#SBATCH --job-name=sp_perf_intel_sw
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out

# NCORE sweep on vasp/cpu_6.4.2-intel (Intel oneAPI + Intel MPI + MKL),
# same 64 ranks/node (1 rank/physical core) as the openmpi sweep, to see
# whether the optimal NCORE differs by compiler/MPI stack.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1

d="runs/Cl-Li_1000K/${CASE_DIR}"
echo "CASE_DIR=$CASE_DIR ($(date))"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE_DIR ($(date))"
