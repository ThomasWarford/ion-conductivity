#!/bin/bash
#SBATCH --job-name=sp_perf
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out

# Generic perf-sweep runner. --nodes/--ntasks-per-node/--cpus-per-task are
# passed at submission time (they vary per sweep point); CASE_DIR names
# which runs/Cl-Li_1000K/<dir> to run vasp_gam in.
#
# Node topology (Swift `standard` partition): 2x AMD EPYC 7532 (Rome),
# 32 physical cores/socket, SMT-2 -> 64 physical cores / 128 logical CPUs
# per node. Per https://vasp.at/wiki/Optimizing_the_parallelization and
# https://vasp.at/wiki/NCORE: NCORE ~ sqrt(ranks), ~100-atom systems favor
# NCORE=4, and VASP does not benefit from SMT oversubscription.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-openmpi
export OMP_NUM_THREADS=1

d="runs/Cl-Li_1000K/${CASE_DIR}"
echo "CASE_DIR=$CASE_DIR nodes=$SLURM_NNODES ntasks=$SLURM_NTASKS ($(date))"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE_DIR ($(date))"
