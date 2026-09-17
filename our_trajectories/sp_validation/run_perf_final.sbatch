#!/bin/bash
#SBATCH --job-name=sp_perf_final
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out

# Finalized recipe (Intel build, NCORE=32, 64 ranks/node, ISPIN=1,
# LREAL=False, PREC=Accurate, EDIFF=1e-5, ENCUT=520) on CASE_DIR, to get
# real per-system throughput instead of extrapolating from old-settings
# ratios.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1

d="runs/${CASE_DIR}/perf_final"
echo "CASE_DIR=$CASE_DIR ($(date))"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE_DIR ($(date))"
