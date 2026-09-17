#!/bin/bash
#SBATCH --job-name=sp_perf_intel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.out

# Isolate compiler/MPI stack: vasp/cpu_6.4.2-intel (Intel oneAPI + Intel MPI
# + MKL, -march=core-avx2) vs the cpu-openmpi (GCC+OpenMPI) build used
# everywhere else, same physics/parallel settings (64 ranks/node, NCORE=8).

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1

d="runs/Cl-Li_1000K/perf_intel_64r_ncore8"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done ($(date))"
