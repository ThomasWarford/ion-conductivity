#!/bin/bash
#SBATCH --job-name=sp_vaspgam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.out

# Isolate vasp_gam (real-arithmetic, Gamma-only build) vs vasp_std on the
# exact same matpes+NCORE8 INCAR/POSCAR/POTCAR/KPOINTS (ablate_baseline_ncore8
# duplicated into ablate_vaspgam), to see if we've been leaving speed on the
# table by running the general-k-point binary for a Gamma-only calculation.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-openmpi
export OMP_NUM_THREADS=1

d="runs/Cl-Li_1000K/ablate_vaspgam"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done ($(date))"
