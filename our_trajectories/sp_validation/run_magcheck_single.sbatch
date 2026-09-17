#!/bin/bash
#SBATCH --job-name=sp_magcheck_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:20:00
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1
d="runs/magcheck/Li-O-P-W/ispin1"
echo "host: $(hostname)"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done ($(date))"
