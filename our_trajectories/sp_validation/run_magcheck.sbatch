#!/bin/bash
#SBATCH --job-name=sp_magcheck
#SBATCH --array=0-33
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:20:00
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.out

# ISPIN=1 vs ISPIN=2 single points across all 17 systems, to check whether
# any of them actually need spin polarization (energy/force differences,
# converged magnetic moments), rather than assuming non-magnetic.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1

mapfile -t CASES < runs/magcheck/cases.txt
CASE="${CASES[$SLURM_ARRAY_TASK_ID]}"
d="runs/magcheck/$CASE"
echo "task $SLURM_ARRAY_TASK_ID -> $CASE ($(date))"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE ($(date))"
