#!/bin/bash
#SBATCH --job-name=sp_magck_cls
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:20:00
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.out

# Extra Cu-Li-S frames (1000K/1500K/2000K, in addition to the 2500K already
# checked) -- confirming non-magnetic wasn't a one-frame fluke.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1

mapfile -t CASES < runs/magcheck_CuLiS/cases.txt
CASE="${CASES[$SLURM_ARRAY_TASK_ID]}"
d="runs/magcheck_CuLiS/$CASE"
echo "task $SLURM_ARRAY_TASK_ID -> $CASE ($(date))"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE ($(date))"
