#!/bin/bash
#SBATCH --job-name=sp_ablate
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:45:00
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.out

# Ablation: isolate ISPIN / LREAL / PREC's individual contribution to the
# matpes-recipe slowdown vs the old exact/tight/gold settings, on Cl-Li_1000K.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-openmpi
VASP=vasp_gam
export OMP_NUM_THREADS=1

mapfile -t CASES < runs/ablate_cases.txt
CASE="${CASES[$SLURM_ARRAY_TASK_ID]}"
d="runs/Cl-Li_1000K/$CASE"
echo "task $SLURM_ARRAY_TASK_ID -> $CASE  ($(date))"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" "$VASP" > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE ($(date))"
