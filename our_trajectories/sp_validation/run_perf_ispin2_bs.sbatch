#!/bin/bash
#SBATCH --job-name=sp_ispin2_bs
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load vasp/cpu_6.4.2-intel
export OMP_NUM_THREADS=1
mapfile -t CASES < runs/perf_ispin2/cases_bils_lissb.txt
CASE="${CASES[$SLURM_ARRAY_TASK_ID]}"
d="runs/perf_ispin2/$CASE"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE ($(date))"
