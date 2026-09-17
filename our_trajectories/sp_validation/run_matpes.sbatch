#!/bin/bash
#SBATCH --job-name=sp_matpes
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.out

# NVT-MD timing test using the MatPESStaticSet-based recipe from
# ../nvt_md_input.py (see prep_matpes.py), for the same 6 frames as
# run.sbatch, to compare per-step wall time against the MPMorphMDSetGenerator
# exact/tight/gold tiers and extrapolate steps/7days.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-openmpi
VASP=vasp_gam
export OMP_NUM_THREADS=1

mapfile -t CASES < runs/cases.txt
CASE="${CASES[$SLURM_ARRAY_TASK_ID]}"
d="runs/$CASE/matpes"
echo "task $SLURM_ARRAY_TASK_ID -> $CASE  ($(date))"
echo "=== $d ==="
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" "$VASP" > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE ($(date))"
