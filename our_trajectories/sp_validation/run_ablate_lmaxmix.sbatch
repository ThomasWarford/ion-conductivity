#!/bin/bash
#SBATCH --job-name=sp_ablate_lm
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.out

# Isolate LMAXMIX (charge-density mixing angular-momentum cutoff) alone,
# on top of the matpes+NCORE8 baseline (ISPIN=2, LREAL=False, PREC=Accurate,
# EDIFF=1e-5), using vasp_gam. Old exact/tight/gold settings used LMAXMIX=2;
# matpes/MatPESStaticSet defaults to 6 (meant for f-electron systems, none
# of which are present in Cl-Li_1000K).

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-openmpi
export OMP_NUM_THREADS=1

CASES=(ablate_lmaxmix4 ablate_lmaxmix2)
CASE="${CASES[$SLURM_ARRAY_TASK_ID]}"
d="runs/Cl-Li_1000K/$CASE"
echo "task $SLURM_ARRAY_TASK_ID -> $CASE  ($(date))"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" vasp_gam > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done $CASE ($(date))"
