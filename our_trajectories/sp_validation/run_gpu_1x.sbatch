#!/bin/bash
#SBATCH --job-name=sp_gpu_1x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --gpu-bind=closest
#SBATCH --account=drxkp
#SBATCH --partition=gpu
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.out

# Single-GPU comparison (no --exclusive, only gres=gpu:1 -- can share the
# node with other users' jobs on the other 3 GPUs), same physics settings
# as gpu_baseline (ISPIN=2, LREAL=False, PREC=Accurate, EDIFF=1e-5, NCORE=1).

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/gpu_6.4.2-nvidia
VASP_BIN=/nopt/nrel/apps/vasp/gpu/6.4.2/vasp_gam

echo "host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi -L
nvidia-smi --query-gpu=index,uuid,compute_mode,memory.used --format=csv

d="runs/Cl-Li_1000K/gpu_1gpu"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  mpirun -npernode 1 "$VASP_BIN" > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done ($(date))"
