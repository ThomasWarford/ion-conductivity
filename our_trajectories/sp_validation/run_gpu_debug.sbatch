#!/bin/bash
#SBATCH --job-name=sp_gpu_dbg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --account=drxkp
#SBATCH --partition=gpu
#SBATCH --time=00:05:00
#SBATCH --output=logs/sp_gpu_dbg.out

echo "START"
cd "$SLURM_SUBMIT_DIR"
echo "PWD: $(pwd)"
module purge
echo "purged"
module load vasp/gpu_6.4.2-nvidia
echo "module loaded, exit=$?"
which vasp_gam vasp_std mpirun
echo "END"
