#!/bin/bash
#SBATCH --job-name=sp_intel_dbg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=00:05:00
#SBATCH --output=logs/sp_intel_dbg.out
module purge
module load vasp/cpu_6.4.2-intel
which vasp_gam vasp_std
echo "---"
mpirun --version 2>&1 | head -3
