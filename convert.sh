#!/bin/bash
#SBATCH --job-name=convert-traj
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --time=0:30:00
#SBATCH --account=matgen
#SBATCH --output=logs/convert-traj_%A_%a.out
#SBATCH --error=logs/convert-traj_%A_%a.err
#SBATCH --array=5-10

if [ -n "$1" ]; then
    conda run -n dft python convert.py $SLURM_ARRAY_TASK_ID "$1"
else
    conda run -n dft python convert.py $SLURM_ARRAY_TASK_ID
fi