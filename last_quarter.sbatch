#!/bin/bash
#SBATCH --job-name=last-quarter
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --qos=normal
#SBATCH --time=1:00:00
#SBATCH --account=drxkp
#SBATCH --exclude=c11-13
#SBATCH --mem=8G
#SBATCH --output=logs/last-quarter_%A_%a.out
#SBATCH --error=logs/last-quarter_%A_%a.err
#SBATCH --array=0-17,19-23,25-33,35,45,59%1

conda run -n dft python last_quarter.py $SLURM_ARRAY_TASK_ID
