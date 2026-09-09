#!/bin/bash
#SBATCH -J ncsd_screen
#SBATCH -A matgen
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 00:25:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.out

module load python
cd "$SLURM_SUBMIT_DIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

conda run -n dft python -u run_screen.py \
    --temps 1000,1500,2000,2500 \
    --procs 128 \
    --out results/screen

echo "DONE $(date)"
