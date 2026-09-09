#!/bin/bash
#SBATCH -J ncsd_onsager
#SBATCH -A matgen
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 00:25:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.out

# Onsager coefficients for the systems in md_targets.txt that do not have them
# yet.  --merge keeps the rows already in results/onsager.csv for every other
# system, so this adds to the analysis rather than replacing it.
#
#   sbatch onsager.sbatch "Li-Nb-S,Li-Rb-S,..."

module load python
cd "$SLURM_SUBMIT_DIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

conda run -n dft python -u onsager.py \
    --systems "${1:?usage: sbatch onsager.sbatch <comma-separated systems>}" \
    --procs 48 \
    --merge

echo "DONE $(date)"
