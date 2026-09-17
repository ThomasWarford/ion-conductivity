#!/bin/bash
#SBATCH --job-name=sp_validation
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=standard
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.out

# Single-point replication of the MP amorphous_diffusivity (NCSD) AIMD settings.
# Inputs built by prep.py; one array task per frame, runs exact/ then tight/.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/cpu_6.4.2-openmpi          # world-readable build, no vasp-group gate
VASP=vasp_gam

# Without this, each of the 128 MPI ranks spawns one OpenMP thread per
# node core (128), oversubscribing the node ~128x and stalling the SCF
# loop to a crawl (observed: 2 electronic steps in 70 min).
export OMP_NUM_THREADS=1

mapfile -t CASES < runs/cases.txt
CASE="${CASES[$SLURM_ARRAY_TASK_ID]}"
echo "task $SLURM_ARRAY_TASK_ID -> $CASE  ($(date))"

for tier in exact tight gold; do
    d="runs/$CASE/$tier"
    echo "=== $d ==="
    ( cd "$d"
      rm -f OUTCAR OSZICAR vasprun.xml
      srun --mpi=pmi2 --cpu-bind=cores -n "$SLURM_NTASKS" "$VASP" > vasp.out 2>&1
      tail -n 3 OSZICAR
    )
done
echo "done $CASE ($(date))"
