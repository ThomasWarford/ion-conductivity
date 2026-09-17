#!/bin/bash
#SBATCH --job-name=sp_gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --gpu-bind=map_gpu:0,1,2,3
#SBATCH --exclusive
#SBATCH --account=drxkp
#SBATCH --partition=gpu
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.out

# GPU (OpenACC) port comparison, same physics settings as the CPU
# ablate_baseline_ncore8 case (ISPIN=2, LREAL=False, PREC=Accurate, EDIFF=1e-5),
# on Cl-Li_1000K. Per https://vasp.at/wiki/GPU_ports_of_VASP :
#   - exactly 1 MPI rank per GPU (4 GPUs/node here -> ntasks-per-node=4)
#   - NCORE must be 1 (VASP forces this anyway; set explicitly)
#   - KPAR should match GPU count only if k-points divide evenly across it;
#     we're Gamma-only (NKPTS=1), so KPAR is left at its default (1)
#   - the gpu build's binaries are group-"vasp"-restricted (unlike the
#     world-readable cpu-openmpi build). `sg vasp -c` grants that, but
#     prompts for a password with no TTY inside a batch job -- so this
#     script must itself be *submitted* as `sg vasp -c "sbatch ..."`,
#     not call sg internally.
#   - the module's PATH puts a dummy /nopt/nrel/apps/vasp/noaccess ahead
#     of the real gpu build, so plain `vasp_gam` silently resolves to a
#     "you're not in the vasp group" stub -- call it by absolute path.
#   - NREL's own tested example launches with `mpirun -npernode N`, not
#     srun/pmi2, so we follow that rather than the cpu-openmpi recipe

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load vasp/gpu_6.4.2-nvidia
VASP_BIN=/nopt/nrel/apps/vasp/gpu/6.4.2/vasp_gam

d="runs/Cl-Li_1000K/gpu_baseline"
( cd "$d"
  rm -f OUTCAR OSZICAR vasprun.xml
  mpirun -npernode 4 "$VASP_BIN" > vasp.out 2>&1
  tail -n 5 OSZICAR
)
echo "done ($(date))"
