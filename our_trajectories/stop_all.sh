#!/bin/bash
# Safely end the AIMD campaign ~1 day before the cluster closes. Two steps, both
# non-corrupting:
#   1. cancel any link-2 tasks that haven't started yet (nothing running to
#      interrupt, so plain scancel is safe)
#   2. gracefully stop everything still running via VASP's own STOPCAR mechanism
#      (LSTOP = .TRUE. finishes the current ionic step and writes a clean,
#      resumable CONTCAR/OUTCAR -- much safer than scancel on a live VASP process)
#
# Run `squeue -u "$USER"` first to see what this will affect.
set -euo pipefail
cd "$(dirname "$0")/.."

scancel --state=PENDING --name=aimd

for d in our_trajectories/runs/*/link*/; do
  if [ -f "$d/OSZICAR" ]; then
    echo "LSTOP = .TRUE." > "$d/STOPCAR"
    echo "stop requested: $d"
  fi
done
