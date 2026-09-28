#!/bin/bash
# Submit an N-link aftercorr chain for a case manifest. Each link is one sbatch
# array submission; link k+1 depends on link k via --dependency=aftercorr.
set -euo pipefail
CASES_FILE="${1:?usage: submit_chain.sh CASES_FILE N_LINKS [PARTITION] [TIME] [QOS]}"
N_LINKS="${2:?usage: submit_chain.sh CASES_FILE N_LINKS [PARTITION] [TIME] [QOS]}"
PARTITION="${3:-standard}"
TIME="${4:-2-00:00:00}"
QOS="${5:-normal}"

N=$(wc -l < "$CASES_FILE")
PREV=""
for ((link=1; link<=N_LINKS; link++)); do
  CMD=(sbatch --partition="$PARTITION" --qos="$QOS" --time="$TIME"
       --array="0-$((N-1))" --export="ALL,CASES_FILE=$CASES_FILE,LINK=$link"
       --parsable)
  if [ -n "$PREV" ]; then
    CMD+=(--dependency="aftercorr:$PREV")
  fi
  CMD+=(our_trajectories/run_aimd.sbatch)
  JOBID=$(sg vasp -c "${CMD[*]}")
  echo "link $link -> job $JOBID"
  PREV="$JOBID"
done
