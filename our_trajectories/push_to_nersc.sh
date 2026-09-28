#!/bin/bash
# Push ~/singlepoint_tars to NERSC. Safe to re-run: finished files are skipped,
# interrupted ones resume where they stopped.
#
#   tmux new -s push
#   ~/push_to_nersc.sh                                   # singlepoint tars
#   SRC=~/aimd_tars \\
#     DEST=/pscratch/sd/t/twarford/ion-conductivity/aimd_tars \\
#     PAR=2 ~/push_to_nersc.sh                           # AIMD tar
#
# Requires a live cert:  ~/bin/sshproxy -u twarford   (24 h lifetime)
#
# NOTE: /pscratch is PURGED at NERSC (files untouched for ~8 weeks are deleted).
# Once Swift goes away this is the only copy, so move it to CFS or HPSS before
# the purge window closes.

set -uo pipefail

# ---- defaults; override from the environment ------------------------------
NERSC_USER="${NERSC_USER:-twarford}"
DEST="${DEST:-/pscratch/sd/t/twarford/ion-conductivity/singlepoint_tars}"
# ---------------------------------------------------------------------------

HOST=dtn01.nersc.gov
SRC="${SRC:-$HOME/singlepoint_tars}"
KEY="$HOME/.ssh/nersc"
PAR="${PAR:-6}"                      # concurrent rsyncs; raise if the link isn't saturated

SSH_CMD="ssh -i $KEY -o Compression=no -c aes128-gcm@openssh.com"
# -a keeps mtimes (without them a re-run re-sends everything instead of skipping)
# --inplace       write straight into the target, no 228 GB temp copy + rename
# --partial       keep what transferred if the link drops
# --append-verify on resume, checksum the existing prefix and append rather than
#                 re-deltaing the whole file -- safe only because these archives
#                 never change; drop it if you ever regenerate one
# no -z: the contents are already gzip/h5, so compression is pure CPU burn
RS=(-a --inplace --partial --append-verify --human-readable)

cd "$SRC" || exit 1
mkdir -p logs

echo "=== preflight $(date) ==="
[ -f "$KEY" ] || { echo "no cert at $KEY -- run: ~/bin/sshproxy -u $NERSC_USER"; exit 1; }
$SSH_CMD "$NERSC_USER@$HOST" "mkdir -p '$DEST' && echo connected: \$(hostname)" || {
  echo "cannot reach $HOST -- is the cert current? (24 h lifetime)"; exit 1; }

echo
echo "=== small files first (manifests, checksums, logs) ==="
# nullglob so an unmatched pattern vanishes instead of being passed to rsync
# literally (aimd_tars has no *.log, which would otherwise abort the run).
shopt -s nullglob
small=(SHA256SUMS *.manifest *.sha256 *.log)
shopt -u nullglob
if [ ${#small[@]} -gt 0 ]; then
  rsync "${RS[@]}" --info=progress2 -e "$SSH_CMD" \
        "${small[@]}" "$NERSC_USER@$HOST:$DEST/" || exit 1
else
  echo "  (no small files to send)"
fi

echo
echo "=== $(ls *.tar | wc -l) archives, $PAR at a time ==="
for f in *.tar; do
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 5; done
  (
    if rsync "${RS[@]}" -e "$SSH_CMD" "$f" "$NERSC_USER@$HOST:$DEST/" \
         > "logs/$f.rsync.log" 2>&1; then
      echo "  done   $f ($(date +%H:%M:%S))"
    else
      echo "  FAILED $f -- see logs/$f.rsync.log"
    fi
  ) &
  echo "  started $f ($(date +%H:%M:%S))"
done
wait

echo
echo "=== remote verification (re-reads every archive at NERSC) ==="
$SSH_CMD "$NERSC_USER@$HOST" "cd '$DEST' && sha256sum -c SHA256SUMS"
rc=$?

echo
if [ $rc -eq 0 ]; then
  echo "ALL ARCHIVES VERIFIED AT NERSC $(date)"
else
  echo "VERIFICATION FAILED -- re-run this script; it will resume the bad files"
fi
exit $rc
