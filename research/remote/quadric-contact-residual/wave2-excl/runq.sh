#!/bin/bash
# Load-gated sequential job queue.  Usage: runq.sh <queue-name> [--after PID] -- cmd1 ;; cmd2 ;; ...
# Waits until 1-min load < LOADMAX (default 140) (and PID gone) before each command; runs commands sequentially in the lane dir.
QN=$1; shift
AFTER=""
if [ "$1" = "--after" ]; then AFTER=$2; shift 2; fi
[ "$1" = "--" ] && shift
LOADMAX=${LOADMAX:-140}
cd /tmp/hodge/results/quadric-contact-residual
LOG=/tmp/hodge/jobs/$QN.log
echo "[runq $QN] start $(date -u +%FT%TZ) loadmax=$LOADMAX after=$AFTER" >> $LOG
IFS=';;' read -ra CMDS <<< "$*"
for c in "$@"; do :; done
# split on ' ;; '
CMDLIST=()
cur=""
for tok in "$@"; do
  if [ "$tok" = ";;" ]; then CMDLIST+=("$cur"); cur=""; else cur="$cur $tok"; fi
done
[ -n "$cur" ] && CMDLIST+=("$cur")
for c in "${CMDLIST[@]}"; do
  while :; do
    L=$(cut -d' ' -f1 /proc/loadavg | cut -d. -f1)
    if [ -n "$AFTER" ] && kill -0 $AFTER 2>/dev/null; then sleep 60; continue; fi
    if [ "$L" -lt "$LOADMAX" ]; then break; fi
    sleep 60
  done
  echo "[runq $QN] $(date -u +%FT%TZ) load=$(cut -d' ' -f1 /proc/loadavg) running:$c" >> $LOG
  OMP_NUM_THREADS=1 FLINT_NUM_THREADS=1 nice -n 5 bash -c "cd /tmp/hodge/results/quadric-contact-residual/wave2-excl && $c" >> $LOG 2>&1
  echo "[runq $QN] $(date -u +%FT%TZ) rc=$? done:$c" >> $LOG
done
echo "[runq $QN] ALL DONE $(date -u +%FT%TZ)" >> $LOG
