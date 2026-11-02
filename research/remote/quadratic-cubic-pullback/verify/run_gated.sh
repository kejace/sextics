#!/bin/bash
# usage: run_gated.sh <name> <cmd...>   -- waits for load < 140 then runs cmd, logging to /tmp/hodge/jobs/<name>.log
name=$1; shift
while [ "$(cut -d. -f1 /proc/loadavg)" -ge 140 ]; do sleep 60; done
echo "[start $(date -u +%FT%TZ) load $(cut -d' ' -f1 /proc/loadavg)]"
"$@"
echo "[end $(date -u +%FT%TZ) rc=$?]"
