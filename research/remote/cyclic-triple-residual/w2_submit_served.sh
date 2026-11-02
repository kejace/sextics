#!/bin/bash
# wait for the aux job to finish (keeps the 4-running-jobs budget), then queue the served-space sweep
cd /home/kejace/devel/sextics-work
L=cyclic-triple-residual
while true; do
  st=$(timeout 120 python3 -m hodge.rjob status --name search-$L-aux 2>/dev/null | grep -E "^(RUNNING|DONE)" | head -1)
  echo "[$(date +%T)] aux: $st"
  [ "$st" = "DONE" ] && break
  sleep 180
done
echo "[$(date +%T)] submitting served"
timeout 7800 python3 -m hodge.rjob submit --lane $L --name search-$L-served --cpus 8 --wait 7200 -- python3 w2_served.py w2-served-cands.jsonl 2147483647 8 ">" w2_served.out "2>" w2_served.err 2>&1 | grep -v "engine v2"
echo "[$(date +%T)] served exit ${PIPESTATUS[0]}"
