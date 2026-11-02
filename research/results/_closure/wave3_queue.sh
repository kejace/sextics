#!/bin/bash
# wave3_queue.sh TAG THREADS QUEUEFILE [MAX_WALL_S]
# Sequential queue runner (ONE Gröbner computation at a time). Queue lines: "<timeout_s> <min_avail_gb> <file.ms> [maxpairs]"
# ('#' comments allowed). For each entry: skip if a receipt with a decided verdict exists; wait (<= 2 h) until the
# box has >= min_avail_gb available memory (so a big F4 step never pushes the box under the 45 GB guard floor);
# then run wave3_run.sh on it (with WAVE3_MAXPAIRS=maxpairs when given). Entries skipped for memory are retried in a
# second and third pass. Stops after MAX_WALL_S seconds (default 12 h).
TAG=$1; TH=$2; Q=$3; MAXW=${4:-43200}
REC=/tmp/hodge/results/_closure/wave3-receipts.jsonl
RUN=/tmp/hodge/results/_closure/wave3_run.sh
T0=$(date +%s)
decided() { grep -F "\"file\":\"$1\"" "$REC" 2>/dev/null | grep -qE 'ONLY-ORIGIN|POSITIVE-DIMENSIONAL'; }
avail() { free -g | awk '/Mem/{print $7}'; }
for PASS in 1 2 3; do
  echo "[queue $TAG] pass $PASS $(date -u +%FT%TZ)"
  while read -r TO MINAV F MP; do
    case "$TO" in ''|\#*) continue;; esac
    [ -f "$F" ] || { echo "[queue $TAG] missing $F"; continue; }
    decided "$F" && { echo "[queue $TAG] already decided: $F"; continue; }
    NOW=$(date +%s); [ $((NOW - T0)) -gt $MAXW ] && { echo "[queue $TAG] max wall reached"; exit 0; }
    W=0
    while [ "$(avail)" -lt "$MINAV" ]; do
      [ $W -ge 7200 ] && break
      [ $W -eq 0 ] && echo "[queue $TAG] waiting for >= $MINAV GB available (now $(avail) GB) before $F"
      sleep 120; W=$((W + 120))
    done
    if [ "$(avail)" -lt "$MINAV" ]; then echo "[queue $TAG] still only $(avail) GB available; deferring $F"; continue; fi
    echo "[queue $TAG] avail $(avail) GB; running $F (timeout $TO s, maxpairs ${MP:-unlimited})"
    if [ -n "$MP" ]; then WAVE3_MAXPAIRS=$MP bash "$RUN" "$TAG" "$TH" "$TO" "$F"; else bash "$RUN" "$TAG" "$TH" "$TO" "$F"; fi
  done < "$Q"
done
echo "[queue $TAG] queue exhausted $(date -u +%FT%TZ)"
