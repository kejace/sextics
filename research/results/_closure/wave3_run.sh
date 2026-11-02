#!/bin/bash
# wave3_run.sh TAG THREADS TIMEOUT_S FILE [FILE ...]
# Runs `msolve -g 1` (leading ideal only) on each homogeneous-ideal input FILE (mod p, '^' syntax), one after the
# other, writing <dir>/wave3/<base>.wave3-TAG.{lead,msolve.log} and appending a JSON receipt line to
# /tmp/hodge/results/_closure/wave3-receipts.jsonl and <dir>/wave3/receipts.jsonl.  Never overwrites wave-2 files.
TAG=$1; TH=$2; TO=$3; shift 3
REC=/tmp/hodge/results/_closure/wave3-receipts.jsonl
VERD=/tmp/hodge/results/_closure/wave3_verdict.py
TIMEV=""; [ -x /usr/bin/time ] && TIMEV="/usr/bin/time -v"
for F in "$@"; do
  D=$(dirname "$F"); B=$(basename "$F")
  LANE=$(echo "$F" | sed -E 's#^/tmp/hodge/results/([^/]+)/.*#\1#')
  OUTD="$D/wave3"; mkdir -p "$OUTD"
  LEAD="$OUTD/$B.wave3-$TAG.lead"; LOG="$OUTD/$B.wave3-$TAG.msolve.log"
  START=$(date -u +%FT%TZ); T0=$(date +%s)
  echo "[wave3 $TAG] $START start msolve -g 1 -t $TH timeout $TO on $F"
  $TIMEV timeout $TO msolve -v 2 -t $TH -g 1 -f "$F" -o "$LEAD" > "$LOG" 2>&1
  RC=$?; T1=$(date +%s)
  MAXRSS=$(grep "Maximum resident" "$LOG" 2>/dev/null | awk '{print $NF}')
  V=$(python3 "$VERD" "$F" "$LEAD" "$RC")
  LINE="{\"ts\":\"$START\",\"tag\":\"$TAG\",\"lane\":\"$LANE\",\"file\":\"$F\",\"lead\":\"$LEAD\",\"log\":\"$LOG\",\"engine\":\"msolve -g 1\",\"threads\":$TH,\"timeout_s\":$TO,\"rc\":$RC,\"wall_s\":$((T1-T0)),\"maxrss_kb\":${MAXRSS:-null},\"host\":\"$(hostname)\",\"result\":$V}"
  echo "$LINE" >> "$REC"; echo "$LINE" >> "$OUTD/receipts.jsonl"
  echo "[wave3 $TAG] done rc=$RC wall=$((T1-T0))s maxrss_kb=${MAXRSS:-?} $V"
done
echo "[wave3 $TAG] ALL DONE"
