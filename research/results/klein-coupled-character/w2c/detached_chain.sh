#!/bin/bash
# Detached confirmatory char-0 runs (survive the agent session); each completion appends a JSON line to the lane log.
cd /tmp/hodge/results/klein-coupled-character/w2c
LOG=/tmp/hodge/results/klein-coupled-character/wave2-construct.log.jsonl
run() { # kind chart variant outfile
  f=$4; [ -s "$f" ] && return
  timeout 14000 msolve -v 0 -t 1 -P 1 -f "$3" -o "$f"; rc=$?
  head=$(head -c 60 "$f" | tr '\n' ' ')
  printf '{"ts":"%s","step":"detached-char0-run","lens":"construct","engine":2,"input":"%s","output":"%s","exit":%d,"head":"%s"}\n' "$(date +%FT%T)" "$3" "$f" "$rc" "$head" >> "$LOG"
}
for k in 1 2 3 4; do run V2 $k sigmaQ_V2_on_c$k.ms sigmaQP_V2_on_c$k.out; done
for k in 0 1 2 3 4; do run V3 $k sigmaQ_V3_on_c$k.ms sigmaQP_V3_on_c$k.out; done
for k in 1 2 4 5; do run X $k sigmaQ_X_off_c$k.ms sigmaQP_X_off_c${k}_det.out; done
printf '{"ts":"%s","step":"detached-char0-chain-finished","lens":"construct"}\n' "$(date +%FT%T)" >> "$LOG"
