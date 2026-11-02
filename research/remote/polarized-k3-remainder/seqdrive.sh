#!/bin/bash
# usage: seqdrive.sh candidates.json outlog [labels...]  -- runs candidates one at a time
cd /tmp/hodge/results/polarized-k3-remainder
CANDS=$1; OUT=$2; shift 2
/usr/bin/python3 - "$CANDS" "$@" << 'PY' > /tmp/hodge/results/polarized-k3-remainder/.cmds.$$
import json, sys, shlex
C = json.load(open(sys.argv[1])); want = sys.argv[2:]
for fam, label, expr, frame in C:
    if want and label not in want: continue
    print(f"timeout 3000 /usr/bin/python3 runD2.py {shlex.quote(fam)} {shlex.quote(label)} {shlex.quote(expr)} {frame}")
PY
while read -r cmd; do echo "### $(date +%T) $cmd" | cut -c1-200 >> "$OUT"; bash -c "$cmd" >> "$OUT" 2>&1; echo "EXIT=$?" >> "$OUT"; done < /tmp/hodge/results/polarized-k3-remainder/.cmds.$$
echo "SEQDRIVE-DONE $(date +%T)" >> "$OUT"
