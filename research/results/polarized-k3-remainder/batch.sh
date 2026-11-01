#!/bin/bash
# usage: batch.sh candidates.json  (runs 2 jobs in parallel, each with a timeout)
cd /tmp/hodge/results/polarized-k3-remainder
/usr/bin/python3 - "$1" << 'PY' | xargs -P 2 -d '\n' -I{} bash -c '{}'
import json, sys, shlex
C = json.load(open(sys.argv[1]))
for fam, label, expr, frame in C:
    safe = label.replace('/', '_').replace(' ', '_').replace('*', 'x').replace('+', 'p').replace('^', '').replace('(', '').replace(')', '')
    print(f"timeout 2400 /usr/bin/python3 runD.py {shlex.quote(fam)} {shlex.quote(label)} {shlex.quote(expr)} {frame} fast target > run-{safe}.out 2>&1; echo EXIT=$? >> run-{safe}.out")
PY
