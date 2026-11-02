#!/bin/bash
# sequential fixed-D runs (UV frame) with exact Farkas certificates; args: candidates file, optional degree filter
W=/tmp/hodge/results/polarized-k3-remainder/wave2-excl
cd $W
python3 - "$@" <<'PY'
import json, sys, subprocess, time
cands = [json.loads(l) for l in open(sys.argv[1])]
degs = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else None
for c in cands:
    if degs and c['d'] not in degs: continue
    name = c['label'].replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '').replace('^', '').replace('+', 'p').replace('.', '')
    out = f"fx_{name}.out"
    print(time.strftime('%H:%M:%S'), 'running', c['label'], '->', out, flush=True)
    with open(out, 'w') as fh:
        subprocess.run(['python3', 'fixedD_uv.py', c['label'], c['family'], c['expr'], '--components', '--frame', c['frame']], stdout=fh, stderr=subprocess.STDOUT)
    print(time.strftime('%H:%M:%S'), 'done', c['label'], flush=True)
PY
