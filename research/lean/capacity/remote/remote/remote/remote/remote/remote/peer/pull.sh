#!/usr/bin/env bash
# pull remote results/logs into local mirror and merge timings
mkdir -p /tmp/hodge/lean/capacity/remote
rsync -az --exclude out --exclude '*.olean' --exclude src -e "ssh -i $HOME/.ssh/google_compute_engine -o StrictHostKeyChecking=accept-new" kejace@35.229.16.59:/tmp/hodge/lean/capacity/ /tmp/hodge/lean/capacity/remote/
python3 - <<'PY'
import os
src='/tmp/hodge/lean/capacity/remote/results/timings.jsonl'; dst='/tmp/hodge/lean/capacity/results/timings.jsonl'
if os.path.exists(src):
    have=set(open(dst).read().splitlines()) if os.path.exists(dst) else set()
    new=[l for l in open(src).read().splitlines() if l and l not in have]
    open(dst,'a').write(''.join(l+'\n' for l in new)); print('merged',len(new),'new timing records')
PY
