#!/bin/bash
# Pull the remote lane dir mirror and merge all per-process w2log-*.jsonl (remote mirror + local) and the legacy shared
# remote wave2-search.log.jsonl into the canonical local wave2-search.log.jsonl (dedupe by line).
cd /home/kejace/devel/sextics-work && timeout 300 python3 -m hodge.rjob pull --lane cyclic-triple-residual 2>&1 | grep -v "engine v2"
python3 - <<'PY'
import os, glob
dst = '/tmp/hodge/results/cyclic-triple-residual/wave2-search.log.jsonl'
srcs = sorted(glob.glob('/tmp/hodge/remote/cyclic-triple-residual/w2log-*.jsonl') + glob.glob('/tmp/hodge/results/cyclic-triple-residual/w2log-*.jsonl')
              + ['/tmp/hodge/remote/cyclic-triple-residual/wave2-search.log.jsonl'])
have = set(open(dst).read().splitlines()) if os.path.exists(dst) else set()
new = []
for src in srcs:
    if os.path.exists(src):
        for l in open(src).read().splitlines():
            if l and l not in have:
                have.add(l); new.append(l)
with open(dst, 'a') as f:
    for l in new: f.write(l + '\n')
print(f'merged {len(new)} new wave2 log lines from {len(srcs)} files (local total {len(have)})')
PY
