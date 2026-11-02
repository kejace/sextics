"""Fetch the remote per-job logs (wave2-excl/excl-*.log.jsonl) over ssh and append new lines to the local master log."""
import subprocess, sys, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import rjob
MASTER = '/tmp/hodge/results/quadric-contact-residual/wave2-excl.log.jsonl'
r = rjob.ssh('cat /tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-*.log.jsonl 2>/dev/null')
remote = [l for l in r.stdout.splitlines() if l.strip()]
have = set(l.strip() for l in open(MASTER)) if os.path.exists(MASTER) else set()
new = [l for l in remote if l.strip() not in have]
with open(MASTER, 'a') as f:
    for l in new: f.write(l.rstrip('\n') + '\n')
print(f'remote lines {len(remote)}, new merged {len(new)}')
for l in new: print(l[:260])
