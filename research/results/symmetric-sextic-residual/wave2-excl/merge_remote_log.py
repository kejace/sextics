"""Merge remote wave2-excl.log.jsonl (after `rjob pull`) into the local one (dedupe by line)."""
import sys, os
loc = '/tmp/hodge/results/symmetric-sextic-residual/wave2-excl.log.jsonl'
rem = '/tmp/hodge/remote/symmetric-sextic-residual/wave2-excl.log.jsonl'
have = set(open(loc).read().splitlines()) if os.path.exists(loc) else set()
added = 0
if os.path.exists(rem):
    with open(loc, 'a') as f:
        for line in open(rem).read().splitlines():
            if line and line not in have:
                f.write(line + '\n'); have.add(line); added += 1
print('merged', added, 'new lines')
