"""merge remote mirror of the lane into the local lane dir: wave2-excl.log.jsonl lines (dedupe) + packet files (no overwrite of newer local)."""
import os, subprocess, sys
src = '/tmp/hodge/remote/cycle-q2'; dst = '/tmp/hodge/results/cycle-q2'
rl, ll = f'{src}/wave2-excl.log.jsonl', f'{dst}/wave2-excl.log.jsonl'
have = set(open(ll).read().splitlines()) if os.path.exists(ll) else set()
added = 0
if os.path.exists(rl):
    with open(ll, 'a') as f:
        for line in open(rl).read().splitlines():
            if line and line not in have:
                f.write(line + '\n'); have.add(line); added += 1
r = subprocess.run(['rsync', '-a', '--update', f'{src}/wave2-excl/', f'{dst}/wave2-excl/'], capture_output=True, text=True)
print(f'merged {added} log lines; rsync rc {r.returncode} {r.stderr[-200:]}')
