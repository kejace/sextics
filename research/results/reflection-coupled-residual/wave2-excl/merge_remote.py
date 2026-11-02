"""Merge remote wave2-excl receipts (mirror) into the local wave2-excl dir and the lane-level wave2-excl log (union of lines)."""
import os, glob
M = '/tmp/hodge/remote/reflection-coupled-residual/wave2-excl'
Lo = '/tmp/hodge/results/reflection-coupled-residual/wave2-excl'
LOG = '/tmp/hodge/results/reflection-coupled-residual/wave2-excl.log.jsonl'
have_log = set(open(LOG).read().splitlines()) if os.path.exists(LOG) else set()
added = {}
for f in glob.glob(f'{M}/*.jsonl'):
    name = os.path.basename(f)
    loc = f'{Lo}/{name}'
    have = set(open(loc).read().splitlines()) if os.path.exists(loc) else set()
    new = [l for l in open(f).read().splitlines() if l and l not in have]
    if new:
        with open(loc, 'a') as g:
            g.write('\n'.join(new) + '\n')
        with open(LOG, 'a') as g:
            for l in new:
                if l not in have_log:
                    g.write(l + '\n'); have_log.add(l)
        added[name] = len(new)
# also copy non-jsonl outputs that do not exist locally
for f in glob.glob(f'{M}/*'):
    name = os.path.basename(f)
    if not name.endswith(('.jsonl', '.out', '.log', '.err', '.pyc')) and not os.path.isdir(f) and not os.path.exists(f'{Lo}/{name}'):
        os.system(f'cp -p "{f}" "{Lo}/{name}"')
print('merged', added)
