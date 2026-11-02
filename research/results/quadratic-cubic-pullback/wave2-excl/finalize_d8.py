"""Merge degree-8 hesstest receipts and update report/draft/packet (idempotent)."""
import json, glob, re
LANE = '/tmp/hodge/results/quadratic-cubic-pullback'
recs = [json.loads(l) for l in open(f'{LANE}/wave2-excl.hesstest.jsonl') if '"d": 8' in l]
master = f'{LANE}/wave2-excl.log.jsonl'
have = set(open(master).read().splitlines()); added = 0
for r in recs:
    s = json.dumps(r)
    if s not in have:
        open(master, 'a').write(s + '\n'); have.add(s); added += 1
excl = [r['name'] for r in recs if r['result'] == 'EXCLUDED-rigorous']
surv = [r['name'] for r in recs if r['result'].startswith('survives')]
mod = [r['name'] for r in recs if r['result'] == 'excluded-modular']
line = (f"Degree 8 (Hessian test, J_16 block at Hilbert rank): {len(excl)} structured candidates excluded rigorously "
        f"({', '.join(excl)}); survivors: {surv or 'none'}; modular-only: {mod or 'none'}.")
p = '/tmp/hodge/notes/wave2-quadratic-cubic-pullback-excl.md'
s = open(p).read()
if 'Degree 8 (Hessian test' not in s:
    s = s.replace("Degree 8 running.", line)
    s = s.replace("degree 8: see hesstest8.", line)
    open(p, 'w').write(s)
q = '/tmp/hodge/notes/wave2-quadratic-cubic-pullback-draftpost.md'
t = open(q).read()
t = t.replace("Degree 8: [[hesstest d=8 results]].", "Degree 8: " + line)
open(q, 'w').write(t)
r_ = f'{LANE}/wave2-excl/README.md'
u = open(r_).read()
if 'Degree 8 (Hessian test' not in u:
    u = u.replace("Individual structured degree-7 polars:", "Degree-8 structured polars: " + line + "\n    Individual structured degree-7 polars:")
    open(r_, 'w').write(u)
import shutil
shutil.copy(master, f'{LANE}/wave2-excl/receipts/')
shutil.copy(f'{LANE}/wave2-excl.hesstest.jsonl', f'{LANE}/wave2-excl/receipts/')
print('merged', added, '|', line)
