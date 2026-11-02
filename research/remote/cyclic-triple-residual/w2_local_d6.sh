#!/bin/bash
cd /tmp/hodge/results/cyclic-triple-residual
nice -n 19 python3 -u w2_macblocks.py 6 0 5 1073741827 restrict 2>&1 | grep -v "engine v2" | head -2
nice -n 19 python3 -u - <<'PY' 2>&1 | grep -v "engine v2"
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from w2core import *
import numpy as np
F, A, q = load()
Sp = Sect(F, A, parity_group())
r = e1_space(Sp, F, A, q, 6, (0,), verbose=False)
fam = Family(Sp, F, A, q, 6, (0,), r['basis'], verbose=False)
t0 = time.time()
cands = fam.search(n_restarts=24, seed=11, relax=True, maxiter=400)
print('d=6 even relax search 24 restarts: rels', sorted(round(c['rel'], 14) for c in cands)[:8], f'({time.time()-t0:.1f}s)')
best = cands[0]
print('best D (float, top terms):', sorted(((round(v, 6), core.to_str({m: 1})) for m, v in fam.D_float(best['a']).items()), key=lambda t: -abs(t[0]))[:12])
full = fam.search(n_restarts=1, init=best['a'], maxiter=600)
print('full E2 refined from the best relax point: rel', full[0]['rel'])
w2log({'batch': 'relax-d6-even-local', 'restarts': 24, 'rels_sorted': sorted(c['rel'] for c in cands), 'best_a': best['a'].tolist(),
       'best_D_float': {core.to_str({m: 1}): v for m, v in fam.D_float(best['a']).items()}, 'full_rel_from_best': full[0]['rel'], 'seconds': round(time.time()-t0, 1)})
PY
echo LOCAL_D6_DONE
