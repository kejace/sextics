"""Reducible degree-7 E1 survivors: c(x)*X2 (c in even d=6 core kernel) and c(x)*Y2 (c in odd d=6 core kernel). 6-variable engine, mod p."""
from h import *
import json, random, sys
sys.path.insert(0, RES)
from threefold import e1_kernel_core
random.seed(17)
p = core.PRIMES[0]
X2 = core.parse('X2'); Y2 = core.parse('Y2')
ev, _ = e1_kernel_core(6, 0, p); od, _ = e1_kernel_core(6, 1, p)
tests = []
for r in range(2):
    c = {}
    for b in ev: c = core.add(c, core.scale(b, random.randrange(1, p), p), p)
    tests.append((core.mul(c, X2, p), f'c(x)*X2, c = random combo #{r} of the 20 even d=6 core-kernel elements'))
    c = {}
    for b in od: c = core.add(c, core.scale(b, random.randrange(1, p), p), p)
    tests.append((core.mul(c, Y2, p), f'c(x)*Y2, c = random combo #{r} of the 13 odd d=6 core-kernel elements'))
for D, note in tests:
    t0 = time.time()
    out = fixed_d.solve_fixed_D_fast(F, A, q, D, p=p, verbose=False, want_exact=False)
    ent = {'D': note, 'd': 7, 'family': 'reducible-d7-products', 'e1': True, 'consistent': bool(out['consistent']),
           'rank': out['rank'], 'nullity': out['nullity'], 'rows': out.get('rows'), 'cols': out.get('cols'),
           'seconds': round(time.time() - t0, 1), 'note': 'support V(c) u {X2=0} (or {Y2=0}); 6-var compressed solver mod p', 'engine': 2, 'p': p}
    log(ent)
    print({k: ent[k] for k in ('D', 'consistent', 'rank', 'nullity', 'rows', 'cols', 'seconds')}, flush=True)
    if out['consistent']: print('!!!!! CONSISTENT', note, flush=True)
