"""Mixed irreducible degree-8 E1 survivors: D = c1(x) X2^2 + c2(x) Y2^2 + c3(x) X2 Y2, c1,c2 in the even d=6 core kernel,
c3 in the odd d=6 core kernel (D Y-even). 6-variable compressed solver mod p."""
from h import *
import json, random, sys
sys.path.insert(0, RES)
from threefold import e1_kernel_core
random.seed(23)
p = core.PRIMES[0]
X22 = core.parse('X2**2'); Y22 = core.parse('Y2**2'); XY = core.parse('X2*Y2')
ev, _ = e1_kernel_core(6, 0, p); od, _ = e1_kernel_core(6, 1, p)
def rnd(basis):
    c = {}
    for b in basis: c = core.add(c, core.scale(b, random.randrange(1, p), p), p)
    return c
tests = []
for r in range(2):
    D = core.add(core.add(core.mul(rnd(ev), X22, p), core.mul(rnd(ev), Y22, p), p), core.mul(rnd(od), XY, p), p)
    tests.append((D, f'c1 X2^2 + c2 Y2^2 + c3 X2Y2, random c_i from the d=6 core kernels (#{r})'))
for D, note in tests:
    t0 = time.time()
    out = fixed_d.solve_fixed_D_fast(F, A, q, D, p=p, verbose=False, want_exact=False)
    ent = {'D': note, 'd': 8, 'family': 'mixed-d8-kernel6xR(B)2', 'e1': True, 'consistent': bool(out['consistent']),
           'rank': out['rank'], 'nullity': out['nullity'], 'rows': out.get('rows'), 'cols': out.get('cols'),
           'seconds': round(time.time() - t0, 1), 'note': 'irreducible mixed denominator; 6-var compressed solver mod p', 'engine': 2, 'p': p}
    log(ent)
    print({k: ent[k] for k in ('D', 'consistent', 'rank', 'nullity', 'rows', 'cols', 'seconds')}, flush=True)
    if out['consistent']: print('!!!!! CONSISTENT', note, flush=True)
