import sys, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, modlu
from hodge import core
F, A, q, rec = core.load_instance('cycle-q2')
D2 = core.mul(core.parse('X0**6 + X1**6 + X2**6'), core.parse('X0**6 + X1**6 + X2**6'))
for p in [8388593, 8388587, 4194287]:
    S = eq.EqSystem(F, A, q, D2, p=p)
    t0 = time.time(); r = eq.solve_fast(S, seed=1); print('p', p, r, f'{time.time()-t0:.1f}s', flush=True)
