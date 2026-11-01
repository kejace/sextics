import sys, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, modlu
from hodge import core
F, A, q, rec = core.load_instance('cycle-q2')
D2 = core.mul(core.parse('X0**6 + X1**6 + X2**6'), core.parse('X0**6 + X1**6 + X2**6'))
S = eq.EqSystem(F, A, q, D2, p=modlu.P22)
for seed in [1, 2]:
    t0 = time.time(); r = eq.solve_fast(S, seed=seed); print('seed', seed, r, f'{time.time()-t0:.1f}s', flush=True)
