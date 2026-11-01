import sys, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, modlu, numpy as np
from hodge import core
from flint import nmod_mat
F, A, q, rec = core.load_instance('cycle-q2')
D = core.parse('X0**6 + X1**6 + X2**6'); D2 = core.mul(D, D)
S = eq.EqSystem(F, A, q, D2, p=modlu.P22, verbose=True)
Mc = eq._compressed_int(S, k_extra=40, seed=1)
print('compressed shape', Mc.shape, flush=True)
t0 = time.time(); c, r, ra = modlu.consistent_rank(Mc, S.p); t1 = time.time()
print('modlu: consistent', c, 'rank', r, 'rank_aug', ra, f'{t1-t0:.1f}s', flush=True)
k, C1 = Mc.shape
Mn = nmod_mat(k, C1, Mc.ravel().tolist(), S.p)
t0 = time.time(); ra_f = Mn.rank(); t1 = time.time()
M0 = nmod_mat(k, C1 - 1, Mc[:, :-1].ravel().tolist(), S.p)
r_f = M0.rank()
print('flint same matrix: rank', r_f, 'rank_aug', ra_f, f'{t1-t0:.1f}s', flush=True)
# also flint on the big prime path
S2 = eq.EqSystem(F, A, q, D2, p=eq.PRIME)
res = S2.solve()
print('flint big prime:', {k: v for k, v in res.items() if k not in ('pivots', 'x')}, flush=True)
