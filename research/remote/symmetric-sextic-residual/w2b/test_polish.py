import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2lib import *
F, A, q = load()
F0 = core.add(p_k(X, 6), p_k(Y, 6))
fam0, info0 = build_sector_family(F0, A, q, 'triv', 6, 0)
c = fam0.search(n_restarts=1, seed=3, maxiter=30)[0]   # deliberately short: unconverged
print('unpolished rel', c['rel'])
a, x, rel, hist = polish(fam0, c['a'], verbose=False)
print('polished rel', rel, 'hist', [f'{h:.1e}' for h in hist])
# real F: polish the triv/d7/p1 minimum
fam, info = build_sector_family(F, A, q, 'triv', 7, 1)
c = fam.search(n_restarts=1, seed=0, maxiter=100)[0]
a, x, rel, hist = polish(fam, c['a'])
print('real triv/7/1: lbfgs rel', c['rel'], '-> polished', rel)
