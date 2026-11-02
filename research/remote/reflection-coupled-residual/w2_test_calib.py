import sys, time, numpy as np
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual')
import w2rcr
from hodge import core
iid = sys.argv[1]
F, A, q, rec = core.load_instance(iid)
d, D, N, Q = core.load_witness(iid)
print(iid, 'd=', d, 'secA', w2rcr.sector_of(A), 'secD', w2rcr.sector_of(D))
# K=1 family
for mode in ('ob', 'full'):
    t = time.time()
    fam = w2rcr.ReducedFamily(F, A, [D], mode=mode, verbose=True)
    f, x, r = fam.value(np.array([1.0]))
    print(mode, 'K=1: f=%.3e bnorm=%.3e rel=%.3e e1rel=%.2e (%.1fs)' % (f, fam.bnorm([1.0]), np.sqrt(f)/fam.bnorm([1.0]), fam.e1_relres, time.time()-t))
# K=3 family: witness + two Jacobian members of same sector
polys, gens = w2rcr.jacobian_family(F, d, w2rcr.sector_of(D))
basis = [D, polys[0], polys[len(polys)//2]]
fam = w2rcr.ReducedFamily(F, A, basis, mode='full', verbose=True)
print('gradcheck', fam.check_grad())
t = time.time()
res = fam.search(n_restarts=6, seed=1, maxiter=200, verbose=False)
for c in res[:3]:
    print('rel=%.2e a=%s spur=%s iters=%d' % (c['rel'], np.round(c['a'], 4), c['spurious'], c['iters']))
print('search time %.1fs nfev=%d' % (time.time()-t, fam._nfev))
