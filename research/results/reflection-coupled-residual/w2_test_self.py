import sys, time, numpy as np
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual')
import w2rcr
from hodge import core
F, A, q = w2rcr.load()
d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
chi = (0, 0)
polys, gens = w2rcr.jacobian_family(F, d, chi)
print('d', d, 'J family size', len(polys))
t = time.time()
fam = w2rcr.ReducedFamily(F, A, polys, mode='full', verbose=True)
print('build %.1fs' % (time.time() - t), 'nK', fam.nK, 'nQ', fam.nQ, 'rows', fam.SB.n)
rng = np.random.default_rng(0)
a = rng.standard_normal(fam.K); a /= np.linalg.norm(a)
D, N = fam._DN(a)
xK = rng.standard_normal(fam.nK)
t = time.time(); M = fam.kos_matrix(D); print('kos_matrix %.2fs' % (time.time() - t))
lhs = M @ xK
rhs = fam.B(D, fam.kos_to_N(xK))
print('Kos identity relerr %.2e' % (np.linalg.norm(lhs - rhs) / np.linalg.norm(rhs)))
# E1 check for family: sum N_i F_i = A D numerically
print('e1 relres', fam.e1_relres)
for mode, fm in (('full', fam), ('ob', w2rcr.ReducedFamily(F, A, polys, mode='ob', verbose=True))):
    t = time.time(); f, g = fm.value_and_grad(fm.P.T @ a); dt = time.time() - t
    print(mode, 'value %.4e  eval time %.3fs' % (f, dt), 'gradcheck (f, analytic, numeric):', fm.check_grad(1), fm.check_grad(2, 1e-5))
t = time.time()
res = fam.search(n_restarts=3, seed=0, maxiter=100)
print('3 restarts %.1fs' % (time.time() - t), [(round(c['rel'], 6), c['iters'], c['rank']) for c in res])
