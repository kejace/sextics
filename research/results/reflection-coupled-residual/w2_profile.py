import sys, time, numpy as np
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual')
import w2rcr
F, A, q = w2rcr.load()
d = int(sys.argv[1]); mode = sys.argv[2]
polys, gens = w2rcr.jacobian_family(F, d, (0, 0))
t = time.time(); fam = w2rcr.ReducedFamily(F, A, polys, mode=mode, verbose=True); print('build', round(time.time() - t, 1), flush=True)
print('gradcheck', fam.check_grad(1), fam.check_grad(2, 1e-5), flush=True)
rng = np.random.default_rng(0)
for k in range(3):
    u = rng.standard_normal(fam.dimP); u /= np.linalg.norm(u)
    t = time.time(); f, g = fam.value_and_grad(u); print('eval %.3fs f=%.3e' % (time.time() - t, f), flush=True)
t = time.time(); res = fam.search(n_restarts=2, seed=0, maxiter=200)
print('2 restarts %.1fs' % (time.time() - t), [(round(c['rel'], 6), c['iters'], c['nfev'], c['rank']) for c in res], flush=True)
