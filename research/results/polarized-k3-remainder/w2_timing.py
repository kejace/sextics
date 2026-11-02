import sys, time
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np
import w2k3 as W
F, A, q = W.load()
d = int(sys.argv[1])
t0 = time.time()
fam = W.Family(F, A, d, kind='jac', verbose=True, sketch_margin=300)
print('build', round(time.time()-t0,1), flush=True)
rng = np.random.default_rng(0)
u = rng.standard_normal(fam.dimP); u /= np.linalg.norm(u)
t0 = time.time(); f, bn, r, c, a = fam.value(u); print('value', round(time.time()-t0,2), 'relres', (f/bn)**0.5, flush=True)
t0 = time.time(); r, J, c, a, SB = fam.resid_jac(u); print('resid_jac', round(time.time()-t0,2), J.shape, flush=True)
t0 = time.time(); out = fam.lm(seed=1, maxiter=8); print('lm 8 iters', round(time.time()-t0,1), out['relres'], out['status'], out['hist'], flush=True)
t0 = time.time(); print('full residual check', fam.full_residual(out['a']), round(time.time()-t0,1), flush=True)
