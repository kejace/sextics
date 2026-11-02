import sys, time
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from w2lib import *
F, A, q = load()
d, pD = int(sys.argv[1]), int(sys.argv[2])
withk = len(sys.argv) > 3
t0 = time.time(); fam, info = build_jac_family(F, A, q, d, pD, with_kernel=withk, verbose=True); print('info', info, f'{time.time()-t0:.1f}s', flush=True)
t0 = time.time(); u = np.random.default_rng(0).standard_normal(fam.dimP); f, g = fam.value_and_grad(u); print('one eval', f, f'{time.time()-t0:.2f}s', flush=True)
# gradient check
eps = 1e-6; v = np.random.default_rng(1).standard_normal(fam.dimP)
f1, _ = fam.value_and_grad(u + eps * v); f0, _ = fam.value_and_grad(u - eps * v)
print('grad check: fd', (f1 - f0) / (2 * eps), 'analytic', float(g @ v))
t0 = time.time(); c = fam.search(n_restarts=3, seed=1, maxiter=200, verbose=True); print('3 restarts', f'{time.time()-t0:.1f}s', 'best', c[0]['rel'])
# Fermat sanity at d=6 (p6X = E_X(F0)/6 is in the w-space)
if d == 6:
    F0 = core.add(p_k(X, 6), p_k(Y, 6))
    fam0, info0 = build_jac_family(F0, A, q, 6, 0, with_kernel=False, verbose=True)
    c0 = fam0.search(n_restarts=3, seed=2, maxiter=300, verbose=True)
    print('Fermat jac best rel', c0[0]['rel'])
    ok, res = try_confirm(fam0, c0[0], {'family': 'test-fermat-jac'})
    print('fermat confirm', ok, res)
