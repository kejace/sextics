import sys, time
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from w2lib import *
F, A, q = load()
chi, d, pD = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
G = GROUPS[chi]()
t0 = time.time(); L = live_reps(d, pD, G); reps = [r for r, _ in L]; basis = [P for _, P in L]; print('live_reps', len(basis), f'{time.time()-t0:.1f}s', flush=True)
t0 = time.time(); dimK, kv, nN1, rows1 = e1_family_modp(F, A, q, basis, G, pD); print('modp kernel', dimK, nN1, rows1, f'{time.time()-t0:.1f}s', flush=True)
t0 = time.time(); Kb, finfo = e1_family_float(F, A, q, basis, G, pD, expect=dimK); print('float kernel', finfo, f'{time.time()-t0:.1f}s', flush=True)
t0 = time.time(); fam = CharFamily(F, A, q, G, basis, reps=reps, verbose=True); print('family', f'{time.time()-t0:.1f}s', flush=True)
fam.set_subspace(Kb)
t0 = time.time(); u = np.random.default_rng(0).standard_normal(fam.dimP); f, g = fam.value_and_grad(u); print('one eval', f, f'{time.time()-t0:.2f}s', flush=True)
t0 = time.time(); c = fam.search(n_restarts=2, seed=1, maxiter=100, verbose=True); print('2 restarts', f'{time.time()-t0:.1f}s', 'nfev', fam._nfev, flush=True)
