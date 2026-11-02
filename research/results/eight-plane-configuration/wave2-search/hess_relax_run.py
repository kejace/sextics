"""Hessian relaxation search: minimize || proj_{Ker_{2d-6}^+ perp}(w^T Hess F w) ||^2 over the unit sphere of the w-space
(Euler cone G z removed).  Usage: hess_relax_run.py d parD nrestarts workers [seed0]"""
import sys, time, json, os, numpy as np
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-search')
from w2ep import *
d, parD, nres, workers = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
seed0 = int(sys.argv[5]) if len(sys.argv) > 5 else 0
F, A, q = load()
ker = None
if 2 * d - 6 >= 10:
    z = np.load(f'{SD_DIR}/kernel_float_d{2*d-6}_par0.npz'); ker = z['D']
obj = HessianRelax(F, A, d, parD, ker_basis_vecs=ker, verbose=True)
name = f'hessrelax_d{d}_par{parD}'
t0 = time.time()
results, secs = run_restarts(obj, name, range(seed0, seed0 + nres), maxiter=500, workers=workers, keep_rel=1e-5)
summ = summarize(results)
best = min(results, key=lambda r: r[0]['rel'])
entry = {'batch': 'hessian-relaxation', 'family': 'w-space (D = w(F))', 'd': d, 'parity': parD, 'K': obj.K, 'dimP': obj.dimP,
         'n_quadrics': obj.ST.n - obj.kdim, 'restarts': nres, 'seed0': seed0, 'seconds': round(secs, 1), **summ,
         'note': 'necessary condition A (w^T H w) in J(F); zero only on the Euler cone would exclude all D in J_d'}
log_batch(entry)
print(json.dumps(entry), flush=True)
# save the best few a-vectors
goods = sorted([r for r in results if r[1] is not None], key=lambda r: r[0]['rel'])[:20]
if goods:
    np.save(f'{SD_DIR}/{name}.best.npy', np.array([r[1] for r in goods]))
    print('saved', len(goods), 'candidates with rel <', 1e-5)
