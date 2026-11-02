import os, sys, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-search')
import numpy as np
from fptlib import *
# 1. control: witness residual, gradient check, timing
iid = 'family-hesse-web'
inst0 = Instance(iid); d, Dw, Nw, Qw = core.load_witness(iid)
sub = symmetry.stabilizer(inst0.G, [inst0.A, Dw]); inst = Instance(iid, gens=sub)
lab = label_of_poly(inst, Dw); fam = e1_family(inst, d, lab)
t0 = time.time(); kv = KoszulVarPro(inst, d, lab, fam['basis'], verbose=True); print('build', round(time.time()-t0,1), flush=True)
aw, prj = kv.coords_of(Dw); aw /= np.linalg.norm(aw)
print('witness rel (koszul)', kv.rel(aw), flush=True)
rng = np.random.default_rng(1); u = rng.standard_normal(kv.dimP); u /= np.linalg.norm(u)
f0, g = kv.value_and_grad(u); h = 1e-6; gf = np.zeros_like(u)
for k in range(len(u)):
    e = np.zeros_like(u); e[k] = h; gf[k] = (kv.value_and_grad(u+e)[0] - kv.value_and_grad(u-e)[0])/(2*h)
print('grad rel err', np.max(np.abs(g-gf))/np.max(np.abs(g)), flush=True)
r, J = kv.residual_and_jac(u); print('grad vs 2J^T r', np.linalg.norm(2*J.T@r - g)/np.linalg.norm(g), flush=True)
t0 = time.time(); res = kv.search(n_restarts=16, seed=9, maxiter=400); print('lbfgs 16 restarts', round(time.time()-t0,1), 's; rels', [f"{r['rel']:.1e}" for r in res[:6]], flush=True)
t0 = time.time(); res2 = kv.search_lm(n_restarts=16, seed=9, max_nfev=200); print('lm 16 restarts', round(time.time()-t0,1), 's; rels', [f"{r['rel']:.1e}" for r in res2[:6]], flush=True)
# 2. our instance: d=6 sector (0,0,0) build + timing
inst2 = Instance(); fam2 = e1_family(inst2, 6, (0,0,0))
t0 = time.time(); kv2 = KoszulVarPro(inst2, 6, (0,0,0), fam2['basis'], verbose=True); print('build fpt d6', round(time.time()-t0,1), flush=True)
t0 = time.time(); res3 = kv2.search(n_restarts=3, seed=1, maxiter=300); print('fpt d6 3 restarts', round(time.time()-t0,1), 's; rels', [f"{r['rel']:.1e}" for r in res3], 'iters', [r['iters'] for r in res3], flush=True)
# compare with the full formulation at the same point
sv2 = SectorVarPro(inst2, 6, (0,0,0), fam2['basis'])
for r_ in res3[:2]:
    print('  full-system rel at koszul minimizer:', f"{sv2.rel(r_['a']):.3e}", 'koszul rel', f"{r_['rel']:.3e}", flush=True)
