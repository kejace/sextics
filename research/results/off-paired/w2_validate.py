import sys, time; sys.path.insert(0,'/tmp/hodge/results/off-paired')
import numpy as np
import w2lib as W
from hodge import core
F,A,q = W.load()
F0 = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
for hess in (True, False):
    fam, S = W.build_family(F0, A, q, 6, (0,0), 'jac', hessian=hess, verbose=True)
    t=time.time(); res = fam.search(n_restarts=6, seed=1, maxiter=300)
    print('Fermat jac d=6 hess=%s: rels'%hess, ['%.1e'%c['rel'] for c in res], 'time %.1f'%(time.time()-t), 'nfev', fam._nfev, flush=True)
    best = res[0]
    if best['rel'] < 1e-8:
        for fr in fam.reconstruct(best['a'], max_den=50):
            print('  recon', W.D_str(fam.D_exact(fr))[:120], flush=True)
            st, data = fam.confirm(fr)
            print('  confirm ->', st, {k:v for k,v in data.items() if k in ('rank','nullity','rows','cols')}, flush=True)
            break
for kind in ('jac','full'):
    fam, S = W.build_family(F, A, q, 6, (0,0), kind, hessian=False, verbose=True)
    t=time.time(); res = fam.search(n_restarts=8, seed=2, maxiter=300)
    print('off-paired %s d=6 (0,0): rels'%kind, ['%.1e'%c['rel'] for c in res], 'iters', [c['iters'] for c in res], 'time %.1f'%(time.time()-t), 'nfev', fam._nfev, flush=True)
# timing probe at d=8 (0,0) full, staged
t=time.time(); S = W.SectorSpaces(F,A,q,8,(0,0)); print('spaces %.1f'%(time.time()-t), flush=True)
t=time.time(); L = W.E1Lifter(F,A,q,8,(0,0)); print('lifter build %.1f'%(time.time()-t), L.M.shape, flush=True)
t=time.time(); Ns=[L.lift(B) for B in S.ker_basis]; print('lifts %.1f'%(time.time()-t), [round(r,12) for _,r in Ns], flush=True)
t=time.time(); kos = W.koszul_generators(F,8,(0,0)); print('koszul %.1f'%(time.time()-t), len(kos), flush=True)
t=time.time(); fam, S = W.build_family(F, A, q, 8, (0,0), 'full', hessian=False, verbose=True, spaces=S); print('family %.1f'%(time.time()-t), flush=True)
t=time.time(); f,x,r = fam.value(np.random.default_rng(0).standard_normal(fam.K)); print('one eval %.3f s'%(time.time()-t), flush=True)
t=time.time(); res = fam.search(n_restarts=2, seed=3, maxiter=300)
print('off-paired full d=8 (0,0): rels', ['%.1e'%c['rel'] for c in res], 'iters', [c['iters'] for c in res], 'time %.1f'%(time.time()-t), 'nfev', fam._nfev, flush=True)
