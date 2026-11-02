import sys, time; sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import w2lib as W
from hodge import core
import numpy as np
F, A, q = W.load()
t0 = time.time()
fam = W.family(F, A, q, 7, 0, 'sgnX', cache=False, verbose=True)
print('d7 sgnX p0:', {k: v for k, v in fam.items() if k not in ('kernel', 'jac')}, f'{time.time()-t0:.1f}s')
print('  kernel[0] terms', len(fam['kernel'][0]), 'first coeffs', list(fam['kernel'][0].items())[:3])
t0 = time.time()
fam8 = W.family(F, A, q, 8, 0, 'triv', cache=False, verbose=True)
print('d8 triv p0:', {k: v for k, v in fam8.items() if k not in ('kernel', 'jac')}, f'{time.time()-t0:.1f}s')
for P in fam8['kernel']:
    dens = max(Fraction(v).denominator for v in P.values()) if False else max(v.denominator for v in P.values())
    print('  kernel elt terms', len(P), 'max den', dens, 'max num', max(abs(v.numerator) for v in P.values()))
# Fermat sanity: F0 = p6X + p6Y, family J_6^triv (contains p6X), expect rel ~ 0
F0 = {}
for i in range(6):
    e = [0]*6; e[i] = 6; F0[tuple(e)] = core.Fraction(1) if hasattr(core, 'Fraction') else __import__('fractions').Fraction(1)
G = W.group('triv')
pairs = W.jac_family_with_lifts(F0, A, q, 6, 0, G)
print('Fermat J6 triv spanning set size', len(pairs))
basis = [Dk for Dk, _ in pairs]
vp = W.CharVarPro(F0, A, q, basis, G, verbose=True)
res = vp.search(n_restarts=6, seed=1, maxiter=200, verbose=True)
print('Fermat CharVarPro best rel', res[0]['rel'], 'a=', np.round(res[0]['a'], 4))
st, data, fr = W.try_reconstruct_confirm(vp, res[0]['a'], verbose=True)
print('Fermat confirm:', st, fr, {k: v for k, v in data.items() if k in ('rank', 'nullity', 'verified', 'F_divides_D')})
kv = W.KoszulVarPro(F0, A, q, pairs, G, koszul=True, verbose=True)
res = kv.search(n_restarts=6, seed=1, maxiter=200, verbose=True)
print('Fermat KoszulVarPro best rel', res[0]['rel'], 'a=', np.round(res[0]['a'], 4))
kh = W.KoszulVarPro(F0, A, q, pairs, G, koszul=False, verbose=True)
res = kh.search(n_restarts=6, seed=1, maxiter=200, verbose=True)
print('Fermat Hessian-relax best rel', res[0]['rel'], 'a=', np.round(res[0]['a'], 4))
# actual F at d=7 sgnX: family = kernel + jac
basis = fam['kernel'] + fam['jac']
vp = W.CharVarPro(F, A, q, basis, G=W.group('sgnX'), verbose=True)
t0 = time.time(); res = vp.search(n_restarts=4, seed=0, maxiter=200, verbose=True)
print('d7 sgnX CharVarPro best rel', res[0]['rel'], f'{time.time()-t0:.1f}s')
# Koszul version at d=7 sgnX: lifts for kernel elt
Gx = W.group('sgnX')
N0 = W.e1_lift_char(F, A, q, fam['kernel'][0], Gx, verbose=True)
print('lift ok:', N0 is not None)
pairs7 = [(fam['kernel'][0], N0)] + W.jac_family_with_lifts(F, A, q, 7, 0, Gx)
kv = W.KoszulVarPro(F, A, q, pairs7, Gx, koszul=True, verbose=True)
t0 = time.time(); res = kv.search(n_restarts=4, seed=0, maxiter=200, verbose=True)
print('d7 sgnX KoszulVarPro best rel', res[0]['rel'], f'{time.time()-t0:.1f}s')
kh = W.KoszulVarPro(F, A, q, pairs7, Gx, koszul=False, verbose=True)
t0 = time.time(); res = kh.search(n_restarts=4, seed=0, maxiter=200, verbose=True)
print('d7 sgnX Hessian-relax best rel', res[0]['rel'], f'{time.time()-t0:.1f}s')
