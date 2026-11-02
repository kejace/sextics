"""E1 (Jacobian-type) tables, XY frame, engine implementation (independent of the wave-1 (U,V) tensor-ring code).
For target A (instance) and A_K3 = [(X0+Y0)^2(X1+Y1)^2(X2+Y2)^2 - (X0-Y0)^2(X1-Y1)^2(X2-Y2)^2]/64 (K3-sector part),
compute ker( [A] : (R/J)_d -> (R/J)_{d+6} ) for both Y-parities, d = 0..dmax, mod p.
A full-rank (injective) modular result implies injectivity over Q (modular rank <= rational rank).
usage: e1table.py p dmax [out.json]"""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from fractions import Fraction as Fr
p = int(sys.argv[1]); dmax = int(sys.argv[2]); out = sys.argv[3] if len(sys.argv) > 3 else 'e1table.json'
F, A, q, rec = core.load_instance('polarized-k3-remainder')
def var(i):
    e = [0]*6; e[i] = 1; return {tuple(e): Fr(1)}
X = [var(i) for i in range(3)]; Y = [var(3+i) for i in range(3)]
U = [core.add(X[i], Y[i]) for i in range(3)]; V = [core.sub(X[i], Y[i]) for i in range(3)]
def prodsq(L):
    out = {(0,)*6: Fr(1)}
    for l in L: out = core.mul(out, core.mul(l, l))
    return out
AK3 = core.scale(core.sub(prodsq(U), prodsq(V)), Fr(1, 64))
assert core.parity_of(AK3) == 1 and core.degree(AK3) == 6
res = {'p': p, 'targets': {}}
for name, T in [('A', A), ('A_K3', AK3)]:
    res['targets'][name] = {}
    for par in (0, 1):
        for d in range(dmax + 1):
            t0 = time.time()
            k = jacobian.mult_kernel(F, T, d, p=p, parity=par)
            row = {'d': d, 'parity': par, 'ring_dim': k['ring_dim'], 'jacobian_dim': k['jacobian_dim'],
                   'kernel_in_ring': k['dim'], 'kernel_dim_total': k['kernel_dim'], 'seconds': round(time.time()-t0, 1)}
            res['targets'][name][f'{d},{par}'] = row
            print(name, json.dumps(row), flush=True)
            json.dump(res, open(out, 'w'))
print('DONE', flush=True)
