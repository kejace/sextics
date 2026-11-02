"""E2-mod-J relaxation for D = u*w (u = E_X(F)/6, w in V_{d-6}): necessary condition A w^2 E_X(v) in J_{2d}, i.e.
A * w^2 * L^2 M^2 (L^2+M^2) in J_{2d}  (two primes, Hilbert check).  d = 7..12 (for d >= 13, R_{2d} = 0: vacuous)."""
import sys, json, random, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave4')
from w4_tlib import *
from w2c_lib import fam_vert, c, L, M, eX, eY, ppow
F, A, q, rec = core.load_instance(IID)
L2, M2 = core.mul(L, L), core.mul(M, M)
G = core.mul(core.mul(L2, M2), core.add(L2, M2))     # E_X(v) / (40c)
random.seed(11)
for d in range(7, 13):
    for par in (0, 1):
        fam = fam_vert(d - 6, par)          # w of degree d-6, parity par  (D = u w has parity par)
        if not fam: continue
        names = [n for n, _ in fam]
        ws = [(n, P) for n, P in fam]
        for t in range(2):   # two random combinations
            cf = {n: random.randint(1, 97) for n in names}
            P = {}
            for n, Q in fam: P = core.add(P, core.scale(Q, cf[n]))
            ws.append(('rand' + str(t) + ':' + ','.join('%s*%d' % (n, cf[n]) for n in names), P))
        res = {}
        t0 = time.time()
        for p in core.PRIMES[:2]:
            polys = [core.mul(core.mul(A, core.mul(w, w)), G) for _, w in ws]
            polys = [reduce_mod(P, p) for P in polys]
            B = TBlock(F, 2 * d, (1 + 0) % 2 if False else core.parity_of(polys[0]), polys, p)
            res[p] = {'in_J2d': [k for k in range(len(ws)) if any(v == {k: 1} for v in B.vecs)], 'kernel_dim': B.info['kernel_dim'],
                      'hilbert_ok': B.info['hilbert_ok'], 'ring_dim': B.info['ring_dim'], 'ring_dim_pred': B.info['ring_dim_pred']}
        line = {'what': 'E2modJ-relaxation-u*w', 'd': d, 'parD': par, 'n_w': len(ws), 'w_names': [n for n, _ in ws][:len(fam)], 'prime': core.PRIMES[:2],
                'kernel_dim': [res[p]['kernel_dim'] for p in core.PRIMES[:2]], 'hilbert_ok': [res[p]['hilbert_ok'] for p in core.PRIMES[:2]],
                'in_J2d': [res[p]['in_J2d'] for p in core.PRIMES[:2]],
                'result': 'A w^2 L^2M^2(L^2+M^2) in J_%d for %d of %d w (basis + 2 random) [membership = e_k in kernel]; kernel_dim %s; hilbert %s' % (
                    2 * d, len(res[core.PRIMES[0]]['in_J2d']), len(ws), [res[p]['kernel_dim'] for p in core.PRIMES[:2]], [res[p]['hilbert_ok'] for p in core.PRIMES[:2]]),
                'seconds': round(time.time() - t0, 1), 'file': 'w4_relax.py'}
        log(line); print(json.dumps({k: line[k] for k in ('d', 'parD', 'n_w', 'kernel_dim', 'hilbert_ok', 'result', 'seconds')}), flush=True)
