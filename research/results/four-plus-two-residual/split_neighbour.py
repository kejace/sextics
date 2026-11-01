"""Control: split neighbour H0 = f(X0,X1) + g(Y0,Y1) (K removed; F0 = H0 + B is in the Chapter 2a family, so a certificate
for (F0, X0^3 Y0^3) exists by theorem). Reduced model with core H0: sweep its core E1-kernel at d=6..8."""
import sys, random
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual')
import threefold as T
from threefold import *
random.seed(13)
p = core.PRIMES[0]
P = core.parse
H0 = add(P('X0**6 + 5*X0**5*X1/7 + 10*X0*X1**5/7 + X1**6'), P('Y0**6 + 5*Y0**5*Y1/7 + 15*Y0*Y1**5/7 + Y1**6'))
# monkeypatch the core H used by e1_kernel_core
T.H = H0
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
nb = int(sys.argv[2]) if len(sys.argv) > 2 else 6
for d in range(6, dmax + 1):
    for par in (0, 1):
        basis, mons = T.e1_kernel_core(d, par, p)
        print(f'[split neighbour H0] d={d} par={par}: core E1-kernel mod J(H0) dim {len(basis)}', flush=True)
        tests = []
        for r in range(3):
            D = {}
            for b in basis: D = add(D, scale(b, random.randrange(1, p), p), p)
            if D: tests.append((D, f'[split neighbour H0=f+g] random combo #{r} of {len(basis)} core E1-kernel elements, d={d} par={par}'))
        for i, b in enumerate(basis[:nb]):
            tests.append((b, f'[split neighbour H0=f+g] core E1-kernel rref basis #{i}, d={d} par={par}, {len(b)} terms'))
        for D, note in tests:
            r = reduced_test(H0, A, CORE, D, p, uidx=2)
            ent = {'D': note, 'd': d, 'family': 'split-neighbour-control', 'e1': True, 'consistent': r['consistent'],
                   'rank': r['rank'], 'nullity': r['nullity'], 'rows': r['rows'], 'cols': r['cols'], 'seconds': r['seconds'],
                   'note': 'CONTROL on H0=f+g (not our F): reduced model U^6+H0, target aU, mod p', 'engine': 'threefold-v1', 'p': p}
            log(ent)
            print({k: ent[k] for k in ('D', 'consistent', 'rank', 'nullity', 'rows', 'cols', 'seconds')}, flush=True)
