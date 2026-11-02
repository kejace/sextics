"""Sector screening: for each character sector and parity, degrees dlist: compute the twisted E1-admissible family
(kernel beyond J + J part), then solve the full twisted E1+E2 system at random points (dense random, sparse random,
and each kernel representative), at the given prime.  Records ranks; any consistent point is lifted exactly."""
import sys, time, json, argparse, random, pickle
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import numpy as np
import w2x_lib as L
import charsolve as CS
from hodge import core
ap = argparse.ArgumentParser()
ap.add_argument('--prime', type=int, default=L.P1); ap.add_argument('--degrees', default='8,9,10,11,12')
ap.add_argument('--chis', default='triv,sgnX,sgnY,sgnsgn'); ap.add_argument('--npoints', type=int, default=4); ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args()
F, A, q = L.load(); p = a.prime
random.seed(a.seed)
for d in [int(x) for x in a.degrees.split(',')]:
    for chi in a.chis.split(','):
        for pD in (0, 1):
            t0 = time.time()
            fam = L.SectorFamily(F, A, d, pD, chi, p, verbose=True)
            rec = {'exp': 'sector-family', 'd': d, 'pD': pD, 'chi': chi, 'prime': p, 'rows': fam.rows, 'nN': fam.nN, 'nD': fam.nD, 'dimK': fam.dimK, 'dimJd': fam.dimJd, 'ker_R': fam.ker_R,
                   'radial': len(fam.radial) if hasattr(fam, 'radial') else 0, 'seconds': round(time.time() - t0, 1), 'job': f'screen_{p}'}
            L.log(rec); print(json.dumps(rec), flush=True)
            if fam.dimK == 0: continue
            pairs, nr = fam.complement_basis()
            if not pairs: continue
            Ds = [D for D, _ in pairs]
            tests = []
            for i, D in enumerate(Ds): tests.append((f'basis{i}', D))
            n = len(Ds)
            for t in range(a.npoints):
                D = {}
                for Dk in Ds: D = core.add(D, core.scale(Dk, random.randrange(1, p), p), p)
                tests.append((f'dense{t}', D))
            for t in range(a.npoints):
                D = {}
                for Dk in random.sample(Ds, min(n, 2)): D = core.add(D, core.scale(Dk, random.randrange(1, p), p), p)
                tests.append((f'sparse2_{t}', D))
            hits = []
            for lab, D in tests:
                t1 = time.time()
                S = CS.CharSystem(F, A, q, core.reduce_mod(D, p), fam.G, p=p)
                r = S.solve()
                rr = {'exp': 'sector-screen', 'd': d, 'pD': pD, 'chi': chi, 'prime': p, 'point': lab, 'result': 'consistent' if r['consistent'] else 'inconsistent',
                      'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows, 'cols': S.ncols, 'seconds': round(time.time() - t1, 1), 'job': f'screen_{p}'}
                if r['consistent']:
                    hits.append((lab, D)); rr['note'] = 'CONSISTENT MOD P -- lift attempted below'
                    st, data = CS.solve_char(F, A, q, {e: (v if v < p//2 else v - p) for e, v in core.reduce_mod(D, p).items()}, fam.G, verbose=True)
                    rr['lift'] = st
                    if st == 'solved':
                        pickle.dump({'D': D, 'N': data['N'], 'Q': data['Q']}, open(f'{L.PK}/WITNESS_{chi}_d{d}_p{pD}_{lab}.pkl', 'wb'))
                L.log(rr); print(json.dumps(rr), flush=True)
            pickle.dump({'chi': chi, 'd': d, 'pD': pD, 'prime': p, 'Ds': Ds, 'radial': nr, 'hits': hits}, open(f'{L.PK}/sector_{chi}_d{d}_p{pD}_{p}.pkl', 'wb'))
