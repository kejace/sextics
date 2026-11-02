"""JOB excl-dp-d9: the degree-9 non-Jacobian E1-kernel (21 lines per parity, sibling lens at p=2^31-1).

1. recompute the kernel at P1 (both parities), save canonical representatives mod P1 (standard monomials);
2. check the (01)-sector split of the kernel (dims per sigma);
3. for each basis element (and the (01)-symmetrised ones) run the relaxation test Psi in J_18 (cheap, rigorous per D)
   and the reduced E2' system at P1 (and P2 for anything passing the relaxation);
4. the 21-dimensional family is NOT exhaustively excluded by this (served(D) is not linear in D); the receipts are a screen
   of the basis directions plus 20 random rational combinations.
Output: d9_<job>.json; log lines exp in {'d9_kernel','d9_relax','d9_e2'}."""
import sys, os, time, json
import numpy as np
import w2lib as W
from hodge import core
from hodge.core import monomials, mul, add, sub, scale, pderiv, reduce_mod, degree

job = os.environ.get('W2JOB', 'd9')
OUT = f'{W.PKT}/d9_{job}.json'
res = json.load(open(OUT)) if os.path.exists(OUT) else {}


def save():
    json.dump(res, open(OUT, 'w'), indent=1, default=str)


p = W.P1
d = 9
bases = {}
for par in (0, 1):
    key = f'kernel_par{par}'
    t0 = time.time()
    k = W.e1_kernel(d, p, par=par, sigma=None, verbose=True)
    res[key] = {'dim': k['dim'], 'ring_dim': k['ring_dim'], 'target_dim': k['target_dim'], 'rank': k['rank'],
                'good_reduction': (k['ring_dim'] == W.koszul_hilbert(d, par) and k['target_dim'] == W.koszul_hilbert(d + 6, 1 - par)),
                'seconds': round(time.time() - t0, 1)}
    bases[par] = k['basis']
    json.dump([core.to_str({e: W.Fraction(v if v <= p // 2 else v - p) for e, v in b.items()}) for b in k['basis']],
              open(f'{W.PKT}/d9_kernel_par{par}_p{p}.json', 'w'))
    W.log('d9_kernel', d=9, prime=p, par=par, result=f"ker={k['dim']}", rank=k['rank'], seconds=res[key]['seconds'],
          note=f"R_9={k['ring_dim']} -> R_15={k['target_dim']}; good reduction {res[key]['good_reduction']}")
    save()
    # (01)-sector split
    for sg in (1, -1):
        ks = W.e1_kernel(d, p, par=par, sigma=sg, verbose=True)
        res[f'kernel_par{par}_sigma{sg}'] = {'dim': ks['dim'], 'ring_dim': ks['ring_dim'], 'target_dim': ks['target_dim']}
        W.log('d9_kernel', d=9, prime=p, par=par, sigma=sg, result=f"ker={ks['dim']}", rank=ks['rank'], note='(01)-sector')
        save()

# relaxation + reduced E2' for each basis element (both parities) and random combinations
reds18 = {par: W.Red(18, p, par, None, verbose=True) for par in (0, 1)}
rng = np.random.default_rng(2026)
for par in (0, 1):
    B = bases[par]
    tests = [(f'par{par}_basis{i}', b) for i, b in enumerate(B)]
    for r in range(10):
        c = rng.integers(1, 100, size=len(B))
        Dc = {}
        for ci, b in zip(c, B):
            Dc = add(Dc, scale(b, int(ci), p), p)
        tests.append((f'par{par}_random{r}', Dc))
    for label, D in tests:
        if label in res:
            continue
        t0 = time.time()
        u = W.solve_E1(D, p)
        if u is None:
            res[label] = {'error': 'E1 failed?!'}; save(); continue
        Ps = W.psi(D, u, p)
        epar = core.parity_of(Ps)
        inJ = reds18[epar].in_J(Ps)
        rec = {'d': 9, 'par': par, 'relax_in_J': bool(inJ), 'seconds_relax': round(time.time() - t0, 1)}
        W.log('d9_relax', D=label, d=9, prime=p, result=('Psi in J (relaxation passes)' if inJ else 'Psi NOT in J: excluded'),
              seconds=rec['seconds_relax'])
        if inJ:
            for pr in (p, W.P2):
                t1 = time.time()
                uu = u if pr == p else W.solve_E1(D, pr)
                S = W.ReducedE2(D, pr, u=uu, verbose=True)
                r2 = S.solve(compressed=True, verbose=True)
                rec[f'e2_p{pr}'] = r2
                W.log('d9_e2', D=label, d=9, prime=pr, result=('consistent' if r2['consistent'] else 'inconsistent'), rank=r2['rank'],
                      nullity=r2['nullity'], seconds=round(time.time() - t1, 1), note=f"reduced E2 {r2['rows']}x{r2['cols']}")
                if r2['consistent']:
                    W.log('d9_e2_ALERT', D=label, d=9, prime=pr, result='CONSISTENT mod p -- lift needed')
        res[label] = rec; save()
        print(label, rec, flush=True)
print('DONE', flush=True)
