"""JOB excl-dp-d9: the degree-9 non-Jacobian E1-kernel (21 lines per parity; forced in the (01)-invariant sector).

1. kernel at P1, P2, P3, P4 (both parities), canonical bases (flint nullspace in reduced form over standard monomials);
   CRT + rational reconstruction of each basis vector -> rational representatives D9_k (verified: A*D9_k in J at a 5th prime);
2. for each rational D9_k: relaxation Psi in J_18 (mod P1, P2) and reduced E2' system (mod P1, P2; (01)-equivariant
   restriction when D9_k is a (01)-eigenvector), plus 10 random rational combinations per parity.
Rigor: for a rational D with A*D in J over Q, 'Psi(D;u_p) notin J_p' at a good prime proves that D (and D + gF) fails E2
(regular-sequence argument, README §1.3); E2' inconsistency mod p alone is a screen unless rank-certified.
Output: d9_<job>.json, d9_kernel_par<par>_rational.json; log lines exp in {'d9_kernel','d9_lift','d9_relax','d9_e2'}."""
import sys, os, time, json
import numpy as np
from fractions import Fraction
import w2lib as W
from hodge import core, fixed_d
from hodge.core import monomials, mul, add, sub, scale, pderiv, reduce_mod, degree

job = os.environ.get('W2JOB', 'd9')
OUT = f'{W.PKT}/d9_{job}.json'
res = json.load(open(OUT)) if os.path.exists(OUT) else {}


def save():
    json.dump(res, open(OUT, 'w'), indent=1, default=str)


d = 9
primes = [W.P1, W.P2, W.P3, 2147483629, 2147483587]
rational = {}
for par in (0, 1):
    vecs = []; used = []; std_ref = None; lifted = None
    for p in primes[:4]:
        t0 = time.time()
        k = W.e1_kernel(d, p, par=par, sigma=None, verbose=True)
        key = f'kernel_par{par}_p{p}'
        res[key] = {'dim': k['dim'], 'ring_dim': k['ring_dim'], 'target_dim': k['target_dim'], 'rank': k['rank'],
                    'seconds': round(time.time() - t0, 1)}
        W.log('d9_kernel', d=9, prime=p, par=par, result=f"ker={k['dim']}", rank=k['rank'], seconds=res[key]['seconds'])
        save()
        # basis vectors over the standard monomials of R_9^{par} at this prime
        std = [m for m in W.Red(d, p, par, None).std_mons]
        if std_ref is None:
            std_ref = std
        if std != std_ref:
            print('standard monomials differ at', p, '-- skipping prime for the lift', flush=True)
            continue
        sidx = {m: i for i, m in enumerate(std)}
        B = np.zeros((k['dim'], len(std)), dtype=np.int64)
        for r, b in enumerate(k['basis']):
            for e, v in b.items():
                B[r, sidx[e]] = v
        # canonical form: rref of the basis matrix (unique given the pivot structure)
        Br, piv = W.rref_np(B, p)
        vecs.append((Br, piv)); used.append(p)
        if len(vecs) >= 2 and all(v[1] == vecs[0][1] for v in vecs):
            lifted = []
            ok = True
            for r in range(Br.shape[0]):
                lv = W.lift_vector([v[0][r] for v in vecs], used)
                if lv is None:
                    ok = False; break
                lifted.append(lv)
            if ok:
                print(f'par {par}: lifted {len(lifted)} vectors with {len(used)} primes ({time.time()-t0:.0f}s)', flush=True)
                break
            lifted = None
    if lifted is None:
        res[f'lift_par{par}'] = {'status': 'failed', 'primes': used}
        W.log('d9_lift', d=9, par=par, prime=used, result='rational lift FAILED', note='more primes needed')
        save(); continue
    polys = []
    for lv in lifted:
        P = {std_ref[j]: v for j, v in enumerate(lv) if v != 0}
        polys.append(P)
    # verify E1 at a fresh prime
    pv = primes[4]
    okE1 = all(fixed_d.e1_solvable(W.F, W.A, 1, P, p=pv)[0] for P in polys)
    maxh = max(max(abs(v.numerator), v.denominator) for P in polys for v in P.values())
    res[f'lift_par{par}'] = {'status': 'ok', 'primes': used, 'n': len(polys), 'E1_at_fresh_prime': okE1, 'max_height_bits': maxh.bit_length(),
                             'terms': [len(P) for P in polys]}
    json.dump([core.to_str(P) for P in polys], open(f'{W.PKT}/d9_kernel_par{par}_rational.json', 'w'))
    W.log('d9_lift', d=9, par=par, prime=used, result=f'{len(polys)} rational kernel representatives, E1 at fresh prime {pv}: {okE1}',
          note=f'max height {maxh.bit_length()} bits')
    save()
    rational[par] = polys

# relaxation + reduced E2' on each rational representative and random combinations
reds18 = {}
rng = np.random.default_rng(2026)
for par, polys in rational.items():
    tests = [(f'par{par}_basis{i}', P) for i, P in enumerate(polys)]
    for r in range(10):
        c = rng.integers(-20, 21, size=len(polys))
        Dc = {}
        for ci, P in zip(c, polys):
            if ci:
                Dc = add(Dc, scale(P, Fraction(int(ci))))
        tests.append((f'par{par}_random{r}', Dc))
    for label, D in tests:
        if label in res:
            continue
        rec = {'d': 9, 'par': par, 'terms': len(D)}
        for p in (W.P1, W.P2):
            t0 = time.time()
            Dp = reduce_mod(D, p)
            u = W.solve_E1(D, p)
            if u is None:
                rec[f'p{p}'] = {'error': 'E1 failed mod p'}; continue
            Ps = W.psi(Dp, u, p)
            epar = core.parity_of(Ps)
            if (18, epar, p) not in reds18:
                reds18[(18, epar, p)] = W.Red(18, p, epar, None, verbose=True)
            inJ = reds18[(18, epar, p)].in_J(Ps)
            rec[f'p{p}'] = {'relax_in_J': bool(inJ), 'seconds_relax': round(time.time() - t0, 1)}
            W.log('d9_relax', D=label, d=9, prime=p, result=('Psi in J (relaxation passes)' if inJ else 'Psi NOT in J: D excluded (and D + gF)'),
                  seconds=rec[f'p{p}']['seconds_relax'])
            if inJ:
                t1 = time.time()
                sg = None
                S = W.ReducedE2(D, p, u=u, verbose=True)
                r2 = S.solve(compressed=True, verbose=True)
                rec[f'p{p}']['e2'] = r2
                W.log('d9_e2', D=label, d=9, prime=p, result=('consistent' if r2['consistent'] else 'inconsistent'), rank=r2['rank'],
                      nullity=r2['nullity'], seconds=round(time.time() - t1, 1), note=f"reduced E2 {r2['rows']}x{r2['cols']}")
                if r2['consistent']:
                    W.log('d9_e2_ALERT', D=label, d=9, prime=p, result='CONSISTENT mod p -- exact lift needed')
            save()
        res[label] = rec; save()
        print(label, rec, flush=True)
print('DONE', flush=True)
