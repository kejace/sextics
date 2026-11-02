#!/usr/bin/env python3
"""Referee's own generator of the degree-6 Hessian relaxation  w^T Hess(F) w in J_6  (w linear, D = w(F)),
sectors even/odd/mixed, at a fresh prime; runs msolve -g 1 and parses with ref_ltparse. Usage: ref_quad6.py prime"""
import subprocess
from ref_common import *
from ref_ltparse import check
p = int(sys.argv[1])
F, A = load(); Fp = modp(F, p)
Fi = [pd(Fp, i, p) for i in range(NV)]
H = [[pd(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
ei = lambda i: tuple(1 if j == i else 0 for j in range(NV))
for sector in ['even', 'odd', 'mixed']:
    t0 = time.time()
    par_of = {'even': lambda i: [pi(i)], 'odd': lambda i: [1 - pi(i)], 'mixed': lambda i: [0, 1]}[sector]
    basis = [(i, m) for i in range(NV) for par in par_of(i) for m in monos(1, par)]
    n = len(basis); bidx = {b: k for k, b in enumerate(basis)}
    gauge = []
    if all((i, ei(i)) in bidx for i in range(NV)):
        v = [0] * n
        for i in range(NV): v[bidx[(i, ei(i))]] = 1
        gauge.append(v)
    tpars = [0] if sector != 'mixed' else [0, 1]
    mons = [m for tp in tpars for m in monos(6, tp)]; idx = {m: j for j, m in enumerate(mons)}
    gens = [shift(Fi[i], m, p) for i in range(NV) for tp in tpars for m in monos(1, (tp - pi(i)) % 2)]
    G = nmod_mat(len(gens), len(mons), p); fill_rows(G, gens, idx)
    Gr, rk = G.rref(); assert rk == len(gens) == 18 * len(tpars), (rk, len(gens))
    piv = pivots_of(Gr, rk); nonpiv = [j for j in range(len(mons)) if j not in set(piv)]
    pairs = [(k, l) for k in range(n) for l in range(k, n)]
    V = nmod_mat(len(pairs), len(mons), p)
    for r, (k, l) in enumerate(pairs):
        (ik, mk), (il, ml) = basis[k], basis[l]
        for mm, c in shift(H[ik][il], tuple(a + b for a, b in zip(mk, ml)), p).items(): V[r, idx[mm]] = c
    VP = nmod_mat(len(pairs), rk, [int(V[r, j]) for r in range(len(pairs)) for j in piv], p)
    NF = V - VP * Gr
    drop = set()
    if gauge:
        Gm = nmod_mat(len(gauge), n, [x for v in gauge for x in v], p); Gmr, grk = Gm.rref(); drop = set(pivots_of(Gmr, grk))
    keep = [k for k in range(n) if k not in drop]
    eqs = []
    for j in nonpiv:
        terms = []
        for r, (k, l) in enumerate(pairs):
            if k in drop or l in drop: continue
            c = int(NF[r, j])
            if k != l: c = 2 * c % p
            if c: terms.append(f'{c}*c{k}^2' if k == l else f'{c}*c{k}*c{l}')
        if terms: eqs.append('+'.join(terms))
    fn = os.path.join(HERE, f'ref_quad6_{sector}_p{p}.ms'); out = fn.replace('.ms', '.g1.out')
    with open(fn, 'w') as f: f.write(','.join(f'c{k}' for k in keep) + '\n' + f'{p}\n' + ',\n'.join(eqs) + '\n')
    r = subprocess.run(['msolve', '-v', '0', '-t', '4', '-g', '1', '-f', fn, '-o', out], capture_output=True, text=True, timeout=1800)
    rec = check(out)
    log(check='referee-quad6-own-generator', sector=sector, prime=p, unknowns=len(keep), gauge_dim=len(gauge), dim_J6_target=rk,
        equations=len(eqs), msolve_exit=r.returncode, **{k: rec[k] for k in rec if k not in ('file', 'prime', 'nvars')},
        verdict=('EMPTY mod p (m-primary leading ideal) -> no degree-6 denominator in this sector over C' if rec.get('m_primary_mod_p') else 'NOT CERTIFIED'),
        seconds=round(time.time() - t0, 1))
