#!/usr/bin/env python3
"""Normalized E2 test for E1-kernel representatives (D, N0) with N0(F) = A D (mod p):
   B_D(N0) + T_D(K) = F*C,  K skew (15 entries of degree d-4, parities forced), C in S_{2d-6}^-;  rows S_{2d}^-.
   T_D(K) = sum_{i<j} ( D[(d_i K_ij) F_j - (d_j K_ij) F_i] - K_ij (D_i F_j - D_j F_i) ).
Every E1 numerator is N0 + K grad F (Koszul), so inconsistency of this linear system means no certificate with this D
(modular statement at prime p; rigorous over C iff rank attains its rational value, which is not certified here).
Usage: python3 wave2_excl_e2reps.py reps.json workers [indices...]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
from multiprocessing import Pool
fn = sys.argv[1]; workers = int(sys.argv[2]); which = [int(x) for x in sys.argv[3:]]
data = json.load(open(fn)); p = data['p']; d = data['d']
F, A, q = load()
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
def parse(dct): return {tuple(int(x) for x in k.split()): v for k, v in dct.items()}
def one(ri):
    rep = data['reps'][ri]
    t0 = time.time()
    D = parse(rep['D']); N0 = [parse(x) for x in rep['N']]
    lhs = {}
    for i in range(NV): lhs = add(lhs, mul(N0[i], Fi[i], p), p)
    assert not sub(lhs, mul(Ap, D, p), p), 'E1 fails for rep'
    Di = [pderiv(D, i, p) for i in range(NV)]
    divN = {}
    for i in range(NV): divN = add(divN, pderiv(N0[i], i, p), p)
    ND = {}
    for i in range(NV): ND = add(ND, mul(N0[i], Di[i], p), p)
    rhs = scale(sub(mul(D, divN, p), ND, p), p - 1, p)
    pars = fixed_d.numerator_parities(D, Ap)
    cols = []
    DFij = {}
    for i in range(NV):
        for j in range(i + 1, NV):
            DFij[(i, j)] = sub(mul(Di[i], Fi[j], p), mul(Di[j], Fi[i], p), p)
    nK = 0
    for i in range(NV):
        for j in range(i + 1, NV):
            for m in monomials(d - 4, (pars[i] - par_i(j)) % 2):
                mm = {m: 1}
                t = mul(D, sub(mul(pderiv(mm, i, p), Fi[j], p), mul(pderiv(mm, j, p), Fi[i], p), p), p)
                cols.append(sub(t, mul(mm, DFij[(i, j)], p), p)); nK += 1
    for m in monomials(2 * d - 6, 1):
        cols.append(scale(shift(Fp, m, 1, p), p - 1, p))
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    M = nmod_mat(R, C + 1, p)
    for j, c in enumerate(cols):
        for e, v in c.items(): M[rows[e], j] = v
    for e, v in rhs.items(): M[rows[e], C] = v
    r_aug = M.rank()
    M0 = nmod_mat(R, C, p)
    for j, c in enumerate(cols):
        for e, v in c.items(): M0[rows[e], j] = v
    r0 = M0.rank()
    rec = {'exp': 'e2-normalized-fixed-D', 'reps_file': fn, 'index': ri, 'D': rep['kind'], 'd': d, 'prime': p, 'parD': data.get('parD'),
           'terms_D': len(D), 'rows': R, 'cols': C, 'nK': nK, 'nC': C - nK, 'rank': r0, 'rank_aug': r_aug,
           'consistent': r0 == r_aug, 'seconds': round(time.time() - t0, 1),
           'result': 'CONSISTENT mod p: candidate certificate -> lift and verify!' if r0 == r_aug else 'inconsistent mod p: no certificate with this D (modular)'}
    log2(**rec); print(json.dumps(rec), flush=True)
    return rec
if __name__ == '__main__':
    idx = which if which else list(range(len(data['reps'])))
    with Pool(workers) as pool: res = pool.map(one, idx)
    with open(fn.replace('.json', '_e2.json'), 'w') as f: json.dump(res, f, indent=1)
