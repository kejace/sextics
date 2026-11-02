"""Koszul-reduced E2 for Jacobian denominators D = w(F) (exact reduction, see report):
   find c_ij in S_{d-4} (parity-restricted), Q in S_{2d-6}^{odd} with
   [D div N0 - N0(D)] + sum_{i<j} [D (F_j d_i c_ij - F_i d_j c_ij) - c_ij (F_j D_i - F_i D_j)] = F Q,   N0 = A w.
Modular (any prime) with optional exact rational solve (flint fmpq_mat) for small systems.
Usage: python3 w2x_koszul.py <prime|exact> <tag> <which>   which in {partials, J6cands, randomJ:d:n:seed}
"""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2x_lib import *
from fractions import Fraction
from flint import fmpq_mat, fmpq, nmod_mat

def BD(D, N, p):
    divN = {}; ND = {}
    for i in range(NV):
        if N[i]:
            divN = add(divN, pderiv(N[i], i, p), p)
            ND = add(ND, mul(N[i], pderiv(D, i, p), p), p)
    return sub(mul(D, divN, p), ND, p)

def koszul_system(w, p):
    """w: list of 6 polys (coefficients of the field, exact or mod p). Returns (cols, colpolys, rhs, d, D)."""
    Fp, Fi = partials(p) if p else (F, [pderiv(F, i) for i in range(NV)])
    Ap = reduce_mod(A, p) if p else A
    wp = [reduce_mod(x, p) if p else x for x in w]
    D = {}
    for i in range(NV):
        if wp[i]: D = add(D, mul(wp[i], Fi[i], p), p)
    d = degree(D); pD = parity_of(D); pA = parity_of(A)
    N0 = [mul(Ap, wp[i], p) if wp[i] else {} for i in range(NV)]
    pars = fixed_d.numerator_parities(D, A)
    cols = []; colpolys = []
    degc = d - 4
    for i in range(NV):
        for j in range(i + 1, NV):
            pc = (pars[i] - (0 if j < 3 else 1)) % 2
            for m in monomials(degc, pc):
                K = [{} for _ in range(NV)]
                K[i] = shift(Fi[j], m, 1, p)
                K[j] = shift(Fi[i], m, (p - 1) if p else -1, p)
                cols.append(('c', i, j, m)); colpolys.append(BD(D, K, p))
    Fneg = scale(Fp, (p - 1) if p else -1, p)
    for m in monomials(2 * d - 6, pA):
        cols.append(('Q', m)); colpolys.append(shift(Fneg, m, 1, p))
    rhs = scale(BD(D, N0, p), (p - 1) if p else -1, p)
    return cols, colpolys, rhs, d, D

def solve_mod(w, p):
    t0 = time.time()
    cols, colpolys, rhs, d, D = koszul_system(w, p)
    rows = {}
    for col in colpolys:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, col in enumerate(colpolys):
        for e, v in col.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p); del flat
    Maug_rank = M.rank()
    # rank of M without rhs
    flat2 = [0] * (R * C)
    for j, col in enumerate(colpolys):
        for e, v in col.items(): flat2[rows[e] * C + j] = v
    M0 = nmod_mat(R, C, flat2, p); del flat2
    r0 = M0.rank()
    return {'d': d, 'p': p, 'rows': R, 'cols': C, 'rank': r0, 'rank_aug': Maug_rank, 'consistent': Maug_rank == r0, 'seconds': round(time.time() - t0, 1)}

def solve_exact(w):
    """exact rational consistency via fmpq_mat rank of [M] and [M|b] (small systems only)."""
    t0 = time.time()
    cols, colpolys, rhs, d, D = koszul_system(w, None)
    rows = {}
    for col in colpolys:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    M = fmpq_mat(R, C + 1)
    for j, col in enumerate(colpolys):
        for e, v in col.items(): M[rows[e], j] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
    for e, v in rhs.items(): M[rows[e], C] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
    ra = M.rank()
    M0 = fmpq_mat(R, C)
    for j, col in enumerate(colpolys):
        for e, v in col.items(): M0[rows[e], j] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
    r0 = M0.rank()
    return {'d': d, 'p': 'exact-Q', 'rows': R, 'cols': C, 'rank': r0, 'rank_aug': ra, 'consistent': ra == r0, 'seconds': round(time.time() - t0, 1)}

def unit_field(i, m=None):
    w = [{} for _ in range(NV)]
    w[i] = {m if m else tuple([0] * NV): 1}
    return w

if __name__ == '__main__':
    mode = sys.argv[1]; tag = sys.argv[2]; which = sys.argv[3]
    p = None if mode == 'exact' else int(mode)
    logfn = f'{PK}/wave2-excl-remote-{tag}.jsonl' if tag != 'local' else None
    X = [tuple(1 if k == i else 0 for k in range(NV)) for i in range(NV)]
    cands = []
    if which == 'partials':
        for i in range(NV): cands.append((f'dF/dZ{i}', unit_field(i)))
    elif which == 'J6cands':
        def vf(pairs):
            w = [{} for _ in range(NV)]
            for (zi, i, c) in pairs: w[i] = add(w[i], {X[zi]: c}, None)
            return w
        cands += [('E_0(F)=X0 F_X0 + Y0 F_Y0', vf([(0,0,1),(3,3,1)])), ('E_1(F)', vf([(1,1,1),(4,4,1)])), ('E_2(F)', vf([(2,2,1),(5,5,1)])),
                  ('eps_X(F)', vf([(0,0,1),(1,1,1),(2,2,1)])), ('X0 F_X0', vf([(0,0,1)])), ('Y0 F_Y0', vf([(3,3,1)])),
                  ('X1 F_X0 - X0 F_X1', vf([(1,0,1),(0,1,-1)])), ('R_0(F) = Y0 F_X0 - X0 F_Y0', vf([(3,0,1),(0,3,-1)])),
                  ('X0 F_Y0', vf([(0,3,1)])), ('Y0 F_X0', vf([(3,0,1)])), ('R_0+R_1', vf([(3,0,1),(0,3,-1),(4,1,1),(1,4,-1)])),
                  ('Y1 F_X0 - X0 F_Y1', vf([(4,0,1),(0,4,-1)])), ('Y0 F_X0 + X0 F_Y0', vf([(3,0,1),(0,3,1)]))]
    elif which.startswith('randomJ'):
        _, d, n, seed = which.split(':'); d, n, seed = int(d), int(n), int(seed)
        random.seed(seed)
        for r in range(n):
            for par in (0, 1):
                w = [{} for _ in range(NV)]
                for i in range(NV):
                    pm = (par + (0 if i < 3 else 1)) % 2
                    for m in monomials(d - 5, pm):
                        c = random.randint(-3, 3)
                        if c: w[i][m] = c
                cands.append((f'randomJ_d{d}_par{par}_s{seed}_r{r}', w))
    for name, w in cands:
        try:
            res = solve_exact(w) if p is None else solve_mod(w, p)
        except Exception as ex:
            res = {'error': f'{type(ex).__name__}: {ex}'}
        e = {'exp': 'a3-koszul-reduced-E2', 'D': name, 'family': which, **res, 'prime': res.get('p'), 'result': ('inconsistent' if res.get('consistent') is False else 'CONSISTENT' if res.get('consistent') else 'error'), 'job': tag, 'note': 'Koszul-reduced E2 (exactly equivalent to the full fixed-D system for D in J, d<=9)'}
        print(json.dumps(e), flush=True)
        wlog(e, fn=logfn)
