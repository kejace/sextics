"""Koszul-reduced E2 solver (mod p).

Given D (degree d, pure parity) and a particular E1 solution N0 (N0(F) = A D, deg N0_i = d+1),
decide whether some N = N0 + K with K = sum_{i<j} c_ij (F_j d_i - F_i d_j) (deg c_ij = d-4)
satisfies  D div N - N(D) = F * R  (deg R = 2d-6).
Unknowns: c_ij monomials (with parities) and R monomials. Rows: monomials of degree 2d (odd parity).
Uses random row compression as in fixed_d.solve_compressed.
"""
import sys, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from lib import *
from hodge.core import monomials, shift, mul, add, sub, scale, pderiv, reduce_mod, yparity, parity_of
from flint import nmod_mat
import numpy as np, scipy.sparse as sp

def BD(D, N, p):
    """D div N - N(D) mod p"""
    divN = {}; ND = {}
    for i in range(6):
        if N[i]:
            divN = add(divN, pderiv(N[i], i, p), p)
            ND = add(ND, mul(N[i], pderiv(D, i, p), p), p)
    return sub(mul(D, divN, p), ND, p)

def e1_particular(D, p):
    """Solve N(F) = A D mod p (least structure); returns N or None."""
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); Dp = reduce_mod(D, p)
    d = core.degree(D); pars = fixed_d.numerator_parities(D, A)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    cols = []; colpolys = []
    for i in range(6):
        for m in monomials(d + 1, pars[i]):
            cols.append((i, m)); colpolys.append(shift(Fi[i], m, 1, p))
    rhs = mul(Ap, Dp, p)
    rows = {}
    for col in colpolys:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, col in enumerate(colpolys):
        for e, v in col.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    Mr, rank = M.rref()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        pivots.append(j)
    if pivots and pivots[-1] == C:
        return None
    N = [{} for _ in range(6)]
    for i, j in enumerate(pivots):
        v = int(Mr[i, C])
        if v:
            ii, m = cols[j]; N[ii][m] = v
    return N

def koszul_e2(D, N0=None, p=core.PRIMES[0], k_extra=40, seed=1, verbose=True):
    t0 = time.time()
    Fp = reduce_mod(F, p); Dp = reduce_mod(D, p)
    d = core.degree(D); pD = parity_of(D); pA = parity_of(A)
    if N0 is None:
        N0 = e1_particular(D, p)
        if N0 is None:
            return {'e1': False}
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    pars = fixed_d.numerator_parities(D, A)
    cols = []; colpolys = []
    degc = d - 4
    # Koszul generators kappa_ij = F_j d_i - F_i d_j ; N_i must have parity pars[i];
    # component i of c*kappa_ij is c F_j: parity(c) + parity(F_j) = pars[i]
    for i in range(6):
        for j in range(i + 1, 6):
            pc = (pars[i] - (0 if j < 3 else 1)) % 2
            assert pc == (pars[j] - (0 if i < 3 else 1)) % 2
            for m in monomials(degc, pc):
                K = [{} for _ in range(6)]
                K[i] = shift(Fi[j], m, 1, p)
                K[j] = shift(Fi[i], m, p - 1, p)
                cols.append(('c', i, j, m)); colpolys.append(BD(Dp, K, p))
    Fneg = scale(Fp, p - 1, p)
    for m in monomials(2 * d - 6, pA):
        cols.append(('R', m)); colpolys.append(shift(Fneg, m, 1, p))
    rhs = scale(BD(Dp, N0, p), p - 1, p)
    rows = {}
    for col in colpolys:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    t1 = time.time()
    if verbose: print(f'[koszul-e2] d={d} rows={R} cols={C} build {t1-t0:.1f}s', file=sys.stderr, flush=True)
    ri, ci, vals = [], [], []
    for j, col in enumerate(colpolys):
        for e, v in col.items(): ri.append(rows[e]); ci.append(j); vals.append(v)
    for e, v in rhs.items(): ri.append(rows[e]); ci.append(C); vals.append(v)
    M = sp.csr_matrix((np.array(vals, dtype=np.int64), (np.array(ri), np.array(ci))), shape=(R, C + 1))
    k = min(R, C + 1 + k_extra)
    Cm = fixed_d._compress(M, p, k, seed)
    Mn = nmod_mat(k, C + 1, Cm.ravel().tolist(), p)
    Mr, rank = Mn.rref()
    t2 = time.time()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == C)
    out = {'e1': True, 'consistent': consistent, 'rank': rank, 'nullity': C - rank, 'rows': R, 'cols': C,
           't_build': round(t1 - t0, 1), 't_solve': round(t2 - t1, 1)}
    if verbose: print(f'[koszul-e2] consistent={consistent} rank={rank} nullity={C-rank} solve {t2-t1:.1f}s', file=sys.stderr, flush=True)
    return out

if __name__ == '__main__':
    # self-test on a d=5 partial (known inconsistent) -- just for timing
    D = core.pderiv(F, 0)
    print(koszul_e2(D))
