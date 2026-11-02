"""Exact (over Q) inconsistency replay of fixed-denominator certificate systems for D = m F_i (m in {1, Z_j}): Koszul form
N = A w + K(p) (w = m e_i), E2 reduced mod 7F (Y0-monic division), p restricted to the sigma/parity sector of D (lossless),
exact ranks with fmpz_mat: inconsistent iff rank [M | b] = rank M + 1.  By the prescribed-support theorem (q = 1) this
excludes every certificate with poles inside V(D)."""
import sys, json, time
from fractions import Fraction
from math import gcd
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl/referee')
import ref_lib as L
from flint import fmpz_mat
F, Fi, H = L.Fpolys(None); A = L.red(L.A, None); pA, kA = 1, 3
Frest = {e: Fraction(-v, 7) for e, v in F.items() if e != (0,0,0,6,0,0)}
memo = {}
def nf_mono(e):
    if e[3] < 6: return {e: Fraction(1)}
    if e in memo: return memo[e]
    base = list(e); base[3] -= 6; base = tuple(base); out = {}
    for e2, v in Frest.items():
        ee = tuple(a + b for a, b in zip(base, e2))
        for e3, v3 in nf_mono(ee).items(): out[e3] = out.get(e3, 0) + v * v3
    memo[e] = {kk: v for kk, v in out.items() if v}; return memo[e]
def nf(P):
    out = {}
    for e, v in P.items():
        for e2, v2 in nf_mono(e).items(): out[e2] = out.get(e2, 0) + v * v2
    return {kk: v for kk, v in out.items() if v}
def test(name, i, m):
    t0 = time.time()
    D = L.shift(Fi[i], m, 1, None); d = L.degree(D)
    par = L.par(next(iter(D))); k = L.wt(next(iter(D)))
    assert all(L.par(e) == par and L.wt(e) == k for e in D)
    Dd = [L.pderiv(D, c, None) for c in range(6)]
    Am = L.shift(A, m, 1, None)                       # N^0 = A w, w = m e_i
    b = nf(L.add(L.mul(D, L.pderiv(Am, i, None), None), L.scale(L.mul(Am, Dd[i], None), -1, None), None))
    cols = []
    for c in range(6):
        for dd in range(c + 1, 6):
            pp = (pA + par + L.PI[c] + L.PI[dd]) % 2; kp = (kA + k + L.W[c] + L.W[dd]) % 7
            for mm in L.monomials(d - 4, pp, kp):
                dm_c, dm_d = L.pderiv({mm: 1}, c, None), L.pderiv({mm: 1}, dd, None)
                col = L.add(L.mul(D, L.add(L.mul(dm_c, Fi[dd], None), L.scale(L.mul(dm_d, Fi[c], None), -1, None), None), None),
                            L.scale(L.shift(L.add(L.mul(Fi[dd], Dd[c], None), L.scale(L.mul(Fi[c], Dd[dd], None), -1, None), None), mm, 1, None), -1, None), None)
                cols.append(nf(col))
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in b: rows.setdefault(e, len(rows))
    den = 1
    for c in cols + [b]:
        for v in c.values(): den = den * v.denominator // gcd(den, v.denominator)
    nr, nc = len(rows), len(cols)
    M = fmpz_mat(nr, nc + 1)
    for c, col in enumerate(cols):
        for e, v in col.items(): M[rows[e], c] = int(v * den)
    for e, v in b.items(): M[rows[e], nc] = int(v * den)
    rAug = M.rank()
    M0 = fmpz_mat(nr, nc)
    for c, col in enumerate(cols):
        for e, v in col.items(): M0[rows[e], c] = int(v * den)
    r0 = M0.rank()
    verdict = 'INCONSISTENT (exact over Q)' if rAug == r0 + 1 else 'consistent'
    print(f'{name}: d={d} par={par} k={k} rows {nr} cols {nc} rank {r0} aug {rAug} -> {verdict} [{time.time()-t0:.1f}s]', flush=True)
    return {'D': name, 'd': d, 'par': par, 'k': k, 'rows': nr, 'cols': nc, 'rank': r0, 'rank_aug': rAug, 'verdict': verdict}
res = []
Z = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
for i in range(6): res.append(test(f'F_{Z[i]}', i, (0,)*6))
for j in range(6):
    for i in range(6): res.append(test(f'{Z[j]}*F_{Z[i]}', i, tuple(1 if t == j else 0 for t in range(6))))
json.dump(res, open('referee_fixedD_exact.json', 'w'), indent=1)
print('FIXEDD_DONE', sum(r['verdict'].startswith('INCONSISTENT') for r in res), 'of', len(res), 'inconsistent exactly over Q')
