"""Exact (mod 32003) bihomogeneous E2 system for a whole sector family, F-reduced Koszul form, for Singular.

  python3 -u w2_exactE2.py d parD k [family J|KJ] [timeout_s] [prime]

Unknowns: a (family coordinates in the chart complementary to the F-multiples: the coordinates a_{(0, x*g)} are set
to 0) and p (Koszul coefficients p_ij, degree d-4, character-restricted).  Equations (one per monomial of degree 2d,
Y-odd, weight 3+2k, not divisible by Y0^6):   sum_{k,l} a_k a_l R(B_{D_k}(N_l))[m] + sum_{k,j} a_k p_j R(L_j(D_k))[m] = 0.
The solution set is a cone containing {a = 0}; we saturate by (a) and ask Singular for the dimension: 0 (only the origin)
means that no sigma-eigen denominator of this (d, parity, weight) has a certificate -- over Q-bar, since the ideal is
homogeneous and mod-p emptiness of a projective scheme implies emptiness of the generic fibre.
"""
import os, sys, time, json, subprocess, pickle
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import w2lib as L
from hodge import core
from hodge.core import NV, shift, mul, add, sub, scale, pderiv
import kcc2
from kcc import W, wt, wmonomials
from kcc2 import Ctx, PI, PAIRS
from fractions import Fraction

F, A, q = L.load()
d, parD, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
family = sys.argv[4] if len(sys.argv) > 4 else 'J'
timeout = int(sys.argv[5]) if len(sys.argv) > 5 else 600
p = int(sys.argv[6]) if len(sys.argv) > 6 else 32003
t0 = time.time()
if family == 'J':
    fam = L.jacobian_family(F, A, d, parD, k)
else:
    from w2_varpro import build_family
    fam = build_family(d, parD, k, family, L.LANE + '/w2_sectors.pkl')
    for _, Nk, _ in fam:
        for Ni in Nk:
            for v in Ni.values():
                assert not isinstance(v, float), 'exact lifts required'
K = len(fam)
ctx = Ctx(F, A, p)
Fp = ctx.F; Fi = ctx.Fi
a_, pA = 3, 1


def rem_F_modp(P):
    t6 = (0, 0, 0, 6, 0, 0)
    rest = {e: v for e, v in Fp.items() if e != t6}
    P = dict(P)
    while True:
        big = [e for e in P if e[3] >= 6 and P[e]]
        if not big:
            break
        e = max(big, key=lambda e: e[3])
        c = P.pop(e)
        m = list(e); m[3] -= 6; m = tuple(m)
        for e2, v2 in rest.items():
            e3 = tuple(x + y for x, y in zip(m, e2))
            P[e3] = (P.get(e3, 0) - c * v2) % p
    return {e: v for e, v in P.items() if v}


Dp = [core.reduce_mod(Dk, p) for Dk, _, _ in fam]
Np = [[core.reduce_mod(Ni, p) for Ni in Nk] for _, Nk, _ in fam]
# p-columns
pcols = []
if d - 4 >= 0:
    for (i, j) in PAIRS:
        par = (pA + parD + PI[i] + PI[j]) % 2
        kp = (a_ + k + W[i] + W[j]) % 7
        for m in wmonomials(d - 4, par, kp):
            pcols.append((i, j, m))
rows = {}
def ridx(e):
    if e not in rows:
        rows[e] = len(rows)
    return rows[e]
# T_k[m][j] coefficients: for each family member k and p-column j the reduced column polynomial
Tcoef = {}   # (k, j) -> dict row -> coeff
for kk in range(K):
    D = Dp[kk]
    Di = [pderiv(D, i, p) for i in range(NV)]
    DF = [mul(D, Fi[i], p) for i in range(NV)]
    G = {(i, j): sub(mul(Fi[j], Di[i], p), mul(Fi[i], Di[j], p), p) for (i, j) in PAIRS}
    for cidx, (i, j, m) in enumerate(pcols):
        col = scale(shift(G[(i, j)], m, 1, p), p - 1, p)
        if m[i]:
            m2 = list(m); m2[i] -= 1
            col = add(col, shift(DF[j], tuple(m2), m[i], p), p)
        if m[j]:
            m2 = list(m); m2[j] -= 1
            col = sub(col, shift(DF[i], tuple(m2), m[j], p), p)
        col = rem_F_modp(col)
        if col:
            Tcoef[(kk, cidx)] = {ridx(e): v for e, v in col.items()}
# B_kl = R(B_{D_k}(N_l))
divN = []
for l in range(K):
    acc = {}
    for i in range(NV):
        if Np[l][i]:
            acc = add(acc, pderiv(Np[l][i], i, p), p)
    divN.append(acc)
Dgrad = [[pderiv(Dp[kk], i, p) for i in range(NV)] for kk in range(K)]
Bcoef = {}
for kk in range(K):
    for l in range(K):
        v = mul(Dp[kk], divN[l], p)
        for i in range(NV):
            if Np[l][i]:
                v = sub(v, mul(Np[l][i], Dgrad[kk][i], p), p)
        v = rem_F_modp(v)
        if v:
            Bcoef[(kk, l)] = {ridx(e): c for e, c in v.items()}
R = len(rows)
# chart: drop a_{(0, x g)}
drop = set()
if d >= 6:
    tagidx = {fam[i][2]: i for i in range(K)}
    for g in wmonomials(d - 6, parD, k):
        xg = tuple(g[j] + (1 if j == 0 else 0) for j in range(6))
        key = ('J', 0, xg)
        if key in tagidx:
            drop.add(tagidx[key])
keep = [i for i in range(K) if i not in drop]
avar = {i: f'a{j}' for j, i in enumerate(keep)}
pvar = [f'p{j}' for j in range(len(pcols))]
eqs = [[] for _ in range(R)]
for (kk, l), vec in Bcoef.items():
    if kk in drop or l in drop:
        continue
    for r, c in vec.items():
        eqs[r].append(f'{c}*{avar[kk]}*{avar[l]}')
for (kk, j), vec in Tcoef.items():
    if kk in drop:
        continue
    for r, c in vec.items():
        eqs[r].append(f'{c}*{avar[kk]}*{pvar[j]}')
polys = ['+'.join(e) for e in eqs if e]
path = f'{L.LANE}/w2-exactE2-d{d}-p{parD}-k{k}-{family}.sing'
with open(path, 'w') as f:
    f.write(f'ring r = {p}, ({",".join(avar[i] for i in keep)},{",".join(pvar)}), dp;\n')
    f.write('option(redSB);\nLIB "elim.lib";\n')
    f.write('ideal I = ' + ',\n'.join(polys) + ';\n')
    f.write('ideal MA = ' + ','.join(avar[i] for i in keep) + ';\n')
    f.write('ideal G = std(I);\n')
    f.write('print("DIM_I " + string(dim(G)) + " GBSIZE " + string(size(G)));\n')
    f.write('list S = sat(G, MA);\n')
    f.write('ideal GS = std(S[1]);\n')
    f.write('print("DIM_SAT " + string(dim(GS)) + " GBSIZE " + string(size(GS)) + " VARS " + string(nvars(r)) + " NEQ " + string(size(I)));\n')
    f.write('if (dim(GS) == 0) { print("VDIM_SAT " + string(vdim(GS))); }\n')
    f.write('quit;\n')
print(f'system: K={K} kept a={len(keep)} p={len(pcols)} rows={R} equations={len(polys)} written {path} [{time.time()-t0:.1f}s]', flush=True)
t1 = time.time()
try:
    r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=timeout)
    out = [l for l in r.stdout.splitlines() if l.startswith(('DIM', 'VDIM'))]
    res = {'lines': out, 'seconds': round(time.time() - t1, 1)}
    for l in out:
        if l.startswith('DIM_SAT'):
            res['dim_sat'] = int(l.split()[1])
        if l.startswith('DIM_I'):
            res['dim_I'] = int(l.split()[1])
except subprocess.TimeoutExpired:
    res = {'lines': ['timeout'], 'seconds': timeout}
print('singular:', res, flush=True)
L.wlog({'family': f'exactE2-{family}', 'd': d, 'parD': parD, 'k': k, 'K': K, 'n_a': len(keep), 'n_p': len(pcols), 'rows': R,
        'n_eq': len(polys), 'prime': p, 'singular': res, 'seconds': round(time.time() - t0, 1),
        'note': 'exact bihomogeneous E2 (F-reduced Koszul form) saturated by (a): dim_sat 0 => no sigma-eigen denominator of this sector has any certificate (over Q-bar)'})
print('EXACTE2_DONE', flush=True)
