"""Annihilator of x0^3 in the ternary Jacobian ring R_f = Q[X]/J(f), and E1 tests for D = alpha(X) * h(Y)."""
import sys, json, itertools
import lib
from lib import *
from flint import nmod_mat
p = core.PRIMES[0]
fp = core.reduce_mod(f, p)
fi = [core.pderiv(fp, i, p) for i in range(3)]
def xmons(deg):
    return [m for m in core.monomials(deg) if m[3] == m[4] == m[5] == 0]
def ann_x03_Rf(k, cube=(3,0,0,0,0,0)):
    """{alpha in S(X)_k : x0^3 alpha in J(f)} modulo J(f)_k, basis mod p (canonical)."""
    rows = xmons(k + 3); idx = {m: j for j, m in enumerate(rows)}
    cols = []
    for i in range(3):
        for m in xmons(k + 3 - 5):
            cols.append(core.shift(fi[i], m, 1, p))
    nJ = len(cols)
    amons = xmons(k)
    for m in amons:
        e = tuple(a + b for a, b in zip(m, cube)); cols.append({e: 1})
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[idx[e] * C + j] = v
    X, nul = nmod_mat(R, C, flat, p).nullspace()
    nA = len(amons)
    if nul == 0: return amons, [], 0
    proj = nmod_mat(nA, nul, [int(X[nJ + j, kk]) for j in range(nA) for kk in range(nul)], p)
    E, rk = proj.transpose().rref()
    # J(f)_k
    colsJ = []
    for i in range(3):
        for m in xmons(k - 5):
            colsJ.append(core.shift(fi[i], m, 1, p))
    idxA = {m: j for j, m in enumerate(amons)}
    if colsJ:
        flatJ = [0] * (len(colsJ) * nA)
        for j, c in enumerate(colsJ):
            for e, v in c.items(): flatJ[j * nA + idxA[e]] = v
        Jm = nmod_mat(len(colsJ), nA, flatJ, p); Jr, rkJ = Jm.rref()
        pivJ = set()
        for i in range(rkJ):
            j = 0
            while j < nA and int(Jr[i, j]) == 0: j += 1
            pivJ.add(j)
        rowsAll = [[int(Jr[i, j]) for j in range(nA)] for i in range(rkJ)] + [[int(E[i, j]) for j in range(nA)] for i in range(rk)]
        Sr, rkS = nmod_mat(len(rowsAll), nA, [v for r in rowsAll for v in r], p).rref()
        ker = []
        for i in range(rkS):
            j = 0
            while j < nA and int(Sr[i, j]) == 0: j += 1
            if j not in pivJ: ker.append({amons[jj]: int(Sr[i, jj]) for jj in range(nA) if int(Sr[i, jj])})
    else:
        rkJ = 0
        ker = [{amons[jj]: int(E[i, jj]) for jj in range(nA) if int(E[i, jj])} for i in range(rk)]
    return amons, ker, rkJ
for k in range(3, 9):
    amons, ker, rkJ = ann_x03_Rf(k)
    print(f'Ann_{{R_f}}(x0^3) in degree {k}: dim mod J(f) = {len(ker)} (J(f)_k dim {rkJ})')
    if k == 5:
        ker5 = ker
    if k == 6:
        ker6 = ker
# E1 tests of D = alpha * h(Y) for alpha in Ann_{R_f}(x0^3)_5 and various h
Mp = core.reduce_mod(M, p)
def e1(D):
    ok, info = fixed_d.e1_solvable(F, A, q, {e: int(v) for e, v in D.items()}, p)
    return ok
print('E1 tests (mod p) for D = alpha(X)*h(Y), alpha in Ann_{R_f}(x0^3)_5 basis:')
for ai, al in enumerate(ker5):
    for hn, h in [('M', Mp), ('Y0^2', {(0,0,0,2,0,0):1}), ('Y0*Y1', {(0,0,0,1,1,0):1}), ('Y1^2', {(0,0,0,0,2,0):1}), ('Y2^2', {(0,0,0,0,0,2):1}),
                  ('L', core.reduce_mod(L, p)), ('X0^2', {(2,0,0,0,0,0):1}), ('1*X0*X1', {(1,1,0,0,0,0):1})]:
        D = core.mul(al, h, p)
        print(f'   alpha_{ai} * {hn:6s}: E1 = {e1(D)}')
print('E1 tests for alpha in Ann_{R_f}(x0^3)_6 basis times M, and alpha_5 * M * Y_j:')
for ai, al in enumerate(ker6):
    print(f'   beta_{ai} * M: E1 = {e1(core.mul(al, Mp, p))}')
for ai, al in enumerate(ker5):
    for j, yn in [(3,'Y0'),(4,'Y1'),(5,'Y2')]:
        D = core.mul(core.mul(al, Mp, p), {tuple(1 if i==j else 0 for i in range(6)): 1}, p)
        print(f'   alpha_{ai} * M * {yn}: E1 = {e1(D)}')
json.dump({'ker5': [core.to_str(a) for a in ker5], 'ker6': [core.to_str(a) for a in ker6]}, open(f'{WORK}/annf_x03.json', 'w'))
