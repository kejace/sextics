"""Check Prop. 4.1 numerically: dim Ann_R(L)_7 (6 variables, t=5/7) vs the t-independent prediction
sum_m dim Ann_{R_f}(L)_m * dim (R_g)_{7-m}; also compute the same for F0 = f+g directly."""
import json, time
import lib
from lib import *
from flint import nmod_mat
p = core.PRIMES[0]
def xmons(k): return [m for m in core.monomials(k) if m[3] == m[4] == m[5] == 0]
def ann_dim(Fpoly, T, d, par, vars6=True, xonly=False):
    """dim {D in S_d (parity par, or X-only) : T*D in J(Fpoly)} mod J(Fpoly)_d."""
    Fp = core.reduce_mod(Fpoly, p); Tp = core.reduce_mod(T, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    dT = core.degree(T); parT = core.parity_of(T)
    mons = lambda k, pr: (xmons(k) if xonly else core.monomials(k, pr))
    rows = mons(d + dT, (par + parT) % 2); idx = {m: k for k, m in enumerate(rows)}
    cols = []
    for i in (range(3) if xonly else range(6)):
        pi = 0 if i < 3 else 1
        for m in mons(d + dT - 5, (par + parT - pi) % 2):
            cols.append(core.shift(Fi[i], m, 1, p))
    nJ = len(cols)
    Dm = mons(d, par)
    for m in Dm: cols.append(core.shift(Tp, m, 1, p))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[idx[e] * C + j] = v
    X, nul = nmod_mat(R, C, flat, p).nullspace()
    nD = len(Dm)
    if nul == 0: return 0
    proj = nmod_mat(nD, nul, [int(X[nJ + j, k]) for j in range(nD) for k in range(nul)], p)
    rkE = proj.transpose().rank()
    # J_d
    colsJ = []
    for i in (range(3) if xonly else range(6)):
        pi = 0 if i < 3 else 1
        for m in mons(d - 5, (par - pi) % 2): colsJ.append(core.shift(Fi[i], m, 1, p))
    rJ = 0
    if colsJ:
        idxD = {m: k for k, m in enumerate(Dm)}
        flatJ = [0] * (len(colsJ) * nD)
        for j, c in enumerate(colsJ):
            for e, v in c.items(): flatJ[j * nD + idxD[e]] = v
        rJ = nmod_mat(len(colsJ), nD, flatJ, p).rank()
    return rkE - rJ
Rg = [1,3,6,10,15,18,19,18,15,10,6,3,1]
annf = {m: ann_dim(f, L, m, 0, xonly=True) for m in range(0, 8)}
print('Ann_{R_f}(L)_m for m=0..7:', annf)
pred7 = sum(annf[m] * Rg[7 - m] for m in range(0, 8))
print('prediction for dim Ann_R(L)_7 (t-independent):', pred7)
F0 = core.add(f, g)
for name, FF in [('F (t=5/7)', F), ('F0 (t=0)', F0)]:
    for d in (5, 6, 7):
        print(f'  {name}: dim Ann(L)_{d} (Y-even) =', ann_dim(FF, L, d, 0), flush=True)
print('cross-check for M: dim Ann_R(M)_7 at t=5/7 =', ann_dim(F, M, 7, 0), ' at t=0 =', ann_dim(F0, M, 7, 0))
