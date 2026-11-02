"""Degree-6 Jacobian-type denominators that are eigenvectors of the diagonal group G:
    D = F_w = sum_a w_a F_a,   w_a = sum_{c in class(a)} w_ac z_c,
classes of variables with equal G-weights: {X0,X1}, {Y0,Y1}, {X2}, {Y2}  ->  10 parameters w_ac
(this contains the 6-parameter 'weighted Euler' family D = sum_a lambda_a z_a F_a).  chi_D = chi_F, so the
sector restriction is lossless.  E1 solution N = A w + (m/6) z + K (K Koszul, quadratic coefficients in the sector),
E2 (m drops out):   F_w w(A) + A (F_w div w - w(F_w)) + F_w div K - K(F_w) = F Q.
HESSIAN RELAXATION over the 10 parameters: S in Sym^2(C^10) (55), T = w (x) k (10 x n_k), Q in the sector;
exact kernel over Q, projection W to Sym^2, rank-one test by Groebner (sympy, 10 variables).
If W has no nonzero rank-one point, NO denominator of this 10-parameter family admits a certificate (over C)."""
import sys, time, json, os, math, itertools
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
from w2lib import log
from hodge import core, symmetry
from hodge.core import monomials, mul, add, sub, shift, pderiv, scale, NV
from fractions import Fraction
from flint import fmpz_mat
import sympy

JOB = sys.argv[1] if len(sys.argv) > 1 else 'excl-d6jac'
F, A = W.F, W.A
Fi = [pderiv(F, i) for i in range(NV)]
Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
Ai = [pderiv(A, i) for i in range(NV)]
Z = [{tuple(1 if j == i else 0 for j in range(NV)): Fraction(1)} for i in range(NV)]
classes = [(0, 1), (3, 4), (2,), (5,)]
PAIRS = [(a, c) for cl in classes for a in cl for c in cl]     # (a, c): w_a contains w_ac z_c
nP = len(PAIRS)

# sector filter from a sample member (all members have the same characters)
w_sample = {}
Dsample = {}
for (a, c) in PAIRS:
    Dsample = add(Dsample, mul(Z[c], Fi[a]))
allowed, sub_ = W.sector_allowed(Dsample)
assert len(sub_) == 24, len(sub_)
d = 6
pars = W.numerator_parities(Dsample)
pA = core.parity_of(A)

t0 = time.time()
cols, info = [], []
# S-columns: s_{PQ} for P=(a,c), Q=(b,e) (P<=Q): coefficient of w_P w_Q in
#   F_w w(A) + A F_w div w - A w(F_w)
#   F_w = sum_P w_P z_c F_a ;  w(A) = sum_Q w_Q z_e A_b ; div w = sum_Q w_Q [c==a] (i.e. d_a(z_c) = delta_{ac})
#   w(F_w) = sum_{b} w_b d_b(F_w) = sum_{Q=(b,e)} w_Q z_e * sum_{P=(a,c)} w_P ( delta_{bc} F_a + z_c F_ab )
def term(P, Q):
    a, c = P; b, e = Q
    t = mul(mul(Z[c], Fi[a]), mul(Z[e], Ai[b]))                      # F_w w(A)
    if b == e: t = add(t, mul(A, mul(Z[c], Fi[a])))                  # A F_w div w
    # - A w(F_w): Q-part w_b = w_Q z_e acting on P-part of F_w
    u = {}
    if b == c: u = add(u, Fi[a])
    u = add(u, mul(Z[c], Fij[a][b]))
    t = sub(t, mul(A, mul(Z[e], u)))
    return t
for i1 in range(nP):
    for i2 in range(i1, nP):
        P, Q = PAIRS[i1], PAIRS[i2]
        t = term(P, Q) if i1 == i2 else add(term(P, Q), term(Q, P))
        cols.append(t); info.append(('S', P, Q))
nS = len(cols)
# T-columns: w_P * k_{ij,mu}:  F_w div K - K(F_w) with K = sum_{i<j} k_ij (F_j d_i - F_i d_j):
#   div K = sum (d_i k_ij) F_j - (d_j k_ij) F_i ;  K(F_w) = sum k_ij (F_j d_i F_w - F_i d_j F_w)
kinfo = []
for i in range(NV):
    for j in range(i + 1, NV):
        park = (pars[i] - (0 if j < 3 else 1)) % 2
        ms = monomials(2, park)
        mFj = next(iter(Fi[j]))
        ms = [m for m in ms if allowed('N', i, tuple(x + y for x, y in zip(m, mFj)))]
        for m in ms: kinfo.append((i, j, m))
nk = len(kinfo)
for P in PAIRS:
    a, c = P
    FwP = mul(Z[c], Fi[a])                      # the w_P-part of F_w
    dFwP = [pderiv(FwP, i) for i in range(NV)]
    for (i, j, m) in kinfo:
        mp = {m: Fraction(1)}
        dmi, dmj = pderiv(mp, i), pderiv(mp, j)
        t = mul(FwP, sub(mul(dmi, Fi[j]), mul(dmj, Fi[i])))
        t = sub(t, mul(mp, sub(mul(Fi[j], dFwP[i]), mul(Fi[i], dFwP[j]))))
        cols.append(t); info.append(('T', P, (i, j, m)))
nT = len(cols) - nS
msQ = symmetry.filter_monomials(monomials(2 * d - 6, pA), 'Q', -1, allowed)
for m in msQ:
    cols.append(scale(shift(F, m, 1), -1)); info.append(('Q', m))
nQ = len(msQ)
M, rows, scales = W.build_exact(cols)
X, nul = M.nullspace()
rank = M.ncols() - nul
Wrows = []
for k in range(nul):
    vec = [Fraction(int(X[j, k]), scales[j]) for j in range(nS)]
    if any(vec): Wrows.append(vec)
rW = 0; quads = []; cone = None; gb = []
if Wrows:
    den = 1
    for vec in Wrows:
        for x in vec: den = den * x.denominator // math.gcd(den, x.denominator)
    Wm = fmpz_mat(len(Wrows), nS, [int(x * den) for vec in Wrows for x in vec])
    R, dd, rW = Wm.rref()
    Wb = fmpz_mat(rW, nS, [int(R[i, j]) for i in range(rW) for j in range(nS)])
    Lam, nlam = Wb.nullspace()
    v = sympy.symbols(f'w0:{nP}')
    Sidx = [(i1, i2) for i1 in range(nP) for i2 in range(i1, nP)]
    for k in range(nlam):
        q = 0
        for jj, (i1, i2) in enumerate(Sidx):
            cc = int(Lam[jj, k])
            if cc: q += cc * v[i1] * v[i2]
        if q != 0: quads.append(sympy.expand(q))
    if quads:
        G = sympy.groebner(quads, *v, order='grevlex', domain='QQ')
        gb = [str(g) for g in G.exprs]
        lead = [sympy.Poly(g, *v).monoms(order='grevlex')[0] for g in G.exprs]
        pure = [any(all(mm[j] == 0 for j in range(nP) if j != i) and mm[i] > 0 for mm in lead) for i in range(nP)]
        cone = 0 if all(pure) else None
    else:
        cone = nP
log(JOB, {'exp': 'd6-jacobian-eigenvector-relaxation', 'prime': None, 'rows': M.nrows(), 'cols': M.ncols(), 'nS': nS, 'nT': nT, 'nQ': nQ,
          'n_params': nP, 'pairs': PAIRS, 'rank': rank, 'nullity': nul, 'dim_W': rW, 'n_quadrics': len(quads), 'cone_dim': cone,
          'groebner_head': gb[:12], 'seconds': round(time.time() - t0, 1),
          'result': ('no nonzero rank-one tensor in W: no denominator of the 10-parameter eigenvector Jacobian family (deg 6) has a certificate over C [EXACT]'
                     if cone == 0 else f'inconclusive: W has rank-one points (cone dim {cone}); Groebner basis saved'),
          'note': 'family D = sum_a w_a F_a, w_a = sum_{c in class(a)} w_ac z_c, classes {X0,X1},{Y0,Y1},{X2},{Y2}; sector-restricted Koszul/Q unknowns (lossless)'})
W.save_json(JOB, 'relax', {'pairs': PAIRS, 'dim_W': rW, 'quadrics': [str(q) for q in quads], 'groebner': gb})
