"""Controls for the degree-5 Hessian relaxation (w2_d5.py), exact over Q.
relax(Atarget) returns dim W (projection of the kernel of [Psi(S)+Phi(T)-FQ] to Sym^2) and the dimension of the
affine cone {v : v(x)v in W} computed with sympy's groebner (independent of Singular).
Controls: (1) Atarget = Y0 * dF/dX0 in J(F): every D serves it, so W must be all of Sym^2 (dim 21) and the cone all of C^6;
          (2) Atarget = A = X0^3 Y0^3: the instance (expect dim W = 8, cone dim 0)."""
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
from w2lib import log
from hodge import core
from hodge.core import monomials, mul, add, sub, shift, pderiv, scale, NV
from fractions import Fraction
from flint import fmpz_mat
import sympy

JOB = sys.argv[1] if len(sys.argv) > 1 else 'excl-d5-control'
F = W.F
Fi = [pderiv(F, i) for i in range(NV)]
Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
Z = [{tuple(1 if j == i else 0 for j in range(NV)): Fraction(1)} for i in range(NV)]
Sidx = [(a, b) for a in range(NV) for b in range(a, NV)]

def relax(At):
    Ai = [pderiv(At, i) for i in range(NV)]
    cols = []
    for a, b in Sidx:
        if a == b: col = sub(mul(Fi[a], Ai[a]), mul(At, Fij[a][a]))
        else: col = sub(add(mul(Fi[a], Ai[b]), mul(Fi[b], Ai[a])), scale(mul(At, Fij[a][b]), 2))
        cols.append(col)
    nS = len(cols)
    for a in range(NV):
        for i in range(NV):
            for j in range(i + 1, NV):
                for c in range(NV):
                    col = {}
                    if c == i: col = add(col, mul(Fi[a], Fi[j]))
                    if c == j: col = sub(col, mul(Fi[a], Fi[i]))
                    col = sub(col, mul(Z[c], sub(mul(Fi[j], Fij[a][i]), mul(Fi[i], Fij[a][j]))))
                    cols.append(col)
    for mq in monomials(4):
        cols.append(scale(shift(F, mq, 1), -1))
    M, rows, scales = W.build_exact(cols)
    X, nul = M.nullspace()
    Wrows = []
    for k in range(nul):
        vec = [Fraction(int(X[j, k]), scales[j]) for j in range(nS)]
        if any(vec): Wrows.append(vec)
    if not Wrows:
        return 0, 0, M.nrows(), M.ncols(), nul, []
    import math
    den = 1
    for vec in Wrows:
        for x in vec: den = den * x.denominator // math.gcd(den, x.denominator)
    Wm = fmpz_mat(len(Wrows), nS, [int(x * den) for vec in Wrows for x in vec])
    R, dd, rW = Wm.rref()
    Wb = fmpz_mat(rW, nS, [int(R[i, j]) for i in range(rW) for j in range(nS)])
    Lam, nlam = Wb.nullspace()
    v = sympy.symbols('v0:6')
    quads = []
    for k in range(nlam):
        qexpr = 0
        for jj, (a, b) in enumerate(Sidx):
            c = int(Lam[jj, k])
            if c: qexpr += c * v[a] * v[b]
        if qexpr != 0: quads.append(sympy.expand(qexpr))
    if not quads:
        return rW, 6, M.nrows(), M.ncols(), nul, []
    G = sympy.groebner(quads, *v, order='grevlex', domain='QQ')
    # Krull dimension of the affine cone: check whether each variable has a pure power as leading term (=> dim 0)
    lead = [sympy.Poly(g, *v).monoms(order='grevlex')[0] for g in G.exprs]
    pure = [any(all(m[j] == 0 for j in range(6) if j != i) and m[i] > 0 for m in lead) for i in range(6)]
    dim0 = all(pure)
    return rW, (0 if dim0 else None), M.nrows(), M.ncols(), nul, [str(g) for g in G.exprs][:20]

for name, At in (('Y0*dF/dX0 (in J(F), served by every D)', mul(core.parse('Y0'), Fi[0])),
                 ('A = X0^3 Y0^3 (instance target)', W.A)):
    t0 = time.time()
    rW, cone, nr, nc, nul, gb = relax(At)
    log(JOB, {'exp': 'd5-relaxation-control', 'target': name, 'rows': nr, 'cols': nc, 'nullity': nul, 'dim_W': rW,
              'cone_dim_sympy': cone, 'groebner_head': gb[:8], 'seconds': round(time.time() - t0, 1),
              'note': 'cone_dim_sympy 0 = only v=0 (pure powers of every variable lead the Groebner basis); 6 = no quadrics (W = Sym^2)'})
