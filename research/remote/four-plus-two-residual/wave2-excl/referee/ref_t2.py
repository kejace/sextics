"""Independent degree-5 Hessian relaxation.  D = F_v.  E2 (m eliminated):  F_v A_v - A F_vv + F_v div K - K(F_v) = F Q.
Relax v(x)v -> S (21), v(x)l -> T (6*15*6 = 540), Q (126 quartics).  Exact nullspace (one global denominator)."""
import sys, time, json, subprocess, math; sys.path.insert(0, '/tmp/hodge/work/referee42')
import reflib as R
from fractions import Fraction
from flint import fmpz_mat
LOG = '/tmp/hodge/work/referee42/out/ref_t2.jsonl'
F, A = R.load_instance()
NV = 6
Fi = [R.pderiv(F, i) for i in range(NV)]; Fij = [[R.pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
def relax(target, tag):
    Ai = [R.pderiv(target, i) for i in range(NV)]
    t0 = time.time(); cols = []; Sidx = []
    for a in range(NV):
        for b in range(a, NV):
            if a == b: col = R.padd(R.pmul(Fi[a], Ai[a]), R.pmul(target, Fij[a][a]), -1)
            else: col = R.padd(R.padd(R.pmul(Fi[a], Ai[b]), R.pmul(Fi[b], Ai[a])), R.pscale(R.pmul(target, Fij[a][b]), Fraction(2)), -1)
            cols.append(col); Sidx.append((a, b))
    nS = len(cols)
    for a in range(NV):
        for i in range(NV):
            for j in range(i + 1, NV):
                for c in range(NV):
                    col = {}
                    if c == i: col = R.padd(col, R.pmul(Fi[a], Fi[j]))
                    if c == j: col = R.padd(col, R.pmul(Fi[a], Fi[i]), -1)
                    col = R.padd(col, R.pmul(R.z(c), R.padd(R.pmul(Fi[j], Fij[a][i]), R.pmul(Fi[i], Fij[a][j]), -1)), -1)
                    cols.append(col)
    nT = len(cols) - nS
    for m in R.monomials(4):
        cols.append(R.pscale(R.pshift(F, m), Fraction(-1)))
    nQ = len(cols) - nS - nT
    M, rows, L = R.build_exact(cols)
    X, nul = M.nullspace()
    rank = M.ncols() - nul
    # project kernel onto S-coordinates (no column scaling was applied, so coordinates are exact)
    vecs = [[int(X[j, k]) for j in range(nS)] for k in range(nul)]
    vecs = [v for v in vecs if any(v)]
    if vecs:
        Wm = fmpz_mat(len(vecs), nS, [x for v in vecs for x in v]); Rr, den, rW = Wm.rref()
        Wb = fmpz_mat(rW, nS, [int(Rr[i, j]) for i in range(rW) for j in range(nS)])
        Lam, nlam = Wb.nullspace()
        quads = []
        for k in range(nlam):
            terms = [f'({int(Lam[jj, k])})*v{a}*v{b}' for jj, (a, b) in enumerate(Sidx) if int(Lam[jj, k])]
            quads.append('+'.join(terms))
        Wsupp = sorted({Sidx[j] for i in range(rW) for j in range(nS) if int(Rr[i, j])})
    else:
        rW, nlam, quads, Wsupp = 0, nS, [f'v{a}*v{b}' for (a, b) in Sidx], []
    sing = 'ring r=0,(v0,v1,v2,v3,v4,v5),dp; ideal I=' + ','.join(quads) + '; ideal G=std(I); print(dim(G)); print(G); quit;'
    path = f'/tmp/hodge/work/referee42/out/ref_t2_{tag}.sing'; open(path, 'w').write(sing)
    out = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=600).stdout
    dim_cone = int(out.strip().splitlines()[0]) if quads and quads[0] else NV
    R.log(LOG, {'exp': 'd5-relaxation-independent', 'target': tag, 'rows': M.nrows(), 'cols': M.ncols(), 'global_denominator': L, 'nS': nS, 'nT': nT, 'nQ': nQ,
                'rank_Q': rank, 'nullity_Q': nul, 'dim_W': rW, 'W_support_pairs': Wsupp, 'n_quadrics': nlam, 'affine_cone_dim_Singular': dim_cone,
                'groebner': out.strip().replace('\n', ' ')[:300], 'seconds': round(time.time() - t0, 1)})
relax(A, 'A')
relax(R.pmul({(0, 0, 0, 1, 0, 0): Fraction(1)}, Fi[0]), 'control_Y0_FX0')
