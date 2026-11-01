"""Span scan: for many invariant denominators D, compute the space of G-invariant Y-odd sextic
targets served by D (mod p) and accumulate their span; check whether [A] enters the span (mod J)."""
import sys, time, json, pickle
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from invsolve import *
import numpy as np
from flint import nmod_mat
F, A, q = load()
p = P0
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
QX, QY = e_k(X, 2), e_k(Y, 2); e1X, e1Y = e_k(X, 1), e_k(Y, 1); e3X, e3Y = e_k(X, 3), e_k(Y, 3)
p2X, p2Y, p3X, p3Y, p4X, p4Y, p6X, p6Y = p_k(X,2), p_k(Y,2), p_k(X,3), p_k(Y,3), p_k(X,4), p_k(Y,4), p_k(X,6), p_k(Y,6)
def dlt(vs): return mul(mul(sub(vs[0], vs[1]), sub(vs[0], vs[2])), sub(vs[1], vs[2]))
DX, DY = dlt(X), dlt(Y)
def radial(F, idx):
    out = {}
    for i in idx: out = add(out, mul(var(i), core.pderiv(F, i)))
    return out
EXF, EYF = radial(F, (0,1,2)), radial(F, (3,4,5))
e2sqX = sub(pw(QX, 2), sc(mul(e1X, e3X), 2)); e2sqY = sub(pw(QY, 2), sc(mul(e1Y, e3Y), 2))
Treps = orbit_reps(6, 1, canon_codes)   # 19 invariant odd sextic orbit representatives
Tidx = {r: i for i, r in enumerate(Treps)}
nT = len(Treps)
def vec_of_inv_poly(P):
    Pp = core.reduce_mod(P, p)
    return [Pp.get(r, 0) for r in Treps]
# J part: invariant odd sextics in J
Jvecs = []
Fi = [core.pderiv(F, i) for i in range(6)]
for slot, m in [(0, (0,0,0,1,0,0)), (3, (1,0,0,0,0,0))]:
    Jvecs.append(vec_of_inv_poly(sym(core.shift(Fi[slot], m, 1))))
Avec = vec_of_inv_poly(A)
def rank_of(vecs):
    if not vecs: return 0
    M = nmod_mat(len(vecs), nT, [v[i] for v in vecs for i in range(nT)], p)
    return M.rank()
def in_span(vecs, v):
    return rank_of(vecs + [v]) == rank_of(vecs)
print('J^G odd sextic dims:', rank_of(Jvecs), ' A in J?', in_span(Jvecs, Avec), ' sector dim:', nT - rank_of(Jvecs))
cands = []
def C(lab, D): cands.append((lab, D))
C('QX', QX); C('QY', QY); C('p2X', p2X); C('p2Y', p2Y); C('e1X^2', pw(e1X,2)); C('e1Y^2', pw(e1Y,2))
for c in (1, -1, 2, -2, Fraction(1,2), Fraction(-1,2), 3, -3, Fraction(1,3), Fraction(-1,3), 5, -5, Fraction(2,3), Fraction(-2,3)):
    C(f'QX+({c})QY', add(QX, sc(QY, Fraction(c))))
for c in (1, -1, 2, -2):
    C(f'p2X+({c})p2Y', add(p2X, sc(p2Y, Fraction(c)))); C(f'QX+({c})p2Y', add(QX, sc(p2Y, Fraction(c)))); C(f'p2X+({c})QY', add(p2X, sc(QY, Fraction(c))))
for c in (0, 1, -1, 2, -2, 3, -3, 6, -6):
    C(f'p3Y+({c})e3Y', add(p3Y, sc(e3Y, Fraction(c)))); C(f'p3X+({c})e3X', add(p3X, sc(e3X, Fraction(c))))
C('e3Y', e3Y); C('e3X', e3X); C('e1Y*QY', mul(e1Y, QY)); C('e1X*QX', mul(e1X, QX)); C('e1Y^3', pw(e1Y,3)); C('e1X^3', pw(e1X,3))
C('e1X*QY', mul(e1X, QY)); C('e1X*p2Y', mul(e1X, p2Y)); C('e1Y*QX', mul(e1Y, QX)); C('e1Y*p2X', mul(e1Y, p2X))
C('e1X*e1Y^2', mul(e1X, pw(e1Y,2))); C('e1Y*e1X^2', mul(e1Y, pw(e1X,2)))
C('DX', DX); C('DY', DY)
C('p4X', p4X); C('p4Y', p4Y); C('p4X+p4Y', add(p4X, p4Y)); C('p4X-p4Y', sub(p4X, p4Y)); C('QX*QY', mul(QX, QY)); C('p2X*p2Y', mul(p2X, p2Y))
C('e2(X^2)', e2sqX); C('e2(Y^2)', e2sqY); C('e2(X^2)+e2(Y^2)', add(e2sqX, e2sqY)); C('e2(X^2)-e2(Y^2)', sub(e2sqX, e2sqY))
C('QX^2-QY^2', sub(pw(QX,2), pw(QY,2))); C('QX^2+QY^2', add(pw(QX,2), pw(QY,2))); C('e3X*e1X', mul(e3X, e1X)); C('e3Y*e1Y', mul(e3Y, e1Y))
C('e3X*e1Y', mul(e3X, e1Y)); C('e3Y*e1X', mul(e3Y, e1X)); C('p3X*e1Y', mul(p3X, e1Y)); C('p3Y*e1X', mul(p3Y, e1X)); C('DX*e1X', mul(DX, e1X)); C('DY*e1Y', mul(DY, e1Y))
C('DX*e1Y', mul(DX, e1Y)); C('DY*e1X', mul(DY, e1X))
C('p6X', p6X); C('p6Y', p6Y); C('p6X+p6Y', add(p6X, p6Y)); C('EX(F)', EXF); C('EY(F)', EYF); C('EX(F)-EY(F)', sub(EXF, EYF)); C('EX(F)+2EY(F)', add(EXF, sc(EYF,2)))
C('DX^2', pw(DX,2)); C('DY^2', pw(DY,2)); C('DX^2+DY^2', add(pw(DX,2), pw(DY,2))); C('DX*DY', mul(DX, DY)); C('e3X*e3Y', mul(e3X, e3Y)); C('p3X*p3Y', mul(p3X, p3Y))
C('DX*e2(X^2)*e1X', mul(mul(DX, e2sqX), e1X)); C('DX*e2(X^2)*e1Y', mul(mul(DX, e2sqX), e1Y))
if __name__ == '__main__':
    span = list(Jvecs)
    results = []
    for lab, D in cands:
        if not is_invariant(D):
            print(f'{lab:22s} not invariant'); continue
        t0 = time.time()
        rk, basis, S = served_invariant_targets(F, q, D)
        vecs = [[b.get(r, 0) % p for r in Treps] for b in basis]
        before = rank_of(span)
        new = rank_of(span + vecs) - before
        served_dim = rank_of(Jvecs + vecs) - rank_of(Jvecs)   # beyond J
        span = span + vecs
        hasA = in_span(span, Avec)
        results.append((lab, core.degree(D), rk, served_dim, new, hasA))
        log(lab, core.degree(D), 'served-scan', None, None, None, None, time.time() - t0,
            f'served invariant targets: total {rk} (beyond J: {served_dim}); new span dirs {new}; cumulative span {rank_of(span)}; A in span: {hasA}; engine=inv-v1')
        print(f'{lab:22s} d={core.degree(D)} served={rk} beyondJ={served_dim} new={new} cum={rank_of(span)} A_in_span={hasA} ({time.time()-t0:.1f}s)', flush=True)
        if hasA:
            print('!!!!! A enters the span', flush=True)
    pickle.dump(results, open(f'{RES}/served_scan_results.pkl', 'wb'))
