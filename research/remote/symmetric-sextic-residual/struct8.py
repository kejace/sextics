"""E1 test of structured invariant degree-8/9 candidates; express kernel reps in e-monomial basis."""
import sys, time, json, pickle, itertools
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from invsolve import *
F, A, q = load()
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
e = {('X', k): e_k(X, k) for k in (1, 2, 3)}; e.update({('Y', k): e_k(Y, k) for k in (1, 2, 3)})
def dlt(vs): return mul(mul(sub(vs[0], vs[1]), sub(vs[0], vs[2])), sub(vs[1], vs[2]))
DX, DY = dlt(X), dlt(Y)
QX, QY = e[('X',2)], e[('Y',2)]
e2sqX = sub(pw(QX, 2), sc(mul(e[('X',1)], e[('X',3)]), 2)); e2sqY = sub(pw(QY, 2), sc(mul(e[('Y',1)], e[('Y',3)]), 2))
def emonomials(deg, pD):
    """all monomials in e1X,e2X,e3X,e1Y,e2Y,e3Y of weighted degree deg with Y-degree parity pD"""
    out = []
    for a in itertools.product(range(deg + 1), repeat=6):
        w = a[0] + 2*a[1] + 3*a[2] + a[3] + 2*a[4] + 3*a[5]
        if w != deg: continue
        ydeg = a[3] + 2*a[4] + 3*a[5]
        if ydeg % 2 != pD: continue
        P = {(0,)*6: Fraction(1)}
        for k, (blk, j) in enumerate([('X',1),('X',2),('X',3),('Y',1),('Y',2),('Y',3)]):
            if a[k]: P = mul(P, pw(e[(blk, j)], a[k]))
        out.append((a, P))
    return out
structured = {
 'DX^2*QX': mul(pw(DX,2), QX), 'DX^2*p2X': mul(pw(DX,2), p_k(X,2)), 'DX^2*QY': mul(pw(DX,2), QY), 'DX^2*p2Y': mul(pw(DX,2), p_k(Y,2)),
 'DY^2*QY': mul(pw(DY,2), QY), 'DY^2*p2Y': mul(pw(DY,2), p_k(Y,2)), 'DY^2*QX': mul(pw(DY,2), QX), 'DY^2*p2X': mul(pw(DY,2), p_k(X,2)),
 'e2(X^2)^2': pw(e2sqX, 2), 'e2(Y^2)^2': pw(e2sqY, 2), 'e2(X^2)*e2(Y^2)': mul(e2sqX, e2sqY), 'e2(X^2)*QX^2': mul(e2sqX, pw(QX,2)),
 'e2(X^2)*p4X': mul(e2sqX, p_k(X,4)), 'e2(X^2)*QY^2': mul(e2sqX, pw(QY,2)), 'e2(X^2)*p4Y': mul(e2sqX, p_k(Y,4)),
 'p8X': p_k(X,8), 'p8Y': p_k(Y,8), 'p4X^2': pw(p_k(X,4),2), 'QX^4': pw(QX,4), 'p8X+p8Y': add(p_k(X,8), p_k(Y,8)),
 'DX^2*e1X^2': mul(pw(DX,2), pw(e[('X',1)],2)), 'DX^2*e1Y^2': mul(pw(DX,2), pw(e[('Y',1)],2)),
 'DX*DY*e1X*e1Y': mul(mul(DX, DY), mul(e[('X',1)], e[('Y',1)])),
}
for lab, D in structured.items():
    ok, info = e1_inv(F, A, q, D) if is_invariant(D) else (None, 'not invariant')
    print(f'{lab:22s} deg={core.degree(D)} E1={ok}', flush=True)
# express degree-8 even kernel reps + J in e-monomial coordinates: solve linear system over F_p
p = P0
info, polys = pickle.load(open(f'{RES}/e1inv_d8_p0.pkl', 'rb'))
basis = emonomials(8, 0)
print('number of e-monomials of degree 8 (even):', len(basis))
Dreps = orbit_reps(8, 0, canon_codes)
import numpy as np
from flint import nmod_mat
def coords(P):
    Pp = core.reduce_mod(P, p)
    return [Pp.get(r, 0) for r in Dreps]
B = [coords(P) for _, P in basis]
M = nmod_mat(len(Dreps), len(B), [B[c][r] for r in range(len(Dreps)) for c in range(len(B))], p)
print('rank of e-monomial basis matrix:', M.rank(), 'of', len(Dreps))
# solve M x = coords(Pk) for each kernel rep -> x = e-monomial coordinates (mod p)
for i, Pk in enumerate(polys):
    rhs = coords(Pk)
    Maug = nmod_mat(len(Dreps), len(B) + 1, [ (B[c][r] if c < len(B) else rhs[r]) for r in range(len(Dreps)) for c in range(len(B) + 1)], p)
    Mr, rk = Maug.rref()
    x = [0]*len(B)
    for r in range(rk):
        j = 0
        while j < len(B) and int(Mr[r, j]) == 0: j += 1
        if j < len(B): x[j] = int(Mr[r, len(B)])
    terms = [(basis[j][0], (v if v < p//2 else v - p)) for j, v in enumerate(x) if v]
    print(f'kernel rep {i} in e-monomial coords ({len(terms)} terms):')
    for a, v in terms: print('    ', v, ' * e1X^%d e2X^%d e3X^%d e1Y^%d e2Y^%d e3Y^%d' % a)
