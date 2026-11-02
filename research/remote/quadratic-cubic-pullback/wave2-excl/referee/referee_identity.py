#!/usr/bin/env python3
"""Exact check (own code) of the E2 reduction used by the packet: for D = w.grad F (w constant or polynomial),
N = A w + H z/6 + c.grad F (c antisymmetric), the certificate defect E = A D^2 - D N(F) + F (D div N - N(D)) equals
F (Phi(w) + L(c,w)) with Phi = D (grad A.w + A div w) - A w(D), L = sum_{i,j} F_j [(d_i c_ij) D - c_ij d_i D];
in particular the H-terms cancel.  Also checks the Euler translation invariance of the Hessian condition:
z^T Hess F w = 5 w(F) and z^T Hess F z = 30 F.  Random small integer data, exact rationals."""
import sys, os, json, random, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from hodge import core
import referee_replay as R
F, A, Fi, Fij = R.F, R.A, R.Fi, R.Fij
random.seed(11)
def rpoly(deg, lo=-3, hi=3):
    return {m: Fraction(random.randint(lo, hi)) for m in R.monos(deg) if random.random() < 0.5} or {R.monos(deg)[0]: Fraction(1)}
def ssum(it):
    out = {}
    for t in it: out = core.add(out, t)
    return out
res = {}
for d in (5, 6, 7):
    w = [rpoly(d - 5) for _ in range(6)]
    H = rpoly(d)
    c = [[{} for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            c[i][j] = rpoly(d - 4); c[j][i] = core.scale(c[i][j], -1)
    D = ssum(core.mul(w[i], Fi[i]) for i in range(6))
    z = [{tuple(1 if j == i else 0 for j in range(6)): Fraction(1)} for i in range(6)]
    N = [core.add(core.add(core.mul(A, w[i]), core.scale(core.mul(H, z[i]), Fraction(1, 6))), ssum(core.mul(c[i][j], Fi[j]) for j in range(6))) for i in range(6)]
    E = core.defect(F, A, D, 1, N)
    Phi = core.sub(core.mul(D, core.add(ssum(core.mul(core.pderiv(A, i), w[i]) for i in range(6)), core.mul(A, ssum(core.pderiv(w[i], i) for i in range(6))))),
                   core.mul(A, ssum(core.mul(w[i], core.pderiv(D, i)) for i in range(6))))
    L = ssum(core.mul(Fi[j], core.sub(core.mul(core.pderiv(c[i][j], i), D), core.mul(c[i][j], core.pderiv(D, i)))) for i in range(6) for j in range(6) if c[i][j])
    ok = (core.sub(E, core.mul(F, core.add(Phi, L))) == {})
    # Hessian form: Phi == -A w^T Hess F w  modulo J?  check the exact identity w(D) = sum_k w_k (sum_l (d_k w_l) F_l) + w^T Hess F w
    wD = ssum(core.mul(w[k], core.pderiv(D, k)) for k in range(6))
    wHw = ssum(core.mul(core.mul(w[k], w[l]), Fij[k][l]) for k in range(6) for l in range(6))
    rest = ssum(core.mul(core.mul(w[k], core.pderiv(w[l], k)), Fi[l]) for k in range(6) for l in range(6))
    ok2 = (core.sub(wD, core.add(wHw, rest)) == {})
    res[d] = dict(E_eq_F_PhiL=ok, wD_eq_wHw_plus_Jpart=ok2, deg_D=core.degree(D), N_degrees=[core.degree(n) for n in N])
# Euler identities
z = [{tuple(1 if j == i else 0 for j in range(6)): Fraction(1)} for i in range(6)]
zHz = ssum(core.mul(core.mul(z[k], z[l]), Fij[k][l]) for k in range(6) for l in range(6))
w = [rpoly(1) for _ in range(6)]
zHw = ssum(core.mul(core.mul(z[k], w[l]), Fij[k][l]) for k in range(6) for l in range(6))
wF = ssum(core.mul(w[l], Fi[l]) for l in range(6))
res['euler'] = dict(zHz_eq_30F=(core.sub(zHz, core.scale(F, 30)) == {}), zHw_eq_5wF=(core.sub(zHw, core.scale(wF, 5)) == {}))
R.log(dict(check='E2_reduction_identity_exact', exact=True, results=res, note='own exact check of T3/T5 and of the Euler translation invariance'))
