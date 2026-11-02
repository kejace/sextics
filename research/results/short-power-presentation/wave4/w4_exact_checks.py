"""Exact (Fraction) checks of the structural identities used in the vertical-families argument."""
import sys, json, random
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave4')
from w4_tlib import *
from w2c_lib import c, L, M, eX, eY, ppow
F, A, q, rec = core.load_instance(IID)
Fi = [pderiv(F, i) for i in range(6)]
X = [{tuple(1 if k == i else 0 for k in range(6)): Fraction(1)} for i in range(6)]
EX = lambda P: core.add(core.add(core.mul(X[0], pderiv(P, 0)), core.mul(X[1], pderiv(P, 1))), core.mul(X[2], pderiv(P, 2)))
u = core.scale(core.add(core.add(core.mul(X[0], Fi[0]), core.mul(X[1], Fi[1])), core.mul(X[2], Fi[2])), Fraction(1, 6))
v = core.scale(core.add(core.add(core.mul(X[3], Fi[3]), core.mul(X[4], Fi[4])), core.mul(X[5], Fi[5])), Fraction(1, 6))
L2, M2 = core.mul(L, L), core.mul(M, M)
hL = core.scale(core.add(core.add(core.mul(L2, core.mul(L2, L2)), core.scale(core.mul(core.mul(L2, L2), M2), 10)), core.scale(core.mul(L2, core.mul(M2, M2)), 5)), 2 * c)
out = {}
out['u == e(X) + 2c(L^6+10L^4M^2+5L^2M^4)'] = core.sub(u, core.add(eX, hL)) == {}
out['u + v == F'] = core.sub(core.add(u, v), F) == {}
EXv = EX(v)
target = core.scale(core.mul(core.mul(L2, M2), core.add(L2, M2)), 40 * c)
out['E_X(v) == 40c L^2M^2(L^2+M^2)'] = core.sub(EXv, target) == {}
# E2 with the particular numerator N0 = (A w/6)(X0,X1,X2,0,0,0) for D = u w:  D div N0 - N0(D) == (A w^2/6) E_X(v)
random.seed(7)
w = core.add(core.add(core.scale(core.mul(L, L), Fraction(3)), core.scale(core.mul(L, M), Fraction(-2))), core.scale(core.mul(M, M), Fraction(5, 7)))
w = core.add(w, core.scale(eX, Fraction(0)))  # degree-2 w (D of degree 8)
D = core.mul(u, w)
N0 = [core.scale(core.mul(core.mul(A, w), X[i]), Fraction(1, 6)) for i in range(3)] + [{}, {}, {}]
divN = {}
for i in range(6): divN = core.add(divN, pderiv(N0[i], i))
ND = {}
for i in range(6): ND = core.add(ND, core.mul(N0[i], pderiv(D, i)))
lhs = core.sub(core.mul(D, divN), ND)
rhs = core.scale(core.mul(core.mul(A, core.mul(w, w)), EXv), Fraction(1, 6))
out['D div N0 - N0(D) == (A w^2/6) E_X(v)  [w = 3L^2-2LM+5/7 M^2]'] = core.sub(lhs, rhs) == {}
# E1 particular solution: sum N0_i F_i == A D
s = {}
for i in range(6): s = core.add(s, core.mul(N0[i], Fi[i]))
out['sum N0_i F_i == A u w'] = core.sub(s, core.mul(A, D)) == {}
# also a degree-7 w including e(X): w = e(X) L + L^7 - 2 L M^6
w2 = core.add(core.add(core.mul(eX, L), ppow(L, 7)), core.scale(core.mul(L, ppow(M, 6)), Fraction(-2)))
D2 = core.mul(u, w2)
N2 = [core.scale(core.mul(core.mul(A, w2), X[i]), Fraction(1, 6)) for i in range(3)] + [{}, {}, {}]
divN2 = {}; ND2 = {}
for i in range(6): divN2 = core.add(divN2, pderiv(N2[i], i)); ND2 = core.add(ND2, core.mul(N2[i], pderiv(D2, i)))
out['identity for w = e(X)L + L^7 - 2LM^6 (deg 13)'] = core.sub(core.sub(core.mul(D2, divN2), ND2), core.scale(core.mul(core.mul(A, core.mul(w2, w2)), EXv), Fraction(1, 6))) == {}
print(json.dumps(out, indent=1))
log({'what': 'exact-identities', 'prime': None, 'result': out, 'file': 'w4_exact_checks.py'})
json.dump({'u': core.to_str(core.scale(u, 7**7)), 'note': '7^7 * u, u = E_X(F)/6 = e(X) + 2c(L^6+10L^4M^2+5L^2M^4)'}, open('u_poly.json', 'w'))
