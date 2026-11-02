"""Structured candidate denominators (UV frame: slots 0..2 = U, 3..5 = V; converted to XY by the sweep)."""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, json, itertools
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import w2lib as L
from hodge import core
from hodge.core import mul, add, sub, scale, power, pderiv
from fractions import Fraction as Fr
f = L.f_ternary()
S = core.to_str
def sw(P): return L.swapUV(P)
def sym(P): return add(P, sw(P))
def anti(P): return sub(P, sw(P))
U = [L.var(i) for i in range(3)]; V = [L.var(3 + i) for i in range(3)]
fU = f; fV = sw(f)
fi = [pderiv(f, i) for i in range(3)]                      # degree 5, in U
fij = [[pderiv(fi[i], j) for j in range(3)] for i in range(3)]   # degree 4
fijk = [[[pderiv(fij[i][j], k) for k in range(3)] for j in range(3)] for i in range(3)]  # degree 3
fiV = [sw(x) for x in fi]; fijV = [[sw(x) for x in r] for r in fij]; fijkV = [[[sw(x) for x in r] for r in rr] for rr in fijk]
def dot(a, b):
    out = {}
    for x, y in zip(a, b): out = add(out, mul(x, y))
    return out
def pw(P, k): return power(P, k)
C = []
def addc(fam, label, P):
    if not P: return
    par = None
    C.append((fam, label, S(P), 'uv'))
# ---- degree 8
# Hessian 2x2 cofactors of f (degree 8 ternary covariants), sym/anti
cof = {}
cof['00'] = sub(mul(fij[1][1], fij[2][2]), mul(fij[1][2], fij[1][2]))
cof['11'] = sub(mul(fij[0][0], fij[2][2]), mul(fij[0][2], fij[0][2]))
cof['22'] = sub(mul(fij[0][0], fij[1][1]), mul(fij[0][1], fij[0][1]))
cof['01'] = sub(mul(fij[0][2], fij[1][2]), mul(fij[0][1], fij[2][2]))
cof['02'] = sub(mul(fij[0][1], fij[1][2]), mul(fij[0][2], fij[1][1]))
cof['12'] = sub(mul(fij[0][1], fij[0][2]), mul(fij[1][2], fij[0][0]))
for k, P in cof.items():
    addc('hess-cofactor', f'cof{k}(U)+cof{k}(V)', sym(P)); addc('hess-cofactor', f'cof{k}(U)-cof{k}(V)', anti(P))
tr_cof = add(add(cof['00'], cof['11']), cof['22'])
addc('hess-cofactor', 'trace-cof(U)+(V)', sym(tr_cof)); addc('hess-cofactor', 'trace-cof(U)-(V)', anti(tr_cof))
# Hessian pairing and Laplacian pairings (degree 8, swap-symmetric)
HP = {}
for i in range(3):
    for j in range(3):
        HP = add(HP, mul(fij[i][j], fijV[i][j]))
addc('hess-pairing', 'sum_ij f_ij(U) f_ij(V)', HP)
lapU = add(add(fij[0][0], fij[1][1]), fij[2][2]); lapV = sw(lapU)
addc('hess-pairing', 'lap(U) lap(V)', mul(lapU, lapV))
addc('hess-pairing', 'lap(U)^2 + lap(V)^2', sym(mul(lapU, lapU)))
addc('hess-pairing', 'lap(U)^2 - lap(V)^2', anti(mul(lapU, lapU)))
addc('hess-pairing', 'sum_i f_ii(U) f_ii(V)', dot([fij[i][i] for i in range(3)], [fijV[i][i] for i in range(3)]))
# third-polar matrix P_ij(U;V) = sum_k f_ijk(U) V_k  (degree 4), pairings (degree 8)
PUV = [[dot([fijk[i][j][k] for k in range(3)], V) for j in range(3)] for i in range(3)]
PVU = [[sw(PUV[i][j]) for j in range(3)] for i in range(3)]
T1 = {}
for i in range(3):
    for j in range(3):
        T1 = add(T1, mul(PUV[i][j], PVU[i][j]))
addc('third-polar', 'sum_ij P_ij(U;V) P_ij(V;U)', T1)
T2 = {}
for i in range(3):
    for j in range(3):
        T2 = add(T2, mul(PUV[i][j], fijV[i][j]))
addc('third-polar', 'sum_ij P_ij(U;V) f_ij(V) + swap', sym(T2)); addc('third-polar', 'sum_ij P_ij(U;V) f_ij(V) - swap', anti(T2))
T3 = {}
for i in range(3):
    for j in range(3):
        T3 = add(T3, mul(PUV[i][j], PUV[i][j]))
addc('third-polar', 'sum_ij P_ij(U;V)^2 + swap', sym(T3)); addc('third-polar', 'sum_ij P_ij(U;V)^2 - swap', anti(T3))
addc('third-polar', 'det P(U;V) ... skipped (deg 12)', {})
# polar octics
pol = dot(fi, [pw(V[i], 3) for i in range(3)])
addc('polar', 'sum_i f_i(U) V_i^3 + swap', sym(pol)); addc('polar', 'sum_i f_i(U) V_i^3 - swap', anti(pol))
polU = dot(fi, [pw(U[i], 3) for i in range(3)])
addc('polar', 'sum_i f_i(U) U_i^3 + swap', sym(polU)); addc('polar', 'sum_i f_i(U) U_i^3 - swap', anti(polU))
polm = dot(fi, [mul(mul(U[i], U[i]), V[i]) for i in range(3)])
addc('polar', 'sum_i f_i(U) U_i^2 V_i + swap', sym(polm)); addc('polar', 'sum_i f_i(U) U_i^2 V_i - swap', anti(polm))
polm2 = dot(fi, [mul(mul(V[i], V[i]), U[i]) for i in range(3)])
addc('polar', 'sum_i f_i(U) V_i^2 U_i + swap', sym(polm2)); addc('polar', 'sum_i f_i(U) V_i^2 U_i - swap', anti(polm2))
# second polar quadric  Q2(U;V) = sum_ij f_ij(U) V_i V_j (deg 6) times quadrics -> products (controls)
Q2 = {}
for i in range(3):
    for j in range(3):
        Q2 = add(Q2, mul(mul(fij[i][j], V[i]), V[j]))
UV = dot(U, V)
addc('control-product', 'Q2(U;V) (U.V) + swap', sym(mul(Q2, UV)))
addc('control-product', 'g * (U.V)  [g = f(U)-f(V)]', mul(sub(fU, fV), UV))
# Fermat / power sums
F8 = add(add(pw(U[0], 8), pw(U[1], 8)), pw(U[2], 8))
addc('fermat', 'sum U^8 + sum V^8', sym(F8)); addc('fermat', 'sum U^8 - sum V^8', anti(F8))
# Pluecker minors w_k = U_i V_j - U_j V_i (swap-antisymmetric)
w = [sub(mul(U[1], V[2]), mul(U[2], V[1])), sub(mul(U[2], V[0]), mul(U[0], V[2])), sub(mul(U[0], V[1]), mul(U[1], V[0]))]
addc('pluecker', 'sum w_k^4', add(add(pw(w[0], 4), pw(w[1], 4)), pw(w[2], 4)))
addc('pluecker', 'sum_k w_k^2 w_{k+1}^2', add(add(mul(pw(w[0], 2), pw(w[1], 2)), mul(pw(w[1], 2), pw(w[2], 2))), mul(pw(w[2], 2), pw(w[0], 2))))
addc('pluecker', 'HP-like: sum_k w_k^2 * lap(U)lap(V)^0 .. (sum w_k^2)*(sum_i f_ii(U)+f_ii(V))', mul(add(add(pw(w[0], 2), pw(w[1], 2)), pw(w[2], 2)), sym(lapU)))
# mixed: f-polar pairing with Pluecker: sum_i f_i(U) f_i(V) is deg 10 (below); here (deg 8) sum_k w_k^2 * HP? no (deg 10).  grad f(U) . (U x V)?? deg 7 -> times (U.V): product.  skip.
# "second polar of f at V in direction U": Q2(V;U) * ... skip.
# ---- degree 9
pol4 = dot(fi, [pw(V[i], 4) for i in range(3)])
addc('polar-d9', 'sum_i f_i(U) V_i^4 + swap', sym(pol4)); addc('polar-d9', 'sum_i f_i(U) V_i^4 - swap', anti(pol4))
pol4U = dot(fi, [pw(U[i], 4) for i in range(3)])
addc('polar-d9', 'sum_i f_i(U) U_i^4 + swap', sym(pol4U)); addc('polar-d9', 'sum_i f_i(U) U_i^4 - swap', anti(pol4U))
pol31 = dot(fi, [mul(pw(U[i], 3), V[i]) for i in range(3)])
addc('polar-d9', 'sum_i f_i(U) U_i^3 V_i + swap', sym(pol31)); addc('polar-d9', 'sum_i f_i(U) U_i^3 V_i - swap', anti(pol31))
# ---- degree 10
PP = dot(fi, fiV)
addc('polar-pairing-d10', 'sum_i f_i(U) f_i(V)', PP)
PP1 = dot(fi, [fiV[(i + 1) % 3] for i in range(3)])
addc('polar-pairing-d10', 'sum_i f_i(U) f_{i+1}(V) + swap', sym(PP1)); addc('polar-pairing-d10', 'sum_i f_i(U) f_{i+1}(V) - swap', anti(PP1))
HPUV = mul(HP, UV)
addc('control-product-d10', 'HP * (U.V)', HPUV)
# grad f(U) x grad f(V) . (U x V)?  = det[grad f(U), grad f(V), U x V]: deg 5+5+2 = 12. skip.  det[grad f(U), grad f(V), U] deg 11 -> sym/anti
detU = {}
M = [fi, fiV, U]
def det3(M):
    out = {}
    for perm, sgn in [((0,1,2),1),((1,2,0),1),((2,0,1),1),((0,2,1),-1),((2,1,0),-1),((1,0,2),-1)]:
        t = mul(mul(M[0][perm[0]], M[1][perm[1]]), M[2][perm[2]])
        out = add(out, scale(t, sgn))
    return out
d11 = det3([fi, fiV, U])
addc('det-d11', 'det[grad f(U), grad f(V), U] + swap', sym(d11)); addc('det-d11', 'det[grad f(U), grad f(V), U] - swap', anti(d11))
# ---- degree 12 (E1 only)
Hf = det3([fij[0], fij[1], fij[2]])
addc('hess-det-d12', 'det Hess f(U) + det Hess f(V)', sym(Hf)); addc('hess-det-d12', 'det Hess f(U) - det Hess f(V)', anti(Hf))
Hmix = det3([[add(fij[i][j], fijV[i][j]) for j in range(3)] for i in range(3)])
addc('hess-det-d12', 'det[Hess f(U) + Hess f(V)]', Hmix)
Hmix2 = det3([[sub(fij[i][j], fijV[i][j]) for j in range(3)] for i in range(3)])
addc('hess-det-d12', 'det[Hess f(U) - Hess f(V)]', Hmix2)
HP3 = det3(PUV)
addc('hess-det-d12', 'det P(U;V) + swap', sym(HP3)); addc('hess-det-d12', 'det P(U;V) - swap', anti(HP3))
C = [c for c in C if c[2] != '0']
json.dump(C, open('w2-candidates.json', 'w'), indent=0)
from collections import Counter
print(len(C), 'candidates;', Counter(core.degree(core.parse(c[2])) for c in C))
for c in C:
    P = L.uv_to_xy(core.parse(c[2]))
    print(core.degree(P), core.parity_of(P), len(P), c[0], '|', c[1])
