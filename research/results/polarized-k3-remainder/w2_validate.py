"""Validate the reduced E2 system against the cleared identity (exact mod p) at small degree."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np
from fractions import Fraction as Fr
import w2lib as L
from hodge import core
from hodge.core import NV, mul, add, sub, scale, shift, pderiv, power
d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
pD = int(sys.argv[2]) if len(sys.argv) > 2 else 0
p = L.P0
F, A, q = L.load()
t0 = time.time()
S = L.ReducedE2.jacobian(d, pD, build_float=True)
print('build', round(time.time()-t0,1), 's; K', S.K, 'dimP', S.dimP, 'nc', S.nc, 'nred', S.nred, 'nQ', S.S_Q.n)
rng = np.random.default_rng(0)
# random rational a and c
a = [Fr(int(x)) for x in rng.integers(-3, 4, size=S.K)]
D = {}
for k, (i0, m, Dk) in enumerate(S.basis):
    if a[k]: D = add(D, scale(Dk, a[k]))
N0 = [{} for _ in range(NV)]
for k, (i0, m, Dk) in enumerate(S.basis):
    if a[k]: N0[i0] = add(N0[i0], scale(mul(A, {m: Fr(1)}), a[k]))
# check dF(N0) = A D
dFN = {}
for i in range(NV): dFN = add(dFN, mul(N0[i], S.Fi[i]))
print('E1 exact (dF(N0) - A D == 0):', not sub(dFN, mul(A, D)))
cvals = {}
for col, (i, j, n) in enumerate(S.cpairs):
    v = int(rng.integers(-2, 3))
    if v: cvals[(i, j, n)] = Fr(v)
N = S.numerators_from_c(N0, cvals)
E = core.defect(F, A, D, 1, N)
# B_D(N)
divN = {}; NdD = {}
for i in range(NV):
    divN = add(divN, pderiv(N[i], i)); NdD = add(NdD, mul(N[i], pderiv(D, i)))
B = sub(mul(D, divN), NdD)
print('E == F * B_D(N):', not sub(E, mul(F, B)))
# compare with assembled exact system: Mc c + b0 == vec(B)  (mod p)
M, R, C = S.exact_system_modp(D, N0, p)
Mc = np.array([[int(M[r, c]) for c in range(S.nc)] for r in range(R)], dtype=object)  # slow but small
cvec = np.array([int(cvals.get(cp, 0)) % p for cp in S.cpairs], dtype=object)
b0 = np.array([(-int(M[r, C])) % p for r in range(R)], dtype=object)
lhs = (Mc.dot(cvec) + b0) % p
Bv = S.S_B.vec(B, dtype=np.int64, p=p)
print('assembled Mc c + b0 == vec(B) mod p:', all(int(lhs[r]) == int(Bv[r]) for r in range(R)))
# F-multiple vectors: D(a_F) = F h ?
SF = S.SF
h0 = L.monomials(d-6, pD)[0]
aF = SF[:, 0]
DF = L.a_to_D(S.basis, aF)
print('F-multiple check:', not sub(DF, shift(F, h0, 1)) if DF else 'empty')
# float residual at a random point vs exact structure: value of f and gauge
af = np.array([float(x) for x in a]); af /= np.linalg.norm(af)
t1 = time.time(); asm = S.assemble(af); res = S.residual(af, asm, with_jac=True); print('float residual', res['f'], 'rel', res['rel'], 'rank', res['rank'], 'time', round(time.time()-t1,2))
# finite-difference check of Kaufman Jacobian (directional)
u = S.P.T @ af; u /= np.linalg.norm(u)
def resf(u):
    a_ = S.P @ u
    asm_ = S.assemble(a_); r_ = S.residual(a_, asm_, with_jac=True)
    return r_['f'], r_['r'], r_['J'] @ S.P
f0, r0, J0 = resf(u)
v = rng.standard_normal(len(u)); v -= (v @ u) * u; v /= np.linalg.norm(v)
eps = 1e-6
f1, r1, _ = resf((u + eps*v)/np.linalg.norm(u+eps*v))
f2, r2, _ = resf((u - eps*v)/np.linalg.norm(u-eps*v))
dfd = (f1 - f2) / (2*eps); dan = 2 * r0 @ (J0 @ v)
print('directional derivative FD vs Kaufman:', dfd, dan)
# exact consistency of this random D (expected inconsistent)
t2 = time.time(); out = S.exact_consistent(D, N0, p=p); print('exact:', {k: out[k] for k in ('consistent','rank','nullity','rows','cols','seconds')})
