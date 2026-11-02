import sys, time, numpy as np
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-search')
import w2ep
from w2ep import *
from fractions import Fraction
F, A, q = load()
d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
parD = int(sys.argv[2]) if len(sys.argv) > 2 else 0
fam = jacobian_family(F, A, d, parD)
t0 = time.time(); obj = ReducedE2(F, A, fam, verbose=True); print('build', round(time.time()-t0,1))
rng = np.random.default_rng(1)
# Test 1: exact residual comparison at d (small ints)
w = rng.integers(-2, 3, size=obj.K).astype(float)
c = rng.integers(-2, 3, size=obj.nc).astype(float)
Fi = [pderiv(F, i) for i in range(6)]
# exact D, N
D = {}; N = [dict() for _ in range(6)]
for k, (i, m) in enumerate(fam['wdesc']):
    if w[k]:
        D = add(D, scale(shift(Fi[i], m, 1), Fraction(int(w[k]))))
        N[i] = add(N[i], scale(shift(A, m, 1), Fraction(int(w[k]))))
# Koszul
for p, (i, j) in enumerate(obj.pairs):
    Sc = obj.cspaces[p]; cp = c[obj.cslices[p]]
    cij = {Sc.mons[t]: Fraction(int(cp[t])) for t in range(Sc.n) if cp[t]}
    if cij:
        N[i] = add(N[i], mul(cij, Fi[j])); N[j] = sub(N[j], mul(cij, Fi[i]))
# check E1 exact
e1 = sub(sum([mul(N[i], Fi[i]) for i in range(6)][1:], mul(N[0], Fi[0])) if False else None, None) if False else None
S = mul(N[0], Fi[0])
for i in range(1, 6): S = add(S, mul(N[i], Fi[i]))
print('E1 exact (N(F) == A D):', sub(S, mul(A, D)) == {})
# B_D(N) exact
divN = {}
for i in range(6): divN = add(divN, pderiv(N[i], i))
B = mul(D, divN)
for i in range(6): B = sub(B, mul(N[i], pderiv(D, i)))
Bv = obj.SR.vec(B)
Dv = obj.D_of(w); Nv = obj.N_of(w)
b = obj.b_vec(Dv, Nv)
# T_D c full (unreduced): recompute T without reduction
Dk = obj.SD.keys; nz = np.nonzero(Dv)[0]
T = np.zeros((obj.SR.n, obj.nc))
dD = []
for i in range(6):
    v = obj.DerD[i] @ Dv; nzi = np.nonzero(v)[0]; dD.append((obj.SdD[i].keys[nzi], v[nzi]))
for p, (i, j) in enumerate(obj.pairs):
    sl = obj.cslices[p]
    T[:, sl] += mult_matrix(Dk[nz], Dv[nz], obj.Sphi[p], obj.SR) @ obj.Phi[p]
    T[:, sl] -= mult_matrix(dD[i][0], dD[i][1], obj.Spsi[p], obj.SR) @ obj.Psi[p]
    T[:, sl] += mult_matrix(dD[j][0], dD[j][1], obj.SpsiP[p], obj.SR) @ obj.PsiP[p]
r_full = T @ c - b
print('Test1 max|r_full - B_exact| =', np.max(np.abs(r_full - Bv)), ' |B| =', np.max(np.abs(Bv)))
# Test 2: reducer kills F-multiples
G = {m: Fraction(int(x)) for m, x in zip(obj.SQ.mons, rng.integers(-3, 4, size=obj.SQ.n)) if x}
FG = obj.SR.vec(mul(F, G))
print('Test2 |Red(F G)| =', np.linalg.norm(obj.red.reduce(FG)), ' |FG| =', np.linalg.norm(FG), ' |Red(random)| =', np.linalg.norm(obj.red.reduce(rng.standard_normal(obj.SR.n))))
# Test 3: gradient check
g, fd, f0 = obj.fd_check(seed=3)
print('Test3 grad check: analytic', g, 'fd', fd, 'f0', f0)
# Test 4: timing
t0 = time.time(); u = rng.standard_normal(obj.dimP); obj.value_and_grad(u); print('one value_and_grad', round(time.time()-t0, 3), 's')
t0 = time.time(); r = obj.minimize_from(0, maxiter=50); print('50-iter minimize', round(time.time()-t0, 2), 's', 'f', r['f'], 'rel', r['rel'], 'iters', r['iters'], 'spurious', r['spurious'])
