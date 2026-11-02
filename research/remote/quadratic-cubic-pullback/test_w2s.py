#!/usr/bin/env python3
"""Self-tests of w2s.py: exact reduction identity, point tables vs exact evaluation, Fermat control."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import numpy as np
from fractions import Fraction
import w2s
from w2s import core, mul, add, sub, scale, shift, pderiv, monomials, NV

rng = np.random.default_rng(1)


def randpoly(deg, par, nterms=8, den=3):
    lst = list(monomials(deg, par))
    out = {}
    for _ in range(nterms):
        e = lst[rng.integers(len(lst))]
        out[e] = out.get(e, 0) + Fraction(int(rng.integers(-5, 6)), den)
    return {e: c for e, c in out.items() if c}


# 1. exact reduction identity on the actual F, d = 6, D in J_6 (w linear), random c (deg 2) and G (deg 6)
d = 6
pN = w2s.numerator_parities()
w = [randpoly(1, 0 if i < 3 else 1, 3) for i in range(NV)]
D = {}
for i in range(NV):
    D = add(D, mul(w[i], w2s.Fi[i]))
u = [mul(w2s.A, w[i]) for i in range(NV)]
cpar = w2s.c_parities()
c = {ij: randpoly(d - 4, par, 4) for ij, par in cpar.items()}
G = randpoly(d, 0, 5)
ok = w2s.reduced_identity_check(D, u, c, G)
print('exact reduction identity (d=6, with G):', ok)
assert ok
ok = w2s.reduced_identity_check(D, u, c, None)
print('exact reduction identity (d=6, no G):', ok)

# 2. point tables vs exact: family J6, random a, random c  -> P(z_s) from tables == exact P evaluated
fam = w2s.jacobian_family(6)
pv = w2s.PointVarPro(fam, S=40, seed=3, verbose=True)
a = rng.standard_normal(fam.K)
# exact D, u from the family
frac = [Fraction(float(x)).limit_denominator(1000) for x in a]
Dex = w2s.D_exact_from_family(fam, frac)
uex = [{} for _ in range(NV)]
for k, (i, m) in enumerate(fam.labels):
    uex[i] = add(uex[i], scale(shift(w2s.A, m, 1), frac[k]))
af = np.array([float(x) for x in frac])
# random c as a flat vector consistent with pv's layout
cvec = rng.standard_normal(pv.nc)
cdict = {}
for b, ((i, j), par, n) in enumerate(pv.cblocks):
    lst, _ = w2s.mons(pv.dc, par)
    o = pv.coff[b]
    cdict[(i, j)] = {lst[t]: Fraction(float(cvec[o + t])).limit_denominator(1000) for t in range(n)}
    for t in range(n):
        cvec[o + t] = float(cdict[(i, j)][lst[t]])
# exact P
N = [dict(uex[i]) for i in range(NV)]
for (i, j), cij in cdict.items():
    N[i] = add(N[i], mul(cij, w2s.Fi[j])); N[j] = sub(N[j], mul(cij, w2s.Fi[i]))
divN = {}
for i in range(NV):
    divN = add(divN, pderiv(N[i], i))
Pex = mul(Dex, divN)
for i in range(NV):
    Pex = sub(Pex, mul(N[i], pderiv(Dex, i)))
Z = pv.Z
E_ = np.array(list(Pex.keys())); C_ = np.array([float(v) for v in Pex.values()])
Pex_at = (np.prod(Z[:, None, :] ** E_[None, :, :], axis=2) * C_[None, :]).sum(axis=1)
Dv, dDv, uv, divuv = pv._state(af)
p0 = Dv * divuv - sum(uv[i] * dDv[i] for i in range(NV))
Em = pv._Ematrix(Dv, dDv)
Ptab = Em @ cvec + p0
print('point table vs exact P: max abs diff', np.max(np.abs(Ptab - Pex_at)), 'scale', np.max(np.abs(Pex_at)))
assert np.max(np.abs(Ptab - Pex_at)) < 1e-8 * max(1, np.max(np.abs(Pex_at)))

# 3. gradient check (finite differences in reduced coords)
pv2 = w2s.PointVarPro(fam, S=200, seed=5)
u0 = rng.standard_normal(pv2.dimP); u0 /= np.linalg.norm(u0)
f0, g0 = pv2.value_and_grad(u0)
eps = 1e-6
errs = []
for t in range(5):
    e = np.zeros_like(u0); e[rng.integers(len(u0))] = 1.0
    fp, _ = pv2.value_and_grad(u0 + eps * e); fm, _ = pv2.value_and_grad(u0 - eps * e)
    fd = (fp - fm) / (2 * eps)
    errs.append(abs(fd - g0 @ e) / max(abs(fd), 1e-12))
print('gradient relative errors (FD):', ['%.2e' % x for x in errs])

# 4. Fermat control: F0 = sum Z_i^6, A = X0^3 Y0^3, J6 family; known solution w = (X0, X1, X2, 0,0,0)
F0 = core.parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
w2s.configure(F0, w2s.A)
famF = w2s.jacobian_family(6)
pvF = w2s.PointVarPro(famF, S=300, seed=7, verbose=True)
a_star = np.zeros(famF.K)
for k, (i, m) in enumerate(famF.labels):
    if i < 3 and m[i] == 1 and sum(m) == 1:
        a_star[k] = 1.0
a_star /= np.linalg.norm(a_star)
rel, f, c = pvF.rel(a_star)
print('Fermat control at the known witness: f=%.3e rel=%.3e' % (f, rel))
assert rel < 1e-24
t0 = time.time()
res = pvF.search(n_restarts=20, seed=11, maxiter=400, verbose=False)
print('Fermat control random restarts: best rel %.3e (spurious %s) in %.1fs; rels:' % (res[0]['rel'], res[0]['spurious'], time.time() - t0),
      ['%.1e' % r['rel'] for r in res[:10]])
# exact confirmation of the best candidate through rational reconstruction + exact reduced solve (Fermat)
best = res[0]
for fr in w2s.reconstruct_vector(best['a']):
    Dc = w2s.D_exact_from_family(famF, fr)
    print('  reconstructed D =', core.to_str(Dc)[:100])
    out = w2s.exact_reduced_solve(Dc, verbose=True, return_solution=True)
    print('  exact reduced solve:', {k: v for k, v in out.items() if k not in ('solution', 'u')})
    if out.get('consistent'):
        break
# exact reduced solve on the known Fermat witness denominator
Dw = core.parse('X0**6 + X1**6 + X2**6')
out = w2s.exact_reduced_solve(Dw, verbose=True, return_solution=True)
print('Fermat D=X0^6+X1^6+X2^6 exact reduced:', {k: v for k, v in out.items() if k not in ('solution', 'u')})
# verify the mod-p solution via core.defect mod p
p = out['p']
N, Q = w2s.numerators_from_reduced(Dw, out['u'], out['solution'], p)
F0p = core.reduce_mod(F0, p); Ap = core.reduce_mod(w2s.A, p); Dp = core.reduce_mod(Dw, p)
E = core.defect(F0p, Ap, Dp, 1, N, p)
F2Q = mul(mul(F0p, F0p, p), Q, p)
print('  mod-p witness identity E == F^2 Q:', sub(E, F2Q, p) == {})
# back to the real instance: polar d=5 D = F_0 (E1 true) -> expect inconsistent; join012 -> E1 false
w2s.configure(*core.load_instance('quadratic-cubic-pullback')[:2])
D5 = w2s.Fi[0]
out = w2s.exact_reduced_solve(D5, verbose=True)
print('ours, D = F_0 (polar, d=5):', out)
EOF_MARKER = None
