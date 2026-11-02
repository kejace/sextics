#!/usr/bin/env python3
"""Numerical (float) Levenberg-Marquardt screen of the degree-8 Hessian relaxation in a tau-pure sector:
   A * w^T Hess(F) w  in  J_16^-   (w cubic vector field, D = w(F) of degree 8, parity(D) = sector),
i.e. the component of A w^T H w orthogonal to J_16^- vanishes (468 complex equations, 156 (even) / 159 (odd) unknowns
after the Euler gauge w -> w + l*Euler (l quadratic) is fixed by dropping the X0-divisible coefficients of w_0).
This is the dossier's item 5.3(b): the first search at degree 8. A scale-free residual < 1e-8 would be a candidate to
promote exactly (then E2 must still be solved for K and C). Generic emptiness of this system is NOT certifiable by F4
(Froeberg count positive through Macaulay degree 8), so only numerical evidence is possible here.
Usage: python3 wave2_excl_lm_d8.py sector nstarts seed    (sector in even|odd)"""
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
import numpy as np
from scipy.optimize import least_squares
sector = sys.argv[1]; nstarts = int(sys.argv[2]); seed = int(sys.argv[3])
rng = np.random.default_rng(seed)
F, A, q = load()
Fs = scale(F, 9973)                      # integer coefficients; J unchanged
As = A
Fi = [pderiv(Fs, i) for i in range(NV)]
H = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
mons16 = monomials(16, 1); idx16 = {m: j for j, m in enumerate(mons16)}; n16 = len(mons16)
cache = f'{OUT}/lm_d8_complement_J16odd.npy'
t0 = time.time()
if os.path.exists(cache):
    C = np.load(cache)
else:
    gens = []
    for i in range(NV):
        for m in monomials(11, (1 - par_i(i)) % 2):
            gens.append(shift(Fi[i], m, 1))
    G = np.zeros((len(gens), n16))
    for r, g in enumerate(gens):
        for mm, v in g.items(): G[r, idx16[mm]] = float(v)
    # orthonormal basis of the orthogonal complement of the row space of G (= J_16^-): right singular vectors for zero singular values
    U, s, Vt = np.linalg.svd(G, full_matrices=True)
    rk = int((s > 1e-9 * s[0]).sum())
    C = np.ascontiguousarray(Vt[rk:])
    print(f'dim J_16^- = {rk} (expect 9684), complement = {C.shape[0]} (expect 468), SVD {time.time()-t0:.0f}s', flush=True)
    np.save(cache, C)
assert C.shape[0] == 468, C.shape
pd = 0 if sector == 'even' else 1
basis = []
for i in range(NV):
    par = (par_i(i) + pd) % 2
    for m in monomials(3, par):
        if i == 0 and m[0] >= 1: continue      # gauge: w_0 has no X0-divisible part
        basis.append((i, m))
nw = len(basis)
T = np.zeros((468, nw, nw))
Af = {m: float(v) for m, v in As.items()}
for k in range(nw):
    ik, mk = basis[k]
    for l in range(k, nw):
        il, ml = basis[l]
        P = shift(H[ik][il], tuple(a + b for a, b in zip(mk, ml)), 1)
        if not P: continue
        P = mul(As, P)
        v = np.zeros(n16)
        for mm, c in P.items(): v[idx16[mm]] = float(c)
        cv = C @ v
        T[:, k, l] += cv
        if l != k: T[:, l, k] += cv
print(f'sector={sector} nw={nw} (gauge removed), tensor built ({time.time()-t0:.0f}s)', flush=True)
L = rng.standard_normal(nw) + 1j * rng.standard_normal(nw)
def resid(x):
    w = x[:nw] + 1j * x[nw:]
    r = np.einsum('mkl,k,l->m', T, w, w)
    lin = L @ w - 1.0
    return np.concatenate([r.real, r.imag, [lin.real, lin.imag]])
def jac(x):
    w = x[:nw] + 1j * x[nw:]
    Jc = 2 * np.einsum('mkl,l->mk', T, w)
    top = np.hstack([Jc.real, -Jc.imag]); bot = np.hstack([Jc.imag, Jc.real])
    lt = np.hstack([L.real, -L.imag]); lb = np.hstack([L.imag, L.real])
    return np.vstack([top, bot, lt[None, :], lb[None, :]])
scaleT = np.sqrt((T ** 2).sum() / 468)
best = []
t1 = time.time()
for s_ in range(nstarts):
    x0 = rng.standard_normal(2 * nw)
    w0 = x0[:nw] + 1j * x0[nw:]; w0 /= (L @ w0); x0 = np.concatenate([w0.real, w0.imag])
    res = least_squares(resid, x0, jac=jac, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=600)
    w = res.x[:nw] + 1j * res.x[nw:]
    val = np.linalg.norm(resid(res.x)[:-2]) / max(1e-300, scaleT * np.linalg.norm(w) ** 2)
    best.append(float(val))
    if val < 1e-8:
        print(f'!!! start {s_}: scale-free residual {val:.3e}', flush=True)
        np.save(f'{OUT}/lm_d8_{sector}_hit_{seed}_{s_}.npy', w)
    if (s_ + 1) % 25 == 0:
        print(f'{s_+1} starts, best so far {min(best):.3e}, median {np.median(best):.3e} ({time.time()-t1:.0f}s)', flush=True)
best.sort()
out = {'exp': 'lm-hessian-relaxation-d8', 'd': 8, 'sector': sector, 'unknowns': nw, 'equations_complex': 468, 'nstarts': nstarts, 'seed': seed,
       'best10_scalefree_residual': best[:10], 'median': float(np.median(best)), 'seconds': round(time.time() - t0, 1),
       'result': 'no near-solution (all scale-free residuals >= 1e-8)' if best[0] >= 1e-8 else 'NEAR-SOLUTION FOUND: promote exactly'}
print(json.dumps(out), flush=True)
log2(**out)
