"""Numerical (float) LM search for the full degree-7 Hessian relaxation  w^T Hess(F) w in J_8 (both parities),
w in C^126 (all 6 x 21 quadratic coefficients, mixed tau-parity allowed), gauge not fixed, normalization |w|=1 via
random complex linear slice L(w)=1.  Reports best residuals.  Usage: python3 lm_mixed7.py nstarts seed [sector]
sector: mixed (default) | even | odd  (even/odd restrict w to the tau-pure sector as a cross-check of msolve)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
import numpy as np
from scipy.optimize import least_squares
nstarts = int(sys.argv[1]); seed = int(sys.argv[2]); sector = sys.argv[3] if len(sys.argv) > 3 else 'mixed'
rng = np.random.default_rng(seed)
F, A, q = load()
Fs = scale(F, 9973)   # integer coefficients; J unchanged
Fi = [pderiv(Fs, i) for i in range(NV)]
H = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
mons8 = monomials(8, 0) + monomials(8, 1); idx8 = {m: j for j, m in enumerate(mons8)}; n8 = len(mons8)
# J_8 generators (all parities)
gens = []
for i in range(NV):
    for m in monomials(3, 0) + monomials(3, 1):
        gens.append(shift(Fi[i], m, 1))
G = np.zeros((len(gens), n8))
for r, g in enumerate(gens):
    for mm, v in g.items(): G[r, idx8[mm]] = float(v)
U, s, Vt = np.linalg.svd(G, full_matrices=True)
rk = int((s > 1e-9 * s[0]).sum())
C = Vt[rk:]            # orthonormal basis of the complement of J_8  (rows), shape (n8-rk, n8)
print(f'dim J_8 = {rk} (expect 336), complement = {C.shape[0]} (expect 951)', flush=True)
# w basis
basis = []
for i in range(NV):
    pars = [0, 1] if sector == 'mixed' else ([(0 if i < 3 else 1)] if sector == 'even' else [(1 if i < 3 else 0)])
    for par in pars:
        for m in monomials(2, par):
            if i == 0 and m[0] >= 1: continue   # gauge complement: drop w_0 coefficients of X0*l
            basis.append((i, m))
nw = len(basis)
T = np.zeros((C.shape[0], nw, nw))
for k in range(nw):
    ik, mk = basis[k]
    for l in range(k, nw):
        il, ml = basis[l]
        P = shift(H[ik][il], tuple(a + b for a, b in zip(mk, ml)), 1)
        if not P: continue
        v = np.zeros(n8)
        for mm, c in P.items(): v[idx8[mm]] = float(c)
        cv = C @ v
        T[:, k, l] += cv
        if l != k: T[:, l, k] += cv
print(f'nw={nw} (gauge removed), tensor built', flush=True)
L = rng.standard_normal(nw) + 1j * rng.standard_normal(nw)
def resid(x):
    w = x[:nw] + 1j * x[nw:]
    r = np.einsum('mkl,k,l->m', T, w, w)
    lin = L @ w - 1.0
    return np.concatenate([r.real, r.imag, [lin.real, lin.imag]])
def jac(x):
    w = x[:nw] + 1j * x[nw:]
    Jc = 2 * np.einsum('mkl,l->mk', T, w)       # d r / d w (complex)
    # real/imag blocks: r = Jc w ; d(Re r)/d u = Re Jc, d(Re r)/d v = -Im Jc ; d(Im r)/d u = Im Jc, d(Im r)/d v = Re Jc
    top = np.hstack([Jc.real, -Jc.imag]); bot = np.hstack([Jc.imag, Jc.real])
    lt = np.hstack([L.real, -L.imag]); lb = np.hstack([L.imag, L.real])
    return np.vstack([top, bot, lt[None, :], lb[None, :]])
best = []
t0 = time.time()
for s_ in range(nstarts):
    x0 = rng.standard_normal(2 * nw)
    w0 = x0[:nw] + 1j * x0[nw:]; w0 /= (L @ w0); x0 = np.concatenate([w0.real, w0.imag])
    res = least_squares(resid, x0, jac=jac, method='lm', xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=400)
    w = res.x[:nw] + 1j * res.x[nw:]
    val = np.linalg.norm(resid(res.x)[:-2]) / max(1e-300, np.linalg.norm(w) ** 2)   # scale-free residual
    best.append(float(val))
    if val < 1e-8:
        print(f'!!! start {s_}: scale-free residual {val:.3e}', flush=True)
        np.save(f'/tmp/hodge/results/dense-rational/lm_{sector}_hit_{seed}_{s_}.npy', w)
    if (s_ + 1) % 50 == 0:
        print(f'{s_+1} starts, best so far {min(best):.3e}, median {np.median(best):.3e} ({time.time()-t0:.0f}s)', flush=True)
best.sort()
out = {'sector': sector, 'nw': nw, 'neq': 2 * C.shape[0], 'nstarts': nstarts, 'seed': seed, 'best10': best[:10],
       'median': float(np.median(best)), 'seconds': round(time.time() - t0, 1)}
print(json.dumps(out), flush=True)
with open('/tmp/hodge/results/dense-rational/lm_mixed7.jsonl', 'a') as f: f.write(json.dumps(out) + '\n')
