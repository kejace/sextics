"""Numerical screening of the homogeneous quadratic obstruction system (from obsys *.eqs.json, coefficients mod p lifted to
symmetric residues -- NOTE: only meaningful if the system was built with exact rational coefficients; here we instead rebuild
the quadrics over the reals from the exact rational data when available, otherwise skip).  This script expects an .eqs.json
produced with p = 0 (exact)  -- see obsys_exact.py.  Gauss-Newton / Levenberg-Marquardt on the unit sphere from random starts;
reports the smallest residuals.  Screening only (no rigor)."""
import sys, json, time
import numpy as np
fn = sys.argv[1]; nstart = int(sys.argv[2]) if len(sys.argv) > 2 else 200; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
J = json.load(open(fn)); n = J['nU']; eqs = J['eqs']
m = len(eqs)
# build dense symmetric matrices per equation as a (m, n, n) tensor (float64)
Qs = np.zeros((m, n, n))
for k, eq in enumerate(eqs):
    for v, a, b in eq:
        v = float(v)
        if a == b: Qs[k, a, a] += v
        else: Qs[k, a, b] += v / 2; Qs[k, b, a] += v / 2
# normalise equations
for k in range(m):
    nrm = np.linalg.norm(Qs[k]);
    if nrm > 0: Qs[k] /= nrm
def resid(w):
    return np.einsum('kab,a,b->k', Qs, w, w)
def jac(w):
    return 2 * np.einsum('kab,b->ka', Qs, w)
rng = np.random.default_rng(seed)
best = []
t0 = time.time()
for s in range(nstart):
    w = rng.standard_normal(n); w /= np.linalg.norm(w)
    lam = 1e-3
    for it in range(200):
        r = resid(w); Jm = jac(w)
        # project out the radial direction (homogeneity) by adding the sphere constraint row
        A = np.vstack([Jm, 2 * w[None, :]]); bb = np.concatenate([-r, [0.0]])
        step, *_ = np.linalg.lstsq(A, bb, rcond=None)
        w2 = w + step; w2 /= np.linalg.norm(w2)
        if np.linalg.norm(resid(w2)) < np.linalg.norm(r): w = w2
        else:
            w = w + 0.3 * step; w /= np.linalg.norm(w)
        if np.linalg.norm(resid(w)) < 1e-13: break
    f = np.linalg.norm(resid(w))
    best.append((f, s))
    if s % 20 == 0:
        print(f'start {s}: resid {f:.3e}  (best so far {min(best)[0]:.3e}) {time.time()-t0:.0f}s', flush=True)
best.sort()
print('best residuals:', [f'{f:.3e}' for f, _ in best[:10]])
json.dump({'file': fn, 'n': n, 'm': m, 'nstart': nstart, 'best': best[:20]}, open(fn.replace('.eqs.json', f'.num{seed}.json'), 'w'))
