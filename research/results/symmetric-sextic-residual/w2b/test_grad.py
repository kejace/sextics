import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2lib import *
from scipy.optimize import minimize
F, A, q = load()
fam, info = build_sector_family(F, A, q, 'triv', 7, 1)
rng = np.random.default_rng(5)
u = rng.standard_normal(fam.dimP); u /= np.linalg.norm(u)
f, g = fam.value_and_grad(u)
for trial in range(3):
    v = rng.standard_normal(fam.dimP); eps = 1e-6
    f1, _ = fam.value_and_grad(u + eps * v); f0, _ = fam.value_and_grad(u - eps * v)
    print(f'grad check: fd {(f1-f0)/(2*eps):.8g} analytic {float(g@v):.8g}')
res = minimize(fam.value_and_grad, u, jac=True, method='L-BFGS-B', options={'maxiter': 400, 'gtol': 1e-12, 'ftol': 1e-18, 'maxcor': 20})
a = fam.P @ res.x; a /= np.linalg.norm(a)
print('L-BFGS:', res.message, 'nit', res.nit, 'nfev', res.nfev, 'f', res.fun, 'rel', fam.rel(a)[0], '|grad|', np.linalg.norm(res.jac))
a2, x2, rel2, hist = polish(fam, a)
print('GN polish from L-BFGS point:', [f'{h:.3e}' for h in hist])
u2 = fam.P.T @ a2
res2 = minimize(fam.value_and_grad, u2, jac=True, method='L-BFGS-B', options={'maxiter': 400, 'gtol': 1e-14, 'ftol': 1e-20, 'maxcor': 20})
a3 = fam.P @ res2.x; a3 /= np.linalg.norm(a3)
print('L-BFGS from GN point:', res2.message, 'nit', res2.nit, 'rel', fam.rel(a3)[0])
# try a different optimizer: trust-region / CG from the original point
res3 = minimize(fam.value_and_grad, u, jac=True, method='CG', options={'maxiter': 2000, 'gtol': 1e-12})
a4 = fam.P @ res3.x; a4 /= np.linalg.norm(a4)
print('CG:', res3.message, 'nit', res3.nit, 'rel', fam.rel(a4)[0])
# objective along the segment from a (lbfgs) to a2 (GN)
for t in (0, 0.25, 0.5, 0.75, 1.0):
    at = (1-t)*a + t*a2; print(f'  t={t}: rel {fam.rel(at)[0]:.4e}')
