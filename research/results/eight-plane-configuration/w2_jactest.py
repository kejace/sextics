import sys, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
from w2lib import _col
fam = JacobianFamily(6, 0, want_C=False, verbose=False)
rng = np.random.default_rng(3)
u0 = rng.standard_normal(fam.dimP); u0 /= np.linalg.norm(u0)
a0, na0 = fam.a_of_u(u0)
o = fam.linear_solve(a0)
c0 = o['c']; Ubk = o['Ubk']; U = fam.U
def r_fix(u):
    a, _ = fam.a_of_u(u)
    D = _col(fam.Dfull @ a); Om = fam.pat.omega(D)
    g = Om @ (_col(fam.Nfull @ a) + fam.Kmat @ c0)
    g = g - U @ (U.T @ g); g = g - Ubk @ (Ubk.T @ g)
    return g
J = fam.jac_u(u0)
h = 1e-6; errs = []
for k in range(5):
    du = np.zeros(fam.dimP); du[k] = h
    fd = (r_fix(u0 + du) - r_fix(u0 - du)) / (2*h)
    errs.append(float(np.linalg.norm(fd - J[:, k]) / max(np.linalg.norm(fd), 1e-300)))
print('fixed-x FD relative errors (should be ~1e-7):', errs)
# residual at u0 equals r_fix(u0)?
print('|r - r_fix|/|r| =', np.linalg.norm(o['r'] - r_fix(u0)) / np.linalg.norm(o['r']))
