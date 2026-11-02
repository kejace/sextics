import sys, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
from scipy.optimize import least_squares
fam = JacobianFamily(6, int(sys.argv[1]) if len(sys.argv) > 1 else 0, want_C=False, verbose=False)
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
u0 = rng.standard_normal(fam.dimP); u0 /= np.linalg.norm(u0)
traj = []
def fun(u):
    r = fam.residual_u(u); traj.append(float(r @ r)); return r
t0 = time.time()
res = least_squares(fun, u0, jac=fam.jac_u, method='lm', max_nfev=400, xtol=1e-15, ftol=1e-15, gtol=1e-15)
a, _ = fam.a_of_u(res.x); f, o = fam.value(a)
print('nfev', res.nfev, 'status', res.status, 'final f', f, 'rel', math.sqrt(f)/o['bnorm'], 'time', round(time.time()-t0,1))
import numpy as np
tr = np.array(traj)
idx = sorted(set([0,1,2,3,5,8,12,20,30,50,80,120,200,300,len(tr)-1]))
print('cost trajectory:', [(i, float(f'{tr[i]:.4g}')) for i in idx if i < len(tr)])
# trf for comparison
traj.clear(); t0 = time.time()
res2 = least_squares(fun, u0, jac=fam.jac_u, method='trf', max_nfev=200, xtol=1e-15, ftol=1e-15, gtol=1e-15)
a2, _ = fam.a_of_u(res2.x); f2, o2 = fam.value(a2)
print('trf: nfev', res2.nfev, 'status', res2.status, 'final f', f2, 'rel', math.sqrt(f2)/o2['bnorm'], 'time', round(time.time()-t0,1))
