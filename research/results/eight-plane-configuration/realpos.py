import sys, numpy as np, time
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
from scipy.optimize import minimize
F, A, q, rec = core.load_instance('eight-plane-configuration')
exps = np.array([e for e in F]); coefs = np.array([float(v) for v in F.values()])
def f(z):
    z = np.asarray(z, float); z = z/np.linalg.norm(z)
    return float(np.sum(coefs * np.prod(z[None,:]**exps, axis=1)))
rng = np.random.default_rng(0)
Z = rng.normal(size=(200000,6)); Z /= np.linalg.norm(Z,axis=1)[:,None]
vals = np.array([np.sum(coefs*np.prod(z[None,:]**exps,axis=1)) for z in Z])
print('random sample min/max on sphere:', vals.min(), vals.max(), flush=True)
best=np.inf; bz=None
for k in range(40):
    z0 = Z[np.argsort(vals)[k]]
    res = minimize(f, z0, method='Powell', options={'xtol':1e-8,'ftol':1e-12,'maxiter':3000})
    v = f(res.x)
    if v < best: best=v; bz=res.x/np.linalg.norm(res.x)
print('min after local minimization from 40 best samples:', best, np.round(bz,4), flush=True)
