import sys, time
sys.path.insert(0,'/tmp/hodge/results/one-reflection-skew')
import w2skew as W
import numpy as np
d = 6
Dm = W.mons5(d, 0, 0)
vp = W.VP5(d, 0, 0, Dm, np.eye(len(Dm)))
rng = np.random.default_rng(5)
u = rng.standard_normal(vp.K)
t0 = time.time(); f0, g = vp.value_and_grad(u); print('one evaluation', round(time.time()-t0, 3), 's', flush=True)
for k in rng.choice(vp.K, 6, replace=False):
    h = 1e-6
    up = u.copy(); up[k] += h; um = u.copy(); um[k] -= h
    fd = (vp.value_and_grad(up)[0] - vp.value_and_grad(um)[0]) / (2 * h)
    print('grad check k=%d analytic=%.6e fd=%.6e relerr=%.2e' % (k, g[k], fd, abs(g[k] - fd) / max(abs(fd), 1e-12)), flush=True)
# 6-var gradient check too, d=7 chi=(0,0,0), mons family (small)
mons = W.mons6_sector(7, (0, 0, 0))
vp6 = W.VP6(7, mons, np.eye(len(mons)), verbose=True)
u = rng.standard_normal(vp6.dimP)
t0 = time.time(); f0, g = vp6.value_and_grad(u); print('vp6 one evaluation', round(time.time()-t0, 3), 's', flush=True)
for k in rng.choice(vp6.dimP, 6, replace=False):
    h = 1e-6
    up = u.copy(); up[k] += h; um = u.copy(); um[k] -= h
    fd = (vp6.value_and_grad(up)[0] - vp6.value_and_grad(um)[0]) / (2 * h)
    print('vp6 grad check k=%d analytic=%.6e fd=%.6e relerr=%.2e' % (k, g[k], fd, abs(g[k] - fd) / max(abs(fd), 1e-12)), flush=True)
print('DONE', flush=True)
