import sys, time, numpy as np
from w2lib import *
from w2blocks import D8, S0, P0
t0 = time.time()
Dg = mul(S0, F)
fam = SectorVarPro([Dg], exclude_gauge=False, name='gauge-test')
f, x, r = fam.value(np.array([1.0]))
print('gauge SectorVarPro value (expect ~0):', f, 'bnorm', fam.bnorm([1.0]), flush=True)
fam_s = SectorVarPro([D8], exclude_gauge=False, name='D8-sector')
fs, xs, rs = fam_s.value(np.array([1.0]))
print('D8* sector f=%.6e rel=%.3e' % (fs, fs / fam_s.bnorm([1.0])**2), flush=True)
kv = KoszulVarPro(6, [CHI0], verbose=True)
a = np.zeros(kv.K)
for k, (i, m) in enumerate(kv.Cbasis):
    if sum(m) == 1 and m[i] == 1: a[k] = 1.0 / 6
fk, xk, rk = kv.value(a)
print('Koszul gauge value (expect ~0):', fk, 'bnorm', kv.bnorm(a), 'excluded_dim', kv.excluded_dim, flush=True)
hr = HessRelax(6, [CHI0], verbose=True)
a = np.zeros(hr.K)
for k, (i, m) in enumerate(hr.Cbasis):
    if sum(m) == 1 and m[i] == 1: a[k] = 1.0 / 6
print('Hess gauge value (expect ~0):', hr.value(a), 'Tnorm', hr.Tnorm, 'excluded', hr.excluded_dim, flush=True)
hr5 = HessRelax(5, ALLCHARS, verbose=True)
res = hr5.search(n_restarts=20, seed=1)
print('Hess d=5 all sectors: best rel', res[0]['rel'], 'worst', res[-1]['rel'], flush=True)
hr6 = HessRelax(6, ALLCHARS, verbose=True)
t1 = time.time()
res = hr6.search(n_restarts=20, seed=1)
print('Hess d=6 all sectors: best rel', res[0]['rel'], 'median', res[len(res)//2]['rel'], 'time/restart', (time.time()-t1)/20, flush=True)
kv7 = KoszulVarPro(7, [CHI0], verbose=True)
t1 = time.time(); res = kv7.search(n_restarts=3, seed=1)
print('Koszul d=7 sector 0: best rel', res[0]['rel'], 'time/restart', (time.time()-t1)/3, 'iters', [c['iters'] for c in res], flush=True)
for d in (10, 12):
    cf = colon_family(d, CHI0, verbose=True)
    basis = cf['jac'] + cf['extras']
    t1 = time.time()
    fam = SectorVarPro(basis, verbose=True, name=f'd{d}')
    t2 = time.time()
    u = np.random.default_rng(0).standard_normal(fam.dimP)
    for _ in range(3):
        fv, g = fam.value_and_grad(u)
    print(f'd={d}: K={fam.K} rows={fam.rows} cols={fam.ncols} build {t2-t1:.1f}s eval {(time.time()-t2)/3:.3f}s f={fv:.3e}', flush=True)
    t1 = time.time(); res = fam.search(n_restarts=2, seed=1)
    print(f'd={d}: 2 restarts best rel {res[0]["rel"]:.3e} iters {[c["iters"] for c in res]} time/restart {(time.time()-t1)/2:.1f}s', flush=True)
print('total', time.time() - t0)
