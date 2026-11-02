import sys, time, numpy as np
import w2lib as W
from w2lib import *
from hodge import varpro
from w2blocks import D8, GS, GP, S0 as S0p, P0 as P0p
t0 = time.time()
# 1. gauge direction test (SectorVarPro, no gauge exclusion): D = S0 * F  (degree 8, sector 0)
Dg = mul(S0p, F)
fam = SectorVarPro([Dg], exclude_gauge=False, verbose=True, name='gauge-test')
f, x, r = fam.value(np.array([1.0]))
print('gauge SectorVarPro value (expect ~0):', f, 'bnorm', fam.bnorm([1.0]))
# 2. compare with engine VarProFamily on D8* (single element): same residual?
fam_e = varpro.VarProFamily(F, A, q, [D8], exclude_F_multiples=False, verbose=True)
fe, xe, re_ = fam_e.value(np.array([1.0]))
fam_s = SectorVarPro([D8], exclude_gauge=False, verbose=True, name='D8-sector')
fs, xs, rs = fam_s.value(np.array([1.0]))
print('D8*: engine f=%.6e rel=%.3e | sector f=%.6e rel=%.3e' % (fe, fe / fam_e.bnorm([1.0])**2, fs, fs / fam_s.bnorm([1.0])**2))
# 3. KoszulVarPro gauge test at d=6 sector 0: a = gauge direction (C = Z/6 -> D = F)
kv = KoszulVarPro(6, [CHI0], verbose=True)
# gauge vector: C_i = Z_i/6 -> basis entries (i, e_i)
a = np.zeros(kv.K)
for k, (i, m) in enumerate(kv.Cbasis):
    if sum(m) == 1 and m[i] == 1:
        a[k] = 1.0 / 6
fk, xk, rk = kv.value(a)
print('Koszul gauge value (expect ~0):', fk, 'bnorm', kv.bnorm(a), 'excluded_dim', kv.excluded_dim)
# 4. HessRelax gauge test and d=5 landscape
hr = HessRelax(6, [CHI0], verbose=True)
a = np.zeros(hr.K)
for k, (i, m) in enumerate(hr.Cbasis):
    if sum(m) == 1 and m[i] == 1:
        a[k] = 1.0 / 6
print('Hess gauge value (expect ~0):', hr.value(a), 'Tnorm', hr.Tnorm, 'excluded', hr.excluded_dim)
hr5 = HessRelax(5, ALLCHARS, verbose=True)
res = hr5.search(n_restarts=20, seed=1)
print('Hess d=5 all sectors: best rel', res[0]['rel'], 'worst', res[-1]['rel'])
hr6 = HessRelax(6, ALLCHARS, verbose=True)
t1 = time.time()
res = hr6.search(n_restarts=20, seed=1)
print('Hess d=6 all sectors: best rel', res[0]['rel'], 'median', res[len(res)//2]['rel'], 'time/restart', (time.time()-t1)/20)
# 5. timing SectorVarPro at d=10 and d=12 one sector
for d in (10, 12):
    cf = colon_family(d, CHI0, verbose=True)
    basis = cf['jac'] + cf['extras']
    t1 = time.time()
    fam = SectorVarPro(basis, verbose=True, name=f'd{d}')
    t2 = time.time()
    u = np.random.default_rng(0).standard_normal(fam.dimP)
    for _ in range(3):
        fv, g = fam.value_and_grad(u)
    print(f'd={d}: K={fam.K} rows={fam.rows} cols={fam.ncols} build {t2-t1:.1f}s eval {(time.time()-t2)/3:.3f}s f={fv:.3e}')
print('total', time.time() - t0)
