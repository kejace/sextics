import sys, time, numpy as np
from w2lib import *
for d in (9, 10, 12):
    cf = colon_family(d, CHI0 if d % 2 == 0 else (1,0,0,0), verbose=True)
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
