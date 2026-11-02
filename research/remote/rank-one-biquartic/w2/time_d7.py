import sys, time
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/w2'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
from w2lib import *
from hodge import core
F, A, q, rec = core.load_instance('rank-one-biquartic')
for d in (7,):
    t0 = time.time(); S = ReducedE2(F, A, d, 0, verbose=True); print('build', time.time()-t0, 'rf_maxabs', S.rf_maxabs, 'RF nnz', S.RF.nnz)
    t0 = time.time(); fam = full_family(F, A, d, 0, verbose=True); print('family', time.time()-t0, fam.K, fam.info)
    print('E1 check', check_family_e1(F, A, fam))
    famJ = j_family(F, A, d, 0, S.Dm, S.Ncols); print('J family K', famJ.K)
    rng = np.random.default_rng(0); u = rng.standard_normal(fam.K)
    for kz in (False, True):
        t0 = time.time(); f, g = S.value_and_grad(u, fam, koszul=kz); print(f'eval koszul={kz}: f={f:.3e} {time.time()-t0:.2f}s')
    t0 = time.time(); f, rel, info = S.value(u, fam); print('rel random', rel, time.time()-t0)
    import resource; print('maxrss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
