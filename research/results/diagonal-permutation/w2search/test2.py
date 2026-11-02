import sys, time
sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation/w2search'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
from hodge import core
import fastvp
F, A, q, rec = core.load_instance('contact46006')
_w = core.load_witness('contact46006'); Dw = _w[1]
labels, Bs, N0s = fastvp.jacobian_family(F, A, 6, 0)
N0w, res = fastvp.e1_particular(F, A, [Dw], verbose=True)
rng = np.random.default_rng(0); pick = rng.choice(len(Bs), size=5, replace=False)
fam = fastvp.FastFamily(F, A, [Dw] + [Bs[i] for i in pick], N0w + [N0s[i] for i in pick], verbose=True, seed=3)
print('fd check', fam.fd_check())
t0 = time.time()
cands = fam.search(n_restarts=3, seed=1, verbose=True, maxiter=400)
c = cands[0]
print('best rel', c['rel'], 'rel_lbfgs', c['rel_lbfgs'], 'a', np.array2string(np.array(c['a']), precision=12), 'time', round(time.time()-t0,1), c['msg'])
for fr in fam.reconstruct(np.array(c['a'])):
    print('reconstructed', fr); st, data = fam.confirm(fr); print('confirm ->', st, {k: v for k, v in data.items() if k in ('rank','nullity','note','verify')}); break
