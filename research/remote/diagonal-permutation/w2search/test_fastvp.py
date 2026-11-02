import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation/w2search'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
from hodge import core, jacobian
import fastvp
iid = sys.argv[1] if len(sys.argv) > 1 else 'contact46006'
F, A, q, rec = core.load_instance(iid)
_w = core.load_witness(iid); Dw, Nw, Qw = _w[1], _w[2], _w[3]
d = core.degree(Dw); print('instance', iid, 'q', q, 'd', d, 'parD', core.parity_of(Dw), 'parA', core.parity_of(A))
assert q == 1
# E1 particular solution for the witness D (float) and a couple of Jacobian family members
t0 = time.time()
labels, Bs, N0s = fastvp.jacobian_family(F, A, d, core.parity_of(Dw))
print('jacobian family size', len(Bs), 'built', round(time.time()-t0,1))
N0w, res = fastvp.e1_particular(F, A, [Dw], verbose=True)
print('e1 particular residual', res)
rng = np.random.default_rng(0)
pick = rng.choice(len(Bs), size=min(5, len(Bs)), replace=False)
fam = fastvp.FastFamily(F, A, [Dw] + [Bs[i] for i in pick], N0w + [N0s[i] for i in pick], verbose=True, seed=3)
a = np.zeros(fam.K); a[0] = 1
f, phin = fam.value(a)
print('witness direction: f=%.3e phi2=%.3e rel=%.3e' % (f, phin, f / phin))
a2 = rng.standard_normal(fam.K)
f2, phin2 = fam.value(a2)
print('random direction: f=%.3e phi2=%.3e rel=%.3e' % (f2, phin2, f2 / phin2))
print('fd check (f, g.v, fd):', fam.fd_check())
t0 = time.time()
cands = fam.search(n_restarts=4, seed=1, verbose=True, maxiter=300)
print('best rel', cands[0]['rel'], 'a', np.round(cands[0]['a'], 4), 'time', round(time.time()-t0,1))
for c in cands[:2]:
    for fr in fam.reconstruct(c['a']):
        print('reconstructed', fr)
        st, data = fam.confirm(fr)
        print('confirm ->', st, {k: v for k, v in data.items() if k in ('rank','nullity','note','verify')})
        break
