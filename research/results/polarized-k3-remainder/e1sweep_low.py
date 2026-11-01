import sys, time, json; sys.path.insert(0,'/tmp/hodge/results/polarized-k3-remainder')
from uvring import *
f = load_f(); R = TernaryJacobianRing(f, maxdeg=16); T = TensorRing(R, maxdeg=16)
Ah = core.parse(json.load(open('uv_forms.json'))['Ahat'], P)
vA = T.vector_of(Ah, 6)
std6 = R.data[6]['std']
print('z0^2z1^2z2^2 standard in R_f,6:', (2,2,2) in std6, 'dim R_f,6 =', R.dim(6), flush=True)
print('[Ahat] nonzero in R_F,6:', any(vA), '; K3-block components:', [(T.bases[6][i], x) for i, x in enumerate(vA) if x and T.bases[6][i][0] in (0, 6)], flush=True)
res = {}
for d in range(0, 10):
    t = time.time(); M = T.mult_matrix(Ah, 6, d); rank = M.rank()
    res[d] = dict(d=d, dimRd=T.dim(d), dimRd6=T.dim(d+6), rank=rank, kernel=T.dim(d)-rank, seconds=round(time.time()-t,1))
    print(res[d], flush=True)
json.dump(res, open('e1-kernel-sweep-low.json','w'), indent=1)
