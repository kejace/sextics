import sys, time, json; sys.path.insert(0,'.')
from uvring import *
from flint import nmod_mat
f = load_f(); R = TernaryJacobianRing(f, maxdeg=24)
T = TensorRing(R, maxdeg=24)
print('dim R_F by degree:', [T.dim(d) for d in range(25)], flush=True)
# Ahat in (U,V) frame: slots 0..2 = U, 3..5 = V
Ah = core.parse(json.load(open('uv_forms.json'))['Ahat'], P)
vA = T.vector_of(Ah, 6)
print('class of Ahat in R_6 nonzero:', any(vA), 'nnz', sum(1 for x in vA if x), flush=True)
# K3 sector = components (6,0) and (0,6)
k3 = [(i, T.bases[6][i]) for i, x in enumerate(vA) if x and T.bases[6][i][0] in (0, 6)]
print('K3-sector components of [Ahat]:', [(T.R.data[a]['std'][i], T.R.data[b]['std'][j], x) for (k, (a, b, i, j)), x in zip(k3, [vA[k] for k, _ in k3])], flush=True)
res = {}
for d in range(0, 13):
    t = time.time()
    M = T.mult_matrix(Ah, 6, d)
    rank = M.rank()
    n = T.dim(d); m = T.dim(d + 6)
    ker = n - rank
    # symmetric / antisymmetric split of the kernel under the swap
    info = dict(d=d, dimRd=n, dimRd6=m, rank=rank, kernel=ker, seconds=round(time.time() - t, 1))
    if ker:
        X, nullity = M.nullspace()
        Sw = swap_matrix(T, d)
        # kernel basis as columns: X[:, :nullity]; swap acts; compute dims of +/- eigenspaces within kernel
        K = nmod_mat(n, nullity, [int(X[i, j]) for i in range(n) for j in range(nullity)], P)
        SK = Sw * K
        plus = nmod_mat(n, 2 * nullity, [int(K[i, j]) if j < nullity else int((K[i, j - nullity] + SK[i, j - nullity]) % P) for i in range(n) for j in range(2 * nullity)], P)
        # dim of (1+swap)K = rank of SK + K columns
        KpS = nmod_mat(n, nullity, [int((K[i, j] + SK[i, j]) % P) for i in range(n) for j in range(nullity)], P)
        KmS = nmod_mat(n, nullity, [int((K[i, j] - SK[i, j]) % P) for i in range(n) for j in range(nullity)], P)
        info['kernel_sym'] = KpS.rank(); info['kernel_anti'] = KmS.rank()
    res[d] = info
    print(info, flush=True)
json.dump(res, open('e1-kernel-sweep.json', 'w'), indent=1)
