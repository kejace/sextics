import sys, time, json; sys.path.insert(0,'.')
from uvring import *
from flint import nmod_mat
f = load_f(); R = TernaryJacobianRing(f, maxdeg=24)
T = TensorRing(R, maxdeg=24)
Ah = core.parse(json.load(open('uv_forms.json'))['Ahat'], P)
res = json.load(open('e1-kernel-sweep.json')) if __import__('os').path.exists('e1-kernel-sweep.json') else {}
for d in [int(a) for a in sys.argv[1:]]:
    t = time.time()
    M = T.mult_matrix(Ah, 6, d)
    n = T.dim(d); m = T.dim(d + 6)
    X, nullity = M.nullspace()
    rank = n - nullity
    info = dict(d=d, dimRd=n, dimRd6=m, rank=rank, kernel=nullity, seconds=round(time.time() - t, 1))
    if nullity:
        Sw = swap_matrix(T, d)
        K = nmod_mat(n, nullity, [int(X[i, j]) for i in range(n) for j in range(nullity)], P)
        SK = Sw * K
        KpS = nmod_mat(n, nullity, [(int(K[i, j]) + int(SK[i, j])) % P for i in range(n) for j in range(nullity)], P)
        KmS = nmod_mat(n, nullity, [(int(K[i, j]) - int(SK[i, j])) % P for i in range(n) for j in range(nullity)], P)
        info['kernel_sym'] = KpS.rank(); info['kernel_anti'] = KmS.rank()
        # how many kernel elements are 'pure' in a single bidegree (a,b)? (bidegree structure of kernel)
        info['seconds'] = round(time.time() - t, 1)
    res[str(d)] = info
    print(info, flush=True)
    json.dump(res, open('e1-kernel-sweep.json', 'w'), indent=1)
