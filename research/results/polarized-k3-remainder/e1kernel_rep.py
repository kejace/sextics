"""Extract sparse swap-symmetric representatives of the E1 kernel at degree d (mod p), lifted to small
integer polynomials when possible: we look for kernel elements supported on few standard-monomial pairs by
rref with a column order that favours sparsity, then rationally reconstruct."""
import sys, time, json; sys.path.insert(0,'.')
from uvring import *
from flint import nmod_mat
d = int(sys.argv[1]); nrep = int(sys.argv[2]) if len(sys.argv) > 2 else 3
f = load_f(); R = TernaryJacobianRing(f, maxdeg=24); T = TensorRing(R, maxdeg=24)
Ah = core.parse(json.load(open('uv_forms.json'))['Ahat'], P)
M = T.mult_matrix(Ah, 6, d)
n = T.dim(d)
# symmetric subspace: we want kernel vectors v with swap(v) = v. Solve [M; Sw - I] v = 0.
Sw = swap_matrix(T, d)
rowsM = M.nrows()
flat = [int(M[i, j]) for i in range(rowsM) for j in range(n)] + [((int(Sw[i, j]) - (1 if i == j else 0)) % P) for i in range(n) for j in range(n)]
Big = nmod_mat(rowsM + n, n, flat, P)
X, nullity = Big.nullspace()
print('symmetric E1-kernel dim at d =', d, ':', nullity, flush=True)
# rref of the transposed kernel basis gives reduced vectors; count supports
K = nmod_mat(nullity, n, [int(X[i, j]) for j in range(nullity) for i in range(n)], P)
Kr, rk = K.rref()
reps = []
for i in range(rk):
    supp = [j for j in range(n) if int(Kr[i, j])]
    reps.append((len(supp), i))
reps.sort()
out = []
for size, i in reps[:nrep]:
    poly = {}
    for j in range(n):
        v = int(Kr[i, j])
        if v:
            a, b, ii, jj = T.bases[d][j]
            mu = R.data[a]['std'][ii]; mv = R.data[b]['std'][jj]
            r = core.ratrecon(v, P)
            poly[(mu[0], mu[1], mu[2], mv[0], mv[1], mv[2])] = r if r is not None else v
    ok = all(not isinstance(c, int) for c in poly.values())
    s = core.to_str({e: c for e, c in poly.items()}) if ok else None
    out.append({'support': size, 'rational_lift_ok': ok, 'expr_uv': s})
    print(size, ok, (s or '')[:300], flush=True)
json.dump(out, open(f'e1kernel-reps-d{d}.json', 'w'), indent=1)
