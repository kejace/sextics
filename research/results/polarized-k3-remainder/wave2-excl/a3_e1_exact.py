"""(a) re-verification 3: exact-over-Q ranks of multiplication by [Ahat]: R_d -> R_{d+6} on R_F = R_f (x) R_f,
with the swap-symmetric / antisymmetric split (swap = block swap U<->V = tau).  No primes involved."""
import sys, time, json
sys.path.insert(0, '.')
from k3lib import *
from logutil import log
f = ternary_f(); R = TernaryRingQ(f, maxdeg=18); T = TensorRingQ(R, maxdeg=18)
Fh, Ah, q = load_uv()
Ah64 = scale(Ah, 64)   # integer coefficients
res = {}
for d in range(0, 13):
    t0 = time.time()
    cols = T.mult_columns(Ah64, 6, d)
    n, m = T.dim(d), T.dim(d + 6)
    M = sparse_cols_to_fmpz(cols, m)
    rank = M.rank()
    # swap action on basis of degree d: (a,b,i,j) -> (b,a,j,i); symmetric subspace basis: e_t + e_{s(t)}, anti: e_t - e_{s(t)}
    perm = [T.index[d][(b, a, j, i)] for (a, b, i, j) in T.bases[d]]
    def rank_on(sign):
        seen = set(); vecs = []
        for t in range(n):
            s = perm[t]
            if t in seen or s in seen: continue
            seen.add(t); seen.add(s)
            if s == t:
                if sign == 1: vecs.append({t: 1})
            else:
                vecs.append({t: 1, s: sign})
        # image vectors = M v
        imgs = []
        for v in vecs:
            col = {}
            for t, c in v.items():
                for r, x in cols[t].items():
                    col[r] = col.get(r, 0) + c * x
            imgs.append({r: x for r, x in col.items() if x})
        return len(vecs), sparse_cols_to_fmpz(imgs, m).rank()
    ns, rs = rank_on(1); na, ra = rank_on(-1)
    res[d] = dict(d=d, dimRd=n, dimRd6=m, rank=rank, kernel=n - rank, sym_dim=ns, sym_rank=rs, sym_kernel=ns - rs,
                  anti_dim=na, anti_rank=ra, anti_kernel=na - ra, seconds=round(time.time() - t0, 1))
    print(res[d], flush=True)
    log(exp='a3_e1_exact', rigor='exact over Q (fmpz_mat.rank), UV tensor ring', **res[d])
json.dump(res, open('a3_e1_exact.json', 'w'), indent=1)
