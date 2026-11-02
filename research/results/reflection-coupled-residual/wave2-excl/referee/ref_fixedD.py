"""Referee fixed-D certificate test (q = 1), written from the contract identity
    E(N,Q) = A D^2 - D sum_i N_i F_i + F (D sum_i d_i N_i - sum_i N_i d_i D) - F^2 Q = 0,
restricted (losslessly, by Klein averaging) to the character sector: N_i in S_{d+1}^{chi_A+chi_D+chi(z_i)}, Q in S_{2d-6}^{chi_A}.
Certification of inconsistency over Q: exact a-priori kernel vectors (zG,0) and (F m e_i, D d_i m - m d_i D) give
nullity_Q >= their rank; a row subset with rank_p = ncols - that rank pins rank_Q(M); a row subset of [M|b] with rank one larger
proves inconsistency over Q (rank_p <= rank_Q, rank_Q([M|b]) <= rank_Q(M)+1).
Usage: python3 ref_fixedD.py <tag> <Dfile.json:key or expr> <p> [nrows_subset]"""
import sys, json, time, random
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee')
import reflib as R
from reflib import NV, monomials, chi, cadd, CHI_Z, CHI_A, shift, mulp, addp
from flint import nmod_mat
tag, src, p = sys.argv[1], sys.argv[2], int(sys.argv[3])
nsub = int(sys.argv[4]) if len(sys.argv) > 4 else 0
OUT = '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee/ref_fixedD.jsonl'
if ':' in src and src.split(':')[0].endswith('.json'):
    fn, key = src.split(':', 1); obj = json.load(open(fn))
    for k in key.split('/'):
        obj = obj[int(k)] if isinstance(obj, list) else obj[k]
    expr = obj
else:
    expr = src
D = R.parse(expr); d = D.degree(); cD = D.chars(); assert len(cD) == 1, cD
cD = cD[0]
S = R.Setup(p); Dp = D.modp(p); Di = [D.deriv(i).modp(p) for i in range(NV)]
t0 = time.time()
# columns
cols = []
for i in range(NV):
    for m in monomials(d + 1, cadd(cadd(CHI_A, cD), CHI_Z[i])): cols.append(('N', i, m))
for m in monomials(2 * d - 6, CHI_A): cols.append(('Q', -1, m))
colidx = {c: j for j, c in enumerate(cols)}
F2 = mulp(S.Fp, S.Fp, p)
colpolys = []
for (kind, i, m) in cols:
    if kind == 'N':
        mono = {m: 1}
        P = {e: (-v) % p for e, v in shift(mulp(Dp, S.Fi[i], p), m, p).items()}          # -D m F_i
        dm = {}
        if m[i]:
            f = list(m); f[i] -= 1; dm = {tuple(f): m[i] % p}
        inner = addp(mulp(Dp, dm, p), mulp(mono, Di[i], p), p, s=-1)                     # D d_i m - m d_i D
        P = addp(P, mulp(S.Fp, inner, p), p)                                              # + F (...)
    else:
        P = {e: (-v) % p for e, v in shift(F2, m, p).items()}                            # -F^2 m
    colpolys.append(P)
b = {e: (-v) % p for e, v in mulp(mulp(S.Ap, Dp, p), Dp, p).items()}                     # RHS: E(N,Q) = -A D^2
rows = {}
for P in colpolys:
    for e in P: rows.setdefault(e, len(rows))
for e in b: rows.setdefault(e, len(rows))
nr, nc = len(rows), len(cols)
print(f'system {nr} x {nc} (d={d}, chi_D={cD}) built in {time.time()-t0:.1f}s', flush=True)
# a-priori exact kernel vectors
vecs = []
for G in monomials(d, cadd(CHI_A, cD)):
    v = {}
    for i in range(NV):
        m = list(G); m[i] += 1; v[colidx[('N', i, tuple(m))]] = 1
    vecs.append(v)
nzG = len(vecs)
if d >= 5:
    for i in range(NV):
        for m in monomials(d - 5, cadd(cadd(CHI_A, cD), CHI_Z[i])):
            v = {}
            for e, c in shift(S.Fp, m, p).items(): v[colidx[('N', i, e)]] = c
            mono = {m: 1}; dm = {}
            if m[i]:
                f = list(m); f[i] -= 1; dm = {tuple(f): m[i] % p}
            Qm = addp(mulp(Dp, dm, p), mulp(mono, Di[i], p), p, s=-1)
            for e, c in Qm.items(): v[colidx[('Q', -1, e)]] = c
            vecs.append(v)
# verify they are exact kernel vectors mod p (M v = 0) using the column polynomials
bad = 0
for v in vecs:
    acc = {}
    for j, c in v.items():
        for e, w in colpolys[j].items(): acc[e] = (acc.get(e, 0) + c * w) % p
    if any(acc.values()): bad += 1
Vm = nmod_mat(len(vecs), nc, [0] * (len(vecs) * nc), p)
for r, v in enumerate(vecs):
    for j, c in v.items(): Vm[r, j] = c % p
r_ap = Vm.rank()
print(f'a-priori vectors {len(vecs)} (zG {nzG}), all in kernel mod p: {bad == 0}, rank mod p = {r_ap}  => rank_Q(M) <= {nc - r_ap}', flush=True)
# row subset ranks
random.seed(7)
allrows = list(range(nr))
if nsub and nsub < nr:
    sub = sorted(random.sample(allrows, nsub))
else:
    sub = allrows
subidx = {r: k for k, r in enumerate(sub)}
flat = [0] * (len(sub) * (nc + 1))
for j, P in enumerate(colpolys):
    for e, w in P.items():
        r = rows[e]
        if r in subidx: flat[subidx[r] * (nc + 1) + j] = w
for e, w in b.items():
    r = rows[e]
    if r in subidx: flat[subidx[r] * (nc + 1) + nc] = w
t1 = time.time()
Maug = nmod_mat(len(sub), nc + 1, flat, p)
rk_aug = Maug.rank()
M0 = nmod_mat(len(sub), nc, [flat[r * (nc + 1) + j] for r in range(len(sub)) for j in range(nc)], p)
rk = M0.rank()
print(f'row subset {len(sub)}/{nr}: rank_p(M_sub) = {rk}, rank_p([M|b]_sub) = {rk_aug}  ({time.time()-t1:.1f}s)', flush=True)
certified = (rk == nc - r_ap) and (rk_aug == rk + 1) and bad == 0
rec = dict(exp='ref-fixedD', tag=tag, d=d, chi_D=list(cD), p=p, rows=nr, cols=nc, apriori_vectors=len(vecs), apriori_in_kernel=(bad == 0),
           apriori_rank_p=r_ap, subset_rows=len(sub), rank_p_subset=rk, rank_p_aug_subset=rk_aug, inconsistent_modp=(rk_aug == rk + 1),
           certified_inconsistent_over_Q=certified, seconds=round(time.time() - t0, 1),
           note='certified: rank_Q(M) = ncols - apriori_rank and rank_Q([M|b]) = rank_Q(M)+1, so no (N,Q) over Q or C' if certified else 'NOT certified (increase subset or check)')
R.receipt(OUT, rec)
