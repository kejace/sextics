"""Hessian relaxation on an E1 family with non-Jacobian kernel classes:
   D = sum_k c_k K_k + g(F),   g in V_chi,   K_k the E1-kernel representatives of character chi (degree d).
Obstruction (mod J):  ob = sum_{k,l} c_k c_l B(K_k,K_l) + sum_k c_k L_k(g) + Q(g), with
   B(K_k,K_l) = [B_{K_k}(N_l)] symmetrised (N_l an E1 solution for K_l),   L_k(g) = [A (K_k div g - 2 g(K_k))],
   Q(g) = -[A g^T H g].   All pieces vanish on T = {zG}, so the system lives on (c, W = V/T).
Usage: python3 hessK.py <d> <p> <chi ab> <reps-json>     -> writes hessK_d<d>_chi<ab>_p<p>.ms/.json
"""
import sys, json, time, random
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
from hodge import core
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv
from flint import nmod_mat

d, p = int(sys.argv[1]), int(sys.argv[2]); c = tuple(int(x) for x in sys.argv[3]); reps_file = sys.argv[4]
reps = [core.reduce_mod(core.parse(s), p) for s in json.load(open(reps_file)) if s]
assert all(all(L.chi(e) == c for e in K) for K in reps)
Ap = core.reduce_mod(L.A, p); Fp = core.reduce_mod(L.F, p); Fi = [pderiv(Fp, i, p) for i in range(NV)]
_, H = L.hessian(p)
psi = L.cadd(L.CHI_A, L.cadd(c, c))        # character of every obstruction term (= chi(A))
R = L.CharReducer(2 * d, psi, p)
t0 = time.time()

# E1 solutions N_k for each K_k (character-restricted E1 system, mod p)
def e1_solve(K):
    cols = [(i, m) for i in range(NV) for m in L.mons_chi(d + 1, L.cadd(L.cadd(L.CHI_A, c), L.CHI_Z[i]))]
    rows = L.mons_chi(d + 6, L.cadd(L.CHI_A, c)); ridx = {m: k for k, m in enumerate(rows)}
    AD = mul(Ap, K, p)
    C = len(cols)
    flat = [0] * (len(rows) * (C + 1))
    for j, (i, m) in enumerate(cols):
        for e, v in shift(Fi[i], m, 1, p).items():
            flat[ridx[e] * (C + 1) + j] = v
    for e, v in AD.items():
        flat[ridx[e] * (C + 1) + C] = v
    Mr, rk = nmod_mat(len(rows), C + 1, flat, p).rref()
    N = [dict() for _ in range(NV)]
    for r in range(rk):
        j = 0
        while int(Mr[r, j]) == 0: j += 1
        assert j < C, 'E1 inconsistent for a kernel rep?!'
        v = int(Mr[r, C])
        if v: N[cols[j][0]][cols[j][1]] = v
    acc = {}
    for i in range(NV): acc = add(acc, mul(N[i], Fi[i], p), p)
    assert sub(acc, AD, p) == {}
    return N
Ns = [e1_solve(K) for K in reps]
print(f'E1 solved for {len(reps)} kernel reps ({time.time()-t0:.0f}s)', flush=True)

def BD(D, N):
    divN = {}
    for i in range(NV): divN = add(divN, pderiv(N[i], i, p), p)
    ND = {}
    for i in range(NV): ND = add(ND, mul(N[i], pderiv(D, i, p), p), p)
    return sub(mul(D, divN, p), ND, p)

# W basis
V = L.gbasis(d, c); T = L.trivial_g(d, c)
piv = set()
if T:
    Tm = nmod_mat(len(T), len(V), [t.get(v, 0) for t in T for v in V], p)
    Tr, rk = Tm.rref()
    for r in range(rk):
        j = 0
        while int(Tr[r, j]) == 0: j += 1
        piv.add(j)
W = [V[k] for k in range(len(V)) if k not in piv]
nK, nW = len(reps), len(W)
n = nK + nW
print(f'chi={c}: kernel reps {nK}, dim V={len(V)}, dim T={len(T)}, dim W={nW}, total unknowns {n}', flush=True)

# coefficient vectors
polys, keys = [], []
for k in range(nK):
    for l in range(k, nK):
        P = BD(reps[k], Ns[l]) if k == l else add(BD(reps[k], Ns[l]), BD(reps[l], Ns[k]), p)
        polys.append(P); keys.append((k, l))
for k in range(nK):
    K = reps[k]
    for a, (i, m) in enumerate(W):
        mono = {m: 1}
        P = mul(Ap, sub(mul(K, pderiv(mono, i, p), p), scale(mul(mono, pderiv(K, i, p), p), 2, p), p), p)
        polys.append(P); keys.append((k, nK + a))
for a in range(nW):
    ia, ma = W[a]
    for b in range(a, nW):
        ib, mb = W[b]
        P = shift(mul(Ap, H[(min(ia, ib), max(ia, ib))], p), tuple(x + y for x, y in zip(ma, mb)), p - 1 if a == b else p - 2, p)
        polys.append(P); keys.append((nK + a, nK + b))
# sanity: L_k and B vanish on T
for t in T[:2]:
    gt = [dict() for _ in range(NV)]
    for (i, m), v in t.items(): gt[i][m] = v
    for K in reps[:1]:
        Pt = {}
        for i in range(NV):
            Pt = add(Pt, mul(Ap, sub(mul(K, pderiv(gt[i], i, p), p), scale(mul(gt[i], pderiv(K, i, p), p), 2, p), p), p), p)
        assert not any(R.cls(Pt)), 'L(T) != 0 ?!'
print('sanity L(T)=0 ok', flush=True)
Km = R.cls_many(polys)
quadrics = []
for kk in range(R.nul):
    qd = {}
    for r, key in enumerate(keys):
        v = int(Km[r, kk])
        if v: qd[key] = v
    if qd: quadrics.append(qd)
obK_nonzero = sum(1 for kk in range(R.nul) if int(Km[0, kk]))
Qm = nmod_mat(len(quadrics), len(keys), [qd.get(key, 0) for qd in quadrics for key in keys], p)
rq = Qm.rank()
tag = f'd{d}_chi{c[0]}{c[1]}'
fn = f'{L.OUT}/hessK_{tag}_p{p}.ms'
with open(fn, 'w') as f:
    f.write(','.join([f'c{k}' for k in range(nK)] + [f'w{a}' for a in range(nW)]) + '\n' + f'{p}\n')
    names = [f'c{k}' for k in range(nK)] + [f'w{a}' for a in range(nW)]
    lines = []
    for qd in quadrics:
        lines.append('+'.join(f'{v}*{names[a]}^2' if a == b else f'{v}*{names[a]}*{names[b]}' for (a, b), v in qd.items()))
    f.write(',\n'.join(lines) + '\n')
json.dump({'d': d, 'p': p, 'chis': [list(c)], 'n': n, 'nK': nK, 'nW': nW, 'n_quadrics': len(quadrics), 'rank_quadrics': rq,
           'ob_K0_nonzero_coords': obK_nonzero, 'good_reduction': R.good, 'Wbasis': [[i, list(m)] for (i, m) in W],
           'reps_file': reps_file}, open(fn.replace('.ms', '.json'), 'w'))
L.receipt('hessK', exp='export-kernel-family', d=d, prime=p, chi=c, nK=nK, dim_W=nW, n_quadrics=len(quadrics), rank_quadrics=rq,
          ob_K0_nonzero_coords=obK_nonzero, dim_RJ_2d_psi=R.nul, good_reduction=R.good, file=fn, seconds=round(time.time() - t0, 1))
print('HESSK DONE', flush=True)
