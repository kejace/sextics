"""Hessian relaxation of E2 on the Jacobian family D = g(F) = sum g_i F_i  (q = 1).

For D = g(F) and the natural E1 solution N = A g, the E2-mod-J obstruction is
    ob(g) = [ A * sum_{i,j} g_i g_j F_ij ]  in (S/J)_{2d}      (F_ij = Hessian of F),
and ob(g) = 0 is necessary for a certificate with denominator g(F).  B(t, .) = 0 for t in T = {z G}.

modes:
  phi  <d> <p>          d=5,6: linear map Phi on symmetric tuples h_ij in S_{2d-10}; kernel per target
                        character vs the trivial kernel {(z_i z_j Q)}.  ker Phi = trivial  =>  only D in (F).
  mac3 <d> <p> <chi>    quadratic system on W_chi = V_chi / T_chi, Macaulay matrix in degree 3 (random row
                        compression); full column rank mod p  =>  I_Q contains all cubics  =>  V(I)={0}.
  export <d> <p> <chis> write msolve input for the quadratic system on the direct sum of W_chi (chis like 00,01).
"""
import sys, json, time, itertools, random
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
from hodge import core
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv
from flint import nmod_mat
import numpy as np
import scipy.sparse as sp

mode = sys.argv[1]; d = int(sys.argv[2]); p = int(sys.argv[3])
Ap = core.reduce_mod(L.A, p)
Fi, H = L.hessian(p)
pairs = [(i, j) for i in range(NV) for j in range(i, NV) if H[(i, j)]]
print('nonzero Hessian pairs:', len(pairs), 'zero pairs:', [(i, j) for i in range(NV) for j in range(i, NV) if not H[(i, j)]], flush=True)


def compress_rank(rows_sparse, ncols, k_extra=16, seed=1):
    """rows_sparse: list of dict col->val (mod p). Returns (rank of random compression, k)."""
    R = len(rows_sparse)
    ri, ci, vals = [], [], []
    for r, row in enumerate(rows_sparse):
        for c, v in row.items():
            ri.append(r); ci.append(c); vals.append(v % p)
    M = sp.csr_matrix((np.array(vals, dtype=np.int64), (np.array(ri), np.array(ci))), shape=(R, ncols))
    k = min(R, ncols + k_extra)
    rng = np.random.default_rng(seed)
    Rm = rng.integers(0, p, size=(k, R), dtype=np.int64)
    MT = M.T.tocsr()
    lo = (MT @ (Rm & 0xFFFF).T) % p
    hi = (MT @ (Rm >> 16).T) % p
    Cm = np.ascontiguousarray(((lo + (hi * (1 << 16)) % p) % p).T)
    t0 = time.time()
    rk = nmod_mat(k, ncols, Cm.ravel().tolist(), p).rank()
    return rk, k, round(time.time() - t0, 1)


if mode == 'phi':
    n = 2 * d
    hm = monomials(2 * d - 10)
    cols = [(i, j, m) for (i, j) in pairs for m in hm]
    blocks = {}
    for idx, (i, j, m) in enumerate(cols):
        psi = L.cadd(L.cadd(L.CHI_A, L.cadd(L.CHI_Z[i], L.CHI_Z[j])), L.chi(m))
        blocks.setdefault(psi, []).append(idx)
    verdict = True
    total_ker, total_triv = 0, 0
    for psi, idxs in sorted(blocks.items()):
        t0 = time.time()
        R = L.CharReducer(n, psi, p)
        polys = []
        for idx in idxs:
            i, j, m = cols[idx]
            polys.append(scale(shift(mul(Ap, H[(i, j)], p), m, 1, p), 1 if i == j else 2, p))
        K = R.cls_many(polys)                     # len(idxs) x nul
        X, nullity = K.transpose().nullspace()    # kernel vectors over the block columns
        # trivial kernel in this block: Q in S_{2d-12} with chi(A)+chi(Q) == psi
        triv = []
        if 2 * d - 12 >= 0:
            pos = {cols[idx]: t for t, idx in enumerate(idxs)}
            for Q in monomials(2 * d - 12):
                if L.cadd(L.CHI_A, L.chi(Q)) != psi:
                    continue
                v = {}
                for (i, j) in pairs:
                    m = list(Q); m[i] += 1; m[j] += 1; m = tuple(m)
                    v[pos[(i, j, m)]] = 1
                triv.append(v)
        # sanity: trivial vectors are in the kernel
        triv_ok = True
        if triv:
            T = nmod_mat(len(triv), len(idxs), [v.get(t, 0) for v in triv for t in range(len(idxs))], p)
            TK = T * K
            triv_ok = all(int(TK[a, b]) == 0 for a in range(len(triv)) for b in range(R.nul))
        ok = (nullity == len(triv)) and triv_ok
        verdict &= ok
        total_ker += nullity; total_triv += len(triv)
        L.receipt('hess', exp='phi-kernel', d=d, prime=p, psi=psi, n_cols=len(idxs), dim_S_2d=R.dimS, rankJ_2d=R.rank,
                  fermat_2d=R.fermat, good_reduction=R.good, dim_RJ_2d=R.nul, rank_phi=len(idxs) - nullity,
                  ker_phi=nullity, trivial_ker=len(triv), trivial_in_kernel=triv_ok, kernel_is_trivial=ok,
                  seconds=round(time.time() - t0, 1))
    L.receipt('hess', exp='phi-verdict', d=d, prime=p, nonzero_pairs=len(pairs), n_cols=len(cols), ker_phi_total=total_ker,
              trivial_total=total_triv, kernel_is_trivial=verdict,
              note='ker Phi = trivial kernel => every g with ob(g)=0 has g_i g_j = z_i z_j Q for all pairs with F_ij != 0 '
                   '=> g = +-G z => D = g(F) in (F): no admissible degree-%d denominator in J_%d' % (d, d) if verdict else
                   'kernel larger than trivial: quadratic system needed')

elif mode in ('mac3', 'export'):
    chis = [tuple(int(x) for x in s) for s in sys.argv[4:]]
    # W basis per character: V_chi basis minus pivots of T_chi (T = {zG}); ob factors through V/T since B(T,.)=0
    Wbasis = []   # list of (chi, i, m)
    Tinfo = {}
    for c in chis:
        V = L.gbasis(d, c)
        vidx = {v: k for k, v in enumerate(V)}
        T = L.trivial_g(d, c)
        piv = set()
        if T:
            Tm = nmod_mat(len(T), len(V), [t.get(v, 0) for t in T for v in V], p)
            Tr, rk = Tm.rref()
            assert rk == len(T)
            for r in range(rk):
                j = 0
                while int(Tr[r, j]) == 0: j += 1
                piv.add(j)
        Tinfo[c] = (len(V), len(T))
        Wbasis += [(c,) + V[k] for k in range(len(V)) if k not in piv]
        print(f'chi={c}: dim V={len(V)} dim T={len(T)} dim W={len(V)-len(T)}', flush=True)
    n = len(Wbasis)
    print('total W dim', n, flush=True)
    # sanity check B(t, e) = 0 for trivial t, random e (a few); reducers cached by target character
    reducers = {}
    def getred(psi):
        if psi not in reducers:
            reducers[psi] = L.CharReducer(2 * d, psi, p)
        return reducers[psi]
    random.seed(1)
    for c in chis:
        for t in L.trivial_g(d, c)[:3]:
            gt = [dict() for _ in range(NV)]
            for (i, m), v in t.items(): gt[i][m] = v
            for _ in range(2):
                (cc, i, m) = random.choice(Wbasis)
                ge = [dict() for _ in range(NV)]; ge[i][m] = 1
                P = L.bil_poly(gt, ge, H, Ap, p)
                psi = L.chi(next(iter(P))) if P else None
                if P:
                    R = getred(psi)
                    assert not any(R.cls(P)), 'B(T, .) != 0 ?!'
    print('sanity B(T,.)=0 ok', flush=True)
    # quadric coefficient vectors by target character
    pair_list = [(a, b) for a in range(n) for b in range(a, n)]
    bypsi = {}
    for (a, b) in pair_list:
        ca, ia, ma = Wbasis[a]; cb, ib, mb = Wbasis[b]
        psi = L.cadd(L.cadd(L.CHI_A, ca), cb)
        bypsi.setdefault(psi, []).append((a, b))
    quadrics = []   # list of dict {(a,b): coeff}  (a<=b), coefficient of w_a w_b
    for psi, plist in sorted(bypsi.items()):
        t0 = time.time()
        R = getred(psi)
        polys = []
        for (a, b) in plist:
            ca, ia, ma = Wbasis[a]; cb, ib, mb = Wbasis[b]
            P = shift(mul(Ap, H[(min(ia, ib), max(ia, ib))], p), tuple(x + y for x, y in zip(ma, mb)), 1 if a == b else 2, p)
            polys.append(P)
        K = R.cls_many(polys)      # len(plist) x nul
        Kl = [[int(K[r, k]) for k in range(R.nul)] for r in range(len(plist))]
        for k in range(R.nul):
            qd = {}
            for r, (a, b) in enumerate(plist):
                v = Kl[r][k]
                if v: qd[(a, b)] = v
            if qd:
                quadrics.append(qd)
        print(f'psi={psi}: pairs={len(plist)} nul={R.nul} good={R.good} quadrics so far {len(quadrics)} ({time.time()-t0:.0f}s)', flush=True)
    # drop linearly dependent quadrics (rank of coefficient matrix)
    pidx = {pr: k for k, pr in enumerate(pair_list)}
    Qm = nmod_mat(len(quadrics), len(pair_list), [qd.get(pr, 0) for qd in quadrics for pr in pair_list], p)
    rq = Qm.rank()
    print(f'quadrics: {len(quadrics)}, independent: {rq}, in {n} variables', flush=True)
    if mode == 'mac3':
        c = chis
        cub = list(itertools.combinations_with_replacement(range(n), 3))
        cidx = {t: k for k, t in enumerate(cub)}
        rows = []
        for qd in quadrics:
            for cvar in range(n):
                row = {}
                for (a, b), v in qd.items():
                    key = tuple(sorted((a, b, cvar)))
                    row[cidx[key]] = (row.get(cidx[key], 0) + v) % p
                rows.append(row)
        print(f'Macaulay-3: rows={len(rows)} cols={len(cub)}', flush=True)
        t0 = time.time()
        # first try: random row subset of size cols+64 (rank of a subset is a lower bound for the rank)
        random.seed(7)
        sub_rows = random.sample(range(len(rows)), min(len(rows), len(cub) + 64))
        flat = [0] * (len(sub_rows) * len(cub))
        for r, ri_ in enumerate(sub_rows):
            base = r * len(cub)
            for cc, v in rows[ri_].items():
                flat[base + cc] = v
        rk = nmod_mat(len(sub_rows), len(cub), flat, p).rank(); k = len(sub_rows); method = 'row-subset'
        del flat
        if rk < len(cub) and len(rows) > len(sub_rows):
            rk, k, trk = compress_rank(rows, len(cub)); method = 'random-compression'
        full = (rk == len(cub))
        L.receipt('hess', exp='macaulay3', d=d, prime=p, chi=c, dim_V=sum(v[0] for v in Tinfo.values()), dim_T=sum(v[1] for v in Tinfo.values()), dim_W=n, method=method,
                  n_quadrics=len(quadrics), rank_quadrics=rq, mac_rows=len(rows), mac_cols=len(cub), compressed_rows=k,
                  rank=rk, full_column_rank=full, good_reduction=all(R.good for R in reducers.values()),
                  seconds=round(time.time() - t0, 1),
                  note=('I_3 = S_3(W) mod p => over Q => V(ob)=T: no degree-%d denominator with character components in %s besides (F)' % (d, c)) if full
                       else 'degree-3 Macaulay matrix not of full column rank; need higher degree / GB')
    else:
        tagc = ''.join(f'{c[0]}{c[1]}' for c in chis)
        fn = f'{L.OUT}/hess_d{d}_W{tagc}_p{p}.ms'
        with open(fn, 'w') as f:
            f.write(','.join(f'w{k}' for k in range(n)) + '\n')
            f.write(f'{p}\n')
            lines = []
            for qd in quadrics:
                terms = []
                for (a, b), v in qd.items():
                    mon = f'w{a}^2' if a == b else f'w{a}*w{b}'
                    terms.append(f'{v}*{mon}')
                lines.append('+'.join(terms))
            f.write(',\n'.join(lines) + '\n')
        json.dump({'d': d, 'p': p, 'chis': chis, 'Wbasis': [[list(c), i, list(m)] for (c, i, m) in Wbasis], 'n': n,
                   'n_quadrics': len(quadrics), 'rank_quadrics': rq, 'Tinfo': {str(k): v for k, v in Tinfo.items()},
                   'good_reduction': all(R.good for R in reducers.values())},
                  open(fn.replace('.ms', '.json'), 'w'))
        L.receipt('hess', exp='export-msolve', d=d, prime=p, chis=chis, dim_W=n, n_quadrics=len(quadrics), rank_quadrics=rq, file=fn)
print('HESS DONE', flush=True)
