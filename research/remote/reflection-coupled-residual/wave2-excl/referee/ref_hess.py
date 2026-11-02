"""Referee: Hessian relaxation of E2 on the Jacobian family D = g(F) = sum_i g_i F_i (q = 1), written from scratch.
ob(g) = -[A * sum_{i,j} g_i g_j F_ij] in (S/J)_{2d}; B(zG, .) = 0 so ob lives on W = V/T.
Modes:
  rank20 <p>                     degree 5: the classes [A F_ij] (F_ij != 0) are independent in (S/J)_10  => g = 0.
  export <d> <p> <chis...>       quadrics on W = (+)_chi W_chi (chis like 00 01 10 11), msolve input ref_hess_d<d>_W<chis>_p<p>.ms
"""
import sys, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee')
import reflib as R
from reflib import NV, monomials, chi, cadd, CHI_Z, CHI_A, shift, mulp, addp
from flint import nmod_mat
OUT = '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee'
mode = sys.argv[1]

def gbasis(d, c):
    """basis of V_c: tuples (i, m), m monomial of degree d-5 with chi(m)+chi(z_i) = c (so that m*F_i has character c)."""
    return [(i, m) for i in range(NV) for m in monomials(d - 5, cadd(c, CHI_Z[i]))]

def trivial(d, c):
    """T_c = {z*G : G in S_{d-6}^c} as vectors over gbasis(d, c)."""
    out = []
    for G in monomials(d - 6, c) if d >= 6 else []:
        v = {}
        for i in range(NV):
            m = list(G); m[i] += 1; v[(i, tuple(m))] = 1
        out.append(v)
    return out

if mode == 'rank20':
    p = int(sys.argv[2]); S = R.Setup(p)
    pairs = [(i, j) for i in range(NV) for j in range(i, NV) if S.Fij[(i, j)]]
    zero_pairs = [(i, j) for i in range(NV) for j in range(i, NV) if not S.Fij[(i, j)]]
    polys = [mulp(S.Ap, S.Fij[(i, j)], p) for (i, j) in pairs]
    bychar = {}
    for P, pr in zip(polys, pairs): bychar.setdefault(chi(next(iter(P))), []).append((P, pr))
    total_rank = 0; good = True; info = {}
    for c, lst in bychar.items():
        Rd = R.Reducer(S, 10, c); good &= Rd.good
        rk = Rd.rank_mod_J([P for P, _ in lst]); total_rank += rk
        info[str(c)] = dict(n=len(lst), rank_mod_J=rk, good=Rd.good)
    R.receipt(f'{OUT}/ref_hess.jsonl', dict(exp='ref-rank20', d=5, p=p, nonzero_pairs=len(pairs), zero_pairs=[[i, j] for i, j in zero_pairs],
              total_rank_mod_J=total_rank, independent=(total_rank == len(pairs)), good_reduction_deg10=good, blocks=info))
    print('rank20 verdict: independent =', total_rank == len(pairs), 'zero pairs', zero_pairs)

elif mode == 'export':
    d = int(sys.argv[2]); p = int(sys.argv[3]); chis = [tuple(int(x) for x in s) for s in sys.argv[4:]]
    S = R.Setup(p); t0 = time.time()
    Wb = []; dims = {}
    for c in chis:
        V = gbasis(d, c); vidx = {v: k for k, v in enumerate(V)}
        T = trivial(d, c); piv = set()
        if T:
            Tm = nmod_mat(len(T), len(V), [t.get(v, 0) for t in T for v in V], p)
            Tr, rk = Tm.rref(); assert rk == len(T)
            for r in range(rk):
                j = 0
                while int(Tr[r, j]) == 0: j += 1
                piv.add(j)
        dims[str(c)] = (len(V), len(T), len(V) - len(T))
        Wb += [(c, V[k][0], V[k][1]) for k in range(len(V)) if k not in piv]
    n = len(Wb); print('dims (V,T,W) per chi:', dims, 'total W', n, flush=True)
    # bilinear polys P_ab = A * (g^a_i g^b_j + g^b_i g^a_j)/2 summed ... as monomial pairs: for basis elements a=(ca,i,m), b=(cb,j,m'):
    # contribution to ob: -A*(m*m')*F_ij * w_a w_b  (ordered pairs; we take a<=b and weight 2 for a<b)
    wpairs = [(a, b) for a in range(n) for b in range(a, n)]
    # group by target character psi = chi_A + chi(m m' F_ij)
    groups = {}
    for (a, b) in wpairs:
        (ca, i, m), (cb, j, m2) = Wb[a], Wb[b]
        Fij = S.Fij[(i, j)]
        if not Fij: continue
        mm = tuple(x + y for x, y in zip(m, m2))
        P = shift(mulp(S.Ap, Fij, p), mm, p)
        if a != b: P = {e: (2 * v) % p for e, v in P.items()}
        P = {e: (-v) % p for e, v in P.items()}
        psi = chi(next(iter(P)))
        groups.setdefault(psi, []).append(((a, b), P))
    quad_rows = []   # each quadric: dict (a,b)->coeff
    good = True; nq_total = 0
    for psi, lst in sorted(groups.items()):
        Rd = R.Reducer(S, 2 * d, psi); good &= Rd.good
        C = Rd.coords([P for _, P in lst])      # len(lst) x dimRJ
        nst = len(Rd.standard)
        for t in range(nst):
            q = {}
            for k, ((a, b), _) in enumerate(lst):
                v = C[k][t]
                if v: q[(a, b)] = v
            if q: quad_rows.append(q)
        nq_total += nst
        print(f'target psi={psi}: pairs {len(lst)} coords {nst} good_reduction={Rd.good}', flush=True)
    # rank of the quadrics
    pidx = {pr: k for k, pr in enumerate(wpairs)}
    M = nmod_mat(len(quad_rows), len(wpairs), [0] * (len(quad_rows) * len(wpairs)), p)
    for r, q in enumerate(quad_rows):
        for pr, v in q.items(): M[r, pidx[pr]] = v
    rank_q = M.rank()
    tag = ''.join(f'{c[0]}{c[1]}' for c in chis)
    fn = f'{OUT}/ref_hess_d{d}_W{tag}_p{p}.ms'
    names = [f'w{k}' for k in range(n)]
    with open(fn, 'w') as f:
        f.write(','.join(names) + '\n' + str(p) + '\n')
        terms = []
        for q in quad_rows:
            s = '+'.join(f'{v}*w{a}*w{b}' if a != b else f'{v}*w{a}^2' for (a, b), v in sorted(q.items()))
            terms.append(s)
        f.write(',\n'.join(terms) + '\n')
    meta = dict(exp='ref-hess-export', d=d, p=p, chis=[list(c) for c in chis], dims=dims, dim_W=n, n_quadrics=len(quad_rows), rank_quadrics=rank_q,
                good_reduction=good, file=fn, seconds=round(time.time() - t0, 1))
    json.dump(meta, open(fn.replace('.ms', '.json'), 'w'))
    R.receipt(f'{OUT}/ref_hess.jsonl', meta)
