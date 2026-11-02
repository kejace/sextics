"""k7exact: EXACT rational computation of W_d = {b in B_d : A b in J_{d+6}} for the base-type subspace B_d (Y-even),
without rational reconstruction: pivot rows/cols of the J-generator matrix found mod p give an invertible rational
submatrix; x(b) = S^{-1} (A b)|_pivot rows is exact; Phi(b) = A b - J x(b) on the non-pivot rows is an exact linear map
Q^{dim B} -> Q^{#nonpivot}; W = ker Phi exactly, with E1 lifts N = x(b) for free.  Then J_d cap B exactly (same method
in degree d), and a complement basis = exact non-J generators.  Usage: python3 k7exact.py <d> [--prime P] [--known FILE]
Output: nonj_d<d>.json  [{"label","D","N"}]  (exact; N(F) = A D verified)."""
import sys, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
ap = argparse.ArgumentParser(); ap.add_argument('d', type=int); ap.add_argument('--prime', type=int, default=1000000007); ap.add_argument('--known', default=''); ap.add_argument('--out', default='')
a = ap.parse_args(); d, p = a.d, a.prime
OUT = a.out or f'{PK}/nonj_d{d}.json'
LOGF = f'{PK}/k7exact_d{d}.log.jsonl'
t0 = time.time()
def xmons(k): return [m for m in monomials(k) if m[3] == m[4] == m[5] == 0]
Bpolys = []; Bnames = []
for b_ in range(0, d // 6 + 1):
    for a_ in range(0, (d - 6 * b_) // 2 + 1):
        k = d - 2 * a_ - 6 * b_
        MG = mul(core.power(M, a_), core.power(g, b_))
        for m in xmons(k):
            Bpolys.append(mul(MG, mono(m))); Bnames.append(f'S{k}.M{a_}.g{b_}:' + core.to_str(mono(m)))
nB = len(Bpolys)
print(f'd={d}: base-type subspace dim {nB}', flush=True)

def exact_colspace_test(gens_exact, targets_exact, p):
    """gens_exact: list of exact polys (columns), targets_exact: list of exact polys.
    Returns (Phi matrix as list of dict row->Fraction for each target restricted to non-pivot rows, x-coeffs per target (dict gen index->Fraction))."""
    rows = {}
    for g_ in gens_exact:
        for e in g_: rows.setdefault(e, len(rows))
    for t_ in targets_exact:
        for e in t_: rows.setdefault(e, len(rows))
    R, C = len(rows), len(gens_exact)
    flat = [0] * (R * C)
    for j, g_ in enumerate(gens_exact):
        gp = core.reduce_mod(g_, p)
        for e, v in gp.items(): flat[rows[e] * C + j] = v
    Mx = nmod_mat(R, C, flat, p)
    Mr, rank = Mx.rref(); piv = rref_pivots(Mr, rank, C)
    MTr, r2 = Mx.transpose().rref(); prow = rref_pivots(MTr, r2, R)
    assert r2 == rank
    del Mx, Mr, MTr
    rowpos = {ri: i for i, ri in enumerate(prow)}
    S = fmpq_mat(rank, rank)
    for b_, cj in enumerate(piv):
        for e, v in gens_exact[cj].items():
            ri = rows[e]
            if ri in rowpos:
                v = Fraction(v); S[rowpos[ri], b_] = fmpq(v.numerator, v.denominator)
    T = fmpq_mat(rank, len(targets_exact))
    for tc, t_ in enumerate(targets_exact):
        for e, v in t_.items():
            ri = rows[e]
            if ri in rowpos:
                v = Fraction(v); T[rowpos[ri], tc] = fmpq(v.numerator, v.denominator)
    t1 = time.time()
    X = S.solve(T)
    print(f'   exact solve {rank}x{rank} with {len(targets_exact)} rhs: {time.time()-t1:.0f}s', flush=True)
    nonpiv = [ri for ri in range(R) if ri not in rowpos]
    Phi = []; xs = []
    for tc, t_ in enumerate(targets_exact):
        x = {}
        for b_, cj in enumerate(piv):
            v = X[b_, tc]
            if v != 0: x[cj] = Fraction(int(v.p), int(v.q))
        # residual  t - sum_j x_j g_j  (exact, sparse)
        res = dict((e, Fraction(v)) for e, v in t_.items())
        for cj, xv in x.items():
            for e, v in gens_exact[cj].items():
                res[e] = res.get(e, 0) - xv * Fraction(v)
        res = {e: v for e, v in res.items() if v}
        # residual must vanish on pivot rows by construction; keep the non-pivot part
        for e in res: assert rows[e] not in rowpos, 'residual on a pivot row (bug)'
        Phi.append({rows[e]: v for e, v in res.items()}); xs.append(x)
    return Phi, xs, rank, R

# 1. W = ker Phi  for targets A*b_k against J_{d+6}^- generators
gens = []; glab = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(d + 1, (1 - pi) % 2):
        gens.append(shift(Fi[i], m, 1)); glab.append((i, m))
targets = [mul(A, b) for b in Bpolys]
Phi, xs, rank, R = exact_colspace_test(gens, targets, p)
print(f'J_{d+6}^- rank {rank} of {len(gens)} gens, rows {R}; Phi built ({time.time()-t0:.0f}s)', flush=True)
# kernel of Phi (R-rank rows x nB cols) over Q
rowset = sorted({ri for col in Phi for ri in col}); ridx = {ri: i for i, ri in enumerate(rowset)}
Pm = fmpq_mat(len(rowset), nB)
for j, col in enumerate(Phi):
    for ri, v in col.items(): Pm[ridx[ri], j] = fmpq(v.numerator, v.denominator)
Pr, prank = Pm.rref()
pivc = []
for i in range(prank):
    j = 0
    while j < nB and Pr[i, j] == 0: j += 1
    pivc.append(j)
free = [j for j in range(nB) if j not in set(pivc)]
# kernel basis: for each free column f, vector with 1 at f and -Pr[i,f] at pivot pivc[i]
Wvecs = []
for fcol in free:
    v = [Fraction(0)] * nB; v[fcol] = Fraction(1)
    for i, pc in enumerate(pivc):
        x = Pr[i, fcol]
        if x != 0: v[pc] = -Fraction(int(x.p), int(x.q))
    Wvecs.append(v)
print(f'W = ker Phi: dim {len(Wvecs)} (Phi rank {prank}); {time.time()-t0:.0f}s', flush=True)
# exact E1 lift for each W vector: N = sum_k v_k x_k (linear), verify
def lift_N(v):
    N = [{} for _ in range(NV)]
    for k, vk in enumerate(v):
        if vk == 0: continue
        for cj, xv in xs[k].items():
            i, m = glab[cj]
            N[i][m] = N[i].get(m, 0) + vk * xv
    return [{m: c for m, c in Nk.items() if c} for Nk in N]
def vec_poly(v):
    P_ = {}
    for vk, b in zip(v, Bpolys):
        if vk: P_ = add(P_, scale(b, vk))
    return P_
# 2. J_d cap B exactly: targets b_k against J_d^+ generators
gens7 = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(d - 5, (0 - pi) % 2): gens7.append(shift(Fi[i], m, 1))
Phi7, xs7, rank7, R7 = exact_colspace_test(gens7, Bpolys, p)
rowset7 = sorted({ri for col in Phi7 for ri in col}); ridx7 = {ri: i for i, ri in enumerate(rowset7)}
# J cap B = ker Phi7; we need W modulo (J cap B): evaluate Phi7 on W vectors -> rank of the images = dim W/(W cap J)
imgs = fmpq_mat(len(Wvecs), max(1, len(rowset7)))
for wi, v in enumerate(Wvecs):
    acc = {}
    for k, vk in enumerate(v):
        if vk == 0: continue
        for ri, x in Phi7[k].items(): acc[ri] = acc.get(ri, 0) + vk * x
    for ri, x in acc.items():
        if x: imgs[wi, ridx7[ri]] = fmpq(x.numerator, x.denominator)
# known lower-degree generators times monomials: also map through Phi7 to quotient them out
known_imgs = []
if a.known:
    for rec_ in json.load(open(a.known)):
        Dk = core.parse(rec_['D'])
        for m in monomials(d - core.degree(Dk), 0):
            Pk = mul(Dk, mono(m))
            # express Pk in B coordinates? Not necessarily in B; instead compute its Phi7 image directly
            ph, _, _, _ = exact_colspace_test(gens7, [Pk], p)
            known_imgs.append(ph[0])
Ir, irank = imgs.rref()
print(f'dim W/(W cap J_{d}) = {irank} (W dim {len(Wvecs)}); {time.time()-t0:.0f}s', flush=True)
# choose W vectors whose Phi7-images are independent (greedy, exact) modulo known images
def rank_of(rowsl):
    if not rowsl: return 0
    cols_ = sorted({ri for r_ in rowsl for ri in r_}); cidx = {c: i for i, c in enumerate(cols_)}
    Mq = fmpq_mat(len(rowsl), max(1, len(cols_)))
    for i, r_ in enumerate(rowsl):
        for ri, x in r_.items(): Mq[i, cidx[ri]] = fmpq(x.numerator, x.denominator)
    return Mq.rank()
cur = list(known_imgs); cur_rank = rank_of(cur); chosen = []
for wi, v in enumerate(Wvecs):
    acc = {}
    for k, vk in enumerate(v):
        if vk == 0: continue
        for ri, x in Phi7[k].items(): acc[ri] = acc.get(ri, 0) + vk * x
    acc = {ri: x for ri, x in acc.items() if x}
    if not acc: continue
    rk = rank_of(cur + [acc])
    if rk > cur_rank: cur.append(acc); cur_rank = rk; chosen.append(wi)
print('chosen generators (independent of J and known):', chosen, flush=True)
out = []
for wi in chosen:
    v = Wvecs[wi]; D = vec_poly(v); N = lift_N(v)
    ok = sub(apply_field(N, F), mul(A, D)) == {}
    hD = max(max(abs(Fraction(c).numerator), Fraction(c).denominator) for c in D.values())
    print(f'  generator W[{wi}]: E1 exact verified={ok}; terms {len(D)}; height digits {len(str(hD))}', flush=True)
    if ok: out.append(dict(label=f'K{d}_{wi}', D=core.to_str(D), N=[core.to_str(Nk) for Nk in N]))
json.dump(out, open(OUT, 'w'))
log(dict(exp='nonJ-generators-exact', d=d, prime=p, dimB=nB, dimW=len(Wvecs), dim_W_mod_J=irank, n_new=len(out), file=OUT, seconds=round(time.time() - t0),
         note='exact W = {b in base-type B_d : A b in J} via pivot-submatrix exact solve (no ratrecon); E1 lifts exact; J_d cap B exact'), path=LOGF)
print('wrote', OUT, len(out), 'generators;', round(time.time() - t0), 's')
