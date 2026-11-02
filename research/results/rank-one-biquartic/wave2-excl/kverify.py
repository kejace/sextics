"""kverify: turn reconstructed rational kernel vectors (klift/k7lift output) into exact (D, N) pairs with N(F) = A D over Q.
Usage: python3 kverify.py <rational.json> <d> <out nonj json> [--known nonj_lower.json]
  * decides which rref basis vectors of W are in J_d (exact: rank test mod p for non-membership + exact polar solve for membership)
  * picks a maximal set independent modulo J_d + (products of lower-degree known generators by linear forms) [mod p]
  * exact E1 lift: pivots mod p, fmpq solve of the pivot submatrix, exact verification N(F) = A D."""
import sys, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
ap = argparse.ArgumentParser(); ap.add_argument('ratjson'); ap.add_argument('d', type=int); ap.add_argument('out'); ap.add_argument('--known', default=''); ap.add_argument('--prime', type=int, default=1000000007)
a = ap.parse_args(); d = a.d; p = a.prime
LOGF = f'{PK}/kverify_d{d}.log.jsonl'
J = json.load(open(a.ratjson))
coords = J['coords']
def coord_poly(s):
    kind, expr = s.split(':', 1)
    return core.parse(expr)
cpolys = [coord_poly(s) for s in coords]
def vec_to_poly(v):
    P_ = {}
    for x, cp in zip(v, cpolys):
        x = Fraction(x)
        if x: P_ = add(P_, scale(cp, x))
    return P_
Wv = [vec_to_poly(v) for v in J['W']]
print(f'{len(Wv)} rational vectors of degree {d}', flush=True)
t0 = time.time()
# J-membership mod p and independence modulo J + known*linear
parD = core.parity_of(Wv[0])
red = jacobian._Reducer(F, d, p, parD)
known = []
if a.known:
    for rec_ in json.load(open(a.known)):
        Dk = core.parse(rec_['D'])
        for m in monomials(d - core.degree(Dk)):
            known.append(mul(Dk, {m: Fraction(1)}))
sm = red.standard; sidx = {m: i for i, m in enumerate(sm)}
def nfvec(P_):
    v = red.normal_form(core.reduce_mod(P_, p)); return [v.get(m, 0) for m in sm]
rows_known = [nfvec(x) for x in known]
base_rank = nmod_mat(len(rows_known), len(sm), [x for r_ in rows_known for x in r_], p).rank() if rows_known else 0
chosen = []
cur = list(rows_known); cur_rank = base_rank
for i, Dv in enumerate(Wv):
    nv = nfvec(Dv)
    if not any(nv):
        continue  # in J mod p
    rk = nmod_mat(len(cur) + 1, len(sm), [x for r_ in cur for x in r_] + nv, p).rank()
    if rk > cur_rank:
        cur.append(nv); cur_rank = rk; chosen.append(i)
print('chosen non-J (mod p) generators independent of J + known:', chosen, '(rank over known', base_rank, '->', cur_rank, ')', flush=True)
# exact E1 lift for each chosen vector
out = []
for i in chosen:
    t1 = time.time()
    D = Wv[i]
    pars = fixed_d.numerator_parities(D, A)
    degN = d + 1
    Fip = [core.reduce_mod(Fi[k], p) for k in range(NV)]
    cols = []; colsQ = []
    for k in range(NV):
        for m in monomials(degN, pars[k]):
            cols.append(shift(Fip[k], m, 1, p)); colsQ.append((k, m))
    rhs = core.reduce_mod(mul(A, D), p)
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * C + j] = v
    Mx = nmod_mat(R, C, flat, p)
    Mr, rank = Mx.rref(); piv = rref_pivots(Mr, rank, C)
    MTr, r2 = Mx.transpose().rref(); prow = rref_pivots(MTr, r2, R)
    assert r2 == rank
    S = fmpq_mat(rank, rank); rowpos = {ri: x for x, ri in enumerate(prow)}
    for b_, cj in enumerate(piv):
        k, m = colsQ[cj]
        for e, v in shift(Fi[k], m, 1).items():
            ri = rows[e]
            if ri in rowpos:
                v = Fraction(v); S[rowpos[ri], b_] = fmpq(v.numerator, v.denominator)
    bq = fmpq_mat(rank, 1)
    AD = mul(A, D)
    for e, v in AD.items():
        ri = rows.get(e)
        if ri is not None and ri in rowpos:
            v = Fraction(v); bq[rowpos[ri], 0] = fmpq(v.numerator, v.denominator)
    x = S.solve(bq)
    N = [{} for _ in range(NV)]
    for b_, cj in enumerate(piv):
        v = x[b_, 0]
        if v != 0:
            k, m = colsQ[cj]; N[k][m] = Fraction(int(v.p), int(v.q))
    ok = sub(apply_field(N, F), AD) == {}
    hD = max(max(abs(Fraction(v).numerator), Fraction(v).denominator) for v in D.values())
    hN = max((max(abs(Fraction(v).numerator), Fraction(v).denominator) for Nk in N for v in Nk.values()), default=1)
    print(f'  generator {i}: exact E1 lift verified={ok} (rank {rank}, {time.time()-t1:.0f}s); heights D {len(str(hD))} digits, N {len(str(hN))} digits; terms D {len(D)}', flush=True)
    log(dict(exp='exact-E1-lift', d=d, prime=p, generator=i, verified=ok, rank=rank, height_D_digits=len(str(hD)), height_N_digits=len(str(hN)), seconds=round(time.time() - t1)), path=LOGF)
    if ok:
        out.append(dict(label=f'K{d}_{i}', D=core.to_str(D), N=[core.to_str(Nk) for Nk in N]))
json.dump(out, open(a.out, 'w'))
log(dict(exp='nonJ-generators-exact', d=d, prime=p, n=len(out), file=a.out, seconds=round(time.time() - t0), note='exact rational non-Jacobian E1-kernel generators with exact E1 lifts (N(F)=A D over Q)'), path=LOGF)
print('wrote', a.out, len(out), 'generators;', round(time.time() - t0), 's')
