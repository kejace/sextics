#!/usr/bin/env python3
"""w2_closure.py -- exact (mod p) closure of Jacobian-type degrees for quadratic-cubic-pullback.

At a degree d with K_d^par = J_d^par (every E1-admissible D is sum_l w_l F_l), E2 implies the stage-1 condition
    q(w) := [ B_D(N0) ] = 0 in R_{2d}^odd = S_{2d}^odd / J(F)_{2d},   D = sum w_l F_l, N0 = A w,
a system of rJ = dim R_{2d}^odd quadratic forms in the w-coordinates (n of them).  The forms are invariant under the
spurious translations w -> w + G z/6 (D -> D + F G), so they live on the quotient V = w-space / spurious, dim n - m.
If the Macaulay matrix of the forms in degree k (rows: monomial * form, columns: monomials of degree k in the quotient
coordinates) has full column rank, the forms have no common nonzero zero (the ideal contains (v)^k), hence
    NO denominator of degree d and this parity admits a certificate (E2 fails for all of K_d).
Degree-2 test = rank of the coefficient matrix (rJ x n'(n'+1)/2).  All exact mod p = 2^31-1 (a full-rank mod-p
Macaulay matrix implies full rank over Q, hence the statement over Q/C).

usage: python3 w2_closure.py --d 6 --parity 0 [--kmax 3]
"""
import sys, os, json, time, argparse, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from flint import nmod_mat
from hodge import core, jacobian
from hodge.core import NV, monomials, mul, add, scale, pderiv, reduce_mod, shift, sub

LANE = 'quadratic-cubic-pullback'
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True); ap.add_argument('--parity', type=int, default=0)
ap.add_argument('--kmax', type=int, default=2); ap.add_argument('--instance', default=LANE)
ap.add_argument('--fermat', action='store_true')
a = ap.parse_args()
d, parD = a.d, a.parity
p = core.PRIMES[0]
F, A, q, _ = core.load_instance(a.instance)
if a.fermat:
    F = {tuple(6 if j == i else 0 for j in range(6)): Fraction(1) for i in range(6)}
Fp, Ap = reduce_mod(F, p), reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
t0 = time.time()
# w-basis: w_l in S_{d-5} with parity so that w_l F_l has parity parD  (F_l even for l<3, odd for l>=3)
wlab = []
for l in range(NV):
    par_w = (parD - (0 if l < 3 else 1)) % 2
    for mu in monomials(d - 5, par_w):
        wlab.append((l, mu))
n = len(wlab)
Dcols = [shift(Fi[l], mu, 1, p) for (l, mu) in wlab]
Ncols = [(l, shift(Ap, mu, 1, p)) for (l, mu) in wlab]     # N0 has only component l nonzero
# spurious subspace: w = G z/6 for G in S_{d-6}^{parD} : coordinates in the w-basis
widx = {lab: i for i, lab in enumerate(wlab)}
spur = []
inv6 = pow(6, -1, p)
for G in (monomials(d - 6, parD) if d >= 6 else []):
    v = [0] * n
    for l in range(NV):
        mu = list(G); mu[l] += 1; mu = tuple(mu)
        v[widx[(l, mu)]] = inv6
    spur.append(v)
m = len(spur)
# quotient coordinates: rref the spurious vectors, drop pivot coordinates
pivots = []
if m:
    Sm = nmod_mat(m, n, [x for row in spur for x in row], p)
    Sr, rk = Sm.rref()
    assert rk == m
    for i in range(rk):
        j = 0
        while int(Sr[i, j]) == 0: j += 1
        pivots.append(j)
free = [j for j in range(n) if j not in set(pivots)]
nq = len(free)
# reducer modulo J_{2d}^odd
red = jacobian._Reducer(Fp, 2 * d, p, 1)
sm = red.standard; sidx = {mo: i for i, mo in enumerate(sm)}; rJ = len(sm)
print(f'[closure d={d} par={parD}] n={n} spurious m={m} quotient dim nq={nq}; rJ = dim R_{2*d}^odd = {rJ} '
      f'(J rank {red.jacobian_dim}); reducer {time.time()-t0:.1f}s', flush=True)

def Bpair(a_, b_):
    """B(e_a, e_b) = D_a div N0_b - N0_b . grad D_a   (mod p)."""
    Da = Dcols[a_]
    lb, Nb = Ncols[b_]
    divNb = pderiv(Nb, lb, p)
    out = mul(Da, divNb, p)
    out = sub(out, mul(Nb, pderiv(Da, lb, p), p), p)
    return out

# quadratic form coefficients on the free coordinates: q_r(v) = sum_{a<=b} C[r, (a,b)] v_a v_b
pairs = [(a_, b_) for a_ in range(nq) for b_ in range(a_, nq)]
pidx = {pr: i for i, pr in enumerate(pairs)}
C = [[0] * len(pairs) for _ in range(rJ)]
t1 = time.time()
for (a_, b_) in pairs:
    ia, ib = free[a_], free[b_]
    P = Bpair(ia, ib)
    if a_ != b_:
        P = add(P, Bpair(ib, ia), p)
    nf = red.normal_form(P)
    j = pidx[(a_, b_)]
    for mo, c in nf.items():
        C[sidx[mo]][j] = c
print(f'[closure] {len(pairs)} pair products + normal forms in {time.time()-t1:.1f}s', flush=True)
Cm = nmod_mat(rJ, len(pairs), [x for row in C for x in row], p)
r2 = Cm.rank()
full2 = (r2 == len(pairs))
rec = {'lane': LANE, 'item': 'closure', 'd': d, 'parity': parD, 'n_w': n, 'spurious': m, 'nq': nq, 'rJ': rJ,
       'deg2_forms_rank': r2, 'deg2_dim': len(pairs), 'closed_at_k': 2 if full2 else None, 'p': p,
       'fermat_control': bool(a.fermat)}
print(json.dumps(rec), flush=True)
if not full2 and a.kmax >= 3:
    # Macaulay matrix in degree 3: rows = (variable j) x (nonzero form rows r), columns = cubic monomials in nq vars
    # forms as sparse rows
    rows_forms = []
    Cr, rk2 = Cm.rref()
    for i in range(rk2):
        rows_forms.append({pairs[j]: int(Cr[i, j]) for j in range(len(pairs)) if int(Cr[i, j])})
    trip = {}
    def tidx(t):
        t = tuple(sorted(t))
        if t not in trip:
            trip[t] = len(trip)
        return trip[t]
    ncub = len(list(itertools.combinations_with_replacement(range(nq), 3)))
    for t in itertools.combinations_with_replacement(range(nq), 3):
        tidx(t)
    R = len(rows_forms) * nq
    print(f'[closure] Macaulay degree 3: {R} x {ncub} (rank needed {ncub})', flush=True)
    flat = [0] * (R * ncub)
    r = 0
    for fr in rows_forms:
        for j in range(nq):
            base = r * ncub
            for (a_, b_), c in fr.items():
                flat[base + tidx((a_, b_, j))] = (flat[base + tidx((a_, b_, j))] + c) % p
            r += 1
    t2 = time.time()
    M3 = nmod_mat(R, ncub, flat, p)
    del flat
    r3 = M3.rank()
    rec.update({'deg3_rank': r3, 'deg3_dim': ncub, 'closed_at_k': 3 if r3 == ncub else None, 't_deg3': round(time.time() - t2, 1)})
    print(json.dumps(rec), flush=True)
rec['seconds'] = round(time.time() - t0, 1); rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S')
with open(f'/tmp/hodge/results/{LANE}/w2_closure.jsonl', 'a') as f:
    f.write(json.dumps(rec) + '\n')
