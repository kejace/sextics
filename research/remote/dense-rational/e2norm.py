"""Normalized E2 test for a degree-d denominator D with a given E1 numerator N0 (N0(F) = A D), mod p:
   B_D(N0) + T_D(K) = F*C  with K skew (15 entries of degree d-4), C in S_{2d-6}^- ;  rows S_{2d}^-.
   T_D(K) = sum_{i<j} ( D[(d_i K_ij) F_j - (d_j K_ij) F_i] - K_ij (D_i F_j - D_j F_i) )  (= B_D(K grad F)).
Inconsistent => no certificate with this D (rigorous over C only in the 'rank attains cap' sense: inconsistency mod p
of a system whose coefficients are reductions of a Q-system is rigorous iff rank_aug(p) > rank(p) = rank(Q)... we report the modular statement).
Usage: python3 e2norm.py reps.json [index]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
fn = sys.argv[1]; which = [int(x) for x in sys.argv[2:]]
data = json.load(open(fn)); p = data['p']; d = data['d']
F, A, q = load()
import os
if os.environ.get('HODGE_F'): F = core.parse(os.environ['HODGE_F'])
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
def parse(dct): return {tuple(int(x) for x in k.split()): v for k, v in dct.items()}
for ri, rep in enumerate(data['reps']):
    if which and ri not in which: continue
    t0 = time.time()
    D = parse(rep['D']); N0 = [parse(x) for x in rep['N']]
    # check E1
    lhs = {}
    for i in range(NV): lhs = add(lhs, mul(N0[i], Fi[i], p), p)
    assert not sub(lhs, mul(Ap, D, p), p), 'E1 fails for rep'
    Di = [pderiv(D, i, p) for i in range(NV)]
    divN = {}
    for i in range(NV): divN = add(divN, pderiv(N0[i], i, p), p)
    ND = {}
    for i in range(NV): ND = add(ND, mul(N0[i], Di[i], p), p)
    rhs = sub(mul(D, divN, p), ND, p)          # B_D(N0); we need B_D(N0) + T_D(K) - F C = 0  -> T_D(K) - F C = -B_D(N0)
    rhs = scale(rhs, p - 1, p)
    pars = fixed_d.numerator_parities(D, Ap)
    cols = []; names = []
    DFij = {}
    for i in range(NV):
        for j in range(i + 1, NV):
            DFij[(i, j)] = sub(mul(Di[i], Fi[j], p), mul(Di[j], Fi[i], p), p)
    for i in range(NV):
        for j in range(i + 1, NV):
            # K_ij parity: N_i gets K_ij F_j -> parity(K_ij) = pars[i] - parity(F_j)
            pj = 0 if j < 3 else 1
            for m in monomials(d - 4, (pars[i] - pj) % 2):
                mm = {m: 1}
                t = mul(D, sub(mul(pderiv(mm, i, p), Fi[j], p), mul(pderiv(mm, j, p), Fi[i], p), p), p)
                t = sub(t, mul(mm, DFij[(i, j)], p), p)
                cols.append(t); names.append(('K', i, j, m))
    for m in monomials(2 * d - 6, 1):
        cols.append(scale(shift(Fp, m, 1, p), p - 1, p)); names.append(('C', m))
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    M = nmod_mat(R, C + 1, p)
    for j, c in enumerate(cols):
        for e, v in c.items(): M[rows[e], j] = v
    for e, v in rhs.items(): M[rows[e], C] = v
    t1 = time.time()
    r_aug = M.rank()
    M0 = nmod_mat(R, C, p)
    for j, c in enumerate(cols):
        for e, v in c.items(): M0[rows[e], j] = v
    r0 = M0.rank()
    rec = {'D': f'd={d} E1-kernel rep: ' + rep['kind'], 'd': d, 'family': 'e1kernel-reps-E2norm', 'e1': True,
           'consistent': r0 == r_aug, 'rank': r0, 'nullity': C - r0, 'seconds': round(time.time() - t0, 1),
           'note': f'normalized E2 (K skew + C) mod {p}: rows {R} cols {C} rank {r0} rank_aug {r_aug}; terms(D)={len(D)}', 'engine': 2}
    print(json.dumps(rec), flush=True)
    with open(LOG, 'a') as f: f.write(json.dumps(rec) + '\n')
