"""Rigorous exclusion at Jacobian-type degrees by block Macaulay ranks (G-graded).

Quadrics q_s(w) = [Bil(w.gradF, A w/q)]_s in R_{2d}^- (s standard monomial), w in the w-space W_d (D = w.gradF, parity par).
G-grading: wchar(w_j) = chi(D_j); q_s has w-character chi(s) - chiA; the degree-k Macaulay matrix (rows m*q_s, deg m = k-2;
columns = degree-k monomials in w) is block diagonal by w-character. Full column rank of every block mod p => the ideal
contains all degree-k forms mod p => (rank over F_p <= rank over Q) also over Q => V(q) = {0}: no E1-admissible D of degree d
and parity par satisfies E2. usage: python3 w2_macblocks.py d par [k=3] [p]"""
import os, sys, time, json, itertools
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from w2core import *
from flint import nmod_mat
from hodge import jacobian

d = int(sys.argv[1]); par = int(sys.argv[2]); kdeg = int(sys.argv[3]) if len(sys.argv) > 3 else 3
p = int(sys.argv[4]) if len(sys.argv) > 4 else 1073741827
t0 = time.time()
F, A, q = load()
G = full_group(F, A); S = Sect(F, A, G)
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
restrict = len(sys.argv) > 5 and sys.argv[5] == 'restrict'
W = [(i, m) for i in range(6) for m in core.monomials(d - 5, (par + (0 if i < 3 else 1)) % 2)]
# spurious locus L = {l Z/6 : l in S_{d-6}} (D = l F); V(q) is invariant under translation by L (gauge), so restricting
# to the complement W' = {coordinates (0, m) with X0 | m vanish} is lossless: V(q) = L + (V(q) cap W').
n_dropped = 0
if restrict and d >= 6:
    W2 = [(i, m) for (i, m) in W if not (i == 0 and m[0] >= 1)]
    n_dropped = len(W) - len(W2); W = W2
nW = len(W)
Dk = [core.shift(Fi[i], m, 1, p) for (i, m) in W]
wchar = [sect.cadd(S.chi(m), tuple((-x) % 6 for x in S.chiE[i])) for (i, m) in W]   # chi(D_j) = chi(m) - chi(e_i)
N0k = []
for (i, m) in W:
    N = [{} for _ in range(6)]; N[i] = core.shift(Ap, m, pow(q, -1, p), p); N0k.append(N)
def bil_p(D, N):
    divN, ND = {}, {}
    for i in range(6):
        if N[i]:
            divN = core.add(divN, pderiv(N[i], i, p), p)
            ND = core.add(ND, core.mul(N[i], pderiv(D, i, p), p), p)
    return core.sub(core.mul(D, divN, p), ND, p)
pairs, polys = [], []
for j in range(nW):
    for k in range(j, nW):
        P = bil_p(Dk[j], N0k[k])
        if k != j: P = core.add(P, bil_p(Dk[k], N0k[j]), p)
        pairs.append((j, k)); polys.append(P)
t1 = time.time()
red = jacobian._Reducer(Fp, 2 * d, p, 1)
mons = red.mons; idx = red.idx; C = len(mons); rank = red.rank; piv = red.pivots
flat = [0] * (len(polys) * C)
for r_, P in enumerate(polys):
    base = r_ * C
    for e, v in P.items(): flat[base + idx[e]] = v
V = nmod_mat(len(polys), C, flat, p)
jflat = [0] * (rank * C)
for i, row in enumerate(red.rows):
    for e, v in row.items(): jflat[i * C + idx[e]] = v
Jr = nmod_mat(rank, C, jflat, p)
VP = nmod_mat(len(polys), rank, [int(V[r_, piv[i]]) for r_ in range(len(polys)) for i in range(rank)], p)
NF = V - VP * Jr
pivset = set(piv)
std_cols = [j for j in range(C) if j not in pivset]
t2 = time.time()
# quadrics: dict s -> {pair index: coeff}
quad = {}
for s in std_cols:
    terms = {}
    for r_ in range(len(pairs)):
        v = int(NF[r_, s])
        if v: terms[r_] = v
    if terms: quad[s] = terms
nq = len(quad)
Qm = nmod_mat(len(quad), len(pairs), [quad[s_].get(r_, 0) for s_ in quad for r_ in range(len(pairs))], p)
quad_rank = Qm.rank() if nq else 0
schar = {s: sect.cadd(S.chi(mons[s]), tuple((-x) % 6 for x in S.chiA)) for s in quad}   # w-character of q_s
# degree-k Macaulay blocks (vectorized): rows m*q_s (deg m = k-2), columns = degree-k monomials in w, block-diagonal by w-character
import numpy as np, random
qs_list = list(quad.keys())
Qmat = np.zeros((len(qs_list), len(pairs)), dtype=np.int64)
for qi, s_ in enumerate(qs_list):
    for r_, v in quad[s_].items():
        Qmat[qi, r_] = v
pair_index = {pr: r_ for r_, pr in enumerate(pairs)}
def wchar_of(mon):
    c = (0,) * len(S.chiA)
    for j in mon: c = sect.cadd(c, wchar[j])
    return c
kmons = list(itertools.combinations_with_replacement(range(nW), kdeg))
kidx_by_char = {}
for mon in kmons:
    kidx_by_char.setdefault(wchar_of(mon), []).append(mon)
mult_mons = list(itertools.combinations_with_replacement(range(nW), kdeg - 2))
qchar_groups = {}
for qi, s_ in enumerate(qs_list):
    qchar_groups.setdefault(schar[s_], []).append(qi)
blocks = {}
random.seed(1)
for c, cols in kidx_by_char.items():
    colidx = {mon: t_ for t_, mon in enumerate(cols)}
    C_ = len(cols)
    row_blocks = []
    for mm in mult_mons:
        target = tuple((a - b) % 6 for a, b in zip(c, wchar_of(mm)))   # need schar(q) + wchar(mm) = c
        qis = qchar_groups.get(target)
        if not qis: continue
        # column index for each pair under multiplication by mm (all pairs of the right character land in this block)
        colmap = np.full(len(pairs), -1, dtype=np.int64)
        for r_, (j, k) in enumerate(pairs):
            mon = tuple(sorted(mm + (j, k)))
            ci = colidx.get(mon)
            if ci is not None: colmap[r_] = ci
        sub = Qmat[qis][:, colmap >= 0]
        cm = colmap[colmap >= 0]
        blk = np.zeros((len(qis), C_), dtype=np.int64)
        blk[:, cm] = sub      # injective for fixed mm
        row_blocks.append(blk)
    if not row_blocks:
        blocks[c] = {'rows': 0, 'cols': C_, 'rank': 0, 'HF': C_}; continue
    Mb = np.vstack(row_blocks) % p
    R_ = Mb.shape[0]
    sel = np.arange(R_)
    if R_ > 4 * C_ + 50:
        sel = np.array(sorted(random.sample(range(R_), 4 * C_ + 50)))
    def rank_of(rows_sel):
        return nmod_mat(len(rows_sel), C_, Mb[rows_sel].ravel().tolist(), p).rank()
    rk = rank_of(sel)
    if rk < C_ and len(sel) < R_:
        rk = rank_of(np.arange(R_))
    blocks[c] = {'rows': int(R_), 'cols': int(C_), 'rank': int(rk), 'HF': int(C_ - rk)}
HFk = sum(b['HF'] for b in blocks.values())
res = {'batch': 'macaulay-blocks', 'd': d, 'parity': par, 'k': kdeg, 'p': p, 'nW': nW, 'restricted_to_complement_of_spurious': restrict, 'n_dropped_coords': n_dropped, 'n_quadrics': nq, 'quadric_span_rank': quad_rank, 'n_pairs': len(pairs), 'dim_R2d_odd': C - rank,
       'blocks': {csec(c): b for c, b in sorted(blocks.items())}, 'HF_k_total': HFk, 'rigorous_exclusion': HFk == 0,
       'note': f'degree-{kdeg} Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-{kdeg} forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = {d}, Y-parity {par}',
       't_pairs': round(t1 - t0, 1), 't_reduce': round(t2 - t1, 1), 'seconds': round(time.time() - t0, 1)}
print(json.dumps({k: v for k, v in res.items() if k not in ('blocks', 'note')}), flush=True)
print('blocks:', {csec(c): (b['rows'], b['cols'], b['HF']) for c, b in sorted(blocks.items())}, flush=True)
w2log(res)
