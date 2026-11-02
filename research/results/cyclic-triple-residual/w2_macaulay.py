"""Rigorous exclusion of Jacobian-type degrees via the mod-J relaxation as a quadric system.

For d with Ann([A])_d = 0 in R = S/J(F) (d <= 7 here), every E1-admissible D is D = w.grad F (w of degree d-5,
parity-compatible, w in the w-space W_d of dimension dim J_d), with the canonical lift N0 = A w / q, and E2 forces
    [Bil(D, N0)] = 0 in R_{2d}^-   (necessary),   Bil(D, N) = D div N - N(D).
This is a system of dim R_{2d}^- quadrics q_s(w) in dim W_d unknowns. If the degree-k Macaulay matrix of the q_s has
full column rank mod p (equivalently the Hilbert function of (q_s) vanishes in degree k mod p), it has full column
rank over Q, so the only solution over C is w = 0: no certificate with deg D = d exists.

usage: python3 w2_macaulay.py d parity(0|1) [p] -> writes w2-quadrics-d{d}-par{par}.ms (msolve input) and a JSON summary.
"""
import os, sys, time, json
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from w2core import *
from flint import nmod_mat
from hodge import jacobian

d = int(sys.argv[1]); par = int(sys.argv[2]); p = int(sys.argv[3]) if len(sys.argv) > 3 else 1073741827
F, A, q = load()
t0 = time.time()
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
# w-space: D = sum_i w_i F_i, w_i of degree d-5 with parity par (i<3) / par+1 (i>=3) so that D has Y-parity par
W = []
for i in range(6):
    pi = (par + (0 if i < 3 else 1)) % 2
    for m in core.monomials(d - 5, pi):
        W.append((i, m))
nW = len(W)
Dk = [core.shift(Fi[i], m, 1, p) for (i, m) in W]
N0k = []
for (i, m) in W:
    N = [{} for _ in range(6)]
    N[i] = core.shift(Ap, m, pow(q, -1, p), p)
    N0k.append(N)
def bil_p(D, N):
    divN, ND = {}, {}
    for i in range(6):
        if N[i]:
            divN = core.add(divN, pderiv(N[i], i, p), p)
            ND = core.add(ND, core.mul(N[i], pderiv(D, i, p), p), p)
    return core.sub(core.mul(D, divN, p), ND, p)
# quadric coefficient polys: pair (j,k): Bil(D_j, N0_k) + Bil(D_k, N0_j) (j<k), Bil(D_j, N0_j) (j==k)
pairs = []; polys = []
for j in range(nW):
    for k in range(j, nW):
        P = bil_p(Dk[j], N0k[k])
        if k != j:
            P = core.add(P, bil_p(Dk[k], N0k[j]), p)
        pairs.append((j, k)); polys.append(P)
t1 = time.time()
# reduce modulo J_{2d}^- (rows = pairs, columns = monomials)
red = jacobian._Reducer(Fp, 2 * d, p, 1)
mons = red.mons; idx = red.idx; C = len(mons)
rank = red.rank; piv = red.pivots
V = nmod_mat(len(polys), C, [0] * (len(polys) * C), p)
flat = [0] * (len(polys) * C)
for r_, P in enumerate(polys):
    base = r_ * C
    for e, v in P.items():
        flat[base + idx[e]] = v
V = nmod_mat(len(polys), C, flat, p)
# Jr: rank x C rref rows; NF(V) = V - V[:, piv] * Jr
Jr = nmod_mat(rank, C, [0] * (rank * C), p)
jflat = [0] * (rank * C)
for i, row in enumerate(red.rows):
    for e, v in row.items():
        jflat[i * C + idx[e]] = v
Jr = nmod_mat(rank, C, jflat, p)
VP = nmod_mat(len(polys), rank, [int(V[r_, piv[i]]) for r_ in range(len(polys)) for i in range(rank)], p)
NF = V - VP * Jr
std_cols = [j for j in range(C) if j not in set(piv)]
assert len(std_cols) == C - rank
t2 = time.time()
# quadrics: for each standard monomial column s: q_s(w) = sum_{pairs} NF[pair, s] w_j w_k
quadrics = []
nz_pairs_total = 0
for s in std_cols:
    terms = []
    for r_, (j, k) in enumerate(pairs):
        v = int(NF[r_, s])
        if v:
            terms.append((v, j, k))
    if terms:
        quadrics.append(terms); nz_pairs_total += len(terms)
# rank of the coefficient matrix (quadrics x pairs) = number of independent quadrics
Qm = nmod_mat(len(std_cols), len(pairs), [int(NF[r_, s]) for s in std_cols for r_ in range(len(pairs))], p)
qrank = Qm.rank()
# msolve input
names = [f'w{j}' for j in range(nW)]
path = f'{WD}/w2-quadrics-d{d}-par{par}.ms'
with open(path, 'w') as f:
    f.write(','.join(names) + '\n' + str(p) + '\n')
    lines = []
    for terms in quadrics:
        lines.append('+'.join(f'{v}*{names[j]}*{names[k]}' if j != k else f'{v}*{names[j]}^2' for v, j, k in terms))
    f.write(',\n'.join(lines) + '\n')
summary = {'batch': 'macaulay-prep', 'd': d, 'parity': par, 'p': p, 'nW': nW, 'n_pairs': len(pairs), 'dim_S2d_odd': C,
           'rank_J2d': rank, 'dim_R2d_odd': C - rank, 'n_quadrics_nonzero': len(quadrics), 'quadric_rank': qrank,
           'avg_terms': round(nz_pairs_total / max(1, len(quadrics)), 1), 'msolve_input': path,
           'seconds': round(time.time() - t0, 1), 't_pairs': round(t1 - t0, 1), 't_reduce': round(t2 - t1, 1)}
print(json.dumps(summary), flush=True)
w2log(summary)
