"""E1-kernel sweep: dim of {D in S_d^even : A*D in J(F)} modulo J(F)_d^even, for d = 1..dmax."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import ctr
from ctr import *
from flint import nmod_mat
p = core.PRIMES[0]
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]
res = {}
for d in range(0, dmax + 1):
    t0 = time.time()
    # columns: A*m for even m of degree d ; F_i * m' for m' of degree d+1 with parity so that F_i m' is Y-odd
    cols = []; kinds = []
    Dmons = core.monomials(d, 0)
    for m in Dmons:
        cols.append(core.shift(Ap, m, 1, p)); kinds.append('A')
    for i in range(6):
        par = 1 if i < 3 else 0   # F_i parity: even for i<3 -> m' odd ; odd for i>=3 -> m' even
        for m in core.monomials(d + 1, par):
            cols.append(core.shift(Fi[i], m, 1, p)); kinds.append('J')
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    nA = len(Dmons)
    # kernel of [A] on S_d^even -> R_{d+6}: projection of nullspace to A-coordinates
    if nullity:
        proj = nmod_mat(nA, nullity, [int(X[j, k]) for j in range(nA) for k in range(nullity)], p)
        kA = proj.rank()
    else:
        kA = 0
    # J(F)_d^even dimension: rank of multiplication columns F_i * S_{d-5}
    if d >= 5:
        cols2 = []
        for i in range(6):
            par = 0 if i < 3 else 1
            for m in core.monomials(d - 5, par):
                cols2.append(core.shift(Fi[i], m, 1, p))
        rows2 = {}
        for c in cols2:
            for e in c: rows2.setdefault(e, len(rows2))
        flat2 = [0] * (len(rows2) * len(cols2))
        for j, c in enumerate(cols2):
            for e, v in c.items(): flat2[rows2[e] * len(cols2) + j] = v
        dimJ = nmod_mat(len(rows2), len(cols2), flat2, p).rank() if cols2 else 0
    else:
        dimJ = 0
    res[d] = {'dim_Sd_even': nA, 'dim_Jd_even': dimJ, 'ker_A_in_Sd': kA, 'ker_A_mod_J': kA - dimJ,
              'rows': R, 'cols': C, 'seconds': round(time.time() - t0, 1)}
    print(d, res[d], flush=True)
    log({'engine': 2, 'D': f'E1-kernel sweep degree {d}', 'd': d, 'family': 'e1-kernel', 'e1': None, 'consistent': None,
         'rank': int(kA), 'nullity': int(kA - dimJ), 'seconds': res[d]['seconds'],
         'note': f"dim ker([A]: S_d^even -> R_(d+6)) = {kA}; dim J_d^even = {dimJ}; kernel mod J = {kA - dimJ}"})
json.dump(res, open(f'e1kernel_{dmax}.json', 'w'), indent=1)
