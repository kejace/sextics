"""Extract explicit E1-kernel elements: D in S_d^{par} with A*D in J, D not in J_d; reduce mod J_d,
rational-reconstruct from two primes, verify A*D in J over a third prime; save to e1vec_d{d}_p{par}.json."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
from flint import nmod_mat
from fractions import Fraction

def jac_cols(F, deg, par, p):
    Fi = [core.pderiv(F, i, p) for i in range(6)]
    cols = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(deg - 5, (par - pi) % 2):
            cols.append(core.shift(Fi[i], m, 1, p))
    return cols

def mat(polys, p, rows):
    R = len(rows); C = len(polys)
    flat = [0] * (R * C)
    for j, P in enumerate(polys):
        for e, v in P.items(): flat[rows[e] * C + j] = v
    return nmod_mat(R, C, flat, p)

def kernel_reps(F, A, d, parD, p):
    """Return list of D (dict mod p) spanning ker([A]: R_d -> R_{d+6}) reduced mod J_d (normal form)."""
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    monsD = core.monomials(d, parD); nD = len(monsD)
    Jcols = jac_cols(Fp, d + 6, (parD + 1) % 2, p)
    rows = {}
    for P in Jcols:
        for e in P: rows.setdefault(e, len(rows))
    AD = [core.shift(Ap, m, 1, p) for m in monsD]
    for P in AD:
        for e in P: rows.setdefault(e, len(rows))
    M = mat(Jcols + AD, p, rows)
    X, nul = M.nullspace()
    nJ = len(Jcols)
    # D-parts of nullspace vectors
    Dvecs = nmod_mat(nD, nul, [int(X[nJ + j, k]) for j in range(nD) for k in range(nul)], p)
    # reduce mod J_d: rref of J_d generators in monomial basis of S_d^{parD}
    JD = jac_cols(Fp, d, parD, p)
    idx = {m: i for i, m in enumerate(monsD)}
    K = len(JD)
    flat = [0] * (K * nD)
    for k, g in enumerate(JD):
        for m, v in g.items(): flat[k * nD + idx[m]] = v
    Rr0, r = nmod_mat(K, nD, flat, p).rref()
    Rr = nmod_mat(r, nD, [int(Rr0[i, j]) for i in range(r) for j in range(nD)], p)
    pivots = []
    for i in range(r):
        j = 0
        while int(Rr[i, j]) == 0: j += 1
        pivots.append(j)
    # normal forms: Dvecs^T - Dvecs^T[:, pivots] * Rr
    DT = Dvecs.transpose()
    Dp = nmod_mat(nul, r, [int(DT[k, pj]) for k in range(nul) for pj in pivots], p)
    NF = DT - Dp * Rr
    # independent normal forms: rref of NF
    NFr0, rk = NF.rref()
    reps = []
    for i in range(rk):
        D = {}
        for j in range(nD):
            v = int(NFr0[i, j])
            if v: D[monsD[j]] = v
        reps.append(D)
    return reps, {'nullity': nul, 'dim J_d': r, 'kernel mod J': rk}

if __name__ == '__main__':
    d = int(sys.argv[1]); parD = int(sys.argv[2])
    F, A, q, rec = load()
    t0 = time.time()
    p1, p2, p3 = core.PRIMES[0], core.PRIMES[1], core.PRIMES[2]
    reps1, info1 = kernel_reps(F, A, d, parD, p1)
    print('p1', info1, f'{time.time()-t0:.0f}s', flush=True)
    reps2, info2 = kernel_reps(F, A, d, parD, p2)
    print('p2', info2, f'{time.time()-t0:.0f}s', flush=True)
    out = []
    for D1, D2 in zip(reps1, reps2):
        keys = set(D1) | set(D2)
        D = {}
        ok = True
        for m in keys:
            a, b = D1.get(m, 0), D2.get(m, 0)
            x, mod = core.crt_pair(a, p1, b, p2)
            v = core.ratrecon(x, mod)
            if v is None: ok = False; break
            if v: D[m] = v
        if not ok:
            print('rational reconstruction failed for a representative', flush=True); continue
        # verify A*D in J mod p3
        okE1, info = fixed_d.e1_solvable(F, A, 1, D, p=p3)
        print('E1 (mod p3):', okE1, info, flush=True)
        out.append({'d': d, 'parD': parD, 'D': core.to_str(D), 'nterms': len(D), 'E1_p3': okE1})
    with open(f'/tmp/hodge/results/short-power-presentation/e1vec_d{d}_p{parD}.json', 'w') as f:
        json.dump(out, f, indent=1)
    for o in out: print(o['nterms'], 'terms:', o['D'][:600], flush=True)
