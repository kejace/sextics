"""E1-kernel sweep: for degree d, compute dim of {D in S_d^+ : A*D in J(F)} modulo J(F)_d^+,
i.e. the kernel of multiplication by [A] : R_d^+ -> R_{d+6}^-  (Jacobian ring, mod p).
Also returns an explicit basis of the kernel (as polys mod p) for small d."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
from flint import nmod_mat
p = core.PRIMES[0]

def jac_cols(F, deg, par, p):
    """Columns m*F_i spanning J(F)_deg with given Y-parity (par of the product)."""
    Fi = [core.pderiv(F, i, p) for i in range(6)]
    cols = []
    for i in range(6):
        pi = 0 if i < 3 else 1   # parity of F_i
        for m in core.monomials(deg - 5, (par - pi) % 2):
            cols.append(core.shift(Fi[i], m, 1, p))
    return cols

def rowspace_rank(polys, p, rows=None):
    if rows is None:
        rows = {}
        for P in polys:
            for e in P: rows.setdefault(e, len(rows))
    R = len(rows); C = len(polys)
    flat = [0] * (R * C)
    for j, P in enumerate(polys):
        for e, v in P.items(): flat[rows[e] * C + j] = v
    return nmod_mat(R, C, flat, p), rows

def e1_kernel(F, A, d, verbose=True):
    t0 = time.time()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    parD = 0  # Y-even D
    monsD = core.monomials(d, parD)
    # target space: degree d+6, parity of A*D = odd
    Jcols = jac_cols(Fp, d + 6, 1, p)
    rows = {}
    for P in Jcols:
        for e in P: rows.setdefault(e, len(rows))
    AD = [core.shift(Ap, m, 1, p) for m in monsD]
    for P in AD:
        for e in P: rows.setdefault(e, len(rows))
    R = len(rows)
    # matrix [J | A*D]; kernel vectors with nonzero D-part give E1-solutions.
    # Compute rank of J and rank of [J | AD]; kernel dim of AD mod J = len(monsD) - (rank[J|AD] - rankJ)
    MJ, _ = rowspace_rank(Jcols, p, rows)
    rJ = MJ.rank()
    Mall, _ = rowspace_rank(Jcols + AD, p, rows)
    rall = Mall.rank()
    kdim = len(monsD) - (rall - rJ)
    # J_d^+ dimension (always in kernel)
    JD = jac_cols(Fp, d, 0, p)
    rJD = rowspace_rank(JD, p)[0].rank() if JD else 0
    out = {'d': d, 'dim S_d^+': len(monsD), 'dim J_d^+': rJD, 'dim R_d^+': len(monsD) - rJD,
           'rank J_{d+6}^-': rJ, 'kernel dim (incl J_d^+)': kdim, 'kernel mod J': kdim - rJD,
           'rows': R, 'cols': len(Jcols) + len(monsD), 'seconds': round(time.time() - t0, 1)}
    if verbose: print(json.dumps(out), flush=True)
    return out, (Jcols, AD, monsD, rows)

def e1_kernel_basis(F, A, d):
    """Explicit basis (mod p) of {D : A D in J} as polys, via nullspace of [J | AD]."""
    out, (Jcols, AD, monsD, rows) = e1_kernel(F, A, d, verbose=False)
    M, _ = rowspace_rank(Jcols + AD, p, rows)
    X, nullity = M.nullspace()
    nJ = len(Jcols)
    basis = []
    for k in range(nullity):
        D = {}
        for j, m in enumerate(monsD):
            v = int(X[nJ + j, k])
            if v: D[m] = v
        if D: basis.append(D)
    # reduce to independent D-parts
    return out, basis, monsD

if __name__ == '__main__':
    F, A, q, rec = load()
    ds = [int(a) for a in sys.argv[1:]] or [2, 4, 6]
    for d in ds:
        out, _ = e1_kernel(F, A, d)
        with open('e1kernel.jsonl', 'a') as f: f.write(json.dumps(out) + '\n')
