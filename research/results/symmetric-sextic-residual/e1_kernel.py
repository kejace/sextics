"""E1 kernel sweep: dim { D of degree d, parity pD : A*D in J(F) } modulo J(F)_d, mod p.
Also records a basis of the kernel (as polys mod p, symmetric lifts) for later testing."""
import sys, time, json, pickle
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
import flint
from flint import nmod_mat

p = core.PRIMES[0]
F, A, q = load()
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]

def jac_cols(deg, tp):
    """columns F_i * m spanning J_deg with parity tp"""
    cols = []
    for i in range(6):
        par_m = tp if i < 3 else (tp + 1) % 2
        for m in core.monomials(deg - 5, par_m):
            cols.append(core.shift(Fi[i], m, 1, p))
    return cols

def mat_from_cols(cols, rows):
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    return nmod_mat(R, C, flat, p)

def kernel(d, pD):
    t0 = time.time()
    tp = (1 + pD) % 2
    rowmons = core.monomials(d + 6, tp)
    rows = {e: i for i, e in enumerate(rowmons)}
    J = jac_cols(d + 6, tp) if d + 6 >= 5 else []
    Dmons = core.monomials(d, pD)
    Acols = [core.shift(Ap, m, 1, p) for m in Dmons]
    nJ, nA = len(J), len(Acols)
    M = mat_from_cols(J + Acols, rows)
    X, nullity = M.nullspace()
    # project nullspace onto A-coordinates
    proj = nmod_mat(nA, nullity, [int(X[nJ + j, k]) for j in range(nA) for k in range(nullity)], p)
    T = proj.transpose()
    Tr, rk = T.rref()
    Kbasis = []
    for i in range(rk):
        poly = {}
        for j in range(nA):
            v = int(Tr[i, j])
            if v:
                poly[Dmons[j]] = v
        Kbasis.append(poly)
    dimK = rk
    # J_d of parity pD
    if d >= 5:
        Jd = jac_cols(d, pD)
        rowsd = {e: i for i, e in enumerate(Dmons)}
        MJ = mat_from_cols(Jd, rowsd)
        dimJd = MJ.rank()
        # complement: rref [Jd | K] to find K-vectors independent of J_d
        MK = mat_from_cols(Jd + Kbasis, rowsd)
        Mr, rk2 = MK.rref()
        pivots = []
        for i in range(rk2):
            j = 0
            while j < MK.ncols() and int(Mr[i, j]) == 0: j += 1
            pivots.append(j)
        extra = [Kbasis[j - len(Jd)] for j in pivots if j >= len(Jd)]
    else:
        dimJd = 0
        extra = Kbasis
    dt = time.time() - t0
    return {'d': d, 'pD': pD, 'rows': len(rowmons), 'nJ': nJ, 'nA': nA, 'nullity': nullity,
            'dimK': dimK, 'dimJd': dimJd, 'ker_R': dimK - dimJd, 'seconds': round(dt, 1)}, extra

if __name__ == '__main__':
    dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    out = {}
    for d in range(1, dmax + 1):
        for pD in (0, 1):
            info, extra = kernel(d, pD)
            print(json.dumps(info), flush=True)
            out[(d, pD)] = (info, extra)
            with open(f'{RES}/e1_kernel_d{d}_p{pD}.pkl', 'wb') as f:
                pickle.dump((info, extra), f)
            for Dk in extra[:40]:
                lift = {e: (v if v < p // 2 else v - p) for e, v in Dk.items()}
                print('   ker elt:', core.to_str(lift)[:200], flush=True)
