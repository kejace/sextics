"""E1-kernel sweep: dimension of {D homogeneous of degree d, given Y-parity : A*D in J(F)} (mod p)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from lib import *
from flint import nmod_mat
from hodge.core import monomials, shift, pderiv, reduce_mod, yparity

def e1_kernel(d, parD=0, p=core.PRIMES[0], want_basis=True):
    t0 = time.time()
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    pA = core.parity_of(A)
    degN = d + 1
    cols = []; colpolys = []
    for m in monomials(d, parD):
        cols.append(('D', m)); colpolys.append(shift(Ap, m, 1, p))
    for i in range(6):
        parN = (pA - parD - (0 if i < 3 else 1)) % 2
        for m in monomials(degN, parN):
            cols.append(('N', i, m)); colpolys.append(shift(Fi[i], m, p - 1, p))
    rows = {}
    for col in colpolys:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(colpolys):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    t1 = time.time()
    X, nullity = M.nullspace()
    t2 = time.time()
    nD = sum(1 for c in cols if c[0] == 'D')
    # project nullspace onto D coordinates
    proj = nmod_mat(nD, nullity, [int(X[j, k]) for j in range(nD) for k in range(nullity)], p) if nullity else None
    rk = proj.rank() if nullity else 0
    # dimension of J_d with parity parD: span of m*F_i, deg m = d-5
    nJ = 0
    if d >= 5:
        jcols = []
        for i in range(6):
            parM = (parD - (0 if i < 3 else 1)) % 2
            for m in monomials(d - 5, parM):
                jcols.append(shift(Fi[i], m, 1, p))
        jrows = {}
        for col in jcols:
            for e in col: jrows.setdefault(e, len(jrows))
        flatJ = [0] * (len(jrows) * len(jcols))
        for j, col in enumerate(jcols):
            for e, v in col.items(): flatJ[jrows[e] * len(jcols) + j] = v
        nJ = nmod_mat(len(jrows), len(jcols), flatJ, p).rank() if jcols else 0
    basis = []
    if want_basis and rk:
        T = proj.transpose(); Tr, r2 = T.rref()
        for i in range(r2):
            poly = {}
            for j in range(nD):
                v = int(Tr[i, j])
                if v: poly[cols[j][1]] = v
            basis.append(poly)
    return {'d': d, 'parD': parD, 'rows': R, 'cols': C, 'nD': nD, 'nullity': nullity, 'kerdim': rk, 'dimJ_d': nJ,
            'kermodJ': rk - nJ, 't_build': round(t1 - t0, 1), 't_null': round(t2 - t1, 1)}, basis

if __name__ == '__main__':
    ds = [int(a) for a in sys.argv[1:]] or [1, 2, 3, 4, 5, 6]
    for d in ds:
        for parD in (0, 1):
            info, basis = e1_kernel(d, parD)
            print(json.dumps(info), flush=True)
            with open(f'{WD}/e1kernel.jsonl', 'a') as f:
                f.write(json.dumps(info) + '\n')
            if basis:
                with open(f'{WD}/e1kernel_basis_d{d}_par{parD}.txt', 'w') as f:
                    for b in basis:
                        f.write(core.to_str({k: (v if v < 2**30 else v - core.PRIMES[0]) for k, v in b.items()}) + '\n')
