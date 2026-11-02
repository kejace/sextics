"""E1 kernel sweep: for each degree d, dimension of {D Y-even, deg d : A*D in J(F)} modulo J(F)_d (mod p).
Returns kernel basis vectors (mod p) as polys when nonzero."""
import sys, time, json
import lib
from lib import *
from flint import nmod_mat
p = core.PRIMES[0]
import os
TARGET = os.environ.get('TARGET', 'X0**3*Y0**3')
A = core.parse(TARGET)
parA = core.parity_of(A)
degA = core.degree(A)
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]

def kernel_dim(d, parD=0):
    """D of Y-parity parD, degree d.  A*D has parity parD^1 (A odd).  Work in S_{d+6} of that parity."""
    t0 = time.time()
    parT = (parD + parA) % 2
    rows_mons = core.monomials(d + degA, parT)
    idx = {m: k for k, m in enumerate(rows_mons)}
    cols = []
    # Jacobian columns: F_i * m, parity(F_i) = 0 for i<3 else 1; need parity(m) = parT - parity(F_i)
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(d + degA - 5, (parT - pi) % 2):
            cols.append(core.shift(Fi[i], m, 1, p))
    nJ = len(cols)
    Dmons = core.monomials(d, parD)
    for m in Dmons:
        cols.append(core.shift(Ap, m, 1, p))
    R, C = len(rows_mons), len(cols)
    flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items():
            flat[idx[e] * C + j] = v
    Mx = nmod_mat(R, C, flat, p)
    t1 = time.time()
    X, nullity = Mx.nullspace()
    t2 = time.time()
    # project nullspace onto D-coordinates (cols nJ..C-1): dimension of {D: AD in J}
    nD = len(Dmons)
    if nullity == 0:
        projrank = 0; basis = []
    else:
        proj = nmod_mat(nD, nullity, [int(X[nJ + j, k]) for j in range(nD) for k in range(nullity)], p)
        T = proj.transpose(); Tr, projrank = T.rref()
        basis = []
        for i in range(projrank):
            poly = {Dmons[jj]: int(Tr[i, jj]) for jj in range(nD) if int(Tr[i, jj])}
            basis.append(poly)
    # dim of J_d with parity parD
    colsJ = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(d - 5, (parD - pi) % 2):
            colsJ.append(core.shift(Fi[i], m, 1, p))
    if colsJ:
        idxD = {m: k for k, m in enumerate(Dmons)}
        flatJ = [0] * (nD * len(colsJ))
        for j, c in enumerate(colsJ):
            for e, v in c.items(): flatJ[idxD[e] * len(colsJ) + j] = v
        rJ = nmod_mat(nD, len(colsJ), flatJ, p).rank()
    else:
        rJ = 0
    ker = projrank - rJ
    info = dict(d=d, parD=parD, rows=R, cols=C, nJcols=nJ, nD=nD, nullity=nullity, dim_E1_space=projrank,
                dim_J_d=rJ, kernel_mod_J=ker, t_build=round(t1 - t0, 1), t_null=round(t2 - t1, 1))
    return info, basis, rJ

if __name__ == '__main__':
    ds = [int(x) for x in sys.argv[1:]] or list(range(1, 7))
    for d in ds:
        for parD in (0, 1):
            info, basis, rJ = kernel_dim(d, parD)
            print(json.dumps(info), flush=True)
            log({'D': f'E1-kernel sweep d={d} parity={parD} target={TARGET}', 'd': d, 'family': 'e1-kernel', 'e1': None,
                 'consistent': None, 'rank': info['dim_E1_space'], 'nullity': info['kernel_mod_J'],
                 'seconds': info['t_build'] + info['t_null'], 'note': json.dumps(info)})
            if info['kernel_mod_J'] > 0:
                # save a basis of the E1 space (mod p) reduced: rows beyond J
                with open(f'{WORK}/e1kernel_d{d}_par{parD}_{TARGET.replace("*","")}.json', 'w') as fh:
                    json.dump({'info': info, 'basis_modp': [core.to_str(b) for b in basis]}, fh)
                print(f'  nonzero kernel at d={d} parity={parD}: {info["kernel_mod_J"]} (E1 space {info["dim_E1_space"]}, J_d {rJ})', flush=True)
