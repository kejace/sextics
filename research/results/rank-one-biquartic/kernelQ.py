"""Extract the E1 kernel mod J at degree d (parity parD) in canonical reduced form over Q via CRT+ratrecon."""
import sys, json
import lib
from lib import *
from flint import nmod_mat

def canon_kernel(d, parD, p):
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    parT = (parD + 1) % 2
    rows_mons = core.monomials(d + 6, parT); idx = {m: k for k, m in enumerate(rows_mons)}
    cols = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(d + 1, (parT - pi) % 2):
            cols.append(core.shift(Fi[i], m, 1, p))
    nJ = len(cols)
    Dmons = core.monomials(d, parD); nD = len(Dmons)
    for m in Dmons: cols.append(core.shift(Ap, m, 1, p))
    R, C = len(rows_mons), len(cols)
    flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[idx[e] * C + j] = v
    X, nullity = nmod_mat(R, C, flat, p).nullspace()
    proj = nmod_mat(nD, nullity, [int(X[nJ + j, k]) for j in range(nD) for k in range(nullity)], p)
    E1, rkE = proj.transpose().rref()   # rows = E1-space basis vectors in D-coords
    # J_d basis (parity parD)
    colsJ = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(d - 5, (parD - pi) % 2):
            colsJ.append(core.shift(Fi[i], m, 1, p))
    idxD = {m: k for k, m in enumerate(Dmons)}
    flatJ = [0] * (len(colsJ) * nD)
    for j, c in enumerate(colsJ):
        for e, v in c.items(): flatJ[j * nD + idxD[e]] = v
    Jm = nmod_mat(len(colsJ), nD, flatJ, p)
    Jr, rkJ = Jm.rref()
    # stack [Jr ; E1] and rref: rows beyond rkJ give canonical complement
    rowsAll = [[int(Jr[i, j]) for j in range(nD)] for i in range(rkJ)] + [[int(E1[i, j]) for j in range(nD)] for i in range(rkE)]
    St = nmod_mat(len(rowsAll), nD, [v for r in rowsAll for v in r], p)
    Sr, rkS = St.rref()
    assert rkS == rkE, (rkS, rkE, rkJ)
    # canonical complement: rref rows of the full E1 space whose pivot is not a pivot of Jr
    pivJ = set()
    for i in range(rkJ):
        j = 0
        while j < nD and int(Jr[i, j]) == 0: j += 1
        pivJ.add(j)
    ker = []
    for i in range(rkS):
        j = 0
        while j < nD and int(Sr[i, j]) == 0: j += 1
        if j not in pivJ:
            ker.append([int(Sr[i, jj]) for jj in range(nD)])
    # reduce kernel rows further modulo J rows (eliminate J pivot columns)
    Jrows = [[int(Jr[i, j]) for j in range(nD)] for i in range(rkJ)]
    pivJlist = sorted(pivJ)
    red = []
    for v in ker:
        v = list(v)
        for i, pj in enumerate(sorted(pivJ)):
            # find the J row with pivot pj
            pass
        red.append(v)
    return Dmons, ker, rkJ, rkE

if __name__ == '__main__':
    d = int(sys.argv[1]); parD = int(sys.argv[2])
    p1, p2 = core.PRIMES[0], core.PRIMES[1]
    Dmons, k1, rkJ1, rkE1 = canon_kernel(d, parD, p1)
    Dmons, k2, rkJ2, rkE2 = canon_kernel(d, parD, p2)
    print('kernel dims', len(k1), len(k2), 'J', rkJ1, rkJ2, 'E1', rkE1, rkE2)
    out = []
    for v1, v2 in zip(k1, k2):
        poly = {}
        for j, m in enumerate(Dmons):
            a, mmod = core.crt_pair(v1[j], p1, v2[j], p2)
            r = core.ratrecon(a, mmod)
            if r is None:
                print('ratrecon failed at', m); r = None
            if r: poly[m] = r
        out.append(poly)
        print('---', core.to_str(poly)[:3000])
    with open(f'{WORK}/kernelQ_d{d}_par{parD}.json', 'w') as fh:
        json.dump([core.to_str(pp) for pp in out], fh)
