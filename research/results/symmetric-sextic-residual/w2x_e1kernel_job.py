"""Full E1-kernel sweep with rank bookkeeping and G-isotypic decomposition of the non-J kernel.
usage: python3 w2x_e1kernel_job.py --prime P --degrees 1,2,...  [--tag name]"""
import sys, time, json, argparse, pickle, itertools
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import numpy as np
from flint import nmod_mat
import w2x_lib as L
from hodge import core
ap = argparse.ArgumentParser()
ap.add_argument('--prime', type=int, default=L.P2); ap.add_argument('--degrees', default='1,2,3,4,5,6,7,8')
ap.add_argument('--tag', default='e1kernel')
a = ap.parse_args()
F, A, q = L.load(); p = a.prime
S3 = list(itertools.permutations(range(3)))
# irreps of S3: characters on (identity, transposition, 3-cycle) classes
def cls(pp):
    fixed = sum(1 for i in range(3) if pp[i] == i)
    return 'e' if fixed == 3 else ('t' if fixed == 1 else 'c')
CH = {'triv': {'e': 1, 't': 1, 'c': 1}, 'sgn': {'e': 1, 't': -1, 'c': -1}, 'std': {'e': 2, 't': 0, 'c': -1}}
DIM = {'triv': 1, 'sgn': 1, 'std': 2}
for d in [int(x) for x in a.degrees.split(',')]:
    for pD in (0, 1):
        t0 = time.time()
        fam = L.E1Family(F, A, d, pD, p, verbose=True)
        r, ci = fam.rank_check()
        rec = {'exp': 'e1-kernel', 'd': d, 'pD': pD, 'prime': p, 'rows': fam.rows, 'nA': fam.nA, 'nJ': fam.nJ, 'nullity': fam.nullity,
               'dimK': fam.dimK, 'dimJd': fam.dimJd, 'ci_dimJd': fam.ci_dimJd, 'ker_R': fam.ker_R,
               'rankJtop_modp': r, 'rankJtop_CI': ci, 'rigorous_zero': (fam.ker_R == 0 and r == ci and fam.dimJd == fam.ci_dimJd), 'job': a.tag}
        extras = []
        if fam.ker_R > 0:
            red = L.Reducer(F, d, pD, p)
            Dpolys = [fam.pair(fam.K[i])[0] for i in range(fam.dimK)]
            C = red.coords(Dpolys)
            # basis of the kernel in R_d: rref of C
            Cr, rk = nmod_mat(C.shape[0], C.shape[1], C.ravel().tolist(), p).rref()
            assert rk == fam.ker_R, (rk, fam.ker_R)
            # lift basis vectors back to polynomials in standard monomials
            std = [red.index.mons[j] for j in red.standard]
            basis = []
            for i in range(rk):
                P = {std[j]: int(Cr[i, j]) for j in range(len(std)) if int(Cr[i, j])}
                basis.append(P)
            extras = basis
            # isotypic decomposition
            iso = {}
            for rX in ('triv', 'sgn', 'std'):
                for rY in ('triv', 'sgn', 'std'):
                    projs = []
                    for P in basis:
                        acc = {}
                        for pX in S3:
                            for pY in S3:
                                perm = tuple(pX) + tuple(3 + j for j in pY)
                                c = CH[rX][cls(pX)] * CH[rY][cls(pY)]
                                if c: acc = core.add(acc, core.scale(L.act(P, perm), c % p, p), p)
                        projs.append(acc)
                    Cp = red.coords(projs)
                    rkp = nmod_mat(Cp.shape[0], Cp.shape[1], Cp.ravel().tolist(), p).rank()
                    if rkp: iso[f'{rX}x{rY}'] = {'dim': rkp, 'mult': rkp // (DIM[rX] * DIM[rY])}
            rec['isotypic'] = iso
            pickle.dump({'info': rec, 'basis_modp': basis}, open(f'{L.PK}/{a.tag}_d{d}_p{pD}_{p}.pkl', 'wb'))
        rec['seconds'] = round(time.time() - t0, 1)
        L.log(rec); print(json.dumps(rec), flush=True)
