import sys, time
from flint import nmod_mat
import numpy as np
from refpoly import *
p = int(sys.argv[1]); kos = [int(x) for x in sys.argv[2].split(',')]
F, A, q = load(p); Fi = [pderiv(F, i, p) for i in range(6)]
def slot_par(s, par): return (par + (1 if s >= 3 else 0)) % 2
def mat(cols_polys, rowmons):
    idx = {m: i for i, m in enumerate(rowmons)}
    M = np.zeros((len(rowmons), len(cols_polys)), dtype=np.int64)
    for j, P in enumerate(cols_polys):
        for e, c in P.items(): M[idx[e], j] = c
    return nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
for d in kos:
    for par in (0, 1):
        n = d + 6; t1 = time.time()
        rowsets = [monomials(n-5, slot_par(s, par)) for s in range(6)]
        offs = np.cumsum([0] + [len(r) for r in rowsets]); idx = [{m: i for i, m in enumerate(r)} for r in rowsets]
        cols = [(i, j, m) for i in range(6) for j in range(i+1, 6) for m in monomials(n-10, (slot_par(j, par) + (1 if i >= 3 else 0)) % 2)]
        M = np.zeros((offs[-1], len(cols)), dtype=np.int64)
        for c, (i, j, m) in enumerate(cols):
            for e, v in mulmono(Fi[i], m, 1, p).items(): M[offs[j] + idx[j][e], c] = v
            for e, v in mulmono(Fi[j], m, 1, p).items(): M[offs[i] + idx[i][e], c] = (-v) % p
        r2 = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p).rank()
        rJ = mat([mulmono(Fi[s], m, 1, p) for s in range(6) for m in monomials(n-5, slot_par(s, par))], monomials(n, par)).rank()
        Z1 = offs[-1] - rJ
        print(f"p={p} Koszul degree n={n} parity {par}: rank d2 = {r2}, dim Z1 = {offs[-1]} - rank J_n {rJ} = {Z1} -> H_1 = 0: {r2 == Z1} (rank J_n CI {dimJ_CI(n,par)}) [{time.time()-t1:.0f}s]", flush=True)
print("total done")
