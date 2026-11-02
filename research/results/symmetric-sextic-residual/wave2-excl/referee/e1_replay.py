# Own replay at a fresh prime: Jacobian ranks vs CI, E1 kernels {D : A D in J}/J_d, A notin J, and Koszul H_1 = 0 in degree d+6.
import sys, time
from flint import nmod_mat
import numpy as np
from refpoly import *
p = int(sys.argv[1]); dlist = [int(x) for x in sys.argv[2].split(',')]; kos = [int(x) for x in sys.argv[3].split(',')]
F, A, q = load(p)
Fi = [pderiv(F, i, p) for i in range(6)]
def slot_par(s, par): return (par + (1 if s >= 3 else 0)) % 2     # parity of m with m*F_s of parity par
def mat(cols_polys, rowmons):
    idx = {m: i for i, m in enumerate(rowmons)}
    M = np.zeros((len(rowmons), len(cols_polys)), dtype=np.int64)
    for j, P in enumerate(cols_polys):
        for e, c in P.items(): M[idx[e], j] = c
    return nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
def Jcols(n, par):
    return [mulmono(Fi[s], m, 1, p) for s in range(6) for m in monomials(n-5, slot_par(s, par))]
t0 = time.time()
# A notin J (degree 6, Y-odd)
rows6 = monomials(6, 1); J6 = Jcols(6, 1)
rJ = mat(J6, rows6).rank(); rJA = mat(J6 + [A], rows6).rank()
print(f"p={p}: rank J_6^odd = {rJ} (CI {dimJ_CI(6,1)}), rank [J_6^odd | A] = {rJA} -> A notin J: {rJA == rJ + 1}", flush=True)
for d in dlist:
    for pD in (0, 1):
        n = d + 6; par = (1 + pD) % 2
        rows = monomials(n, par); Jc = Jcols(n, par)
        Ac = [mulmono(A, m, 1, p) for m in monomials(d, pD)]
        t1 = time.time()
        rJ = mat(Jc, rows).rank(); rAJ = mat(Ac + Jc, rows).rank()
        dimK = len(Ac) - (rAJ - rJ)
        dimJd = mat([mulmono(Fi[s], m, 1, p) for s in range(6) for m in monomials(d-5, slot_par(s, pD))], monomials(d, pD)).rank() if d >= 5 else 0
        print(f"d={d} pD={pD}: rank J_{n}^{par} = {rJ} (CI {dimJ_CI(n,par)}) {'OK' if rJ==dimJ_CI(n,par) else 'MISMATCH'}; dimK={dimK} dimJ_d={dimJd} (CI {dimJ_CI(d,pD) if d>=5 else 0}) ker_R={dimK-dimJd}  [{time.time()-t1:.0f}s]", flush=True)
for d in kos:
    for par in (0, 1):
        n = d + 6
        # Koszul d2: columns (i<j, m): slot j gets m F_i, slot i gets -m F_j ; rows: slot s monomials of degree n-5 parity slot_par(s,par)
        rowsets = [monomials(n-5, slot_par(s, par)) for s in range(6)]
        offs = np.cumsum([0] + [len(r) for r in rowsets]); idx = [{m: i for i, m in enumerate(r)} for r in rowsets]
        cols = []
        for i in range(6):
            for j in range(i+1, 6):
                mpar = (slot_par(j, par) + (1 if i >= 3 else 0)) % 2
                for m in monomials(n-10, mpar):
                    cols.append((i, j, m))
        M = np.zeros((offs[-1], len(cols)), dtype=np.int64)
        for c, (i, j, m) in enumerate(cols):
            for e, v in mulmono(Fi[i], m, 1, p).items(): M[offs[j] + idx[j][e], c] = v
            for e, v in mulmono(Fi[j], m, 1, p).items(): M[offs[i] + idx[i][e], c] = (-v) % p
        r2 = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p).rank()
        rJ = mat(Jcols(n, par), monomials(n, par)).rank()
        Z1 = offs[-1] - rJ
        print(f"Koszul degree n={n} parity {par}: rank d2 = {r2}, dim Z1 = 6-slot chains {offs[-1]} - rank J_n {rJ} = {Z1} -> H_1 = 0: {r2 == Z1} (rank J_n CI {dimJ_CI(n,par)})", flush=True)
print(f"total {time.time()-t0:.0f}s")
