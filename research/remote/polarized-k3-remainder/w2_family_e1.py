"""Exact (mod p) dimensions of  W ∩ J(F)_d  for structured linear families W of denominators (cones etc.).
At d <= 9 every E1-admissible D lies in J(F)_d, so dim(W ∩ J_d) = 0 excludes the whole family W at once."""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np
import w2lib as L
from hodge import core
from hodge.core import monomials, mul, add, sub, shift, pderiv
from fractions import Fraction as Fr
from flint import nmod_mat
p = L.P0
F, A, q = L.load()
d = int(sys.argv[1]) if len(sys.argv) > 1 else 8
out = f'w2-family-e1-d{d}.jsonl'
def vecs_of(polys, S):
    return [S.vec(P, dtype=np.int64, p=p) for P in polys]
def inter_dim(Wpolys, pD):
    S = L.space(d, pD)
    Jb = L.jac_basis(d, pD)
    JV = vecs_of([Dk for _, _, Dk in Jb], S)
    WV = vecs_of(Wpolys, S)
    WV = [v for v in WV if v.any()]
    if not WV: return 0, 0, len(JV)
    M = np.array(WV + JV, dtype=np.int64)            # rows = vectors
    r_all = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p).rank()
    r_W = nmod_mat(len(WV), M.shape[1], np.array(WV).ravel().tolist(), p).rank()
    r_J = len(JV)
    return r_W + r_J - r_all, r_W, r_J
def Xonly(vars_, deg, pD):
    """all monomials of degree deg in the given variable indices with Y-parity pD"""
    out = []
    for m in monomials(deg, pD):
        if all(m[i] == 0 for i in range(6) if i not in vars_):
            out.append({m: Fr(1)})
    return out
fams = []
fams.append(('cone vertex plane {X=0}: D in C[X0,X1,X2]', 0, Xonly([0, 1, 2], d, 0)))
fams.append(('cone vertex plane {Y=0}: D in C[Y0,Y1,Y2] (Y-even part)', 0, Xonly([3, 4, 5], d, 0)))
fams.append(('cone vertex plane {Y=0}: D in C[Y0,Y1,Y2] (Y-odd part)', 1, Xonly([3, 4, 5], d, 1)))
for vs in [(0, 1, 3, 4), (0, 1, 2, 3), (0, 3, 4, 5), (0, 2, 3, 5), (0, 1, 3, 5)]:
    names = ','.join(core.VARS[i] for i in vs)
    for pD in (0, 1):
        fams.append((f'cone vertex line: D in C[{names}] (parity {pD})', pD, Xonly(list(vs), d, pD)))
# UV-frame cones: h(U) +/- h(V) over plane curves; and C[U0,U1,V0,V1] sym/anti
def uvmon(e):  # e: 6-tuple in UV slots
    return L.uv_to_xy({e: Fr(1)})
Umons = [m for m in monomials(d, None) if m[3] == m[4] == m[5] == 0]
hU = [uvmon(m) for m in Umons]; hV = [L.uv_to_xy({(0, 0, 0) + m[:3]: Fr(1)}) for m in Umons]
fams.append(('UV cone: h(U) + h(V), h ternary form (Y-even)', 0, [add(a, b) for a, b in zip(hU, hV)]))
fams.append(('UV cone: h(U) - h(V), h ternary form (Y-odd)', 1, [sub(a, b) for a, b in zip(hU, hV)]))
UV4 = [m for m in monomials(d, None) if m[2] == 0 and m[5] == 0]
sym4 = [add(uvmon(m), uvmon((m[3], m[4], m[5], m[0], m[1], m[2]))) for m in UV4]
anti4 = [sub(uvmon(m), uvmon((m[3], m[4], m[5], m[0], m[1], m[2]))) for m in UV4]
fams.append(('UV cone vertex line: D in C[U0,U1,V0,V1] swap-symmetric', 0, sym4))
fams.append(('UV cone vertex line: D in C[U0,U1,V0,V1] swap-antisymmetric', 1, anti4))
# pure (X,Y)-bidegree pieces of R_d (Y-degree k)
for k in range(0, d + 1):
    pD = k % 2
    fams.append((f'(X,Y)-bidegree ({d-k},{k})', pD, [{m: Fr(1)} for m in monomials(d, pD) if m[3] + m[4] + m[5] == k]))
# (U,V)-bidegree sym/anti pieces
for a in range(0, d // 2 + 1):
    b = d - a
    ms = [m for m in monomials(d, None) if m[0] + m[1] + m[2] == a]
    sy = [add(uvmon(m), uvmon((m[3], m[4], m[5], m[0], m[1], m[2]))) for m in ms]
    an = [sub(uvmon(m), uvmon((m[3], m[4], m[5], m[0], m[1], m[2]))) for m in ms]
    fams.append((f'(U,V)-bidegree ({a},{b})+({b},{a}) swap-symmetric', 0, sy))
    fams.append((f'(U,V)-bidegree ({a},{b})-({b},{a}) swap-antisymmetric', 1, an))
for name, pD, W in fams:
    t0 = time.time()
    try:
        k, rW, rJ = inter_dim(W, pD)
        rec = {'experiment': 'family-cap-J', 'd': d, 'family': name, 'pD': pD, 'dim_W': rW, 'dim_J_d': rJ, 'dim_W_cap_J': k, 'seconds': round(time.time() - t0, 1)}
    except Exception as ex:
        rec = {'experiment': 'family-cap-J', 'd': d, 'family': name, 'pD': pD, 'error': str(ex)}
    L.log_line(out, rec); print(json.dumps(rec), flush=True)
