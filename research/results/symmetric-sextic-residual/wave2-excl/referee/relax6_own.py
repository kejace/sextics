# Own reconstruction of the degree-6 Hessian relaxation (no lane-library code): D = w.grad F, w linear (36 params),
# N0 = A w, R0 = D div N0 - N0.grad D must lie in J(F)_12.  Complement of the radial line w = z: a_00 = 0 (35 params).
# Quadrics = coordinates of R0 in R_12 = S_12/J_12 (both parities). Writes an msolve file (power syntax ^).
import sys, time
from flint import nmod_mat
import numpy as np
from refpoly import *
p = int(sys.argv[1]); out = sys.argv[2]
F, A, q = load(p)
Fi = [pderiv(F, i, p) for i in range(6)]
def z(j): e=[0]*6; e[j]=1; return tuple(e)
params = [(i, j) for i in range(6) for j in range(6) if (i, j) != (0, 0)]   # complement of the radial direction
# building blocks
zF = {(i, j): mulmono(Fi[i], z(j), 1, p) for i in range(6) for j in range(6)}          # D-term of a_ij : z_j F_i
Az = {j: mulmono(A, z(j), 1, p) for j in range(6)}                                     # N0 component: a_kl -> A z_l in slot k
dAz = {(k, l): pderiv(Az[l], k, p) for k in range(6) for l in range(6)}                # div(A z_l e_k) = d_k(A z_l)
dzF = {(k, i, j): pderiv(zF[(i, j)], k, p) for k in range(6) for i in range(6) for j in range(6)}
t0 = time.time()
def R(u, v):   # R0 contribution a_u a_v with u=(i,j) the D-index, v=(k,l) the N0-index
    (i, j), (k, l) = u, v
    return add(mul(zF[(i, j)], dAz[(k, l)], p), mul(Az[l], dzF[(k, i, j)], p), p, s=-1)
quads = {}   # key (u_index<=v_index) -> polynomial
n = len(params)
for a in range(n):
    for b in range(a, n):
        P = R(params[a], params[b]) if a == b else add(R(params[a], params[b]), R(params[b], params[a]), p)
        if P: quads[(a, b)] = P
print(f"{len(quads)} nonzero bilinear polynomials of degree 12 built in {time.time()-t0:.0f}s", flush=True)
# reducers mod J_12 per parity
def slot_par(s, par): return (par + (1 if s >= 3 else 0)) % 2
coords_all = {}
for par in (0, 1):
    rows = monomials(12, par); idx = {m: i for i, m in enumerate(rows)}
    gens = [mulmono(Fi[s], m, 1, p) for s in range(6) for m in monomials(7, slot_par(s, par))]
    G = np.zeros((len(gens), len(rows)), dtype=np.int64)
    for r, P in enumerate(gens):
        for e, c in P.items(): G[r, idx[e]] = c
    Gm = nmod_mat(G.shape[0], G.shape[1], G.ravel().tolist(), p)
    Gr, rk = Gm.rref()
    assert rk == dimJ_CI(12, par), (rk, dimJ_CI(12, par))
    Grr = np.array([[int(Gr[i, j]) for j in range(len(rows))] for i in range(rk)], dtype=np.int64)
    piv = [int(np.nonzero(Grr[i])[0][0]) for i in range(rk)]
    nonpiv = [j for j in range(len(rows)) if j not in set(piv)]
    print(f"parity {par}: rank J_12 = {rk} = CI; dim R_12 = {len(nonpiv)} (CI {dimR(12, par)})", flush=True)
    keys = list(quads.keys())
    V = np.zeros((len(keys), len(rows)), dtype=np.int64)
    for r, k in enumerate(keys):
        for e, c in quads[k].items():
            if ypar(e) == par: V[r, idx[e]] = c
    Vm = nmod_mat(V.shape[0], V.shape[1], V.ravel().tolist(), p)
    Pm = nmod_mat(V.shape[0], rk, V[:, piv].ravel().tolist(), p)
    Red = Vm - Pm * nmod_mat(rk, len(rows), Grr.ravel().tolist(), p)
    for r, k in enumerate(keys):
        for cidx, j in enumerate(nonpiv):
            v = int(Red[r, j])
            if v: coords_all.setdefault((par, cidx), {})[k] = v
    # sanity: pivot columns are zero after reduction
    assert all(int(Red[r, j]) == 0 for r in range(min(5, len(keys))) for j in piv)
print(f"{len(coords_all)} nonzero quadrics in {n} unknowns ({time.time()-t0:.0f}s)", flush=True)
with open(out, 'w') as f:
    f.write(','.join(f'a{i+1}' for i in range(n)) + '\n' + str(p) + '\n')
    lines = []
    for key, Qd in coords_all.items():
        terms = []
        for (a, b), c in Qd.items():
            terms.append(f"{c}*a{a+1}^2" if a == b else f"{c}*a{a+1}*a{b+1}")
        lines.append('+'.join(terms))
    f.write(',\n'.join(lines) + '\n')
print("wrote", out, flush=True)
