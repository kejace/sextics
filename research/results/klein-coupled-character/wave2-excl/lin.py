"""(b) Full linearization (Shor relaxation) of the bilinear E2 system over D in J_d (E1-kernel zero), by blocks.

D = sum_a u_a D_a, D_a = m F_i (basis V of J_d via g), N0 = sum_a u_a A m e_i.  Every E1 lift is N0 + K(c), c Koszul
coefficients (deg d-4).  E2:  sum_{a<=b} T_ab C_ab + sum_{a,k} Y_ak L_ak - sum_m Q_m F m = 0  with the rank-one
substitution T = u(x)u (C_aa = H_aa, C_ab = 2 H_ab), Y = u(x)c.  Relaxation: treat (T, Y, Q) as independent unknowns
of a LINEAR system (block diagonal in (parity, sigma-character)); W_T = projection of its solution space onto T.
Gauge tensors (hZ) o v (h in S_{d-6}) lie in W_T.  Decisive iff dim W_T = dim(G o V) in every block.
(Projection dimensions mod p are evidence for the Q statement unless rank_Q M_YQ = rank_p M_YQ is certified.)
usage: python3 lin.py <p> <d> <parity 0|1>
"""
import sys, json, time, math
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
import kx
from hodge import core
from hodge.core import NV, mul, add, sub, scale, shift, pderiv, parse
from flint import nmod_mat
p = int(sys.argv[1]); d = int(sys.argv[2]); parD = int(sys.argv[3])
F, A, q = kx.load()
F7 = scale(F, 7); Fi = [pderiv(F7, i) for i in range(NV)]; Ai = [pderiv(A, i) for i in range(NV)]
Z = [parse(v) for v in ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2')]
t0 = time.time()
V = [(i, m) for i in range(NV) for m in kx.wmonomials(d - 5, (parD + kx.PI[i]) % 2)]
n = len(V)
Dv = [shift(Fi[i], m, 1) for (i, m) in V]           # D_a
Dvi = [[pderiv(D, j) for j in range(NV)] for D in Dv]
def B(D, Di, N):
    """B_D(N) = D div N - N(D)."""
    acc = {}
    for i in range(NV):
        if N[i]:
            acc = add(acc, mul(D, pderiv(N[i], i)))
            acc = sub(acc, mul(N[i], Di[i]))
    return acc
def N0(a):
    i, m = V[a]
    N = [{} for _ in range(NV)]; N[i] = shift(A, m, 1); return N
# C_ab: H_aa = B(D_a, N_a); C_ab = B(D_a,N_b) + B(D_b,N_a)  (= 2 H_ab)
Ccols = {}
Nv = [N0(a) for a in range(n)]
for a in range(n):
    Ccols[(a, a)] = B(Dv[a], Dvi[a], Nv[a])
    for b in range(a + 1, n):
        Ccols[(a, b)] = add(B(Dv[a], Dvi[a], Nv[b]), B(Dv[b], Dvi[b], Nv[a]))
print(f'd={d} parD={parD} n={n}: {len(Ccols)} T-columns built [{round(time.time()-t0,1)}s]', flush=True)
# Koszul basis
PAIRS = [(i, j) for i in range(NV) for j in range(i + 1, NV)]
K = [(i, j, c) for (i, j) in PAIRS for c in kx.wmonomials(d - 4, (1 + parD + kx.PI[i] + kx.PI[j]) % 2)] if d >= 4 else []
Lcols = {}
for a in range(n):
    D, Di = Dv[a], Dvi[a]
    for kidx, (i, j, c) in enumerate(K):
        # L = D[(d_i c)F_j - (d_j c)F_i] - c(D_i F_j - D_j F_i)
        G = sub(mul(Di[i], Fi[j]), mul(Di[j], Fi[i]))
        col = scale(shift(G, c, 1), -1)
        if c[i]:
            c2 = list(c); c2[i] -= 1; col = add(col, shift(mul(D, Fi[j]), tuple(c2), c[i]))
        if c[j]:
            c2 = list(c); c2[j] -= 1; col = sub(col, shift(mul(D, Fi[i]), tuple(c2), c[j]))
        if col: Lcols[(a, kidx)] = col
print(f'{len(K)} Koszul monomials, {len(Lcols)} nonzero Y-columns [{round(time.time()-t0,1)}s]', flush=True)
Qm = kx.wmonomials(2 * d - 6, 1) if 2 * d - 6 >= 0 else []
Qcols = {m: scale(shift(F7, m, 1), -1) for m in Qm}
# gauge tensors (correct scaling: T_aa = u_a v_a, T_ab = (u_a v_b + u_b v_a)/2) in the C-column convention
idx = {(i, m): a for a, (i, m) in enumerate(V)}
gauge_vecs = []
for h in core.monomials(d - 6) if d >= 6 else []:
    gv = {}
    ok = True
    for i in range(NV):
        for m, c in shift(Z[i], h, 1).items():
            if (i, m) not in idx: ok = False; break
            gv[idx[(i, m)]] = kx.modp(c, p)
        if not ok: break
    if ok: gauge_vecs.append(gv)
inv2 = pow(2, -1, p)
def gauge_tensor(u, v):
    T = {}
    for a, ca in u.items():
        for b, cb in v.items():
            if a == b: T[(a, a)] = (T.get((a, a), 0) + ca * cb) % p
            else:
                key = (min(a, b), max(a, b)); T[key] = (T.get(key, 0) + ca * cb * inv2) % p
    return T
gauge_T = [gauge_tensor(gv, {a: 1}) for gv in gauge_vecs for a in range(n)]
print(f'gauge vectors {len(gauge_vecs)}, gauge tensors {len(gauge_T)} [{round(time.time()-t0,1)}s]', flush=True)
# blocks by (parity, character) of column polynomials
def blk(P):
    return (core.parity_of(P), kx.char_of(P))
blocks = {}
for key, P in Ccols.items():
    if P: blocks.setdefault(blk(P), {'T': [], 'Y': [], 'Q': []})['T'].append(key)
for key, P in Lcols.items():
    blocks.setdefault(blk(P), {'T': [], 'Y': [], 'Q': []})['Y'].append(key)
for key, P in Qcols.items():
    blocks.setdefault(blk(P), {'T': [], 'Y': [], 'Q': []})['Q'].append(key)
zeroT = [key for key, P in Ccols.items() if not P]
total_W = 0; total_G = 0; details = {}
allpairs = list(Ccols.keys())
for (par, ch), bl in sorted(blocks.items()):
    tcols = bl['T']; ycols = bl['Y']; qcols = bl['Q']
    cols = [Ccols[k] for k in tcols] + [Lcols[k] for k in ycols] + [Qcols[k] for k in qcols]
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    M = kx.dense(cols, rows, p)
    nT = len(tcols)
    X, nul = M.nullspace()
    rank = len(cols) - nul
    # projection onto T coordinates
    if nul and nT:
        proj = nmod_mat(nul, nT, [int(X[r, c]) for c in range(nul) for r in range(nT)], p)
        dimWT = proj.rank()
    else:
        dimWT = 0
    # rank of M_YQ alone (for the rigor caveat)
    MYQ = kx.dense(cols[nT:], rows, p) if len(cols) > nT else None
    rYQ = MYQ.rank() if MYQ is not None else 0
    # gauge tensors in this block
    tset = {k: i for i, k in enumerate(tcols)}
    gT = []
    for T in gauge_T:
        keys = [k for k in T if k in tset]
        if keys:
            assert all(k in tset for k in T if T[k]), 'gauge tensor straddles blocks'
            gT.append(T)
    dimG = 0
    if gT:
        flat = [0] * (len(gT) * nT)
        for r, T in enumerate(gT):
            for k, v in T.items():
                if k in tset: flat[r * nT + tset[k]] = v
        Gm = nmod_mat(len(gT), nT, flat, p); dimG = Gm.rank()
        # containment check: gauge tensors must lie in W_T (span of proj rows)
        if nul:
            comb = nmod_mat(nul + len(gT), nT, [int(X[r, c]) for c in range(nul) for r in range(nT)] + flat, p)
            assert comb.rank() == dimWT, ('gauge tensor not in W_T', comb.rank(), dimWT)
    total_W += dimWT; total_G += dimG
    details[f'{par},{ch}'] = {'rows': len(rows), 'nT': nT, 'nY': len(ycols), 'nQ': len(qcols), 'rank': rank, 'nullity': nul,
                              'dimWT': dimWT, 'dimGauge': dimG, 'rank_YQ': rYQ, 'nullity_YQ': len(cols) - nT - rYQ}
    print(f'  block par={par} ch={ch}: rows {len(rows)}, T {nT}, Y {len(ycols)}, Q {len(qcols)}, rank {rank}, dim W_T {dimWT}, gauge {dimG}, rank_YQ {rYQ} (nullity_YQ {len(cols)-nT-rYQ})  [{round(time.time()-t0,1)}s]', flush=True)
decisive = (total_W == total_G)
print(f'RESULT lin d={d} parD={parD}: dim W_T = {total_W}, dim gauge = {total_G}, zero T-columns {len(zeroT)} -> {"DECISIVE (mod-p projection; see rigor caveat)" if decisive else "NOT decisive"}', flush=True)
res = {'d': d, 'parD': parD, 'p': p, 'n': n, 'nK': len(K), 'nQ': len(Qm), 'dimWT': total_W, 'dimGauge': total_G, 'decisive': decisive,
       'zero_T_columns': len(zeroT), 'blocks': details, 'seconds': round(time.time() - t0, 1)}
json.dump(res, open(f'lin_d{d}_par{parD}_p{p}.json', 'w'), indent=1)
kx.log(exp='full-linearization', **{k: v for k, v in res.items() if k != 'blocks'})
