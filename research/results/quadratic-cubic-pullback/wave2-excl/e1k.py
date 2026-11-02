"""E1 annihilator kernel K_d^{par} = {D in S_d^{par} : A*D in J(F)} mod p, for quadratic-cubic-pullback.
Usage: python3 e1k.py d parity prime [--basis]
Receipt: e1k_d{d}_par{par}_p{p}.json ; basis (if nontrivial mod J and --basis): e1k_d{d}_par{par}_p{p}_basis.json
Rigor: rank_p(J-block) = Hilbert value  =>  rank_Q(J-block) = that value;  dimK_p >= dimK_Q >= dim J_d (always);
       so dimK_p == dim J_d  =>  K_Q = J_d exactly."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl')
from qx import *
d = int(sys.argv[1]); parD = int(sys.argv[2]); p = int(sys.argv[3]); want_basis = '--basis' in sys.argv
F, A, q = load(); S, R = hilbert()
t0 = time.time()
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
parT = (parD + 1) % 2                      # A is Y-odd
Dmons = monomials(d, parD); nD = len(Dmons)
cols = [shift(Ap, m, 1, p) for m in Dmons]
Jinfo = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(d + 1, (parT - pi) % 2):
        cols.append(shift(Fi[i], m, 1, p)); Jinfo.append((i, m))
rowsl = monomials(d + 6, parT)
rows = {e: j for j, e in enumerate(rowsl)}
Rn, C = len(rowsl), len(cols)
flat = [0] * (Rn * C)
nnz = 0
for j, col in enumerate(cols):
    for e, v in col.items():
        flat[rows[e] * C + j] = v; nnz += 1
M = nmod_mat(Rn, C, flat, p); del flat
t1 = time.time()
print(f'[d={d} par={parD} p={p}] matrix {Rn} x {C} (nD={nD}, nJ={C-nD}) nnz={nnz} built {t1-t0:.1f}s', flush=True)
# J-block rank (must equal Hilbert value dim J_{d+6}^{parT})
MJ = nmod_mat(Rn, C - nD, [0] * (Rn * (C - nD)), p)
for j in range(nD, C):
    pass
# cheaper: build J-block separately
flatJ = [0] * (Rn * (C - nD))
for j, col in enumerate(cols[nD:]):
    for e, v in col.items():
        flatJ[rows[e] * (C - nD) + j] = v
MJ = nmod_mat(Rn, C - nD, flatJ, p); del flatJ
rJ = MJ.rank(); del MJ
hilbJ = S[d + 6][parT] - R[d + 6][parT]
t2 = time.time()
print(f'rank J-block = {rJ} (Hilbert dim J_{d+6}^{parT} = {hilbJ}) {t2-t1:.1f}s', flush=True)
X, nullity = M.nullspace()
rank = C - nullity
t3 = time.time()
print(f'rank stacked = {rank} nullity = {nullity} ({t3-t2:.1f}s)', flush=True)
dimK = nD - (rank - rJ)
Jd = S[d][parD] - R[d][parD]
# projection of the nullspace to the D block: basis of K_full
basis = []
if nullity > 0 and want_basis:
    proj = nmod_mat(nullity, nD, [int(X[j, k]) for k in range(nullity) for j in range(nD)], p)
    Pr, rk = proj.rref()
    assert rk == dimK, (rk, dimK)
    for i in range(rk):
        poly = {}
        for j in range(nD):
            v = int(Pr[i, j])
            if v: poly[Dmons[j]] = v
        basis.append(poly)
rec = {'exp': 'e1kernel', 'd': d, 'parD': parD, 'p': p, 'rows': Rn, 'cols': C, 'nD': nD, 'nnz': nnz,
       'rank_J_block': rJ, 'hilbert_J_target': hilbJ, 'J_rank_attained': rJ == hilbJ,
       'rank_stacked': rank, 'nullity': nullity, 'dimK_p': dimK, 'dim_J_d': Jd, 'kernel_mod_J_p': dimK - Jd,
       'gorenstein_lower_bound': max(0, R[d][parD] - R[d + 6][parT]),
       'rigorous': ('K_Q = J_d exactly' if (rJ == hilbJ and dimK == Jd) else
                    ('dim(K_Q/J_d) = %d exactly (upper bound mod p meets Gorenstein lower bound)' % (dimK - Jd)
                     if (rJ == hilbJ and dimK - Jd == max(0, R[d][parD] - R[d + 6][parT])) else
                     'dim(K_Q/J_d) <= %d (mod-p upper bound only)' % (dimK - Jd))),
       'seconds': round(time.time() - t0, 1)}
print(json.dumps(rec), flush=True)
log(rec)
receipt(f'e1k_d{d}_par{parD}_p{p}.json', rec)
if basis:
    # also save E1 numerator representatives (full nullspace vectors) for non-J elements: reduce K basis mod J
    red = jacobian._Reducer(F, d, p, parD)
    nonJ = []
    for b in basis:
        nf = red.normal_form(b)
        if nf: nonJ.append(b)
    with open(f'{PK}/e1k_d{d}_par{parD}_p{p}_basis.json', 'w') as f:
        json.dump({'p': p, 'd': d, 'parD': parD, 'dimK': dimK, 'dimJ': Jd,
                   'basis': [{' '.join(map(str, k)): v for k, v in b.items()} for b in basis]}, f)
    # nullspace columns give (D, N) pairs with N(F) = A D exactly (mod p): save all columns compactly
    Nm = Jinfo
    vecs = []
    for k in range(nullity):
        Dv = {Dmons[j]: int(X[j, k]) for j in range(nD) if int(X[j, k])}
        if not Dv: continue
        Nv = [{} for _ in range(NV)]
        for t, (i, m) in enumerate(Nm):
            c = int(X[nD + t, k])
            if c: Nv[i][m] = c
        vecs.append({'D': {' '.join(map(str, kk)): v for kk, v in Dv.items()},
                     'N': [{' '.join(map(str, kk)): v for kk, v in Ni.items()} for Ni in Nv]})
    with open(f'{PK}/e1k_d{d}_par{parD}_p{p}_pairs.json', 'w') as f:
        json.dump({'p': p, 'd': d, 'parD': parD, 'pairs': vecs}, f)
    print(f'saved basis ({len(basis)}; {len(nonJ)} not in J) and {len(vecs)} (D,N) nullspace pairs', flush=True)
