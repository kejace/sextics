#!/usr/bin/env python3
"""E1-kernel (direct Macaulay) with a chosen prime: K_full = {D in S_d^{parD} : A D in J_{d+6}} via the nullspace of
[A*S_d^{parD} | S_{d+1}*grad F] (rows S_{d+6}^{parD+1}); kernel mod J = dim K_full - dim J_d^{parD}.
Rigor: K_full(Q) reduces injectively into K_full(p) (saturated lattice), and dim J_d(p) = dim J_d(Q) when the modular
count attains the Koszul value; so kernel_mod_J(p) = 0 proves (J:A)_d = J_d over Q (and C). For d >= 10 the kernel is
at least dim R_d - dim R_{d+6} over Q, so a modular value equal to that bound is exact too.
Usage: python3 wave2_excl_e1kernel.py d parD p [--basis]   -> wave2-excl/e1kernel_d{d}_par{parD}_p{p}.json (+ basis/reps)
"""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
d = int(sys.argv[1]); parD = int(sys.argv[2]); p = int(sys.argv[3]); want_basis = '--basis' in sys.argv
F, A, q = load(); S, Rh = hilbert_counts()
t0 = time.time()
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
parT = (parD + 1) % 2
Dmons = monomials(d, parD)
cols = [shift(Ap, m, 1, p) for m in Dmons]
nD = len(cols)
Nmons = []
for i in range(NV):
    for m in monomials(d + 1, (parT - par_i(i)) % 2):
        cols.append(shift(Fi[i], m, 1, p)); Nmons.append((i, m))
rows = {e: j for j, e in enumerate(monomials(d + 6, parT))}
R, C = len(rows), len(cols)
M = nmod_mat(R, C, p)
for j, col in enumerate(cols):
    for e, v in col.items(): M[rows[e], j] = v
print(f'[d={d} par={parD} p={p}] matrix {R} x {C} built in {time.time()-t0:.1f}s', flush=True)
X, nullity = M.nullspace()
rank = C - nullity
proj = nmod_mat(nullity, nD, [int(X[j, k]) for k in range(nullity) for j in range(nD)], p)
Tr, rk = proj.rref()
Jd = S[d][parD] - Rh[d][parD]          # exact dim J_d^{parD} (Koszul count via the Jacobian-ring Hilbert series)
Jgens = []
for i in range(NV):
    for m in monomials(d - 5, (parD - par_i(i)) % 2): Jgens.append(shift(Fi[i], m, 1, p))
Jidx = {m: j for j, m in enumerate(Dmons)}
JM = nmod_mat(max(len(Jgens), 1), nD, p)
for r, c in enumerate(Jgens):
    for mm, v in c.items(): JM[r, Jidx[mm]] = v
rJ = JM.rank() if Jgens else 0
forced = max(0, Rh[d][parD] - Rh[d + 6][parT])
rec = {'exp': 'e1-kernel-direct', 'd': d, 'parD': parD, 'prime': p, 'rows': R, 'cols': C, 'nD': nD, 'rank': rank,
       'nullity': nullity, 'dim_Kfull_mod_p': rk, 'dim_J_d_exact': Jd, 'dim_J_d_mod_p': rJ, 'kernel_mod_J': rk - Jd,
       'forced_kernel_lower_bound_Q': forced, 'dim_R_d': Rh[d][parD], 'dim_R_d6': Rh[d + 6][parT],
       'rigorous': (rJ == Jd) and (rk - Jd == forced), 'seconds': round(time.time() - t0, 1),
       'result': ('(J:A)_d = J_d over Q' if rk - Jd == 0 and rJ == Jd else f'kernel mod J = {rk-Jd} (exact over Q iff equals forced bound {forced})')}
print(json.dumps(rec), flush=True)
log2(**rec)
with open(f'{OUT}/e1kernel_d{d}_par{parD}_p{p}.json', 'w') as f: json.dump(rec, f)
if want_basis and rk - Jd > 0:
    basis = []
    for i in range(rk):
        poly = {}
        for j in range(nD):
            v = int(Tr[i, j])
            if v: poly[Dmons[j]] = v
        basis.append(poly)
    with open(f'{OUT}/e1kernel_d{d}_par{parD}_p{p}_basis.json', 'w') as f:
        json.dump({'p': p, 'd': d, 'parD': parD, 'dim_Kfull': rk, 'dim_J': Jd,
                   'basis': [{' '.join(map(str, k)): v for k, v in b.items()} for b in basis]}, f)
    # representatives (D, N) with N(F) = A D from nullspace columns whose D-block is not in J
    def notinJ(Dv):
        M2 = nmod_mat(len(Jgens) + 1, nD, p)
        for r, c in enumerate(Jgens):
            for mm, v in c.items(): M2[r, Jidx[mm]] = v
        for j in range(nD): M2[len(Jgens), j] = Dv[j]
        return M2.rank() > rJ
    random.seed(d * 100 + parD)
    reps = []
    cols_nz = [k for k in range(nullity) if any(int(X[j, k]) for j in range(nD))]
    random.shuffle(cols_nz)
    for k in cols_nz:
        v = [int(X[j, k]) for j in range(C)]
        if notinJ(v[:nD]):
            D = {Dmons[j]: v[j] for j in range(nD) if v[j]}
            N = [{} for _ in range(NV)]
            for t, (i, m) in enumerate(Nmons):
                if v[nD + t]: N[i][m] = v[nD + t]
            reps.append({'kind': f'nullspace col {k}', 'D': D, 'N': N})
        if len(reps) >= 6: break
    for t in range(4):
        coef = [random.randrange(p) for _ in range(nullity)]
        v = [sum(coef[k] * int(X[j, k]) for k in range(nullity)) % p for j in range(C)]
        D = {Dmons[j]: v[j] for j in range(nD) if v[j]}
        N = [{} for _ in range(NV)]
        for tt, (i, m) in enumerate(Nmons):
            if v[nD + tt]: N[i][m] = v[nD + tt]
        reps.append({'kind': f'random combo {t}', 'D': D, 'N': N})
    with open(f'{OUT}/e1kernel_d{d}_par{parD}_p{p}_reps.json', 'w') as f:
        json.dump({'p': p, 'd': d, 'parD': parD, 'reps': [{'kind': r['kind'], 'D': {' '.join(map(str, k)): v for k, v in r['D'].items()},
                   'N': [{' '.join(map(str, k)): v for k, v in Ni.items()} for Ni in r['N']]} for r in reps]}, f)
    print(f'saved basis ({rk} vectors) and {len(reps)} E1 representatives', flush=True)
