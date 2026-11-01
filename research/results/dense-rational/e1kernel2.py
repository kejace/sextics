"""Memory-light E1-kernel direct computation (entrywise nmod_mat fill).  Usage: python3 e1kernel2.py d parD [p]
Computes K_full = {D in S_d^{parD} : A D in J}, its dim, kernel mod J, saves a basis (mod p) to e1kernel_d{d}_par{parD}_basis.json"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
d = int(sys.argv[1]); parD = int(sys.argv[2]); p = int(sys.argv[3]) if len(sys.argv) > 3 else P
F, A, q = load()
S, Rh = hilbert_counts()
t0 = time.time()
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
parT = (parD + 1) % 2
Dmons = monomials(d, parD)
cols = [shift(Ap, m, 1, p) for m in Dmons]
nD = len(cols)
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(d + 1, (parT - pi) % 2):
        cols.append(shift(Fi[i], m, 1, p))
rows = {e: j for j, e in enumerate(monomials(d + 6, parT))}
R, C = len(rows), len(cols)
M = nmod_mat(R, C, p)
nnz = 0
for j, col in enumerate(cols):
    for e, v in col.items():
        M[rows[e], j] = v; nnz += 1
print(f'[d={d} par={parD}] matrix {R} x {C}, nnz={nnz}, built in {time.time()-t0:.1f}s', flush=True)
X, nullity = M.nullspace()
rank = C - nullity
t1 = time.time()
print(f'rank={rank} nullity={nullity} ({t1-t0:.1f}s)', flush=True)
proj = nmod_mat(nullity, nD, [int(X[j, k]) for k in range(nullity) for j in range(nD)], p)
Tr, rk = proj.rref()
Jd = S[d][parD] - Rh[d][parD]
basis = []
for i in range(rk):
    poly = {}
    for j in range(nD):
        v = int(Tr[i, j])
        if v: poly[Dmons[j]] = v
    basis.append(poly)
rec = {'mode': 'direct2', 'd': d, 'parD': parD, 'p': p, 'rows': R, 'cols': C, 'rank': rank, 'nullity': nullity,
       'Kfull': rk, 'dimJ': Jd, 'kernel_mod_J': rk - Jd, 'seconds': round(time.time() - t0, 1)}
print(json.dumps(rec), flush=True)
with open('/tmp/hodge/results/dense-rational/e1kernel_sweep.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
with open(f'/tmp/hodge/results/dense-rational/e1kernel_d{d}_par{parD}_basis.json', 'w') as f:
    json.dump({'p': p, 'd': d, 'parD': parD, 'basis': [{' '.join(map(str, k)): v for k, v in b.items()} for b in basis]}, f)
# representatives with E1 numerators: nullspace columns whose D-block is NOT in J (test: D-block not in span of J_d) + random combos
import random
random.seed(d)
Nmons = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(d + 1, (parT - pi) % 2): Nmons.append((i, m))
def vec_to_DN(v):
    D = {}; N = [{} for _ in range(NV)]
    for j in range(nD):
        if v[j]: D[Dmons[j]] = v[j]
    for t, (i, m) in enumerate(Nmons):
        c = v[nD + t]
        if c: N[i][m] = c
    return D, N
Jcols = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(d - 5, (parD - pi) % 2): Jcols.append(shift(Fi[i], m, 1, p))
Jidx = {m: j for j, m in enumerate(Dmons)}
JM = nmod_mat(len(Jcols), nD, p)
for r, c in enumerate(Jcols):
    for mm, v in c.items(): JM[r, Jidx[mm]] = v
rJ = JM.rank()
def notinJ(Dv):
    M2 = nmod_mat(len(Jcols) + 1, nD, p)
    for r, c in enumerate(Jcols):
        for mm, v in c.items(): M2[r, Jidx[mm]] = v
    for j in range(nD): M2[len(Jcols), j] = Dv[j]
    return M2.rank() > rJ
reps = []
cols_nz = [k for k in range(nullity) if any(int(X[j, k]) for j in range(nD))]
random.shuffle(cols_nz)
for k in cols_nz:
    v = [int(X[j, k]) for j in range(C)]
    if notinJ(v[:nD]):
        D, N = vec_to_DN(v); reps.append({'kind': f'nullspace col {k}', 'D': D, 'N': N})
    if len(reps) >= 3: break
for t in range(2):
    coef = [random.randrange(p) for _ in range(nullity)]
    v = [sum(coef[k] * int(X[j, k]) for k in range(nullity)) % p for j in range(C)]
    D, N = vec_to_DN(v); reps.append({'kind': f'random combo {t}', 'D': D, 'N': N})
with open(f'/tmp/hodge/results/dense-rational/e1kernel_d{d}_par{parD}_reps.json', 'w') as f:
    json.dump({'p': p, 'd': d, 'parD': parD, 'reps': [{'kind': r['kind'], 'D': {' '.join(map(str, k)): v for k, v in r['D'].items()},
               'N': [{' '.join(map(str, k)): v for k, v in Ni.items()} for Ni in r['N']]} for r in reps]}, f)
print(f'saved {len(reps)} representatives with E1 numerators', flush=True)
