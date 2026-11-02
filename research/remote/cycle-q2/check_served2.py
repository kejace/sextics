import sys, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, numpy as np
from hodge import core
from flint import nmod_mat
F, A, q, rec = core.load_instance('cycle-q2')
p = 2147483647
Gw = eq.symmetry_group(F); gens = eq.generators(Gw); ng = gens.shape[0]
kA = eq.key_of_poly(A, gens)
kx = [int(eq.keys_of(np.eye(6, dtype=np.int64)[i:i+1], gens)[0]) for i in range(6)]
Fp = core.reduce_mod(F, p); Fi = [core.pderiv(Fp, i, p) for i in range(6)]
def sector(deg, kappa):
    allm = eq.monomial_exps(deg); km = eq.keys_of(allm, gens)
    return [tuple(int(x) for x in e) for e in allm[km == kappa]]
def rank_of(cols):
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    M = np.zeros((len(rows), len(cols)), dtype=np.int64)
    for j, c in enumerate(cols):
        for e, v in c.items(): M[rows[e], j] = v % p
    return nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p).rank(), len(rows)
# B components: B_i has key kA + kx[i]
Bcols = []; meta = []
for i in range(6):
    for m in sector(7, eq.key_add(kA, kx[i], ng)):
        Bcols.append((i, m))
print('B dims', len(Bcols))
J12 = [core.mul({m: 1}, Fi[i], p) for i, m in Bcols]
rJ, nrows = rank_of(J12); print('rank J_12^k =', rJ, 'of S_12^k rows', nrows)
# div map to R_6^k: columns div(B) and J_6^k generators
J6 = []
for i in range(6):
    for m in sector(1, eq.key_add(kA, kx[i], ng)):
        J6.append(core.mul({m: 1}, Fi[i], p))
rJ6, n6 = rank_of(J6); print('J_6^k rank', rJ6, 'S_6^k', n6)
divB = [core.pderiv({m: 1}, i, p) for i, m in Bcols]
r_all, _ = rank_of(J6 + divB)
print('rank [J6 | divB] =', r_all, ' => rank of B->R_6 is', r_all - rJ6)
# exact space dim: {A'=B.gradF : divB in J6} = dim(B-space mod Koszul with divB in J6)
# compute directly: nullspace approach: matrix [divB cols | J6 cols] ; B-vectors with divB + J6-combination = 0
cols = divB + J6
rows = {}
for c in cols:
    for e in c: rows.setdefault(e, len(rows))
M = np.zeros((len(rows), len(cols)), dtype=np.int64)
for j, c in enumerate(cols):
    for e, v in c.items(): M[rows[e], j] = v % p
X, nul = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p).nullspace()
nB = len(divB)
# images A' = B.gradF of nullspace B-parts
imgs = []
for kcol in range(nul):
    P = {}
    for j in range(nB):
        v = int(X[j, kcol])
        if v:
            i, m = Bcols[j]
            P = core.add(P, core.scale(core.mul({m: 1}, Fi[i], p), v, p), p)
    if P: imgs.append(P)
r_img, _ = rank_of(imgs) if imgs else (0, 0)
print('nullity', nul, 'dim of exact targets {B.gradF : divB in J6} =', r_img)

# --- explicit verification of one exact target: A' = B.gradF, N = B/2 + F*C  (mod p)
inv2 = pow(2, p - 2, p)
J6meta = []
for i in range(6):
    for m in sector(1, eq.key_add(kA, kx[i], ng)):
        J6meta.append((i, m))
for kcol in range(nul):
    Bv = [int(X[j, kcol]) for j in range(nB)]
    if not any(Bv): continue
    Cv = [int(X[nB + j, kcol]) for j in range(len(J6meta))]
    # divB + sum Cv_j x_m F_i = 0  => divB/2 = -(1/2) sum Cv_j (x_m e_i).gradF  => C = -(1/2) sum Cv_j x_m e_i
    N = [{} for _ in range(6)]
    for j, v in enumerate(Bv):
        if v:
            i, m = Bcols[j]; N[i] = core.add(N[i], {m: v * inv2 % p}, p)
    for j, v in enumerate(Cv):
        if v:
            i, m = J6meta[j]; N[i] = core.add(N[i], core.scale(core.mul(Fp, {m: 1}, p), (p - v) * inv2 % p, p), p)
    Ap = {}
    for j, v in enumerate(Bv):
        if v:
            i, m = Bcols[j]; Ap = core.add(Ap, core.scale(core.mul({m: 1}, Fi[i], p), v, p), p)
    # E = A' - 2 N.gradF + F div N
    E = dict(Ap)
    for i in range(6):
        E = core.sub(E, core.scale(core.mul(N[i], Fi[i], p), 2, p), p)
        E = core.add(E, core.mul(Fp, core.pderiv(N[i], i, p), p), p)
    print('example exact target: |A\'| terms', len(Ap), ' E == 0 ?', len(E) == 0, ' (nonzero terms:', len(E), ')')
    # check also divB alone
    dB = {}
    for j, v in enumerate(Bv):
        if v:
            i, m = Bcols[j]; dB = core.add(dB, core.scale(core.pderiv({m: 1}, i, p), v, p), p)
    print('   deg divB', core.degree(dB) if dB else None, 'terms', len(dB))
    break
