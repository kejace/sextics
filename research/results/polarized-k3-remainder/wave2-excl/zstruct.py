"""Cheap structural tests of the zero set Z of a relaxation (constraint matrix C over F_p):
 (a) tangent space of Z at random points of the known family W = {a = (U ell/6, 0)} (gauge-fixed coordinates):
     dim T_c Z = nparams - rank(Jac(c)); if it equals dim W (=6 affine) at a random point, Z is smooth of dim W there,
     so W is an irreducible component and no other component passes through a generic point of W.
 (b) linear forms ell with ell * S_1 ⊆ I_2  (=> ell in I^sat, Z ⊆ V(ell)): kernel of a linear system.
 (c) rank-1 / rank-2 search: quadrics in I_2 that are products of linear forms? (only the trivial square test here)
usage: zstruct.py <relax json prefix> [p]"""
import sys, json, random, time
import numpy as np
from flint import nmod_mat
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge.core import monomials
prefix = sys.argv[1]
J = json.load(open(prefix + '.json')); Z = np.load(prefix + '.npz')
C = Z['C'] % J['p']; pairs = [tuple(x) for x in Z['pairs']]; p = J['p']
pars = [(i, tuple(m)) for i, m in J['parameter_basis']]; n = len(pars); d = J['d']
pidx = {pr: k for k, pr in enumerate(pars)}
# reduce constraint rows to a basis (rref) for speed
Cn = nmod_mat(C.shape[0], C.shape[1], C.ravel().tolist(), p)
R, rank = Cn.rref()
rows = [[int(R[i, j]) for j in range(C.shape[1])] for i in range(rank)]
print('constraint rank', rank, 'params', n, 'pairs', len(pairs), flush=True)
# (a) family vectors: a_k = (U_k * ell / 42) for k<3 => parameter (k, m) with m = exponent of U_k*ell-monomial
def family_vector(ell_coeffs):
    v = [0] * n
    for k in range(3):
        for j, c in enumerate(ell_coeffs):
            if not c: continue
            m = [0] * 6; m[k] += 1; m[j] += 1
            key = (k, tuple(m))
            if key in pidx: v[pidx[key]] = (v[pidx[key]] + c) % p
            else: raise KeyError(key)
    return v
inv42 = pow(42, -1, p)
random.seed(11)
def jac_rank_at(c):
    # Jacobian of the quadrics q_r = sum_{(a,b)} R[r,(a,b)] c_a c_b : d q_r / d c_i = sum_b R[r,(i,b)] c_b (+ same for (b,i)); diag: 2 R c_i
    Jm = [[0] * n for _ in range(rank)]
    for r, row in enumerate(rows):
        for j, (a, b) in enumerate(pairs):
            v = row[j]
            if not v: continue
            if a == b: Jm[r][a] = (Jm[r][a] + 2 * v * c[a]) % p
            else:
                Jm[r][a] = (Jm[r][a] + v * c[b]) % p; Jm[r][b] = (Jm[r][b] + v * c[a]) % p
    return nmod_mat(rank, n, [x for rr in Jm for x in rr], p).rank()
def value_at(c):
    vals = []
    for row in rows:
        s = 0
        for j, (a, b) in enumerate(pairs):
            if row[j]: s = (s + row[j] * c[a] * c[b]) % p
        vals.append(s)
    return vals
res = {'prefix': prefix, 'p': p, 'constraint_rank': rank, 'params': n}
for t in range(2):
    ell = [random.randrange(1, p) for _ in range(6)]
    c = [(x * inv42) % p for x in family_vector(ell)]
    vals = value_at(c)
    jr = jac_rank_at(c)
    res[f'family_point_{t}'] = {'all_quadrics_vanish': all(v == 0 for v in vals), 'jac_rank': jr, 'tangent_dim_affine': n - jr, 'family_dim_affine': 6}
    print(json.dumps(res[f'family_point_{t}']), flush=True)
# random non-family point: should not vanish
c = [random.randrange(p) for _ in range(n)]
res['random_point_vanishes'] = all(v == 0 for v in value_at(c))
# (b) linear forms ell with ell*c_i in I_2 for all i.  Work in the quadric space: coordinates = pairs; I_2 = row space of R.
# ell*c_i has coordinates: pair (min(k,i),max(k,i)) gets ell_k.  Condition: ell*c_i in rowspace(R) for each i.
# Compute a basis of the orthogonal complement N of rowspace(R) (nullspace of R, dim = npairs - rank); condition: N^T (ell*c_i) = 0.
Rm = nmod_mat(rank, len(pairs), [x for rr in rows for x in rr], p)
N, nul = Rm.nullspace()           # npairs x nul matrix: columns span {x : R x = 0}
print('quadric-space complement dim', nul, flush=True)
pair_index = {pr: j for j, pr in enumerate(pairs)}
# unknown ell (n coords); for each i and each complement vector nu: sum_k ell_k * nu[pair(k,i)] = 0
cond = []
Nl = [[int(N[j, k]) for k in range(nul)] for j in range(len(pairs))]
for i in range(n):
    for k2 in range(nul):
        row = [0] * n
        for k in range(n):
            pr = (min(k, i), max(k, i))
            row[k] = Nl[pair_index[pr]][k2]
        cond.append(row)
Cond = nmod_mat(len(cond), n, [x for rr in cond for x in rr], p)
K, kdim = Cond.nullspace()
res['linear_forms_in_saturation_deg2'] = kdim
print('linear forms ell with ell*S_1 in I_2:', kdim, flush=True)
json.dump(res, open(prefix + '-zstruct.json', 'w'), indent=1)
