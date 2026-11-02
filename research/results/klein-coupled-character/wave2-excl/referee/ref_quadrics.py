"""Second stage of the degree-6 exclusion, rebuilt from scratch mod p1 (unrestricted Koszul form, Y0-monic reduction mod 7F):
for D_lambda = sum lambda_i Z_i (7F)_i, N = A w + K(p) (w_i = lambda_i Z_i, p antisymmetric of degree 2), the certificate condition
E2: D div N - N(D) = 0 mod F is  M(lambda) p = b(lambda), M linear and b quadratic in lambda.  Constant left-kernel vectors y
(y^T M_i = 0 for all i) give necessary quadrics y^T b(lambda) = 0.  We compute the space Q of all such quadrics, compare with the
lens's 8 exact quadrics (local_d6/euler_quadrics_exact.json), and check with Singular that V(Q) + {l0 = 0} is empty (dim 0)."""
import sys, json, time, subprocess
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl/referee')
import ref_lib as L
from fractions import Fraction
from flint import nmod_mat
p = L.P1; t0 = time.time()
F, Fi, H = L.Fpolys(p); A = L.red(L.A, p)
inv7 = pow(7, -1, p)
Frest = {e: (-v * inv7) % p for e, v in F.items() if e != (0,0,0,6,0,0)}   # Y0^6 = Frest mod F
memo = {}
def nf_mono(e):
    """normal form of the monomial e modulo 7F (Y0-monic division): dict"""
    if e[3] < 6: return {e: 1}
    if e in memo: return memo[e]
    base = list(e); base[3] -= 6; base = tuple(base)
    out = {}
    for e2, v in Frest.items():
        ee = tuple(a + b for a, b in zip(base, e2))
        for e3, v3 in nf_mono(ee).items():
            out[e3] = (out.get(e3, 0) + v * v3) % p
    memo[e] = {k: v for k, v in out.items() if v}; return memo[e]
def nf(P):
    out = {}
    for e, v in P.items():
        for e2, v2 in nf_mono(e).items(): out[e2] = (out.get(e2, 0) + v * v2) % p
    return {k: v for k, v in out.items() if v}
Z = [tuple(1 if t == i else 0 for t in range(6)) for i in range(6)]
D = [L.shift(Fi[i], Z[i], 1, p) for i in range(6)]            # D_i = Z_i F_i
AZ = [L.shift(A, Z[i], 1, p) for i in range(6)]                 # N^0 component i for lambda_i
def T(i, j):   # D_i * d_j(A Z_j) - (A Z_j) * d_j(D_i)
    return L.add(L.mul(D[i], L.pderiv(AZ[j], j, p), p), L.scale(L.mul(AZ[j], L.pderiv(D[i], j, p), p), p - 1, p), p)
pairs = [(i, j) for i in range(6) for j in range(i, 6)]
B = {}
for (i, j) in pairs:
    B[(i, j)] = nf(T(i, i)) if i == j else nf(L.add(T(i, j), T(j, i), p))
print('b_ij built', f'[{time.time()-t0:.0f}s]', flush=True)
mons2 = L.monomials(2)
cols = []
Di = [[L.pderiv(D[i], a, p) for a in range(6)] for i in range(6)]
for i in range(6):
    for a in range(6):
        for b in range(a + 1, 6):
            for m in mons2:
                dm_a, dm_b = L.pderiv({m: 1}, a, p), L.pderiv({m: 1}, b, p)
                col = L.add(L.mul(D[i], L.add(L.mul(dm_a, Fi[b], p), L.scale(L.mul(dm_b, Fi[a], p), p - 1, p), p), p),
                            L.scale(L.shift(L.add(L.mul(Fi[b], Di[i][a], p), L.scale(L.mul(Fi[a], Di[i][b], p), p - 1, p), p), m, 1, p), p - 1, p), p)
                cols.append(nf(col))
print('columns built', len(cols), f'[{time.time()-t0:.0f}s]', flush=True)
rows = {}
for c in cols:
    for e in c: rows.setdefault(e, len(rows))
for b in B.values():
    for e in b: rows.setdefault(e, len(rows))
nr = len(rows); nc = len(cols); nb = len(pairs)
flat = [0] * (nr * (nc + nb))
for c, col in enumerate(cols):
    for e, v in col.items(): flat[rows[e] * (nc + nb) + c] = v
for k, (i, j) in enumerate(pairs):
    for e, v in B[(i, j)].items(): flat[rows[e] * (nc + nb) + nc + k] = v
M = nmod_mat(nr, nc + nb, flat, p)
print('matrix', nr, 'x', nc + nb, f'[{time.time()-t0:.0f}s]', flush=True)
X, nul = M.nullspace()
print('nullity of [M | B]:', nul, 'rank', nc + nb - nul, f'[{time.time()-t0:.0f}s]', flush=True)
# projection of the kernel onto the b-coordinates -> Q^perp
K = nmod_mat(nul, nb, [int(X[nc + k, v]) for v in range(nul) for k in range(nb)], p)
rK = K.rank(); print('dim Q^perp =', rK, '=> dim Q (necessary quadrics) =', nb - rK, flush=True)
X2, nul2 = K.nullspace()
Qbasis = [[int(X2[k, v]) for k in range(nb)] for v in range(nul2)]
# compare with the lens's exact quadrics (keys 'i,j' i<=j, coefficient fractions)
lens = json.load(open('../local_d6/euler_quadrics_exact.json'))['quadrics']
pidx = {pr: k for k, pr in enumerate(pairs)}
ok = True
for qd in lens:
    vec = [0] * nb
    for key, val in qd.items():
        i, j = (int(x) for x in key.split(',')); fr = Fraction(val)
        vec[pidx[(min(i, j), max(i, j))]] = (fr.numerator * pow(fr.denominator, -1, p)) % p
    Kv = K * nmod_mat(nb, 1, vec, p)
    inQ = all(int(Kv[r, 0]) == 0 for r in range(nul))
    ok &= inQ
print('all 8 lens quadrics (reduced mod p1) lie in my quadric space:', ok)
# rank of lens quadrics mod p
LM = nmod_mat(len(lens), nb, [((Fraction(qd.get(f'{i},{j}', '0')).numerator * pow(Fraction(qd.get(f'{i},{j}', '0')).denominator, -1, p)) % p) for qd in lens for (i, j) in pairs], p)
print('rank of the 8 lens quadrics mod p1:', LM.rank(), '| so the spaces coincide:', ok and LM.rank() == nb - rK)
# Singular: my quadric space mod p1 -> dimension of V(Q) and of V(Q, l0)
names = ','.join(f'l{i}' for i in range(6))
def qstr(vec): return '+'.join(f'({v})*l{i}*l{j}' for (i, j), v in zip(pairs, vec) if v)
sg = f'ring R = {p}, ({names}), dp;\nideal Q = ' + ',\n'.join(qstr(v) for v in Qbasis) + ';\nideal SQ = std(Q);\n"dim V(Q) (affine)", dim(SQ), "mult", mult(SQ);\nideal Q0 = Q, l0;\nideal SQ0 = std(Q0);\n"with l0 = 0: dim", dim(SQ0), "vdim", vdim(SQ0);\nideal G = l0-l5, l1-l5, l2-l5, l3-l5, l4-l5;\n"gauge line in V(Q):", size(reduce(Q, std(G))) == 0;\nquit;\n'
open(f'quadrics_p{p}.sing', 'w').write(sg)
r = subprocess.run(['Singular', '-q', f'quadrics_p{p}.sing'], capture_output=True, text=True, timeout=600)
print(r.stdout.strip(), r.stderr.strip()[:200])
# exact replay of the lens's Singular step over Q with its 8 exact quadrics
def qstrQ(qd): return '+'.join(f'({val})*l{key.split(",")[0]}*l{key.split(",")[1]}' for key, val in qd.items())
sgq = f'LIB "primdec.lib";\nring R = 0, ({names}), dp;\nideal Q = ' + ',\n'.join(qstrQ(qd) for qd in lens) + ';\nideal SQ = std(Q);\n"exact: dim V(Q) (affine)", dim(SQ);\nlist Lm = minAssGTZ(Q);\n"n_components", size(Lm);\nint i; for (i = 1; i <= size(Lm); i++) {{ "component", i; print(Lm[i]); }}\nideal SQ0 = std(ideal(Q, l0));\n"exact with l0 = 0: dim", dim(SQ0), "vdim", vdim(SQ0);\nquit;\n'
open('quadrics_exactQ.sing', 'w').write(sgq)
r = subprocess.run(['Singular', '-q', 'quadrics_exactQ.sing'], capture_output=True, text=True, timeout=600)
print(r.stdout.strip(), r.stderr.strip()[:300])
json.dump({'prime': p, 'rows': nr, 'cols_M': nc, 'nullity_MB': nul, 'dim_Qperp': rK, 'dim_Q': nb - rK, 'lens_quadrics_in_Q': ok, 'lens_rank': LM.rank(), 'Q_basis_pairs': pairs, 'Q_basis': Qbasis}, open(f'referee_quadrics_p{p}.json', 'w'))
print('QUADRICS_DONE', f'[{time.time()-t0:.0f}s]')
