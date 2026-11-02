"""From-scratch build of the cleared certificate system for a fixed D (q=1) and consistency test mod p
(dense flint rank of M and [M|b]).  Independent of hodge.*.  Usage: fullsys_v.py <Djson> <p> [tag]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/verify')
from vlib import *
from flint import nmod_mat, ctx
ctx.threads = 4
path, p = sys.argv[1], int(sys.argv[2]); tag = sys.argv[3] if len(sys.argv) > 3 else 'D'
Fs, As, q = sympy_F()
F = poly_from_sympy(Fs); A = poly_from_sympy(As)
D = poly_from_sympy(sp.sympify(json.load(open(path))[0]['D'], locals=dict(zip([str(v) for v in VARS], VARS))))
d = degree(D); parD = parity_of(D); assert parD is not None
Fp, Ap, Dp = reduce_mod(F, p), reduce_mod(A, p), reduce_mod(D, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]; Di = [pderiv(Dp, i, p) for i in range(6)]
degN, degQ = d + 1, 2 * d - 6
# identity: A D^2 - D N(F) + F (D div N - N(D)) - F^2 Q = 0.  Column for N_i = m:  -D m F_i + F (D d_i(m) - m D_i)
t0 = time.time()
cols = []
FD = mul(Fp, Dp, p); F2 = mul(Fp, Fp, p)
G = [add(scale(mul(Dp, Fi[i], p), p - 1, p), scale(mul(Fp, Di[i], p), p - 1, p), p) for i in range(6)]  # -D F_i - F D_i
for i in range(6):
    pi = 0 if i < 3 else 1
    for m in monomials(degN, (1 - pi - parD) % 2):   # N_i F_i must be Y-odd given A D^2 Y-odd: par(N_i) = 1 - par(F_i) - 2parD = 1 - pi
        col = shift(G[i], m, 1, p)
        if m[i]:
            m2 = list(m); m2[i] -= 1
            col = add(col, shift(FD, tuple(m2), m[i], p), p)   # F D d_i(m) = m_i * F D * m/Z_i
        cols.append(col)
nN = len(cols)
for m in monomials(degQ, 1):   # Q Y-odd (A D^2 is Y-odd, F^2 even)
    cols.append(shift(F2, m, p - 1, p))
b = scale(mul(Ap, mul(Dp, Dp, p), p), p - 1, p)   # move A D^2 to rhs with sign: sum cols*x = -A D^2  => we test [M | -AD^2]... sign irrelevant for consistency
rows = {}
for cpoly in cols:
    for e in cpoly: rows.setdefault(e, len(rows))
for e in b: rows.setdefault(e, len(rows))
R, C = len(rows), len(cols)
print(f'[{tag}] d={d} parD={parD} rows={R} cols={C} (N {nN}, Q {C-nN}) built {time.time()-t0:.0f}s', flush=True)
flat = [0] * (R * (C + 1))
for j, cpoly in enumerate(cols):
    for e, v in cpoly.items(): flat[rows[e] * (C + 1) + j] = v
for e, v in b.items(): flat[rows[e] * (C + 1) + C] = v
Maug = nmod_mat(R, C + 1, flat, p); del flat
t1 = time.time(); rA = Maug.rank(); print(f'rank [M|b] = {rA} ({time.time()-t1:.0f}s)', flush=True)
flat = [0] * (R * C)
for j, cpoly in enumerate(cols):
    for e, v in cpoly.items(): flat[rows[e] * C + j] = v
del Maug
M = nmod_mat(R, C, flat, p); del flat
t1 = time.time(); rM = M.rank(); print(f'rank M = {rM} ({time.time()-t1:.0f}s)', flush=True)
cons = (rA == rM)
print(f'[{tag}] consistent mod {p}: {cons}; rank={rM} nullity={C-rM}', flush=True)
logjson(claim=f'full fixed-D system for {tag} (independent builder), consistency mod p', verdict='DATA', p=p, d=d, parD=parD, rows=R, cols=C, rank=rM, rank_aug=rA, consistent=cons, nullity=C - rM)
