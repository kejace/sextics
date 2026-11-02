"""Exact rational E1 lift for a given D: find N (deg d+1) with sum_i N_i F_i = A D, sector-restricted, exact fmpq solve.
Usage: python3 e1_lift.py name 'D expression' out.json   (appends/creates a JSON list)"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from lib import *
name, Dsrc, outfn = sys.argv[1], sys.argv[2], sys.argv[3]
F, A, q, Fz, Az = load()
D = core.parse(Dsrc); Dz = integral(D)
G = semi_invariant_subgroup(group_elements(Fz), [Az, Dz]); gens = generators(G); chi = make_chi(gens)
chiA, chiD = chi_of(Az, chi), chi_of(Dz, chi)
d = core.degree(Dz)
dF = [ideriv(Fz, i) for i in range(6)]
target = cadd(chiA, chiD)
rows = {}; cols = []; meta = []
for i in range(6):
    ci = cadd(target, chi(E6[i]))
    for m in mons(d + 1):
        if chi(m) != ci: continue
        col = {}
        for e, v in dF[i].items():
            k = tuple(a+b for a, b in zip(m, e)); col[k] = col.get(k, 0) + v
        for k in col: rows.setdefault(k, len(rows))
        cols.append(col); meta.append((i, m))
rhs = imul(Az, Dz)
for k in rhs: rows.setdefault(k, len(rows))
R, C = len(rows), len(cols)
flat = [0] * (R * (C + 1))
for j, col in enumerate(cols):
    for k, v in col.items(): flat[rows[k]*(C+1) + j] = v
for k, v in rhs.items(): flat[rows[k]*(C+1) + C] = v
t0 = time.time()
M = fmpq_mat(fmpz_mat(R, C + 1, flat)); Mr, rank = M.rref()
pivots = []
for i in range(rank):
    j = 0
    while j <= C and Mr[i, j] == 0: j += 1
    pivots.append(j)
assert pivots[-1] != C, 'E1 inconsistent: A D not in J'
x = [Fraction(0)] * C
for i, j in enumerate(pivots): x[j] = Fraction(int(Mr[i, C].p), int(Mr[i, C].q))
# N for the integral system: sum N_i dFz_i = Az Dz. Convert to original F, A, D: F = Fz/LF etc.
LF = Fraction(Fz[next(iter(Fz))]) / Fraction(F[next(iter(Fz))]); LA = Fraction(Az[next(iter(Az))]) / Fraction(A[next(iter(Az))]); LD = Fraction(Dz[next(iter(Dz))]) / Fraction(D[next(iter(Dz))])
N = [{} for _ in range(6)]
for (i, m), v in zip(meta, x):
    if v: N[i][m] = v * LF / (LA * LD)
lhs = {}
for i in range(6): lhs = core.add(lhs, core.mul(N[i], core.pderiv(F, i)))
assert core.sub(lhs, core.mul(A, D)) == {}, 'lift check failed'
print(f'{name}: d={d} rows={R} cols={C} rank={rank} nullity={C-rank} exact lift found ({round(time.time()-t0,1)}s); N terms {sum(len(n) for n in N)}')
try: L = json.load(open(outfn))
except Exception: L = []
L = [e for e in L if e['name'] != name]
L.append({'name': name, 'd': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'sector': list(chiD), 'group_order': len(G)})
json.dump(L, open(outfn, 'w'), indent=0)
