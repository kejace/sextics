"""(a) adversarial re-verification: E1-kernel dimensions by sector at new primes; c.i. Hilbert-function check of dim J;
exact (over Q) facts about D8: D8 not in J_8, A*D8 in J_14 with explicit P (saved for the d=8 relaxation)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
from sect import e1_kernel
S = Setup()
F, A = S.F, S.A
primes = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [1000003, 32003]
dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 11
# Hilbert function of a c.i. of six quintics
poly = [1]
for _ in range(6):
    new = [0] * (len(poly) + 4)
    for i, c in enumerate(poly):
        for k in range(5): new[i + k] += c
    poly = new
from math import comb
def dimJ_ci(n): return comb(n + 5, 5) - (poly[n] if n < len(poly) else 0)
# --- exact facts about D8 (sector (0,4)) ---
D8 = core.parse(json.load(open(WD + '/kernel_d8_0_4_0.json'))['D'])
D8z = integral(D8)
c8 = S.chi_poly(D8z)
t0 = time.time()
# D8 not in J_8 exactly: rank over Q of [J_8^c | D8]
gens8 = jac_gens(S, 8, c8)
m8 = [m for m in mons(8) if S.chi(m) == c8]; i8 = {m: i for i, m in enumerate(m8)}
rows = []
for P_ in gens8 + [D8z]:
    v = [0] * len(m8)
    for e, v_ in P_.items(): v[i8[e]] += v_
    rows.append(v)
M = fmpz_mat(rows); rJ = fmpz_mat(rows[:-1]).rank(); rJD = M.rank()
print(f'exact: rank J_8^c = {rJ}, rank [J_8^c | D8] = {rJD}  (D8 in J iff equal)  {time.time()-t0:.1f}s', flush=True)
log({'exp': 'D8-not-in-J-exact', 'D': 'D8', 'd': 8, 'prime': None, 'result': f'rank J_8^sector={rJ}, with D8: {rJD}', 'rank': rJD,
     'seconds': round(time.time() - t0, 1), 'job': 'excl-ctr-kerverify', 'note': 'exact fmpz_mat ranks over Q; D8 not in J(F)_8 iff ranks differ'})
# A*D8 in J_14 exactly: solve sum P_i F_i = A D8 over Q (sector-restricted), fmpq rref
t0 = time.time()
AD8 = imul(S.Az, D8z)
cT = S.chi_poly(AD8)
m14 = [m for m in mons(14) if S.chi(m) == cT]; i14 = {m: i for i, m in enumerate(m14)}
cols = []; meta = []
for i in range(6):
    ci = cadd(cT, S.chiE[i])
    for m in mons(9):
        if S.chi(m) != ci: continue
        cols.append(imul({m: 1}, S.dF[i])); meta.append((i, m))
R, C = len(m14), len(cols)
flat = [0] * (R * (C + 1))
for j, P_ in enumerate(cols):
    for e, v in P_.items(): flat[i14[e] * (C + 1) + j] += v
for e, v in AD8.items(): flat[i14[e] * (C + 1) + C] += v
Maug = fmpq_mat(fmpz_mat(R, C + 1, flat))
Mr, rk = Maug.rref()
piv = []
for i in range(rk):
    j = 0
    while j <= C and Mr[i, j] == 0: j += 1
    piv.append(j)
consistent = not (piv and piv[-1] == C)
print(f'exact: A*D8 in J_14 (Fz-version): consistent={consistent} rows={R} cols={C} rank={rk} {time.time()-t0:.1f}s', flush=True)
Pz = [dict() for _ in range(6)]
if consistent:
    xs = [Fraction(0)] * C
    for i, j in enumerate(piv):
        v = Mr[i, C]; xs[j] = Fraction(int(v.p), int(v.q))
    # clear denominators: scale D8z and P together
    L = 1
    for x in xs: L = L * x.denominator // math.gcd(L, x.denominator)
    for (i, m), x in zip(meta, xs):
        if x: Pz[i][m] = int(x * L)
    D8zz = {e: v * L for e, v in D8z.items()}
    # exact verification: sum_i Pz_i dFz_i == Az * D8zz
    lhs = {}
    for i in range(6): lhs = core.add(lhs, imul(Pz[i], S.dF[i]))
    ok = (lhs == imul(S.Az, D8zz))
    print('exact verification of P(Fz) = Az*D8z:', ok, 'scale L =', L, flush=True)
    json.dump({'L': L, 'D8z': {','.join(map(str, e)): v for e, v in D8zz.items()},
               'Pz': [{','.join(map(str, e)): v for e, v in Pi.items()} for Pi in Pz], 'verified': ok,
               'note': 'Fz = 7F integer; sum_i Pz_i d_i Fz = Az * D8z exactly over Z'}, open(PK + '/e1_d8_exact.json', 'w'))
log({'exp': 'A*D8-in-J-exact', 'D': 'D8', 'd': 8, 'prime': None, 'result': f'consistent={consistent}, verified={consistent and ok}', 'rank': rk,
     'seconds': round(time.time() - t0, 1), 'job': 'excl-ctr-kerverify', 'note': f'exact fmpq rref, sector-restricted E1 system rows={R} cols={C}; P saved to e1_d8_exact.json'})
# --- E1 kernel dims by sector at new primes ---
for p in primes:
    for d in range(5, dmax + 1):
        t0 = time.time()
        res, gens = e1_kernel(F, A, d, p=p, want_vectors=False)
        kerA = sum(r['kerA_in_S'] for r in res.values()); dimJ = sum(r['dimJ'] for r in res.values())
        kmJ = sum(r['ker_mod_J'] for r in res.values()); dimS = sum(r['dim_S'] for r in res.values())
        secs = round(time.time() - t0, 1)
        ok = (dimJ == dimJ_ci(d))
        sectors = {str(c): r['ker_mod_J'] for c, r in res.items() if r['ker_mod_J']}
        print(f'p={p} d={d}: kerA={kerA} dimJ={dimJ} (c.i. {dimJ_ci(d)} ok={ok}) ker mod J={kmJ} sectors={sectors} {secs}s', flush=True)
        log({'exp': 'e1-kernel-dims-reverify', 'family': 'all D of degree d', 'd': d, 'prime': p,
             'result': f'dim ker([A]: S_d -> R_(d+6)) = {kerA}; dim J_d = {dimJ} (c.i. value {dimJ_ci(d)}, match={ok}); kernel mod J = {kmJ}',
             'rank': kerA, 'seconds': secs, 'job': 'excl-ctr-kerverify', 'sectors': sectors,
             'note': 'sect.e1_kernel by G-sector; modular kernel dims are upper bounds for Q when dim J matches the c.i. Hilbert function'})
