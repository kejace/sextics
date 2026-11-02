"""Exact verification over Q that the rational degree-7 even kernel element D7 satisfies A*D7 in J(F):
solve N(F) = A*D7 (deg N_i = 8) mod several primes, CRT + rational reconstruction, exact check with Fractions.
Also: structural facts about D7 (depends on Y only through M; factorization)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/verify')
from vlib import *
from flint import nmod_mat, ctx
ctx.threads = 2
Fs, As, q = sympy_F()
F = poly_from_sympy(Fs); A = poly_from_sympy(As)
D7s = sp.sympify(json.load(open('/tmp/hodge/results/short-power-presentation/e1vec_d7_p0.json'))[0]['D'],
                 locals=dict(zip([str(v) for v in VARS], VARS)))
D7 = poly_from_sympy(D7s)
print('deg D7 =', degree(D7), 'parity', parity_of(D7), 'terms', len(D7), flush=True)
AD = mul(A, D7)                       # exact, degree 13, Y-odd
# structure: D7 as polynomial in X0,X1,X2,M ?
X0, X1, X2, Y0, Y1, Y2 = VARS
Mv = sp.symbols('M')
# substitute Y0 -> M - 3 Y1 - 5 Y2 and check independence of Y1, Y2
sub = sp.expand(D7s.subs(Y0, Mv - 3*Y1 - 5*Y2))
print('D7 depends on Y only via M:', sub.diff(Y1) == 0 and sub.diff(Y2) == 0, flush=True)
print('D7(X,M) factor:', sp.factor(sub), flush=True)
# does it depend on X only through L?  check
Lv = sp.symbols('L')
sub2 = sp.expand(sub.subs(X0, Lv - 2*X1 - 3*X2))
print('D7 depends on X only via L:', sub2.diff(X1) == 0 and sub2.diff(X2) == 0, flush=True)
logjson(claim='D7^+ (e1vec_d7_p0.json) depends on Y only through M', verdict='CONFIRMED' if (sub.diff(Y1) == 0 and sub.diff(Y2) == 0) else 'REFUTED', factor=str(sp.factor(sub))[:500])
# modular solves
primes = [1000000007, 1000000009, 998244353, 754974721]
sols = []
meta = None
for p in primes:
    t0 = time.time()
    Fp = reduce_mod(F, p); ADp = reduce_mod(AD, p)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    cols, meta = [], []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(8, (1 - pi) % 2):     # N_i parity so that N_i F_i is Y-odd
            cols.append(shift(Fi[i], m, 1, p)); meta.append((i, m))
    rows = {}
    for cpoly in cols:
        for e in cpoly: rows.setdefault(e, len(rows))
    for e in ADp: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, cpoly in enumerate(cols):
        for e, v in cpoly.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in ADp.items(): flat[rows[e] * (C + 1) + C] = v
    Mx = nmod_mat(R, C + 1, flat, p)
    Mr, rank = Mx.rref()
    piv = []
    sol = {}
    consistent = True
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        if j == C: consistent = False; break
        piv.append(j); sol[j] = int(Mr[i, C])
    print(f'p={p}: rank={rank} consistent={consistent} pivots={len(piv)} ({time.time()-t0:.0f}s)', flush=True)
    sols.append((p, sol, tuple(piv)))
assert all(s[2] == sols[0][2] for s in sols), 'pivot sets differ between primes'
# CRT + rational reconstruction
from functools import reduce
def crt(residues, mods):
    x, m = 0, 1
    for r, mi in zip(residues, mods):
        t = ((r - x) * pow(m, -1, mi)) % mi
        x += m * t; m *= mi
    return x % m, m
def ratrecon(a, m):
    # standard: find r/s with |r|,s <= sqrt(m/2)
    import math
    bound = math.isqrt(m // 2)
    r0, r1 = m, a % m
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)
N = [{} for _ in range(6)]
fail = 0
for j in sols[0][2]:
    vals = [s[1][j] for s in sols]
    x, m = crt(vals, [s[0] for s in sols])
    v = ratrecon(x, m)
    if v is None: fail += 1; continue
    if v:
        i, mon = meta[j]; N[i][mon] = v
print('reconstruction failures:', fail, flush=True)
# exact check
t0 = time.time()
Fi = [pderiv(F, i) for i in range(6)]
tot = {}
for i in range(6): tot = add(tot, mul(N[i], Fi[i]))
diff = add(tot, scale(AD, -1))
ok = (diff == {})
print('EXACT: sum N_i F_i == A*D7 over Q:', ok, f'({time.time()-t0:.0f}s)', flush=True)
logjson(claim='A*D7^+ in J(F) over Q (exact rational verification with reconstructed N, deg 8)', verdict='CONFIRMED' if ok else ('UNVERIFIED' if fail else 'REFUTED'), reconstruction_failures=fail, primes=primes)
if ok:
    json.dump({'N': [{str(k): str(v) for k, v in Ni.items()} for Ni in N]}, open('d7_N_exact.json', 'w'))
