"""D8 exclusion re-verification by independent methods:
 (i) unrestricted (no sector reduction) fixed-D system at a new prime: rank / rank_aug;
 (ii) sector system: modular nullspace at two primes -> CRT -> rational reconstruction -> exact verification of the
      null vectors over Q (seals rank_Q = rank_p), plus rank_aug mod p  => inconsistency over Q by the modular-rank argument;
 (iii) exact fmpq rref of the sector augmented system (independent of wave-1's fmpz rank)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
from hodge import fixed_d
import sect
S = Setup(); F, A = S.F, S.A
D8 = core.parse(json.load(open(WD + '/kernel_d8_0_4_0.json'))['D'])
p1, p2 = 1000003, 2147483587
# (i)
t0 = time.time()
sysF = fixed_d.FixedDSystem(F, A, 1, D8, p1, verbose=True)
M = sysF.dense_augmented()
ra = M.rank()
C = sysF.ncols
M0 = nmod_mat(sysF.nrows, C, [int(M[i, j]) for i in range(sysF.nrows) for j in range(C)], p1)
r0 = M0.rank()
print(f'(i) unrestricted p={p1}: rows={sysF.nrows} cols={C} rank={r0} rank_aug={ra} ({time.time()-t0:.1f}s)', flush=True)
log({'exp': 'D8-unrestricted-modular', 'D': 'D8', 'd': 8, 'prime': p1, 'result': f'rank={r0} rank_aug={ra} consistent={r0 == ra}', 'rank': r0,
     'seconds': round(time.time() - t0, 1), 'job': 'excl-ctr-d8full', 'note': f'fixed_d.FixedDSystem without sector reduction, rows={sysF.nrows} cols={C}'})
# (ii)
t0 = time.time()
Sy = sect.Sys(F, A, 1, D8, verbose=True)
def null_mod(p):
    M = Sy.dense(p, aug=False)
    X, nullity = M.nullspace()
    return [[int(X[j, k]) for j in range(Sy.ncols)] for k in range(nullity)], nullity
v1, n1 = null_mod(p1); v2, n2 = null_mod(p2)
ra1 = Sy.dense(p1, aug=True).rank(); ra2 = Sy.dense(p2, aug=True).rank()
print(f'(ii) sector system rows={Sy.nrows} cols={Sy.ncols}: nullity p1={n1} p2={n2}; rank_aug p1={ra1} p2={ra2}', flush=True)
# reduce both nullspaces to rref so the bases correspond, then CRT + ratrecon
def rref_rows(vs, p):
    if not vs: return []
    Mv = nmod_mat(len(vs), Sy.ncols, [x for v in vs for x in v], p)
    R, rk = Mv.rref()
    return [[int(R[i, j]) for j in range(Sy.ncols)] for i in range(rk)]
r1 = rref_rows(v1, p1); r2 = rref_rows(v2, p2)
sealed = 0; verified = []
if len(r1) == len(r2):
    cols = Sy.cols; rhs = Sy.rhs
    for a, b in zip(r1, r2):
        xs = [core.crt_pair(x, p1, y, p2)[0] for x, y in zip(a, b)]
        rat = [core.ratrecon(x, p1 * p2) for x in xs]
        if any(r is None for r in rat): verified.append(False); continue
        # exact check: M * rat == 0 (integer system; scale by lcm)
        L = 1
        for r in rat: L = L * r.denominator // math.gcd(L, r.denominator)
        xi = [int(r * L) for r in rat]
        acc = {}
        for j, col in enumerate(cols):
            if xi[j] == 0: continue
            for r_, v in col.items(): acc[r_] = acc.get(r_, 0) + v * xi[j]
        ok = all(v == 0 for v in acc.values())
        verified.append(ok); sealed += ok
print(f'(ii) exact null vectors verified: {sealed}/{len(r1)}  => rank_Q <= {Sy.ncols - sealed}; rank_p >= {Sy.ncols - n1} ({time.time()-t0:.1f}s)', flush=True)
concl = (sealed == n1 and ra1 == Sy.ncols - n1 + 1)
log({'exp': 'D8-sector-modular-rank-sealed', 'D': 'D8', 'd': 8, 'prime': f'{p1},{p2}', 'result': f'nullity mod p = {n1},{n2}; {sealed} exact null vectors verified over Q; rank_aug mod p = {ra1},{ra2}; '
     + ('INCONSISTENT OVER Q (sealed)' if concl else 'not sealed'), 'rank': Sy.ncols - n1, 'seconds': round(time.time() - t0, 1), 'job': 'excl-ctr-d8full',
     'note': f'sector system rows={Sy.nrows} cols={Sy.ncols}; rank_Q = cols - sealed nullity = {Sy.ncols - sealed} exactly, rank_Q(aug) >= rank_p(aug)'})
# (iii) exact fmpq rref of augmented system
t0 = time.time()
Maug = fmpq_mat(Sy.dense(None, aug=True)); Mr, rk = Maug.rref()
piv = []
for i in range(rk):
    j = 0
    while j <= Sy.ncols and Mr[i, j] == 0: j += 1
    piv.append(j)
cons = not (piv and piv[-1] == Sy.ncols)
print(f'(iii) exact fmpq rref: rank_aug={rk} consistent={cons} ({time.time()-t0:.1f}s)', flush=True)
log({'exp': 'D8-sector-exact-fmpq-rref', 'D': 'D8', 'd': 8, 'prime': None, 'result': f'rank_aug over Q = {rk}, consistent={cons}', 'rank': rk,
     'seconds': round(time.time() - t0, 1), 'job': 'excl-ctr-d8full', 'note': 'exact rational rref of the augmented sector system (independent of wave-1 fmpz rank)'})
