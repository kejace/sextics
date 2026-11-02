"""Exact, LLL-reduced integer basis of K4 = ker(eta : R_{f,4} -> R_{f,10}) (quartics g with z0^2 z1^2 z2^2 * g in J(f)_10).
These are the building blocks of the non-Jacobian E1-kernel of the K3-sector target at degree 8 (K4 (x) K4 ⊂ R_4 (x) R_4)."""
import sys, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction as Fr
from flint import fmpz_mat
from hodge import core
uv = json.load(open('/tmp/hodge/results/polarized-k3-remainder/uv_forms.json'))
f = core.parse(uv['f'])
tern = lambda P: {(e[0], e[1], e[2]): v for e, v in P.items()}
fi3 = [tern(core.pderiv(f, i)) for i in range(3)]
def tmons(k): return [(a, b, k-a-b) for a in range(k+1) for b in range(k+1-a)]
def tmul(P, Q):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = (e1[0]+e2[0], e1[1]+e2[1], e1[2]+e2[2]); out[e] = out.get(e, 0) + c1*c2
    return {e: v for e, v in out.items() if v}
m10 = tmons(10); i10 = {m: i for i, m in enumerate(m10)}; m4 = tmons(4)
gens = [tmul({m: Fr(1)}, fi3[i]) for i in range(3) for m in tmons(5)]
cols = gens + [tmul({m: Fr(1)}, {(2,2,2): Fr(1)}) for m in m4]
from flint import fmpq_mat, fmpq
Mq = fmpq_mat(len(m10), len(cols))
for j, c in enumerate(cols):
    for e, v in c.items(): Mq[i10[e], j] = fmpq(int(Fr(v).numerator * 7), int(Fr(v).denominator))
R, rank = Mq.rref()
ng = len(gens); ncol = len(cols)
pivots = []
for i in range(rank):
    j = 0
    while j < ncol and R[i, j] == 0: j += 1
    pivots.append(j)
free = [j for j in range(ncol) if j not in set(pivots)]
import math
vecs = []
for j in free:
    v = [fmpq(0)] * ncol; v[j] = fmpq(1)
    for i, pj in enumerate(pivots): v[pj] = -R[i, j]
    g = v[ng:]
    if any(x != 0 for x in g):
        den = 1
        for x in g: den = den * int(x.q) // math.gcd(den, int(x.q))
        vecs.append([int(x * den) for x in g])
indep = []
for r in vecs:
    if fmpz_mat(indep + [r]).rank() > len(indep): indep.append(r)
L = fmpz_mat(indep).lll()
basis = []
for i in range(L.nrows()):
    g = {m4[j]: int(L[i, j]) for j in range(len(m4)) if int(L[i, j])}
    gg = 0
    for v in g.values(): gg = math.gcd(gg, abs(v))
    g = {m: v // gg for m, v in g.items()}
    basis.append(g)
names = ['z0', 'z1', 'z2']
def gstr(g): return ' + '.join(f"{v}*" + '*'.join(f'{names[i]}^{e}' if e > 1 else names[i] for i, e in enumerate(m) if e) for m, v in sorted(g.items(), reverse=True))
out = {'dim_K4': len(basis), 'basis': [gstr(g) for g in basis], 'exponents': [[list(m), v] for g in basis for m, v in g.items()],
       'check': 'each g satisfies z0^2 z1^2 z2^2 g in J(f)_10 (exact)'}
# exact verification of membership for each basis element
Jz = fmpz_mat(len(m10), ng)
for j, c in enumerate(gens):
    for e, v in c.items(): Jz[i10[e], j] = int(Fr(v) * 7)
rJ = Jz.rank()
ok = True
for g in basis:
    prod = tmul({m: Fr(v) for m, v in g.items()}, {(2,2,2): Fr(1)})
    A = fmpz_mat(len(m10), ng + 1)
    for j, c in enumerate(gens):
        for e, v in c.items(): A[i10[e], j] = int(Fr(v) * 7)
    for e, v in prod.items(): A[i10[e], ng] = int(v)
    ok &= (A.rank() == rJ)
out['membership_verified_exact'] = ok
json.dump(out, open('/tmp/hodge/results/polarized-k3-remainder/wave2-excl/K4-basis-lll.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'exponents'}, indent=1))
