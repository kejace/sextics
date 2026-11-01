"""Search the affine E1 family D = K8 + sum g_i F_i (mod p) for representatives of prescribed shape
(a shape = set of monomials forced to vanish).  Linear algebra in S_8^even only (cheap)."""
import sys, json, time, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, shift, pderiv
from flint import nmod_mat
P = core.PRIMES[0]
F, A, q, rec = core.load_instance('reflection-coupled-residual')
K8 = core.reduce_mod(core.parse(json.load(open('kernel_d8_X03Y03.json'))[0]), P)
Fp = core.reduce_mod(F, P); Fi = [pderiv(Fp, i, P) for i in range(NV)]
mons = monomials(8, 0); idx = {m: i for i, m in enumerate(mons)}; C = len(mons)
par = lambda i: 0 if i < 3 else 1
T = []
for i in range(NV):
    for m in monomials(3, par(i)):
        T.append(shift(Fi[i], m, 1, P))
ng = len(T)
def vec(poly):
    v = [0] * C
    for e, c in poly.items(): v[idx[e]] = c % P
    return v
Tv = [vec(t) for t in T]; Kv = vec(K8)
def test(name, forbid):
    """forbid: predicate on exponent tuple -> monomial must vanish in D'"""
    rows = [j for j, m in enumerate(mons) if forbid(m)]
    if not rows: return
    flat = []
    for j in rows:
        flat.extend([Tv[k][j] for k in range(ng)] + [(-Kv[j]) % P])
    M = nmod_mat(len(rows), ng + 1, flat, P)
    Mr, r = M.rref()
    r0 = nmod_mat(len(rows), ng, [Tv[k][j] for j in rows for k in range(ng)], P).rank()
    ok = (r == r0)
    print(f'{name:45s} constraints={len(rows):4d} rankT={r0:3d} solvable={ok} family dim={ng-r0}', flush=True)
    return ok
V = ['X0','X1','X2','Y0','Y1','Y2']
test('divisible by Y2^2 (D = u*G)', lambda m: m[5] < 2)
test('divisible by Y2^4', lambda m: m[5] < 4)
test('no Y2 (D in C[x])', lambda m: m[5] > 0)
for i in range(6):
    test(f'divisible by {V[i]}', lambda m, i=i: m[i] == 0)
    test(f'divisible by {V[i]}^2', lambda m, i=i: m[i] < 2)
for i, j in itertools.combinations(range(6), 2):
    test(f'independent of {V[i]},{V[j]} (cone, vertex line)', lambda m, i=i, j=j: m[i] > 0 or m[j] > 0)
for i, j in [(0,1),(3,4),(0,3),(0,5),(3,5),(1,4)]:
    test(f'in ideal ({V[i]},{V[j]})^2', lambda m, i=i, j=j: m[i] + m[j] < 2)
    test(f'divisible by {V[i]}*{V[j]}', lambda m, i=i, j=j: m[i] == 0 or m[j] == 0)
test('Y-degree <= 2', lambda m: m[3] + m[4] + m[5] > 2)
test('Y-degree <= 4', lambda m: m[3] + m[4] + m[5] > 4)
test('X-degree <= 2', lambda m: m[0] + m[1] + m[2] > 2)
test('X-degree <= 4', lambda m: m[0] + m[1] + m[2] > 4)
test('only even exponents (poly in squares)', lambda m: any(x % 2 for x in m))
test('in ideal (X0,X1,X2)^4', lambda m: m[0]+m[1]+m[2] < 4)
test('in ideal (Y0,Y1,Y2)^4', lambda m: m[3]+m[4]+m[5] < 4)
test('in ideal (X0,Y0)^4', lambda m: m[0]+m[3] < 4)
test('in ideal (X0,Y0)^6', lambda m: m[0]+m[3] < 6)
test('in (X0,X1,Y0,Y1)^6', lambda m: m[0]+m[1]+m[3]+m[4] < 6)
test('in (X2,Y2)^2', lambda m: m[2]+m[5] < 2)
test('in (X2,Y2)^4', lambda m: m[2]+m[5] < 4)
