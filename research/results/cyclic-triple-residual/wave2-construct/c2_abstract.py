# C2 (abstract): G108 = <d,e,sigma> x <r> x <g3> on the Galois closure Zt; elements (a,b,eps,rho,c) = d^a e^b sigma^eps r^rho g3^c,
# sigma d sigma = e.  Class irrep rho_full: v10 (d->1,e->z), v01 (d->z,e->1), sigma swaps, r->1, g3->z.
# Pseudo-reflections on Zt: d^{+-1}, e^{+-1} (hyperplanes Xb=0, X2=0; stabiliser order 3), d^k e^-k sigma (X2 = z^k Xb; order 2).
# Riemann-Hurwitz on Zt ((3,6) c.i. in P^6, K = O(2)):  deg K_{Zt/K} = 2 - 2[d in K] - 2[e in K] - #{k : d^k e^-k sigma in K}.
import itertools, json, cmath, random
import numpy as np
from collections import Counter
z = cmath.exp(2j*cmath.pi/3)
def mul(x, y):
    a,b,e,r,c = x; a2,b2,e2,r2,c2 = y
    if e: a2,b2 = b2,a2
    return ((a+a2)%3, (b+b2)%3, (e+e2)%2, (r+r2)%2, (c+c2)%3)
E = (0,0,0,0,0)
G = [(a,b,e,r,c) for a in range(3) for b in range(3) for e in range(2) for r in range(2) for c in range(3)]
random.seed(0)
for _ in range(2000):
    x,y,w = (random.choice(G) for _ in range(3)); assert mul(mul(x,y),w)==mul(x,mul(y,w))
def rep(x):
    a,b,e,r,c = x
    D = [[z**b*z**c, 0],[0, z**a*z**c]]
    if e: D = [[D[0][1], D[0][0]],[D[1][1], D[1][0]]]
    return D
def m2(A,B): return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
bad = 0
for _ in range(3000):
    x,y = random.choice(G), random.choice(G)
    P = m2(rep(x), rep(y)); Q = rep(mul(x,y))
    if max(abs(P[i][j]-Q[i][j]) for i in range(2) for j in range(2)) > 1e-9: bad += 1
print('homomorphism violations:', bad)
def closure(gens):
    S = {E}; frontier = [E]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = mul(x,g)
                if y not in S: S.add(y); new.append(y)
        frontier = new
    return frozenset(S)
subs = set()
for g in G: subs.add(closure([g]))
for g,h in itertools.combinations(G,2): subs.add(closure([g,h]))
changed = True
while changed:
    changed = False
    for S in list(subs):
        for g in G:
            if g not in S:
                T = closure(list(S) + [g])
                if T not in subs: subs.add(T); changed = True
print('number of subgroups:', len(subs))
def inv_dim(S):
    P = np.zeros((2,2), dtype=complex)
    for x in S: P += np.array(rep(x))/len(S)
    return int(round(np.linalg.matrix_rank(P, tol=1e-8)))
def can_deg(S):
    return 2 - 2*((1,0,0,0,0) in S) - 2*((0,1,0,0,0) in S) - sum(1 for k in range(3) if (k,(-k)%3,1,0,0) in S)
rows = [(len(S), inv_dim(S), can_deg(S)) for S in subs]
serving = [t for t in rows if t[1] > 0]
print('subgroups serving the class:', len(serving), '-> (order, inv_dim, can_deg):', sorted(Counter(serving).items()))
print('serving subgroups with can_deg < 0:', sum(1 for t in serving if t[2] < 0))
print('minimal canonical degree among serving quotients:', min(t[2] for t in serving))
print('all subgroups (order, inv, deg):', sorted(Counter(rows).items()))
json.dump({'n_subgroups': len(subs), 'serving_multiset': {str(k): v for k,v in Counter(serving).items()},
           'min_deg_serving': min(t[2] for t in serving), 'n_serving_negative': sum(1 for t in serving if t[2] < 0)},
          open('/tmp/hodge/results/cyclic-triple-residual/wave2-construct/c2_abstract.json','w'), indent=1)
