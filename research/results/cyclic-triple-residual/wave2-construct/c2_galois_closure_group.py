# C2: averaging-ideal test on the natural Galois closure Zt = {X2^3+Xb^3+tC=0, X2^3 Xb^3 = H} in P^6,
# full symmetry group G108 = <d,e,sigma> x <r,g3~> acting linearly on C^7 = (X0,X1,X2,Xb,Y0,Y1,Y2).
# For every subgroup K: (i) invariants of the class's isotypic irrep rho_full (2-dim: rho of Gamma tensor the
# characters r->1, g3->zeta3 of the class); (ii) pseudo-reflections in K and the O(1)-degree of K_{Zt/K}
# via Riemann-Hurwitz  deg K_{Zt/K} = 2 - sum_R (e_R - 1) deg R  (Zt is a (3,6) c.i.: K_Zt = O(2)).
import itertools, json, cmath
z = cmath.exp(2j*cmath.pi/3)
w6 = cmath.exp(2j*cmath.pi/6)
def diag(*d): return [[d[i] if i==j else 0 for j in range(7)] for i in range(7)]
def mul(A,B): return [[sum(A[i][k]*B[k][j] for k in range(7)) for j in range(7)] for i in range(7)]
I7 = diag(1,1,1,1,1,1,1)
d = diag(1,1,1,z,1,1,1)            # Xb -> z Xb   (deck of Zt -> X)
e = diag(1,1,z,1,1,1,1)            # X2 -> z X2   (lift of mu3)
sig = [[0]*7 for _ in range(7)]    # swap X2 <-> Xb
for i in range(7): sig[i][i] = 1
sig[2][2]=sig[3][3]=0; sig[2][3]=sig[3][2]=1
r = diag(-1,-1,-1,-1,1,1,1)        # tau: X -> -X  (== Y -> -Y projectively)
# g3 = (2,4,0,4,4,0) on (X0,X1,X2,Y0,Y1,Y2) in zeta6-exponents -> (z, z^2, 1, z^2, z^2, 1); lift fixing Xb
g3 = diag(w6**2, w6**4, 1, 1, w6**4, w6**4, 1)
def key(M): return tuple(complex(round(M[i][j].real,6), round(M[i][j].imag,6)) for i in range(7) for j in range(7))
def closure(gens):
    elems = {key(I7): I7}; frontier = [I7]
    while frontier:
        new = []
        for A in frontier:
            for g in gens:
                B = mul(A,g); k = key(B)
                if k not in elems: elems[k] = B; new.append(B)
        frontier = new
    return elems
G = closure([d,e,sig,r,g3]); print('order of G108 =', len(G))
elems = list(G.values()); keys = list(G.keys())
idx = {k:i for i,k in enumerate(keys)}
# multiplication table
def ki(M): return idx[key(M)]
# projective scalars: the group acts on P^6; two matrices differing by a scalar act identically.  All our generators
# are 'honest' (no scalars introduced since r is not a scalar); check that no non-identity element is scalar:
scal = [k for k in keys if all(abs(G[k][i][j] - (G[k][0][0] if i==j else 0))<1e-6 for i in range(7) for j in range(7))]
print('scalar elements:', len(scal))
# pseudo-reflections: eigenvalue-1 eigenspace of dimension 6
import numpy as np
def fix_dim(M):
    A = np.array(M, dtype=complex) - np.eye(7)
    return 7 - np.linalg.matrix_rank(A, tol=1e-8)
refl = {k: fix_dim(G[k]) for k in keys}
nrefl = sum(1 for k in keys if refl[k]==6)
print('pseudo-reflections (fixed space dim 6):', nrefl)
# the class's 2-dim irrep rho_full: basis v10 (e->z, d->1), v01 (e->1, d->z); sigma swaps; r -> +1; g3 -> z (class char zeta3)
# represent each generator as 2x2
rho = {key(d): [[1,0],[0,z]], key(e): [[z,0],[0,1]], key(sig): [[0,1],[1,0]], key(r): [[1,0],[0,1]], key(g3): [[z,0],[0,z]]}
def m2(A,B): return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
# extend rho to all elements by BFS words
rho_all = {key(I7): [[1,0],[0,1]]}
frontier=[I7]
gens=[d,e,sig,r,g3]
while frontier:
    new=[]
    for A in frontier:
        for g in gens:
            B=mul(A,g); k=key(B)
            if k not in rho_all:
                rho_all[k]=m2(rho_all[key(A)], rho[key(g)]); new.append(B)
    frontier=new
# consistency: rho is a homomorphism (check on random products)
import random
random.seed(1)
bad=0
for _ in range(300):
    a=random.choice(keys); b=random.choice(keys)
    P=m2(rho_all[a], rho_all[b]); Q=rho_all[key(mul(G[a],G[b]))]
    if max(abs(P[i][j]-Q[i][j]) for i in range(2) for j in range(2))>1e-6: bad+=1
print('rho homomorphism violations:', bad)
# enumerate subgroups: closures of <=2 generators, then extend by a third generator
subs = {}
def sub_key(S): return tuple(sorted(idx[k] for k in S))
def closure_keys(gens_k):
    return set(closure([G[k] for k in gens_k]).keys())
one = {}
for k in keys:
    S = closure_keys([k]); one[sub_key(S)] = S
subs.update(one)
two = {}
for a in keys:
    for b in keys:
        if idx[a] < idx[b]:
            S = closure_keys([a,b]); sk = sub_key(S)
            if sk not in subs: two[sk]=S
subs.update(two)
three = {}
for sk, S in list(two.items()):
    for c in keys:
        if c not in S:
            T = closure_keys(list(S)[:0] + [keys[i] for i in sk[:2]] + [c]) if False else closure_keys([keys[i] for i in sk] + [c])
            tk = sub_key(T)
            if tk not in subs: three[tk]=T
subs.update(three)
print('subgroups found (<=3 generators):', len(subs))
# for each subgroup: invariants of rho_full, reflection data, canonical degree
def inv_dim(S):
    P = [[0,0],[0,0]]
    for k in S:
        M = rho_all[k]
        for i in range(2):
            for j in range(2): P[i][j] += M[i][j]/len(S)
    return round(np.linalg.matrix_rank(np.array(P,dtype=complex), tol=1e-6))
def can_deg(S):
    # Riemann-Hurwitz: for each reflection hyperplane fixed pointwise by a nontrivial element of K, the stabiliser
    # (cyclic of order e_R) contributes (e_R - 1) * deg(R) with deg R = 1 (hyperplane section of Zt).
    hyps = {}
    for k in S:
        if refl[k]==6 and k != key(I7):
            A = np.array(G[k],dtype=complex) - np.eye(7)
            # fixed hyperplane = kernel of (M - I): its normal vector = row space generator
            u,s,vh = np.linalg.svd(A); n = vh[-1].conj() if False else vh[0]  # A has rank 1: row space spanned by vh[0]
            n = n/ n[np.argmax(abs(n))]
            hk = tuple(complex(round(c.real,5),round(c.imag,5)) for c in n)
            hyps.setdefault(hk, set()).add(k)
    deg = 2
    for hk, ks in hyps.items():
        e_R = len(ks) + 1   # stabiliser = {identity} u {reflections fixing this hyperplane}
        deg -= (e_R - 1)
    return deg, {hk: len(ks)+1 for hk,ks in hyps.items()}
rows=[]
for sk,S in subs.items():
    iv = inv_dim(S)
    dg, hy = can_deg(S)
    rows.append({'order': len(S), 'inv_dim': iv, 'can_deg': dg, 'n_refl_hyperplanes': len(hy), 'gens_idx': list(sk)[:4]})
serving = [r_ for r_ in rows if r_['inv_dim']>0]
print('subgroups serving the class (rho_full^K != 0):', len(serving))
from collections import Counter
print('their (order, can_deg) multiset:', sorted(Counter((r_['order'], r_['can_deg']) for r_ in serving).items()))
neg = [r_ for r_ in serving if r_['can_deg'] < 0]
print('serving subgroups with negative canonical degree (uniruled candidates):', len(neg))
print('all subgroups: (order, can_deg, inv_dim) multiset:', sorted(Counter((r_['order'], r_['can_deg'], r_['inv_dim']) for r_ in rows).items()))
json.dump({'order_G': len(G), 'n_pseudo_reflections': nrefl, 'n_subgroups': len(subs), 'serving': serving,
           'n_serving_negative': len(neg)}, open('/tmp/hodge/results/cyclic-triple-residual/wave2-construct/c2_galois_closure_group.json','w'), indent=1)
