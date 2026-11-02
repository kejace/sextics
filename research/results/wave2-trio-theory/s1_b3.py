"""S1: hyperoctahedral group B3 = Z/2 wr S3 acting on the Galois closure Ztilde = CI(2,4,6) in P^7
(coordinates x0..x4, y1,y2,y3 with y_i^2 = the three roots of the cubic in the fibre coordinate).
For every subgroup K up to conjugacy: order, number of reflections n_K (elements fixing a hyperplane
section of Ztilde pointwise), pullback of K_{Z/K} = O(4 - n_K), and dim of K-invariants in the irreps
relevant to X = Z/K_X: (3,0) trivial, (21,0) [std of S3, pulled back: the r-even companion sector of W],
(2,1) [the sign-vector rep: the r-odd target sector]."""
import itertools, json
from fractions import Fraction
# elements: (perm, signs) acting on (y1,y2,y3): y_i -> s_i * y_{perm(i)}  (as signed permutation matrices)
def mat(perm, signs):
    M = [[0]*3 for _ in range(3)]
    for i in range(3):
        M[perm[i]][i] = signs[i]
    return tuple(tuple(r) for r in M)
G = []
for perm in itertools.permutations(range(3)):
    for signs in itertools.product([1,-1], repeat=3):
        G.append(mat(perm, signs))
G = sorted(set(G))
assert len(G) == 48
def mul(A,B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
I = tuple(tuple(1 if i==j else 0 for j in range(3)) for i in range(3))
def inv(A):
    return tuple(tuple(A[j][i] for j in range(3)) for i in range(3))  # orthogonal
def order(A):
    k, B = 1, A
    while B != I:
        B = mul(B, A); k += 1
    return k
# reflections: det = -1 and fixed space of dimension 2, i.e. A has eigenvalue -1 with multiplicity 1 -> trace = 1
def det(A):
    return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))
def trace(A): return A[0][0]+A[1][1]+A[2][2]
refl = [A for A in G if order(A)==2 and det(A)==-1 and trace(A)==1]
assert len(refl) == 9
# subgroups: generate from all pairs/triples of generators (B3 subgroups need <=3 generators? use closure of random subsets up to size 3)
def closure(gens):
    S = {I}
    frontier = [I]
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                b = mul(a, g)
                if b not in S:
                    S.add(b); new.append(b)
        frontier = new
    return frozenset(S)
subs = set()
for r in (1,2,3):
    for gens in itertools.combinations(G, r):
        subs.add(closure(gens))
# conjugacy classes of subgroups
def conj(H, g):
    gi = inv(g)
    return frozenset(mul(mul(g,h),gi) for h in H)
classes = []
seen = set()
for H in sorted(subs, key=lambda H:(len(H), sorted(H))):
    if H in seen: continue
    cl = {conj(H,g) for g in G}
    seen |= cl
    classes.append((H, len(cl)))
print("number of subgroups:", len(subs), " conjugacy classes:", len(classes))
# irreps as explicit matrices:
# (2,1): sign-vector rep: basis v1,v2,v3; A acts by the signed permutation matrix itself (3-dim). 
# (3,0): trivial. (21,0): std of S3 via perm part. (1,2): A -> det-twisted? we only need (2,1) and (21,0).
def inv_dim(H, rep):
    # dim of invariants = (1/|H|) sum trace(rep(h))
    s = sum(trace(rep(h)) for h in H)
    return Fraction(s, len(H))
def rep_21(A):  # (2,1): the natural signed-permutation matrix
    return A
def perm_part(A):
    return tuple(tuple(abs(x) for x in r) for r in A)
def rep_std(A):  # (21,0): permutation rep minus trivial: trace = (#fixed points of perm) - 1 ; we fake a matrix with that trace
    P = perm_part(A)
    t = trace(P) - 1
    return ((t,0,0),(0,0,0),(0,0,0))
def rep_sgnperm(A): # (111,0): sign of permutation
    P = perm_part(A)
    # sign of perm
    perm = [next(j for j in range(3) if P[j][i]) for i in range(3)]
    s = 1
    for i in range(3):
        for j in range(i+1,3):
            if perm[i] > perm[j]: s = -s
    return ((s,0,0),(0,0,0),(0,0,0))
def rep_prodsign(A): # (0,3): product of the signs
    return ((det(A)*trace(rep_sgnperm(A)[0:1][0][0:1]) if False else 0,0,0),(0,0,0),(0,0,0))
rows = []
KX = closure([mat((0,2,1),(1,1,1)), mat((0,1,2),(1,-1,1)), mat((0,1,2),(1,1,-1))])  # <(23), e2, e3>
for H, ncl in classes:
    nref = sum(1 for h in H if h in refl)
    d21 = inv_dim(H, rep_21); dstd = inv_dim(H, rep_std)
    rows.append(dict(order=len(H), n_conj=ncl, n_refl=nref, K_pullback_degree=4-nref,
                     inv_21=int(d21), inv_std=int(dstd), is_KX=(H==KX), is_KW=False))
rows.sort(key=lambda r:(r['order'], -r['n_refl']))
for r in rows: print(r)
best = [r for r in rows if r['inv_21']>0]
print("max n_refl among subgroups with (2,1)-invariants:", max(r['n_refl'] for r in best))
print("subgroups with n_refl>=5 (uniruled quotients) and inv_21>0:", [r for r in rows if r['n_refl']>=5 and r['inv_21']>0])
json.dump(rows, open('s1_b3.json','w'), indent=1)
