"""S3: Galois closure of X -> P^4 for cyclic-triple: Gamma = Z/3 wr Z/2 = (Z/3)^2 x| Z/2 acting on
Ztilde = {X2^3 + Xb2^3 + tC = 0, X2^3 Xb2^3 = H} (CI (3,6) in P^6, K = O(2)).
Elements (a,b,s): X2 -> z^a X2, Xb2 -> z^b Xb2, then swap if s=1.  d=(0,1,0), e=(1,0,0), sigma=(0,0,1).
Divisorial fixed loci: <d>-type (fixes {Xb2=0}: inertia order 3 -> 2H), <e>-type ({X2=0}), involutions
sigma d^a e^{-a} (fix {Xb2 = z^a X2}: 1H).  Pullback of K_{Z/K} = O(2 - r_K).
Target irrep: rho = Ind (1,0) = span of the (Z/3)^2-characters (1,0),(0,1); V_1 (mu3-char zeta on X=Z/<d>)
is the (0,1)-line [d -> 1, e -> zeta]."""
import itertools, cmath, json
from fractions import Fraction
z = cmath.exp(2j*cmath.pi/3)
def mul(g,h):  # g=(a,b,s): first h then g? define action as substitution matrices on (X2,Xb2): g acts by M_g = D(a,b)*P^s
    a,b,s = g; c,d,t = h
    # compose: apply h then g.  h: (X2,Xb2) -> coordinates; represent as 2x2 monomial matrices.
    def M(a,b,s):
        D = [[z**a,0],[0,z**b]]
        P = [[0,1],[1,0]] if s else [[1,0],[0,1]]
        return [[sum(D[i][k]*P[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    A, B = M(a,b,s), M(c,d,t)
    C = [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    # decode C back to (a,b,s)
    if abs(C[0][1]) < 1e-9:
        s2 = 0; wa, wb = C[0][0], C[1][1]
    else:
        s2 = 1; wa, wb = C[0][1], C[1][0]
    def lg(w):
        for k in range(3):
            if abs(w - z**k) < 1e-9: return k
        raise ValueError
    return (lg(wa), lg(wb), s2)
G = [(a,b,s) for a in range(3) for b in range(3) for s in range(2)]
I = (0,0,0)
def closure(gens):
    S={I}; fr=[I]
    while fr:
        new=[]
        for x in fr:
            for g in gens:
                y=mul(x,g)
                if y not in S: S.add(y); new.append(y)
        fr=new
    return frozenset(S)
def inv(g):
    for h in G:
        if mul(g,h)==I: return h
subs=set()
for r in (1,2):
    for gens in itertools.combinations(G,r): subs.add(closure(gens))
def conj(H,g):
    gi=inv(g); return frozenset(mul(mul(g,h),gi) for h in H)
classes=[]; seen=set()
for H in sorted(subs,key=lambda H:(len(H),sorted(H))):
    if H in seen: continue
    cl={conj(H,g) for g in G}; seen|=cl; classes.append((H,len(cl)))
# fixed divisors: element g fixes a divisor iff it is an involution (s=1, a+b=0) or a pure d^k / e^k power
def ram(H):
    r=0
    # involutions: each fixes its own hyperplane; inertia group order 2 -> 1H each
    r += sum(1 for g in H if g[2]==1)
    # {Xb2=0}: inertia = {(0,b,0) in H}; {X2=0}: inertia {(a,0,0) in H}
    r += len([g for g in H if g[2]==0 and g[0]==0]) - 1
    r += len([g for g in H if g[2]==0 and g[1]==0]) - 1
    return r
# character of rho = Ind_{(Z/3)^2}^{Gamma} chi_(1,0): chi(a,b,s) = 0 if s=1 else z^a + z^b
def chi_rho(g):
    a,b,s=g
    return 0 if s else z**a + z**b
def chi_rho2(g):
    a,b,s=g
    return 0 if s else z**(2*a) + z**(2*b)
def chi_sgn(g): return -1 if g[2] else 1
rows=[]
names={}
for H,ncl in classes:
    gens=sorted(H)
    inv_rho = sum(chi_rho(g) for g in H)/len(H)
    inv_rho2 = sum(chi_rho2(g) for g in H)/len(H)
    inv_sgn = sum(chi_sgn(g) for g in H)/len(H)
    rows.append(dict(order=len(H), n_conj=ncl, elems=gens, r_K=ram(H), K_pullback_degree=2-ram(H),
                     inv_rho=round(inv_rho.real), inv_rho_conj=round(inv_rho2.real), inv_sgn=round(inv_sgn)))
for r in rows: print(r)
json.dump(rows,open('s3_cyclic.json','w'),indent=1)
print("subgroups with rho-invariants and K_pullback_degree<0 (uniruled):", [r for r in rows if r['inv_rho']>0 and r['K_pullback_degree']<0])
