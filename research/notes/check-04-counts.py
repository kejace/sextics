# Independent re-check of the integer claims in docs/04-known-locus-and-frontier.md
from itertools import product
from math import comb
from collections import Counter

box = list(product(range(5), repeat=6))          # Fermat Jacobian ring monomial basis
def deg(e): return sum(e)
def ydeg(e): return sum(e[3:])

# 4.1 dim S6^+ / S6^-
plus  = sum(comb(8-j,2)*comb(j+2,2) for j in (0,2,4,6))
minus = sum(comb(8-j,2)*comb(j+2,2) for j in (1,3,5))
print("dim S6^+ =", plus, " dim S6^- =", minus)          # 236, 226

# Hilbert/trace series of Fermat Jacobian ring
def poly_pow(p, n):
    r = {0:1}
    for _ in range(n):
        s = {}
        for a,ca in r.items():
            for b,cb in p.items():
                s[a+b] = s.get(a+b,0)+ca*cb
        r = s
    return r
H = poly_pow({i:1 for i in range(5)}, 6)
Th = poly_pow({2*i:1 for i in range(5)}, 3)
for d in (6,12,18):
    print(f"d={d}: [H]={H[d]} [Theta]={Th.get(d,0)}  R^+={(H[d]+Th.get(d,0))//2}  R^-={(H[d]-Th.get(d,0))//2}")
# direct monomial count in Fermat ring, Y-odd
print("Fermat Y-odd counts d=6,12,18:", [sum(1 for e in box if deg(e)==d and ydeg(e)%2==1) for d in (6,12,18)])
print("odd Jacobian directions 226-18 =", minus-18)

# 4.1 table coefficient counts
coordwise_even = comb(8,3)                                   # cubics in six squares: 56
threeY = 28+45+36+10
print("coordinatewise-even coeffs:", coordwise_even, " three-Y-reflection coeffs:", threeY)
def count(cond):
    return sum(1 for e in product(range(7), repeat=6) if deg(e)==6 and ydeg(e)%2==0 and cond(e))
print("check threeY by enumeration:", count(lambda e: all(x%2==0 for x in e[3:])))
print("coordwise-even by enumeration:", count(lambda e: all(x%2==0 for x in e)))
print("one-reflection (even in Y2) coeffs:", count(lambda e: e[5]%2==0))          # 158
print("two-reflection (even in X0,Y0) coeffs:", count(lambda e: e[0]%2==0 and e[3]%2==0))  # 104
print("paired-sign (each X_i+Y_i even) coeffs:", count(lambda e: all((e[i]+e[i+3])%2==0 for i in range(3))))  # 74
print("paired-sign closed form:", sum((2*k0+1)*(2*k1+1)*(2*k2+1)+1 for k0 in range(4) for k1 in range(4) for k2 in range(4) if k0+k1+k2==3)//2)
# Gamma-invariant (paired signs + weights (2,1,0,3,2,1) mod 6) : expect 26
w6=(2,1,0,3,2,1)
print("Gamma-invariant sextic coeffs:", count(lambda e: all((e[i]+e[i+3])%2==0 for i in range(3)) and sum(a*b for a,b in zip(e,w6))%6==0))
# cyclic-7 family: sigma weights (1,2,4,0,0,0), Y-even
w7=(1,2,4,0,0,0)
print("cyclic-7 invariant Y-even sextic coeffs:", count(lambda e: sum(a*b for a,b in zip(e,w7))%7==0))   # 50
# 4.2 two-reflection raw q=1 targets: Y-odd, odd in X0 and Y0
raw = sum(1 for e in product(range(7), repeat=6) if deg(e)==6 and ydeg(e)%2==1 and e[0]%2==1 and e[3]%2==1)
print("two-reflection raw q=1 targets:", raw, " primitive:", raw-2)   # 34, 32
# conic construction weighted sextic in P(2,2,2,2,1,1)
print("weighted sextic coeffs P(2,2,2,2,1,1):", sum(comb((6-m)//2+3,3)*(m+1) for m in (0,2,4,6)))  # 77

# 4.2 embedded snippet (22,74,22)
def invariant(e):
    return (all((e[i]+e[i+3]) % 2 == 0 for i in range(3)) and sum(a*b for a, b in zip(e, w6)) % 6 == 3)
print("paired-sign sector Hodge dims:", [sum(sum(e) == d and invariant(e) for e in box) for d in (6, 12, 18)])

# 4.3 cyclic-7 sector dims via Fermat-ring count (sigma-trivial numerator, Y-odd)
print("cyclic-7 Fermat-ring Y-odd sigma-trivial dims:", [sum(1 for e in box if deg(e)==d and ydeg(e)%2==1 and sum(a*b for a,b in zip(e,w7))%7==0) for d in (6,12,18)])
# 4.3 equivariant Koszul series: coefficient of t^{6q} u^0 v^1 in prod (1 - t^5 u^{-w_i} v^{e_i})/(1 - t u^{w_i} v^{e_i}), u^7=v^2=1
ex=(0,0,0,1,1,1)
def mul(A,B,tmax=12):
    C=Counter()
    for (ta,ua,va),ca in A.items():
        for (tb,ub,vb),cb in B.items():
            if ta+tb<=tmax: C[(ta+tb,(ua+ub)%7,(va+vb)%2)]+=ca*cb
    return C
ser=Counter({(0,0,0):1})
for i in range(6):
    geo=Counter({(t,(t*w7[i])%7,(t*ex[i])%2):1 for t in range(13)})
    num=Counter({(0,0,0):1,(5,(-w7[i])%7,ex[i]):-1})
    ser=mul(mul(ser,geo),num)
print("cyclic-7 Koszul series [t^6 u^0 v^1], [t^12 u^0 v^1]:", ser[(6,0,1)], ser[(12,0,1)])   # 19, 128
# 4.3 exponent matrix det
M=((2,1,-1),(-1,2,1),(1,-1,2))
det=(M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
print("exponent matrix det:", det)
# 4.3 the 19 numerator directions by enumeration (degree 6, Y-odd, sigma weight 0)
dirs=[e for e in product(range(7),repeat=6) if deg(e)==6 and ydeg(e)%2==1 and sum(a*b for a,b in zip(e,w7))%7==0]
print("cyclic-7 numerator directions:", len(dirs), Counter(ydeg(e) for e in dirs))

# 4.4 embedded snippet
def check(extra, anchor, expected):
    def parity(u): return sum(u[3:]) % 2
    G = [u for u in box if sum(u)==6 and parity(u)==0 and extra(u)]
    V = [u for u in box if sum(u)==6 and parity(u)==1 and extra(u)]
    W = {u for u in box if sum(u)==12 and parity(u)==1 and extra(u)}
    signatures, images = [], []
    for u in V:
        image = {tuple(x+y for x,y in zip(u,a)) for a in G if all(x+y<=4 for x,y in zip(u,a))}
        images.append(image)
        signatures.append(tuple(all(x+y<=4 for x,y in zip(u,a)) for a in G))
    assert (len(G),len(V),len(W)) == expected, (len(G),len(V),len(W))
    assert len(set(signatures)) == len(V)
    anchor_image = images[V.index(anchor)]
    assert all(image & anchor_image for image in images)
    assert set().union(*images) == W
    return expected
print("4.4 full check:", check(lambda u: True, (1,1,1,1,1,1), (218,208,866)))
print("4.4 one-reflection check:", check(lambda u: u[5]%2==0, (1,1,1,1,2,0), (144,126,492)))

# 4.7 surface numerics
print("(6,3) surface: p_g = h^0(O(3,1)) =", comb(5,2)*2, "; (3,3) surface: p_g = h^0(O(0,1)) =", 2)
chiO=3; K2=0; e=12*chiO-K2
print("(3,3) surface: chi(O)=3, e =", e, " b2 =", e-2, "; free Z/3 quotient: e =", e//3, " b2 =", e//3-2)
