"""Quartic plane-pair covers X -> Pi_eps x Pi_{-eps}: residual binary quartic q, Galois group of a
Q(i)-specialisation (Frobenius cycle types at split primes), discriminant reducedness on a random
line (mod p), j-invariant variation. Exact over Q(i) where cheap, else mod p = 1 mod 4."""
import sys, json, random, time
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction
from sympy import Poly, symbols, factor_list, GF, gcd, diff
F,A,q,rec = core.load_instance('eight-plane-configuration')
random.seed(20261101)
x=symbols('x')
def sqrt_m1(p):
    for g in range(2,200):
        r=pow(g,(p-1)//4,p)
        if (r*r)%p==p-1: return r
def Fmod(p): return {e:(c.numerator*pow(c.denominator,-1,p))%p for e,c in F.items()}
def binmul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        if u:
            for j,v in enumerate(b): c[i+j]=(c[i+j]+u*v)%p
    return c
def binpow(a,k,p):
    r=[1]
    for _ in range(k): r=binmul(r,a,p)
    return r
def restrict_to_line(Fp,P,Pp,p):
    """coefficients (by t-degree) of F(s P + t P'), homogeneous degree 6."""
    tot=[0]*7
    for e,c in Fp.items():
        term=[c]
        for j in range(6):
            if e[j]: term=binmul(term,binpow([P[j]%p,Pp[j]%p],e[j],p),p)
        for k,v in enumerate(term): tot[k]=(tot[k]+v)%p
    return tot
def plane_points(a,b,eps,i,p):
    P=[a[0],a[1],a[2]]+[(eps[k]*i*a[k])%p for k in range(3)]
    Pp=[b[0],b[1],b[2]]+[(-eps[k]*i*b[k])%p for k in range(3)]
    return P,Pp
res={}
# --- Galois group of one specialisation: fixed integer a,b, many split primes
a=[3,-5,7]; b=[2,11,-4]; eps=(1,1,1)
patterns={}
primes=[1000000009,1000000021,1000000033,1000000087,1000000093,1000000097,1000000103,1000000123,1000000181,1000000207,1000000223,1000000241,1000000271,1000000289,1000000297,1000000321,1000000349,1000000363,1000000403,1000000409,1000000411,1000000427,1000000433,1000000439,1000000453,1000000459,1000000483,1000000513,1000000531,1000000579]
from sympy import isprime
primes=[p for p in primes if isprime(p) and p%4==1]
for p in primes:
    i=sqrt_m1(p); Fp=Fmod(p); P,Pp=plane_points(a,b,eps,i,p)
    tot=restrict_to_line(Fp,P,Pp,p)
    assert tot[0]==0 and tot[6]==0, 'not divisible by st'
    qc=tot[1:6]  # q(s,t)=sum qc[k] s^{4-k} t^k ; q(x,1): coefficient of x^{4-k} is qc[k]
    pol=Poly(sum(qc[k]*x**(4-k) for k in range(5)),x,modulus=p)
    if gcd(pol,pol.diff(x)).degree()>0: patterns.setdefault('ramified',0); patterns['ramified']+=1; continue
    fl=factor_list(pol.as_expr(),modulus=p)[1]
    pat=tuple(sorted([Poly(f,x).degree() for f,m in fl for _ in range(m)],reverse=True))
    patterns[pat]=patterns.get(pat,0)+1
res['galois_specialisation']={'a':a,'b':b,'eps':eps,'n_primes':len(primes),'cycle_types':{str(k):v for k,v in patterns.items()}}
# --- discriminant on a random line in Pi x Pi' (mod one prime): degree and squarefreeness
p=1000000009; i=sqrt_m1(p); Fp=Fmod(p)
lam=symbols('lam')
a0=[random.randrange(p) for _ in range(3)]; a1=[random.randrange(p) for _ in range(3)]
b0=[random.randrange(p) for _ in range(3)]; b1=[random.randrange(p) for _ in range(3)]
# coefficients of q as polynomials in lam of degree 6: evaluate at 7 points and interpolate via sympy
def lagrange_mod(xs,ys,p):
    n=len(xs); coeffs=[0]*n
    for j in range(n):
        num=[1]; den=1
        for m in range(n):
            if m==j: continue
            num=binmul(num,[(-xs[m])%p,1],p); den=den*(xs[j]-xs[m])%p
        f=ys[j]*pow(den,-1,p)%p
        for k in range(n): coeffs[k]=(coeffs[k]+f*num[k])%p
    return coeffs
pts=list(range(1,8)); coeff_polys=[]
vals=[]
for lv in pts:
    aa=[(a0[k]+lv*a1[k])%p for k in range(3)]; bb=[(b0[k]+lv*b1[k])%p for k in range(3)]
    P,Pp=plane_points(aa,bb,eps,i,p); tot=restrict_to_line(Fp,P,Pp,p); vals.append(tot[1:6])
for k in range(5):
    cf=lagrange_mod(pts,[v[k] for v in vals],p)
    coeff_polys.append(Poly(sum(int(cf[d])*lam**d for d in range(len(cf))),lam,modulus=p))
# discriminant of binary quartic sum c_k s^{4-k} t^k with c=[c0..c4] (a=c0,b=c1,c=c2,d=c3,e=c4)
c0,c1,c2,c3,c4=coeff_polys
D=(256*c0**3*c4**3-192*c0**2*c1*c3*c4**2-128*c0**2*c2**2*c4**2+144*c0**2*c2*c3**2*c4-27*c0**2*c3**4
   +144*c0*c1**2*c2*c4**2-6*c0*c1**2*c3**2*c4-80*c0*c1*c2**2*c3*c4+18*c0*c1*c2*c3**3+16*c0*c2**4*c4
   -4*c0*c2**3*c3**2-27*c1**4*c4**2+18*c1**3*c2*c3*c4-4*c1**3*c3**3-4*c1**2*c2**3*c4+c1**2*c2**2*c3**2)
D=Poly(D.as_expr(),lam,modulus=p)
g=gcd(D,D.diff(lam))
res['discriminant_on_line']={'p':p,'degree':D.degree(),'gcd_with_derivative_degree':g.degree(),'squarefree':g.degree()==0}
# --- j-invariant variation (I^3/J^2) at two random points
def IJ(c):
    a,b,cc,d,e=c
    I=(12*a*e-3*b*d+cc*cc)%p
    J=(72*a*cc*e+9*b*cc*d-27*a*d*d-27*e*b*b-2*cc**3)%p
    return I,J
js=[]
for _ in range(3):
    aa=[random.randrange(p) for _ in range(3)]; bb=[random.randrange(p) for _ in range(3)]
    P,Pp=plane_points(aa,bb,eps,i,p); tot=restrict_to_line(Fp,P,Pp,p); I,J=IJ(tot[1:6])
    js.append((pow(I,3,p)*pow(J*J%p,-1,p))%p)
res['I3_over_J2_three_points']=js
# --- twisted diagonal symmetry check: b=a gives q(s,t)=q(t,s)
aa=[random.randrange(p) for _ in range(3)]
P,Pp=plane_points(aa,aa,eps,i,p); tot=restrict_to_line(Fp,P,Pp,p)
res['twisted_diagonal_symmetric']=(tot[1:6]==tot[1:6][::-1])
print(json.dumps(res))
json.dump(res,open('/tmp/hodge/results/eight-plane-configuration/wave2-construct/quartic_cover.json','w'))
