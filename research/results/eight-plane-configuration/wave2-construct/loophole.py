"""Loophole check for (H5): (1) factorisation of Delta restricted to the twisted diagonal {(p,tau p)}
along a random line in Pi (degree 36 univariate): expected F60(X(p))^2 * F06(Y(p))^2 * R^2;
(2) multiplicity of Delta at a generic point of the curve C_X' = {X(p) in C_X} along a random direction of B
(expected 1: Delta smooth there). All mod p = 10^9+9."""
import sys, json, random
sys.path.insert(0,'/home/kejace/devel/sextics-work')
sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration/wave2-construct')
from hodge import core
from sympy import Poly, symbols, gcd, factor_list, GF
exec(open('/tmp/hodge/results/eight-plane-configuration/wave2-construct/quartic_cover.py').read().split("res={}")[0])  # reuse helpers
p=1000000009; i=sqrt_m1(p); Fp=Fmod(p); random.seed(7)
lam=symbols('lam')
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
def disc4(c):
    a,b,cc,d,e=[x%p for x in c]
    return (256*a**3*e**3-192*a**2*b*d*e**2-128*a**2*cc**2*e**2+144*a**2*cc*d**2*e-27*a**2*d**4
   +144*a*b**2*cc*e**2-6*a*b**2*d**2*e-80*a*b*cc**2*d*e+18*a*b*cc*d**3+16*a*cc**4*e
   -4*a*cc**3*d**2-27*b**4*e**2+18*b**3*cc*d*e-4*b**3*d**3-4*b**2*cc**3*e+b**2*cc**2*d**2)%p
def Delta_at(aa,bb):
    P,Pp=plane_points(aa,bb,(1,1,1),i,p); tot=restrict_to_line(Fp,P,Pp,p); return disc4(tot[1:6])
def topoly(cf): return Poly(sum(int(c)*lam**k for k,c in enumerate(cf)),lam,modulus=p)
res={}
# (1) Delta on twisted diagonal along the line p(l) = a0 + l*a1 in Pi
a0=[random.randrange(p) for _ in range(3)]; a1=[random.randrange(p) for _ in range(3)]
xs=list(range(1,40)); ys=[Delta_at([(a0[k]+l*a1[k])%p for k in range(3)],[(a0[k]+l*a1[k])%p for k in range(3)]) for l in xs]
D=topoly(lagrange_mod(xs,ys,p))
F60={e:c for e,c in Fp.items() if e[3]==e[4]==e[5]==0}; F06={e:c for e,c in Fp.items() if e[0]==e[1]==e[2]==0}
def restrict_curve(G,pt_of_l,idx):   # G(pt) as univariate in l, degree 6
    ys=[]
    for l in range(1,8):
        pt=pt_of_l(l); v=0
        for e,c in G.items():
            t=c
            for k in range(3): t=t*pow(pt[k],e[idx[k]],p)%p
            v=(v+t)%p
        ys.append(v)
    return topoly(lagrange_mod(list(range(1,8)),ys,p))
f60=restrict_curve(F60,lambda l:[(a0[k]+l*a1[k])%p for k in range(3)],(0,1,2))
f06=restrict_curve(F06,lambda l:[(i*(a0[k]+l*a1[k]))%p for k in range(3)],(3,4,5))
res['deg_Delta_on_twisted_diag']=D.degree()
q1,r1=D.div(f60**2); q2,r2=q1.div(f06**2)
res['divisible_by_F60^2']=(r1.is_zero); res['divisible_by_F06^2']=(r2.is_zero)
g=gcd(q2,q2.diff(lam)); res['remaining_factor_degree']=q2.degree(); res['remaining_is_square']=(2*g.degree()==q2.degree())
res['remaining_sqrt_squarefree']=(gcd(g,g.diff(lam)).degree()==0)
# (2) multiplicity of Delta at a point of C_X' along a random direction of B
# find p0 in Pi with F60(X(p0))=0: solve univariate in X2
x=symbols('x'); found=None
for trial in range(30):
    c0=random.randrange(p); c1=random.randrange(p)
    ys=[]; 
    for l in range(1,8):
        v=0
        for e,c in F60.items(): v=(v+c*pow(c0,e[0],p)*pow(c1,e[1],p)*pow(l,e[2],p))%p
        ys.append(v)
    pol=topoly(lagrange_mod(list(range(1,8)),ys,p))
    roots=[r for r,_ in Poly(pol.as_expr(),lam,modulus=p).ground_roots().items()] if hasattr(pol,'ground_roots') else []
    if not roots:
        # fall back: factor and take linear factors
        fl=factor_list(pol.as_expr(),modulus=p)[1]
        lin=[f for f,m in fl if Poly(f,lam).degree()==1]
        if lin:
            f=Poly(lin[0],lam,modulus=p); r=(-f.coeffs()[1]*pow(f.coeffs()[0],-1,p))%p; roots=[r]
    if roots:
        found=[c0,c1,int(roots[0])%p]; break
res['point_on_CX']=found
if found:
    a=found; dirs=[[random.randrange(p) for _ in range(3)] for _ in range(2)]  # (xi, eta): b(l) = (a + l xi, a_conj + l eta)
    # tau p = plane_points with b=a gives (a, -i a); generic direction in B: move p and p' independently
    xs=list(range(1,40)); ys=[Delta_at([(a[k]+l*dirs[0][k])%p for k in range(3)],[(a[k]+l*dirs[1][k])%p for k in range(3)]) for l in xs]
    Dl=topoly(lagrange_mod(xs,ys,p))
    # valuation at l=0
    cf=Dl.all_coeffs()[::-1]; v=0
    while v<len(cf) and cf[v]%p==0: v+=1
    res['Delta_at_point']=int(Delta_at(a,a)); res['multiplicity_generic_direction']=v; res['deg_Dl']=Dl.degree()
    # also along the twisted diagonal direction only
    ys=[Delta_at([(a[k]+l*dirs[0][k])%p for k in range(3)],[(a[k]+l*dirs[0][k])%p for k in range(3)]) for l in xs]
    Dl2=topoly(lagrange_mod(xs,ys,p)); cf=Dl2.all_coeffs()[::-1]; v2=0
    while v2<len(cf) and cf[v2]%p==0: v2+=1
    res['multiplicity_along_twisted_diagonal']=v2
print(json.dumps(res))
json.dump(res,open('/tmp/hodge/results/eight-plane-configuration/wave2-construct/loophole.json','w'))
