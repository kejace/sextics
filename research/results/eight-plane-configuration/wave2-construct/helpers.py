import sys
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
from sympy import Poly, symbols
F,A,q,rec = core.load_instance('eight-plane-configuration')
lam=symbols('lam')
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
def disc4(c,p):
    a,b,cc,d,e=[x%p for x in c]
    return (256*a**3*e**3-192*a**2*b*d*e**2-128*a**2*cc**2*e**2+144*a**2*cc*d**2*e-27*a**2*d**4
   +144*a*b**2*cc*e**2-6*a*b**2*d**2*e-80*a*b*cc**2*d*e+18*a*b*cc*d**3+16*a*cc**4*e
   -4*a*cc**3*d**2-27*b**4*e**2+18*b**3*cc*d*e-4*b**3*d**3-4*b**2*cc**3*e+b**2*cc**2*d**2)%p
class Ctx:
    def __init__(self,p):
        self.p=p; self.i=sqrt_m1(p); self.Fp=Fmod(p)
        self.F60={e:c for e,c in self.Fp.items() if e[3]==e[4]==e[5]==0}
        self.F06={e:c for e,c in self.Fp.items() if e[0]==e[1]==e[2]==0}
    def qcoeffs(self,aa,bb,eps=(1,1,1)):
        P,Pp=plane_points(aa,bb,eps,self.i,self.p); tot=restrict_to_line(self.Fp,P,Pp,self.p)
        assert tot[0]==0 and tot[6]==0
        return tot[1:6]
    def Delta_at(self,aa,bb,eps=(1,1,1)): return disc4(self.qcoeffs(aa,bb,eps),self.p)
    def topoly(self,cf): return Poly(sum(int(c)*lam**k for k,c in enumerate(cf)),lam,modulus=self.p)
    def restrict_curve(self,G,pt_of_l,idx):
        p=self.p; ys=[]
        for l in range(1,8):
            pt=pt_of_l(l); v=0
            for e,c in G.items():
                t=c
                for k in range(3): t=t*pow(pt[k],e[idx[k]],p)%p
                v=(v+t)%p
            ys.append(v)
        return self.topoly(lagrange_mod(list(range(1,8)),ys,p))
