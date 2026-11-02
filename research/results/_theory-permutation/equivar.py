"""Equivariant Hodge numbers of primitive H^4 sectors via the Fermat Jacobian ring.
Residue representation on numerators R_d = det(g) * polynomial action.
H^{3,1}=R_6, H^{2,2}=R_12, H^{1,3}=R_18 (Fermat ring C[Z]/(Z_i^5))."""
import itertools
from fractions import Fraction as Fr

NV=6
_MC={}
def monos(d):
    if d not in _MC:
        _MC[d]=[e for e in itertools.product(range(5),repeat=NV) if sum(e)==d]
    return _MC[d]

def det_signed(perm, sign):
    seen=[False]*NV; par=1
    for i in range(NV):
        if not seen[i]:
            j=i;L=0
            while not seen[j]: seen[j]=True;j=perm[j];L+=1
            if L%2==0: par=-par
    s=1
    for x in sign: s*=x
    return par*s

def tr_poly(d, perm, sign):
    tr=0
    for e in monos(d):
        ep=[0]*NV
        for i in range(NV): ep[perm[i]]=e[i]
        if tuple(ep)==e:
            s=1
            for i in range(NV): s*=sign[i]**e[i]
            tr+=s
    return tr

def cyc_type(p,n):
    seen=[False]*n;cyc=[]
    for i in range(n):
        if not seen[i]:
            j=i;L=0
            while not seen[j]: seen[j]=True;j=p[j];L+=1
            cyc.append(L)
    return tuple(sorted(cyc))
def chiS3(p,name):
    ct=cyc_type(p,3)
    if name=='triv':return 1
    if name=='sgn': return {(1,1,1):1,(1,2):-1,(3,):1}[ct]
    if name=='std': return {(1,1,1):2,(1,2):0,(3,):-1}[ct]

def SSR_elems():
    S3=list(itertools.permutations(range(3)))
    E=[]
    for pX in S3:
        for pY in S3:
            for tau in (0,1):
                perm=[0]*6; sign=[1]*6
                for i in range(3): perm[i]=pX[i]
                for i in range(3): perm[3+i]=3+pY[i]
                for i in range(3,6): sign[i]=-1 if tau else 1
                E.append((pX,pY,tau,perm,sign))
    return E

def DP_elems():
    S3=list(itertools.permutations(range(3)))
    E=[]
    for p in S3:
        for tau in (0,1):
            perm=[0]*6; sign=[1]*6
            for i in range(3): perm[i]=p[i]; perm[3+i]=3+p[i]
            for i in range(3,6): sign[i]=-1 if tau else 1
            E.append((p,tau,perm,sign))
    return E

def SSR_table():
    E=SSR_elems(); G=72
    res={}
    for dl,d in [('h31',6),('h22',12),('h13',18)]:
        trres={}
        for (pX,pY,tau,perm,sign) in E:
            trres[(pX,pY,tau)]=det_signed(perm,sign)*tr_poly(d,perm,sign)
        for rx in ['triv','sgn','std']:
            for ry in ['triv','sgn','std']:
                for et in (+1,-1):
                    m=Fr(0)
                    for (pX,pY,tau,perm,sign) in E:
                        chi=chiS3(pX,rx)*chiS3(pY,ry)*(et**tau)
                        m+=chi*trres[(pX,pY,tau)]
                    m=m/G
                    assert m.denominator==1, (dl,rx,ry,et,m)
                    res[(dl,rx,ry,et)]=int(m)
    return res

def DP_table():
    E=DP_elems(); G=12
    res={}
    for dl,d in [('h31',6),('h22',12),('h13',18)]:
        trres={}
        for (p,tau,perm,sign) in E:
            trres[(p,tau)]=det_signed(perm,sign)*tr_poly(d,perm,sign)
        for r in ['triv','sgn','std']:
            for et in (+1,-1):
                m=Fr(0)
                for (p,tau,perm,sign) in E:
                    chi=chiS3(p,r)*(et**tau)
                    m+=chi*trres[(p,tau)]
                m=m/G
                assert m.denominator==1,(dl,r,et,m)
                res[(dl,r,et)]=int(m)
    return res

if __name__=='__main__':
    print("=== SSR  G=S3(X)xS3(Y)x<tau>, residue sector = (rhoX,rhoY, tau-residue-parity) ===")
    T=SSR_table()
    print("sector (rhoX,rhoY,tauRes)  (h31,h22,h13)   [tauRes=+ means residue tau-invariant = target type]")
    for rx in ['triv','sgn','std']:
        for ry in ['triv','sgn','std']:
            for et in (+1,-1):
                tup=(T[('h31',rx,ry,et)],T[('h22',rx,ry,et)],T[('h13',rx,ry,et)])
                if any(tup):
                    mult={'triv':1,'sgn':1,'std':2}
                    irdim=mult[rx]*mult[ry]
                    print(f"  ({rx:4},{ry:4},{'+' if et>0 else '-'})  {tup}   irrep dim/block={irdim}")
    print("\n=== DP  G=S3(diag)x<tau>, residue sector = (rho, tau-residue-parity) ===")
    D=DP_table()
    for r in ['triv','sgn','std']:
        for et in (+1,-1):
            tup=(D[('h31',r,et)],D[('h22',r,et)],D[('h13',r,et)])
            if any(tup):
                print(f"  ({r:4},{'+' if et>0 else '-'})  {tup}")

def verify_A_character():
    # A = X0X1X2Y0Y1Y2 = monomial (1,1,1,1,1,1). Residue action of g=(pX,pY,tau):
    # det(g)*poly(g on A). poly on A: A invariant under any perm; sign = prod sign^1 = (-1)^(3tau).
    # det = sgn(pX)sgn(pY)*(-1)^(3tau). product = sgn(pX)sgn(pY). Independent of tau => char (sgn,sgn,+).
    import itertools
    for pX in itertools.permutations(range(3)):
        for pY in itertools.permutations(range(3)):
            for tau in (0,1):
                perm=[0]*6;sign=[1]*6
                for i in range(3):perm[i]=pX[i]
                for i in range(3):perm[3+i]=3+pY[i]
                for i in range(3,6):sign[i]=-1 if tau else 1
                det=det_signed(perm,sign)
                polyA=(-1)**(3*tau)  # sign factor on A
                resA=det*polyA
                want=chiS3(pX,'sgn')*chiS3(pY,'sgn')  # (sgn,sgn), tau-trivial
                assert resA==want,(pX,pY,tau,resA,want)
    print("VERIFIED: residue character of A = X0X1X2Y0Y1Y2 is exactly (sgn_X (x) sgn_Y, tau:+).")
    print("=> target lives in sector (sgn,sgn,+) = (17,52,17).")
    print("   Delta_X*Delta_Y has residue char (triv,triv,+) = (1,10,1) [served by symmetric Fano].")
