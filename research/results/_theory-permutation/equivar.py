"""Equivariant Hodge numbers of primitive H^4 sectors via the Fermat Jacobian ring.
Valid because the equivariant Koszul Hilbert series is constant on the smooth G-invariant locus.
Residue action of g on a numerator monomial = det(g) * (polynomial action).
Primitive H^{3,1}=R_6, H^{2,2}=R_12, H^{1,3}=R_18 (standard-monomial pieces of C[Z]/(Z_i^5)).
"""
import itertools, numpy as np
from fractions import Fraction

NV=6
def monos(d):  # Fermat ring basis: exponents 0..4, sum d
    return [e for e in itertools.product(range(5),repeat=NV) if sum(e)==d]

def perm_sign_det(perm, sign):
    # det of signed permutation matrix M where M e_i = sign[i] e_{perm[i]}
    # permutation parity:
    seen=[False]*NV; par=1
    for i in range(NV):
        if not seen[i]:
            j=i; L=0
            while not seen[j]:
                seen[j]=True; j=perm[j]; L+=1
            if L%2==0: par=-par
    s=1
    for x in sign: s*=x
    return par*s

def trace_on_Rd(d, perm, sign, twist):
    # trace of (twist * polynomial action) on R_d. polynomial action: g.Z^e = (prod sign[i]^{e[i]}) Z^{e'},
    # e'[perm[i]]=e[i]. Fixed monomials: e with permuted exponent == e.
    tr=0
    for e in monos(d):
        ep=[0]*NV
        for i in range(NV): ep[perm[i]]=e[i]
        if tuple(ep)==e:
            s=1
            for i in range(NV): s*= sign[i]**e[i]
            tr+=s
    return tr*twist

def group_S3xS3_tau():
    # elements: (permX in S3 on {0,1,2}, permY in S3 on {3,4,5}, tauflag)
    import itertools as it
    S3=list(it.permutations(range(3)))
    elems=[]
    for pX in S3:
        for pY in S3:
            for tau in (0,1):
                perm=[0]*6
                for i in range(3): perm[i]=pX[i]
                for i in range(3): perm[3+i]=3+pY[i]
                sign=[1,1,1]+[(-1 if tau else 1)]*3
                elems.append(((pX,pY,tau),perm,sign))
    return elems

def group_diag_S3_tau():
    # DP: simultaneous S3 on pairs (Xi,Yi), plus tau
    import itertools as it
    S3=list(it.permutations(range(3)))
    elems=[]
    for p in S3:
        for tau in (0,1):
            perm=[0]*6
            for i in range(3):
                perm[i]=p[i]; perm[3+i]=3+p[i]
            sign=[1,1,1]+[(-1 if tau else 1)]*3
            elems.append(((p,tau),perm,sign))
    return elems

# characters of S3 (as perms of 3 letters): triv, sgn, std
def chi_S3(p, name):
    # cycle type
    seen=[False]*3; cyc=[]
    for i in range(3):
        if not seen[i]:
            j=i;L=0
            while not seen[j]: seen[j]=True;j=p[j];L+=1
            cyc.append(L)
    cyc=tuple(sorted(cyc))
    if name=='triv': return 1
    if name=='sgn':  return 1 if cyc in [(1,1,1),(3,)] else -1
    if name=='std':  # 2 - fix? std char: id->2, transposition->0, 3cyc->-1
        return {(1,1,1):2,(1,2):0,(3,):-1}[cyc]

def report_SSR():
    elems=group_S3xS3_tau()
    G=len(elems)
    # sectors: (rhoX,rhoY) in {triv,sgn,std}^2 ; tau parity handled by Y-odd restriction
    # We want Y-odd (tau-invariant residue). Build projector to Y-odd by averaging over tau with residue twist.
    names=['triv','sgn','std']
    dims={'triv':1,'sgn':1,'std':2}
    out={}
    for dlabel,d in [('h31',6),('h22',12),('h13',18)]:
        # precompute trace of residue action for each element on R_d
        trR={}
        for key,perm,sign in elems:
            det=perm_sign_det(perm,sign)
            trR[key]=trace_on_Rd(d,perm,sign,det)
        for rx in names:
            for ry in names:
                # multiplicity of (rx (x) ry) in Y-odd part:
                # Y-odd projector = (1/2)(1 - tau_res). tau_res trace included by summing over tau with the irrep extended trivially on tau? 
                # The residue is tau-invariant for Y-odd; we want the S3xS3 isotypic within Y-odd R_d.
                # mult = (1/|S3xS3|) sum_{pX,pY} chiRX(pX)chiRY(pY) * [ (1/2)(trR[(pX,pY,0)] - trR[(pX,pY,1)]) ]
                m=Fraction(0)
                for (key,perm,sign) in elems:
                    (pX,pY,tau)=key
                    yodd = Fraction(1,2)*( (1 if tau==0 else -1) )  # coefficient to pick Y-odd (res invariant): (1 - tau_res)/2 with tau_res value = trR already has det; but we built trR per element; combine
                # redo cleanly below
                out[(dlabel,rx,ry)]=None
    return out
print("loaded")
