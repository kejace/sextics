#!/usr/bin/env python3
"""Gamma-tower for cycle-q2: Gamma = G x| <R> (216), induced rep rho = Ind chi_A, averaging-ideal
test on every subgroup, equivariant Hodge numbers of H^4_prim via the equivariant Hilbert series
Prod_i (1 - (c/lam_i) t^5)/(1 - lam_i t) (partials transform by c/lam_i), Hodge trace
tr(h|H^{5-k,k-1}) = det(h) c^{-k} tr(h|R_{6k-6}).  Output JSON + text summary."""
import itertools, json, sys, cmath
import numpy as np
Z = cmath.exp(2j*cmath.pi/6)
# F as dict of exponent tuples -> coefficient
F = {(6,0,0,0,0,0):1,(0,6,0,0,0,0):1,(0,0,6,0,0,0):1,(0,0,0,6,0,0):1,(0,0,0,0,6,0):1,(0,0,0,0,0,6):1,
     (3,1,0,1,1,0):0.5,(0,3,1,0,1,1):0.5,(1,0,3,1,0,1):0.5}
A = (2,2,3,2,2,1)
perm_R = [1,2,0,4,5,3]   # z_i -> z_{R(i)} index map: (X0,X1,X2,Y0,Y1,Y2) -> (X1,X2,X0,Y1,Y2,Y0)
def pow_perm(p,a):
    q=list(range(6))
    for _ in range(a): q=[p[i] for i in q]
    return q
def apply(h, mono):
    """F(hz): h(z)_i = Z^{w_i} z_{p(i)}; monomial prod z_i^{e_i} -> Z^{sum e_i w_i} prod z_{p(i)}^{e_i}"""
    p,w = h
    e2=[0]*6
    for i,e in enumerate(mono): e2[p[i]]+=e
    return tuple(e2), sum(e*wi for e,wi in zip(mono,w))%6
def preserves(h):
    c=None
    for m,co in F.items():
        m2,ph = apply(h,m)
        if m2 not in F: return None
        ratio = (co/F[m2], ph)  # F∘h has coefficient co*Z^ph on m2; need = c*F[m2]
        if c is None: c=(ratio[0],ph)
        elif c!=(ratio[0],ph): return None
    return c[1]  # phase exponent of c (coefficient ratio must be 1)
def norm(h):
    p,w=h; s=w[5]; return (tuple(p), tuple((x-s)%6 for x in w))
def mul(h1,h2):
    p1,w1=h1; p2,w2=h2
    p=[p2[p1[i]] for i in range(6)]
    w=[(w1[i]+w2[p1[i]])%6 for i in range(6)]
    return norm((p,w))
els=[]
for a in range(3):
    p=pow_perm(perm_R,a)
    for w in itertools.product(range(6),repeat=5):
        h=norm((p,tuple(w)+(0,)))
        if preserves(h) is not None: els.append(h)
idx={h:i for i,h in enumerate(els)}
n=len(els); print("Gamma order",n, "diagonal",sum(1 for h in els if h[0]==list(range(6)) or h[0]==tuple(range(6))))
mt=[[idx[mul(els[i],els[j])] for j in range(n)] for i in range(n)]
e=idx[norm((list(range(6)),(0,)*6))]
inv=[next(j for j in range(n) if mt[i][j]==e) for i in range(n)]
def chi_mono(mono,h):
    """character exponent of monomial mono under h (only meaningful if h maps mono to itself up to scalar)"""
    m2,ph=apply(h,mono); return ph if m2==mono else None
def res_char(P,q,h):
    """residue character exponent of P Omega/F^{q+1}: chi_P + chi_vol - (q+1) chi_F; None if P moved"""
    cp=chi_mono(P,h)
    if cp is None: return None
    cv=sum(h[1])%6; cF=preserves(h)
    return (cp+cv-(q+1)*cF)%6
RA=apply(norm((perm_R,(0,)*6)),A)[0]; R2A=apply(norm((pow_perm(perm_R,2),(0,)*6)),A)[0]
print("A,RA,R2A",A,RA,R2A)
# induced character of rho = Ind_G^Gamma chi_A: on G: sum of the three residue chars; off G: 0
def chi_rho(h):
    if h[0]!=tuple(range(6)): return 0
    return sum(Z**res_char(P,2,h) for P in (A,RA,R2A))
Gidx=[i for i,h in enumerate(els) if h[0]==tuple(range(6))]
K=[i for i in Gidx if res_char(A,2,els[i])==0]
K1=[i for i in Gidx if res_char(RA,2,els[i])==0]; K2=[i for i in Gidx if res_char(R2A,2,els[i])==0]
Kcap=sorted(set(K)&set(K1)&set(K2))
print("|G|",len(Gidx),"|K|",len(K),"|K1|",len(K1),"|K2|",len(K2),"|Kcap|",len(Kcap))
# fibrewise group of the (c1,c2)-fibration: w3-w0 = w4-w1 = w5-w2
Gfib=[i for i in Gidx if (els[i][1][3]-els[i][1][0])%6==(els[i][1][4]-els[i][1][1])%6==(els[i][1][5]-els[i][1][2])%6]
print("|Gfib|",len(Gfib),"|Gfib&K|",len(set(Gfib)&set(K)))
# ---- subgroups: closure enumeration
def closure(S):
    S=set(S); S.add(e); frontier=list(S)
    while frontier:
        new=[]
        for a in list(S):
            for b in frontier:
                for c in (mt[a][b],mt[b][a]):
                    if c not in S: S.add(c); new.append(c)
        frontier=new
    return frozenset(S)
cyc=set(closure([i]) for i in range(n))
subs=set(cyc); frontier=list(subs)
while frontier:
    new=[]
    for H in frontier:
        for C in cyc:
            if C<=H: continue
            J=closure(H|C)
            if J not in subs: subs.add(J); new.append(J)
    frontier=new
subs=sorted(subs,key=lambda H:(len(H),sorted(H)))
print("subgroups",len(subs))
# ---- equivariant Hodge numbers
def mat(h):
    p,w=h; M=np.zeros((6,6),complex)
    for i in range(6): M[i,p[i]]=Z**w[i]
    return M
def hodge_trace(h):
    M=mat(h); lam=np.linalg.eigvals(M); det=np.linalg.det(M); c=Z**preserves(h)
    N=25
    ser=np.zeros(N,complex); ser[0]=1
    for l in lam:
        geo=np.array([l**k for k in range(N)]); ser=np.convolve(ser,geo)[:N]
        num=np.zeros(N,complex); num[0]=1; num[5]=-(c/l); ser=np.convolve(ser,num)[:N]
    return [det*c**(-k)*ser[6*k-6] for k in range(1,6)]
HT=[hodge_trace(h) for h in els]
def inv_dims(H):
    s=np.zeros(5,complex)
    for i in H: s+=np.array(HT[i])
    s/=len(H); return [int(round(x.real)) for x in s]
def isotypic_dims(H,chi):  # multiplicity-weighted dims of chi-isotypic part under subgroup H (chi: dict idx->value)
    s=np.zeros(5,complex)
    for i in H: s+=np.conj(chi[i])*np.array(HT[i])
    s/=len(H); return [int(round(x.real)) for x in s]
print("all:",inv_dims([e]))
print("G-invariant:",inv_dims(Gidx)," K-invariant:",inv_dims(K)," tau:",inv_dims(closure([idx[norm((list(range(6)),(0,0,0,3,3,3))) ]])))
chiA={i:(Z**res_char(A,2,els[i]) if i in Gidx else 0) for i in range(n)}
print("G chi_A sector:",isotypic_dims(Gidx,chiA))
rho={i:chi_rho(els[i]) for i in range(n)}
print("Gamma rho-isotypic (dim rho=3, so divide by 3 for multiplicity):",[x for x in isotypic_dims(range(n),rho)])
# averaging-ideal test on every subgroup
rows=[]
for H in subs:
    m=sum(rho[i] for i in H)/len(H)
    m=int(round(m.real))
    inG = all(els[i][0]==tuple(range(6)) for i in H)
    rows.append({"order":len(H),"in_G":inG,"rho_inv":m,"H4_inv":inv_dims(H),"contains_R":any(els[i][0]==tuple(perm_R) for i in H),
                 "subset_K":set(H)<=set(K),"subset_Kcap":set(H)<=set(Kcap),"gens":sorted(H)[:4]})
good=[r for r in rows if r["rho_inv"]>0]
print("subgroups with rho^H != 0:",len(good),"of",len(subs))
maxg=sorted(good,key=lambda r:-r["order"])[:12]
for r in maxg: print(r["order"],"inG" if r["in_G"] else "nonG","R" if r["contains_R"] else "-","rho^H",r["rho_inv"],"H4^H",r["H4_inv"],"⊆K",r["subset_K"])
json.dump({"gamma_order":n,"K":K,"Kcap":Kcap,"Gfib":Gfib,"subgroups":rows,"els":[(list(h[0]),list(h[1])) for h in els]},open("w2c_gamma.json","w"))
