"""C4: Griffiths-Steenbrink Hodge numbers of H^4(X)_prim^K for every subgroup K of Aut(X) = D6 = <h,tau>,
and h^{4,0}, K3-sector survival for every subgroup of G = (mu6 x mu6) x| Z/2 acting on S6 x S6.
Exact character bookkeeping from Hilb(R_f) = (1,3,6,10,15,18,19,18,15,10,6,3,1) (R_F = R_f (x) R_f in the (U,V) frame).
Residue action: h^k: U^al V^be Omega/F^(p+1) has character eta^{k(|be|+3)}; h^k tau: maps (a,b)-block to (b,a)-block with scalar -eta^{k(|be|+3)}."""
import json, cmath, itertools
hil=[1,3,6,10,15,18,19,18,15,10,6,3,1]
eta=cmath.exp(2j*cmath.pi/6)
def tr_hk(k,p):   # trace of h^k on R_{F,6p} residue classes
    return sum(hil[a]*hil[b]*eta**(k*(b+3)) for a in range(13) for b in range(13) if a+b==6*p)
def tr_hktau(k,p):  # only diagonal block a=b=3p
    n=hil[3*p]; return -n*eta**(k*(3*p+3))
# elements of D6: ('h',k) or ('t',k) meaning h^k tau
elems=[('h',k) for k in range(6)]+[('t',k) for k in range(6)]
def mul(x,y):
    # h^a h^b = h^{a+b}; h^a (h^b tau) = h^{a+b} tau; (h^a tau) h^b = h^{a-b} tau; (h^a tau)(h^b tau) = h^{a-b}
    (s,a),(t,b)=x,y
    if s=='h' and t=='h': return ('h',(a+b)%6)
    if s=='h' and t=='t': return ('t',(a+b)%6)
    if s=='t' and t=='h': return ('t',(a-b)%6)
    return ('h',(a-b)%6)
def closure(gens):
    S={('h',0)}|set(gens); changed=True
    while changed:
        changed=False
        for x in list(S):
            for y in list(S):
                z=mul(x,y)
                if z not in S: S.add(z); changed=True
    return frozenset(S)
subs=set()
for g1 in elems:
    for g2 in elems:
        subs.add(closure([g1,g2]))
subs=sorted(subs,key=lambda K:(len(K),sorted(K)))
def tr(g,p):
    return tr_hk(g[1],p) if g[0]=='h' else tr_hktau(g[1],p)
def name(K):
    gens=sorted(K); return "{"+",".join(("h^%d"%k if s=='h' else "h^%dtau"%k) for s,k in gens)+"}"
rows=[]
print("subgroup K of D6 | order | h^{4,0} h^{3,1} h^{2,2} h^{1,3} h^{0,4} of X/K (invariant primitive part) | K3 (3,1) classes retained | junior elt (uniruled) ")
for K in subs:
    hodge=[]
    for p in range(5):
        v=sum(tr(g,p) for g in K)/len(K); hodge.append(round(v.real))
    # K3 sector retained iff K in ker(rho), rho(h)=-1, rho(tau)=+1  -> rho(h^k tau)=(-1)^k, rho(h^k)=(-1)^k
    k3=all(((k%2)==0) for s,k in K)
    junior=any(s=='h' and k in (1,5) for s,k in K)
    rows.append({"K":name(K),"order":len(K),"hodge":hodge,"K3_retained":k3,"contains_junior":junior})
    print(f"{name(K):40s} | {len(K):2d} | {hodge} | {'yes(19)' if k3 else 'no'} | {'yes -> uniruled (rational)' if junior else 'no -> K-trivial canonical, not uniruled'}")
# ---- cover S6 x S6, group G = (mu6 x mu6) x| Z2 ; H^2(S6)_prim = sum H_k, (h20,h11,h02): 
H=[None,(6,15,0),(3,18,0),(1,19,1),(0,18,3),(0,15,6)]
# elements (a,b,c); product: (a,b,c)(a',b',c') = (a + a'' , b + b'', c+c') where flip conjugates: flip (a',b') flip = (b',a')
def gmul(x,y):
    a,b,c=x; a2,b2,c2=y
    if c==1: a2,b2=b2,a2
    return ((a+a2)%6,(b+b2)%6,(c+c2)%2)
G=[(a,b,c) for a in range(6) for b in range(6) for c in range(2)]
def gclosure(gens):
    S={(0,0,0)}|set(gens); changed=True
    while changed:
        changed=False
        for x in list(S):
            for y in list(S):
                z=gmul(x,y)
                if z not in S: S.add(z); changed=True
    return frozenset(S)
gsubs=set()
for g1 in G:
    for g2 in G:
        gsubs.add(gclosure([g1,g2]))
# add 3-generator subgroups
more=set()
for K in list(gsubs):
    for g3 in G:
        if g3 not in K: more.add(gclosure(list(K)+[g3]))
gsubs|=more
def tr40(g):  # trace on H^{2,0}(S6)_prim (x) H^{2,0}(S6)_prim = sum_{k,l} h20_k h20_l eta^{ka+lb} (c=0) ; c=1: only k=l: h20_k eta^{k(a+b)}
    a,b,c=g
    if c==0: return sum(H[k][0]*H[l][0]*eta**(k*a+l*b) for k in range(1,6) for l in range(1,6))
    return sum(H[k][0]*eta**(k*(a+b)) for k in range(1,6))
def rho(g): a,b,c=g; return (-1)**(a+b+c)
summary={}
for K in gsubs:
    h40=round((sum(tr40(g) for g in K)/len(K)).real)
    k3=all(rho(g)==1 for g in K)
    hasflip=any(g[2]==1 for g in K)
    key=(len(K),h40,k3,hasflip)
    summary[key]=summary.get(key,0)+1
print("\nSubgroups of G=(mu6 x mu6) x| Z2 (order 72) acting on S6 x S6: (order, h^{4,0} of quotient, K3 sector retained, contains flip-type element) -> count")
for key in sorted(summary): print("  ",key,"->",summary[key])
print("Subgroups retaining the K3 sector with h^{4,0}=0:", [k for k in summary if k[2] and k[1]==0])
print("Total subgroups enumerated:",len(gsubs))
json.dump({"D6_rows":rows,"G72_summary":{str(k):v for k,v in summary.items()}},open('/tmp/hodge/results/polarized-k3-remainder/construct/c4_quotients.json','w'),indent=1)
