"""Exact check over Q(i): which block-monomial maps g (permutation of the six coordinates u_j=X_j+iY_j,
v_j=X_j-iY_j preserving the pairing {u_j,v_j}, composed with a diagonal scaling) satisfy F o g = lambda F.
If the 8 planes are the only planes on X, Aut(X) is contained in this family (Part 1, section 4)."""
import sys, json, itertools
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction as Fr
F,A,q,rec=core.load_instance('eight-plane-configuration')
# Gaussian rationals as (re,im)
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gadd(a,b): return (a[0]+b[0],a[1]+b[1])
def ginv(a):
    n=a[0]*a[0]+a[1]*a[1]; return (a[0]/n,-a[1]/n)
def gpow(a,k):
    r=(Fr(1),Fr(0))
    for _ in range(k): r=gmul(r,a)
    return r
ZERO=(Fr(0),Fr(0))
# F in (u,v): X_j=(u_j+v_j)/2, Y_j=(u_j-v_j)/(2i) = -i(u_j-v_j)/2 : coefficients of u_j,v_j in X_j: (1/2,1/2); in Y_j: (-i/2, i/2)
def pmul(a,b):
    c={}
    for ea,ca in a.items():
        for eb,cb in b.items():
            e=tuple(x+y for x,y in zip(ea,eb)); c[e]=gadd(c.get(e,ZERO),gmul(ca,cb))
    return {e:v for e,v in c.items() if v!=ZERO}
def lin(j,coef_u,coef_v):
    e1=[0]*6; e1[j]=1; e2=[0]*6; e2[3+j]=1   # order: u0,u1,u2,v0,v1,v2
    return {tuple(e1):coef_u, tuple(e2):coef_v}
half=(Fr(1,2),Fr(0)); mih=(Fr(0),Fr(-1,2)); pih=(Fr(0),Fr(1,2))
Xl=[lin(j,half,half) for j in range(3)]; Yl=[lin(j,mih,pih) for j in range(3)]
pw={}
def power(L,k,key):
    if (key,k) in pw: return pw[(key,k)]
    r={tuple([0]*6):(Fr(1),Fr(0))}
    for _ in range(k): r=pmul(r,L)
    pw[(key,k)]=r; return r
Ft={}
for e,c in F.items():
    term={tuple([0]*6):(Fr(c),Fr(0))}
    for j in range(3):
        if e[j]: term=pmul(term,power(Xl[j],e[j],('X',j)))
        if e[3+j]: term=pmul(term,power(Yl[j],e[3+j],('Y',j)))
    for ee,v in term.items(): Ft[ee]=gadd(Ft.get(ee,ZERO),v)
Ft={e:v for e,v in Ft.items() if v!=ZERO}
supp=set(Ft)
out={'n_monomials_uv':len(Ft)}
# block-monomial permutations: sigma in S3 on blocks, swap bits per block
perms=[]
for sigma in itertools.permutations(range(3)):
    for bits in itertools.product([0,1],repeat=3):
        pi=[0]*6
        for j in range(3):
            # coordinate u_j -> u_{sigma j} or v_{sigma j}; v_j -> the other
            pi[j]=sigma[j] if bits[j]==0 else 3+sigma[j]
            pi[3+j]=3+sigma[j] if bits[j]==0 else sigma[j]
        perms.append((sigma,bits,tuple(pi)))
results=[]
for sigma,bits,pi in perms:
    # (g z)_{pi(k)} = d_k z_k ; coefficient of z^{m'} in F(gz) is a_m prod_k d_k^{m'_k} with m = m' o pi^{-1}... derive: F(gz)=sum_m a_m prod_j (gz)_j^{m_j}, (gz)_j = d_{pi^{-1} j} z_{pi^{-1} j}
    # => exponent of z_k is m_{pi(k)} =: m'_k, coefficient a_m prod_k d_k^{m'_k}. Condition: a_m d^{m'} = lam a_{m'} for all m', m_j = m'_{pi^{-1}(j)}.
    pinv=[pi.index(k) for k in range(6)]
    def m_of(mp): return tuple(mp[pinv[j]] for j in range(6))
    if any(m_of(mp) not in supp for mp in supp): results.append({'sigma':sigma,'bits':bits,'ok':False,'why':'support not stable'}); continue
    rho={mp:gmul(Ft[mp],ginv(Ft[m_of(mp)])) for mp in supp}   # d^{m'} = lam * rho(m') ... careful: a_m d^{m'} = lam a_{m'} => d^{m'} = lam a_{m'}/a_m
    # determine d_k/d_0 from pairs (m', m'+e_0-e_k)
    ratios=[None]*6; ratios[0]=(Fr(1),Fr(0))
    for k in range(1,6):
        for mp in supp:
            if mp[k]>=1:
                m2=list(mp); m2[k]-=1; m2[0]+=1; m2=tuple(m2)
                if m2 in supp:
                    ratios[k]=gmul(rho[mp],ginv(rho[m2])); break   # d^{m'}/d^{m2} = d_k/d_0 = rho(m')/rho(m2)
    if any(r is None for r in ratios): results.append({'sigma':sigma,'bits':bits,'ok':False,'why':'undetermined'}); continue
    # lambda from m'=(6,0,..) if present else any: lam = d^{m'}/rho(m') with d_0=1
    ok=True; lam=None
    for mp in supp:
        dm=(Fr(1),Fr(0))
        for k in range(6):
            if mp[k]: dm=gmul(dm,gpow(ratios[k],mp[k]))
        l=gmul(dm,ginv(rho[mp]))
        if lam is None: lam=l
        elif l!=lam: ok=False; break
    results.append({'sigma':sigma,'bits':bits,'ok':ok,'d_ratios':[str(r) for r in ratios] if ok else None,'lambda':str(lam) if ok else None})
out['automorphisms']=[r for r in results if r['ok']]; out['all']=results
out['n_candidates']=len(results); out['n_ok']=len(out['automorphisms'])
print(json.dumps({'n_monomials_uv':out['n_monomials_uv'],'n_ok':out['n_ok'],'ok':[(r['sigma'],r['bits'],r['d_ratios'],r['lambda']) for r in out['automorphisms']]}))
json.dump(out,open('/tmp/hodge/results/eight-plane-configuration/wave2-construct/aut_blockmonomial.json','w'))
