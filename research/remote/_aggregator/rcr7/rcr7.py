#!/usr/bin/env python3
"""Aggregator re-derivation of the reflection-coupled-residual degree-7 Hessian relaxation per Klein character.
D in J_7 of character chi  =>  D = w(F), w quadratic vector field with chi(w_i) = chi - chi(z_i);  the relaxation
A*Hess_F(w,w) in J_14  <=>  Hess_F(w,w) in J_8 + C*K8  (the E1 kernel of [A] on R_8^(0,0) is exactly C*K8).
Euler gauge: w -> w + h*E (h linear of character chi) changes Hess_F(w,w) by elements of J, so the coefficient of
Y2^2 / Y0^2 / X0^2 in w_5 / w_3 / w_0 (one per Euler direction) is set to zero.
usage: python3 rcr7.py <chi as two bits, e.g. 10> [p] [msolve timeout s]"""
import sys, json, time, itertools, subprocess, re
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from hodge.core import NV
from flint import nmod_mat
chi_s=sys.argv[1]; target=(int(chi_s[0]),int(chi_s[1]))
p=int(sys.argv[2]) if len(sys.argv)>2 else 998244353
tmo=int(sys.argv[3]) if len(sys.argv)>3 else 600
F,A,q,rec=core.load_instance('reflection-coupled-residual')
Fp=core.reduce_mod(F,p)
Fi=[core.pderiv(Fp,i,p) for i in range(NV)]
Fij=[[core.pderiv(Fi[i],j,p) for j in range(NV)] for i in range(NV)]
txt=open('/tmp/hodge/results/reflection-coupled-residual/k8_expr.txt').read().strip()
txt=txt.split('=',1)[1].strip() if txt.startswith('K8') else txt
K8=core.reduce_mod(core.parse(txt),p)
def chi(e): return (e[5]%2, (e[3]+e[4])%2)
def chi_add(a,b): return ((a[0]+b[0])%2,(a[1]+b[1])%2)
unit=lambda i: tuple(1 if k==i else 0 for k in range(NV))
zchi=[chi(unit(i)) for i in range(NV)]
quad=[e for e in itertools.product(range(3),repeat=NV) if sum(e)==2]
vars_=[(i,m) for i in range(NV) for m in quad if chi(m)==chi_add(target,zchi[i])]
# Euler directions: h linear of character chi: h = z_k with chi(z_k) = chi ; w = h*E has w_i = z_k z_i ; gauge slice: drop (i=k, m=z_k^2)
euler=[k for k in range(NV) if zchi[k]==target]
gauge=[(k, tuple(2 if t==k else 0 for t in range(NV))) for k in euler]
keep=[v for v in vars_ if v not in gauge]
dimV=len(vars_); nv=len(keep)
red=jacobian._Reducer(F,8,p,0)
std=red.standard
nfK=red.normal_form(K8); assert nfK
piv=max(nfK, key=lambda e: abs(nfK[e])); inv=pow(nfK[piv],-1,p)
def quot(nf):
    if piv in nf:
        lam=nf[piv]*inv%p
        nf={e:(nf.get(e,0)-lam*nfK.get(e,0))%p for e in set(nf)|set(nfK)}
        nf={e:c for e,c in nf.items() if c}
    return nf
names=[f'c{t}' for t in range(nv)]
pairs={}
for a in range(nv):
    for b in range(a,nv):
        (i,ma),(j,mb)=keep[a],keep[b]
        P=core.mul(core.mul({ma:1},{mb:1},p),Fij[i][j],p)
        if a!=b: P=core.add(P,P,p)
        pairs[(a,b)]=quot(red.normal_form(P))
coords=sorted({e for nf in pairs.values() for e in nf}); cidx={e:t for t,e in enumerate(coords)}
keys=list(pairs)
M=nmod_mat(len(coords),len(pairs),p)
for t,(a,b) in enumerate(keys):
    for e,c in pairs[(a,b)].items(): M[cidx[e],t]=c
rk=M.rank()
polys=[]
for e in coords:
    terms=[(f'{pairs[(a,b)][e]}*{names[a]}^2' if a==b else f'{pairs[(a,b)][e]}*{names[a]}*{names[b]}') for (a,b) in keys if pairs[(a,b)].get(e,0)]
    if terms: polys.append('+'.join(terms))
polys=sorted(set(polys))
ms=f'rcr7_{chi_s}.ms'; open(ms,'w').write(','.join(names)+'\n'+str(p)+'\n'+',\n'.join(polys)+'\n')
t0=time.time(); status='ok'
try:
    r=subprocess.run(['nice','-n','10','msolve','-g','1','-t','2','-f',ms,'-o',f'rcr7_{chi_s}.g1.out'],capture_output=True,text=True,timeout=tmo)
    if r.returncode!=0: status=f'msolve rc {r.returncode}'
except subprocess.TimeoutExpired:
    status='msolve timeout'
el=round(time.time()-t0,1)
res={'chi':chi_s,'dim_V':dimV,'dim_T':len(gauge),'n_vars':nv,'independent_quadrics':rk,'msolve':status,'seconds':el}
if status=='ok':
    out=open(f'rcr7_{chi_s}.g1.out').read(); body=out[out.find('['):]
    monos=[m.strip().replace(' ','') for m in body.strip().strip('[]:').split(',') if m.strip()]
    gens=[]; pure={}
    for m in monos:
        e={}
        for f in m.split('*'):
            if f=='1': continue
            mm=re.fullmatch(r'([A-Za-z0-9_]+)(?:\^(\d+))?',f)
            if mm: e[mm.group(1)]=e.get(mm.group(1),0)+int(mm.group(2) or 1)
        if e: gens.append(e)
        if len(e)==1: n=list(e)[0]; pure[n]=min(pure.get(n,99),e[n])
    missing=[n for n in names if n not in pure]
    G=[tuple(g.get(n,0) for n in names) for g in gens]
    def divisible(m): return any(all(m[i]>=g[i] for i in range(nv)) for g in G)
    from collections import deque
    start=tuple([0]*nv); seen={start}; dq=deque([start]); cap=20000
    while dq and len(seen)<=cap:
        m=dq.popleft()
        for i in range(nv):
            m2=m[:i]+(m[i]+1,)+m[i+1:]
            if m2 in seen or divisible(m2): continue
            seen.add(m2); dq.append(m2)
    res.update({'n_leading':len(gens),'missing_pure_power_for':missing,'m_primary':not missing,'colength':len(seen) if len(seen)<=cap else None,'max_pure_power':max(pure.values()) if pure else None})
print(json.dumps(res), flush=True)
json.dump(res, open(f'rcr7_{chi_s}.json','w'))
