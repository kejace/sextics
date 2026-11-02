#!/usr/bin/env python3
"""Aggregator re-derivation of quadric-contact-residual degree-7 parity-class-pair exclusions (excl items 3-4).
Classes s = (b0,b2,y) = (e0+e3, e2+e5, e3+e4+e5) mod 2 of the involutions sigma0, sigma2, tau_Y (h1-characters mix freely).
D in J_7 (E1 kernel 0 at d=7) of class support {s,t}:  D = w(F), w quadratic with w_i of class s+class(z_i) or t+class(z_i).
Relaxation: A*Hess_F(w,w) in J_14  <=>  Hess_F(w,w) in J_8 (+ C*D8* in class 0, D8* the unique non-J E1 direction at d=8).
Euler gauge: coefficient of z_k^2 in w_k dropped for every z_k of class s or t.   usage: qcr7.py s t   (e.g. 000 111)"""
import sys, json, time, itertools, subprocess, re
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from hodge.core import NV
from flint import nmod_mat
s_=tuple(int(c) for c in sys.argv[1]); t_=tuple(int(c) for c in sys.argv[2])
p=int(sys.argv[3]) if len(sys.argv)>3 else 998244353
tmo=int(sys.argv[4]) if len(sys.argv)>4 else 300
F,A,q,rec=core.load_instance('quadric-contact-residual')
Fp=core.reduce_mod(F,p)
Fi=[core.pderiv(Fp,i,p) for i in range(NV)]
Fij=[[core.pderiv(Fi[i],j,p) for j in range(NV)] for i in range(NV)]
def cls(e): return ((e[0]+e[3])%2,(e[2]+e[5])%2,(e[3]+e[4]+e[5])%2)
def add(a,b): return tuple((x+y)%2 for x,y in zip(a,b))
unit=lambda i: tuple(1 if k==i else 0 for k in range(NV))
zc=[cls(unit(i)) for i in range(NV)]
# D8* = (F_X0^2 - F_Y0^2)/(X0^2 - Y0^2): solve linear system for the degree-8 quotient
P=core.add(core.mul(Fi[0],Fi[0],p), {e:(-c)%p for e,c in core.mul(Fi[3],Fi[3],p).items()}, p)
m8=[e for e in itertools.product(range(9),repeat=NV) if sum(e)==8]
m10=[e for e in itertools.product(range(11),repeat=NV) if sum(e)==10]
i10={e:t for t,e in enumerate(m10)}
M=nmod_mat(len(m10),len(m8)+1,p)
for t,e in enumerate(m8):
    e1=list(e); e1[0]+=2; M[i10[tuple(e1)],t]=1
    e2=list(e); e2[3]+=2; M[i10[tuple(e2)],t]=(p-1)
for e,c in P.items(): M[i10[e],len(m8)]=c%p
R=M.rref()[0]
# read solution: each pivot row gives var = rhs
D8={}
rows=R.nrows(); cols=R.ncols()
for r in range(rows):
    pivot=None
    for c in range(cols-1):
        if R[r,c]!=0: pivot=c; break
    if pivot is None:
        assert R[r,cols-1]==0, 'inconsistent division'
        continue
    v=int(R[r,cols-1])
    if v: D8[m8[pivot]]=v
chk=core.add(core.mul(D8,{(2,0,0,0,0,0):1,(0,0,0,2,0,0):p-1},p), {e:(-c)%p for e,c in P.items()}, p)
assert not chk, 'D8* division check failed'
red=jacobian._Reducer(F,8,p,0)
nfD=red.normal_form(D8); assert nfD, 'D8* in J?!'
assert cls(next(iter(D8)))==(0,0,0)
piv=max(nfD, key=lambda e: abs(nfD[e])); inv=pow(nfD[piv],-1,p)
def quot(nf):
    if piv in nf:
        lam=nf[piv]*inv%p
        nf={e:(nf.get(e,0)-lam*nfD.get(e,0))%p for e in set(nf)|set(nfD)}
        nf={e:c for e,c in nf.items() if c}
    return nf
quad=[e for e in itertools.product(range(3),repeat=NV) if sum(e)==2]
vars_=[]
for c in (s_,t_):
    for i in range(NV):
        for m in quad:
            if cls(m)==add(c,zc[i]): vars_.append((i,m,c))
gauge={(k, tuple(2 if t==k else 0 for t in range(NV))) for k in range(NV) if zc[k] in (s_,t_)}
keep=[v for v in vars_ if (v[0],v[1]) not in gauge]
nv=len(keep); names=[f'c{t}' for t in range(nv)]
redo=jacobian._Reducer(F,8,p,1)   # odd-parity reducer for cross terms of odd Y-parity
def nf_any(P):
    ev={e:c for e,c in P.items() if (e[3]+e[4]+e[5])%2==0}; od={e:c for e,c in P.items() if (e[3]+e[4]+e[5])%2==1}
    out={}
    if ev: out.update(quot(red.normal_form(ev)))
    if od: out.update(redo.normal_form(od))
    return out
pairs={}
for a in range(nv):
    for b in range(a,nv):
        (i,ma,_),(j,mb,_)=keep[a],keep[b]
        Pab=core.mul(core.mul({ma:1},{mb:1},p),Fij[i][j],p)
        if a!=b: Pab=core.add(Pab,Pab,p)
        pairs[(a,b)]=nf_any(Pab)
coords=sorted({e for nf in pairs.values() for e in nf}); cidx={e:t for t,e in enumerate(coords)}
keys=list(pairs)
Mq=nmod_mat(len(coords),len(pairs),p)
for t,(a,b) in enumerate(keys):
    for e,c in pairs[(a,b)].items(): Mq[cidx[e],t]=c
rk=Mq.rank()
polys=[]
for e in coords:
    terms=[(f'{pairs[(a,b)][e]}*{names[a]}^2' if a==b else f'{pairs[(a,b)][e]}*{names[a]}*{names[b]}') for (a,b) in keys if pairs[(a,b)].get(e,0)]
    if terms: polys.append('+'.join(terms))
polys=sorted(set(polys))
tag=sys.argv[1]+'_'+sys.argv[2]
ms=f'qcr7_{tag}.ms'; open(ms,'w').write(','.join(names)+'\n'+str(p)+'\n'+',\n'.join(polys)+'\n')
t0=time.time(); status='ok'
try:
    r=subprocess.run(['nice','-n','10','msolve','-g','1','-t','2','-f',ms,'-o',f'qcr7_{tag}.g1.out'],capture_output=True,text=True,timeout=tmo)
    if r.returncode!=0: status=f'msolve rc {r.returncode}'
except subprocess.TimeoutExpired: status='msolve timeout'
res={'pair':tag,'dim_V':len(vars_),'n_gauge':len(gauge),'n_vars':nv,'independent_quadrics':rk,'msolve':status,'seconds':round(time.time()-t0,1),'D8_terms':len(D8)}
if status=='ok':
    out=open(f'qcr7_{tag}.g1.out').read(); body=out[out.find('['):]
    monos=[m.strip().replace(' ','') for m in body.strip().strip('[]:').split(',') if m.strip()]
    pure={}; gens=[]
    for m in monos:
        e={}
        for f in m.split('*'):
            if f=='1': continue
            mm=re.fullmatch(r'([A-Za-z0-9_]+)(?:\^(\d+))?',f)
            if mm: e[mm.group(1)]=e.get(mm.group(1),0)+int(mm.group(2) or 1)
        if e: gens.append(e)
        if len(e)==1: n=list(e)[0]; pure[n]=min(pure.get(n,99),e[n])
    missing=[n for n in names if n not in pure]
    res.update({'n_leading':len(gens),'missing_pure_power_for':missing,'m_primary':not missing,'max_pure_power':max(pure.values()) if pure else None})
print(json.dumps(res), flush=True)
open('qcr7_results.jsonl','a').write(json.dumps(res)+'\n')
