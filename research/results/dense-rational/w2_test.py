import sys, time; sys.path.insert(0,'/tmp/hodge/results/dense-rational')
import numpy as np
from w2lib import *
d=int(sys.argv[1]); pD={'even':0,'odd':1,'mixed':None}[sys.argv[2]]
T=time.time()
def lap(msg):
    global T; print(f'{msg}: {time.time()-T:.1f}s', flush=True); T=time.time()
F,A,q = load(); lap('load')
sp = spaces(d,pD); lap('spaces')
fam = jacobian_family(F,A,d,pD,sp); lap('family')
S = NormE2(F,A,q,fam,sp,verbose=True); lap('NormE2 build')
rng=np.random.default_rng(1)
w = rng.standard_normal(S.nw)
st=S.evaluate(w); lap(f'evaluate f={st["f"]:.3e} rel={S.relres(st):.3e} rank={st["rank"]}')
st=S.evaluate(w, need_jac=True); lap('evaluate+jac')
J=S.reduced_jac(st); lap('reduced_jac')
if len(sys.argv)>3 and sys.argv[3]=='check':
    Fi=[pderiv(F,i) for i in range(NV)]
    Dp={}; N0=[{} for _ in range(NV)]
    for k,(i,m) in enumerate(fam.names):
        Dp=add(Dp, scale(shift(Fi[i],m,1), float(w[k]))); N0[i]=add(N0[i], shift(A,m,float(w[k])))
    x=st['x']; N=[dict(N0[i]) for i in range(NV)]
    for c,(i,j,m) in enumerate(S.Kcols):
        if x[c]!=0:
            N[i]=add(N[i], shift(Fi[j],m,float(x[c]))); N[j]=add(N[j], shift(Fi[i],m,-float(x[c])))
    e1={}
    for i in range(NV): e1=add(e1, mul(N[i],Fi[i]))
    e1=sub(e1, mul(A,Dp)); print('E1 residual max', max([abs(v) for v in e1.values()]+[0]), flush=True)
    divN={}; ND={}
    for i in range(NV):
        divN=add(divN,pderiv(N[i],i)); ND=add(ND, mul(N[i],pderiv(Dp,i)))
    BD=sub(mul(Dp,divN),ND); vec=S.sp['rows'].vec(BD)
    # assembled B_D(N): r + bp... recompute b + MK x directly
    Dc=S.Dmat.T@w; N0v=[S.N0mats[i].T@w for i in range(NV)]
    dD=[apply_deriv(S.dD_maps[i],Dc,S.ndDi[i]) for i in range(NV)]; divN0=S.divN0mat.T@w
    S._fill(S.SD,S.FLAT_SD,Dc,S.nT)
    for i in range(NV): S._fill(S.Si[i],S.FLAT_Si[i],dD[i],S.nTi[i])
    b=S.SD@divN0; MK=S.SD@S.V
    for i in range(NV): b-=S.Si[i]@N0v[i]; MK-=S.Si[i]@S.W[i]
    mine=b+MK@x
    print('B_D(N) direct vs assembled: max abs diff', np.max(np.abs(vec-mine)), 'scale', np.max(np.abs(vec)), flush=True)
    lap('polynomial cross-check')
print('fd check (grad.v, fd, raw):', S.fd_check(seed=3), flush=True); lap('fd_check')
res=S.lm(rng.standard_normal(S.dimU), maxiter=30, verbose=True); lap(f'LM rel={res["relres"]:.3e} iters={res["iters"]}')
