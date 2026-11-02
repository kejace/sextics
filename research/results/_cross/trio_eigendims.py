# cheap cross-check of the cyclic-character decomposition of H^{3,1} (R_6), H^{2,2} (R_12), H^{1,3} (R_18) for the three lanes
import re, itertools, sympy as sp
from flint import nmod_mat
p = 1000003
X = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
def readF(lane):
    txt = open(f'/home/kejace/devel/sextics/data/{lane}/instance.rec').read()
    P = sp.Poly(sp.expand(7*sp.sympify(re.search(r'^F: (.*(?:\n  .*)*)', txt, re.M).group(1).replace('\n',' '))), *X)
    return {tuple(e): int(c) for e, c in zip(P.monoms(), P.coeffs())}
def mons(k):
    out=[]
    for c in itertools.combinations_with_replacement(range(6), k):
        e=[0]*6
        for i in c: e[i]+=1
        out.append(tuple(e))
    return out
def pderiv(P, i):
    out={}
    for e,v in P.items():
        if e[i]:
            f=list(e); f[i]-=1; out[tuple(f)]=out.get(tuple(f),0)+v*e[i]
    return out
def sector_dims(F, charfn, k, chars):
    dF=[pderiv(F,i) for i in range(6)]
    res={}
    for ch in chars:
        S=[e for e in mons(k) if charfn(e)==ch]; idx={e:i for i,e in enumerate(S)}
        gens=[]
        for i in range(6):
            for m in mons(k-5):
                g={tuple(a+b for a,b in zip(e,m)): v for e,v in dF[i].items()}
                if g and charfn(next(iter(g)))==ch: gens.append(g)
        M=nmod_mat(len(S), len(gens), p)
        for j,g in enumerate(gens):
            for e,v in g.items(): M[idx[e],j]=v%p
        res[ch]=len(S)-M.rank()
    return res
# ctr: residue character of X^e Omega/F^k under X2 -> zeta X2 is zeta^(e2+1); report by e2 mod 3
F=readF('cyclic-triple-residual')
for k in (6,12):
    d=sector_dims(F, lambda e: e[2]%3, k, [0,1,2])
    print('ctr degree',k,'by e2 mod 3 (residue char zeta^(e2+1)):', d)
# rcr: residue character under Y2 -> -Y2 is (-1)^(e5+1); report by e5 mod 2
F=readF('reflection-coupled-residual')
for k in (6,12):
    d=sector_dims(F, lambda e: e[5]%2, k, [0,1])
    print('rcr degree',k,'by e5 mod 2 (residue char (-1)^(e5+1)):', d)
