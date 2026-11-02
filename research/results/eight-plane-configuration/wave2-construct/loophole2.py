"""Delta|_{Delta_tau} = F60(X(p)) F06(Y(p)) S(p)^2 ? ; tangent cone of Delta at a point of {S=0}."""
import sys, json, random, itertools
sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration/wave2-construct')
from helpers import *
from sympy import gcd, factor_list
p=1000000009; C=Ctx(p); i=C.i; random.seed(11); res={}
def line_pts(a0,a1,l): return [(a0[k]+l*a1[k])%p for k in range(3)]
for attempt in range(12):
    a0=[random.randrange(p) for _ in range(3)]; a1=[random.randrange(p) for _ in range(3)]
    xs=list(range(1,40)); ys=[C.Delta_at(line_pts(a0,a1,l),line_pts(a0,a1,l)) for l in xs]
    D=C.topoly(lagrange_mod(xs,ys,p))
    f60=C.restrict_curve(C.F60,lambda l:line_pts(a0,a1,l),(0,1,2))
    f06=C.restrict_curve(C.F06,lambda l:[(i*c)%p for c in line_pts(a0,a1,l)],(3,4,5))
    q1,r1=D.div(f60); q2,r2=q1.div(f06)
    res['deg_D']=D.degree(); res['divisible_by_F60']=r1.is_zero; res['divisible_by_F06']=r2.is_zero; res['cofactor_degree']=q2.degree()
    S=gcd(q2,q2.diff(lam)); res['S_degree']=S.degree()
    res['cofactor_is_square_of_squarefree']=bool(2*S.degree()==q2.degree() and gcd(S,S.diff(lam)).degree()==0)
    lin=[f for f,m in factor_list(S.as_expr(),modulus=p)[1] if Poly(f,lam).degree()==1]
    if lin:
        f=Poly(lin[0],lam,modulus=p); cs=[int(c) for c in f.all_coeffs()]; r=(-cs[1]*pow(cs[0],-1,p))%p
        res['point_on_S']=line_pts(a0,a1,r); break
if 'point_on_S' in res:
    a=res['point_on_S']; qc=C.qcoeffs(a,a)
    x=symbols('x'); pol=Poly(sum(int(qc[k])*x**(4-k) for k in range(5)),x,modulus=p)
    res['q_factorisation_at_point']=str(factor_list(pol.as_expr(),modulus=p)[1])
    pairs=list(itertools.combinations_with_replacement(range(6),2)); rows=[]; rhs=[]; mults=[]
    for _ in range(40):
        dr=[random.randrange(p) for _ in range(6)]
        xs=list(range(1,40)); ys=[C.Delta_at([(a[k]+l*dr[k])%p for k in range(3)],[(a[k]+l*dr[3+k])%p for k in range(3)]) for l in xs]
        cf=lagrange_mod(xs,ys,p); v=0
        while v<len(cf) and cf[v]==0: v+=1
        mults.append(v); rows.append([dr[u]*dr[w]%p for (u,w) in pairs]); rhs.append(cf[2] if len(cf)>2 else 0)
    res['multiplicities_random_dirs']=sorted(set(mults))
    M=[r[:]+[b] for r,b in zip(rows,rhs)]; n=len(pairs); rk=0; piv=[]
    for c in range(n):
        pr=next((r for r in range(rk,len(M)) if M[r][c]%p),None)
        if pr is None: continue
        M[rk],M[pr]=M[pr],M[rk]; inv=pow(M[rk][c],-1,p); M[rk]=[v*inv%p for v in M[rk]]
        for r in range(len(M)):
            if r!=rk and M[r][c]%p:
                f=M[r][c]; M[r]=[(v-f*w)%p for v,w in zip(M[r],M[rk])]
        piv.append(c); rk+=1
    res['tangent_cone_fit_consistent']=all(all(v%p==0 for v in row) for row in M[rk:])
    coef=[0]*n
    for r,c in enumerate(piv): coef[c]=M[r][n]
    G=[[0]*6 for _ in range(6)]
    for k,(u,w) in enumerate(pairs):
        if u==w: G[u][u]=coef[k]%p
        else: G[u][w]=G[w][u]=coef[k]*pow(2,-1,p)%p
    Am=[row[:] for row in G]; rank=0
    for c in range(6):
        pr=next((r for r in range(rank,6) if Am[r][c]%p),None)
        if pr is None: continue
        Am[rank],Am[pr]=Am[pr],Am[rank]; inv=pow(Am[rank][c],-1,p); Am[rank]=[v*inv%p for v in Am[rank]]
        for r in range(6):
            if r!=rank and Am[r][c]%p:
                f=Am[r][c]; Am[r]=[(v-f*w)%p for v,w in zip(Am[r],Am[rank])]
        rank+=1
    res['tangent_cone_rank']=rank
print(json.dumps(res)); json.dump(res,open('loophole2.json','w'))
