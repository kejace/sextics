# Independent re-verification (sympy) of (i) the chart-reduction divergence criterion for V-degree b (b=2,4) and
# (ii) the (3,2)-weight scalar conversion identity on W_{23,f} (wave-2 B1, B4), with generic symbols.
import sympy as sp, json, time
LOG='/tmp/hodge/results/polarized-k3-remainder/wave4/log.jsonl'
def log(**kw):
    kw['ts']=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()); open(LOG,'a').write(json.dumps(kw, default=str)+'\n')
# ---------- (i) chart reduction ----------
x,a,r,z,lam,t=sp.symbols('x a r z lambda t')
p=sp.Function('p')(x,a); q=sp.Function('q')(r,z)
out={}
for b in (2,4):
    s=sp.Function('s')(x,a,r,z)
    Wx,Wa,Wr,Wz=[sp.Function(n)(x,a,r,z) for n in ('Wx','Wa','Wr','Wz')]
    G=p+lam**6*q
    Wp=Wx*sp.diff(p,x)+Wa*sp.diff(p,a); Wq=Wr*sp.diff(q,r)+Wz*sp.diff(q,z)
    wl=(2*lam**b*s-lam**b*Wp-lam**(b+6)*Wq)/(6*lam**5*q)
    w=[lam**b*Wx,lam**b*Wa,lam**b*Wr,lam**b*Wz,wl]
    # w(G) = 2 lam^b s ?
    wG=sp.simplify(w[0]*sp.diff(G,x)+w[1]*sp.diff(G,a)+w[2]*sp.diff(G,r)+w[3]*sp.diff(G,z)+w[4]*sp.diff(G,lam)-2*lam**b*s)
    # divergence with chart term 2 w_lam/lam
    div=sp.diff(w[0],x)+sp.diff(w[1],a)+sp.diff(w[2],r)+sp.diff(w[3],z)+sp.diff(w[4],lam)+2*w[4]/lam
    divW=sp.diff(Wx,x)+sp.diff(Wa,a)+sp.diff(Wr,r)+sp.diff(Wz,z)
    crit=divW+sp.Rational(b-3,6)*Wp/p-sp.Rational(b+3,6)*Wq/q-sp.Rational(b-3,3)*s/p
    # claim: div = lam^b * [ crit + (something)*G/(p q lam^?) ] ; check div - lam^b*crit is divisible by G: evaluate on G=0 (lam^6 -> -p/q)
    diff=sp.simplify((div-lam**b*crit))
    # substitute lam^6 = -p/q by rewriting in terms of L = lam^6 : check numerator vanishes when lam^6 = -p/q
    num=sp.numer(sp.together(diff)); num=sp.expand(num)
    # replace lam^k with lam^(k mod 6) * (-p/q)^(k//6)
    numsub=sp.expand(num.subs(lam, lam).replace(lambda e: e.is_Pow and e.base==lam and e.exp.is_Integer, lambda e: lam**(e.exp%6)*(-p/q)**(e.exp//6)))
    numsub=sp.simplify(numsub)
    out[b]={'w(G)-2lam^b s':str(wG),'(div - lam^b crit) vanishes on G=0':str(numsub)}
    print('b=%d:'%b,out[b])
log(what='chart-reduction-reverified', result=out)
# ---------- (ii) (3,2) scalar conversion on W_{23,f} ----------
# cover: u^6 = p(x,a), S = u^3 alpha, T = u^2 beta, E_23: alpha^2 + beta^3 = -1, delta beta = alpha, delta alpha = -(3/2) beta^2
# Companion: H(x,a,S,T) = p + S^2 + T^3 ; vector field V=(Vx,Va,VS,VT) with V(H) = B, div V = 0 (functions of x,a,S,T).
X,A,SS,TT,u,al,be=sp.symbols('x a S T u alpha beta')
pp=sp.Function('p')(X,A); B=sp.Function('B')(X,A,SS,TT)
Vx,Va,VS,VT=[sp.Function(n)(X,A,SS,TT) for n in ('Vx','Va','VS','VT')]
Hc=pp+SS**2+TT**3
# derivations on the cover field C(x,a,alpha,beta)[u]/(u^6-p): D_x = d/dx + (p_x/(6u^5)) d/du acting on expressions in x,a,u,alpha,beta
def Dx(e): return sp.diff(e,X)+sp.diff(pp,X)/(6*u**5)*sp.diff(e,u)
def Da(e): return sp.diff(e,A)+sp.diff(pp,A)/(6*u**5)*sp.diff(e,u)
def dl(e): return sp.diff(e,be)*al+sp.diff(e,al)*(-sp.Rational(3,2)*be**2)   # delta (invariant derivation of E_23)
sub={SS:u**3*al, TT:u**2*be}
def ev(e): return e.subs(sub)
f_=6*al*ev(VT)/u**3-4*be*ev(VS)/u**4
Cx=-6*ev(Vx)/u; Ca=-6*ev(Va)/u
lhs=dl(f_)+Dx(Cx)+Da(Ca)
rhs=ev(B)/u**7-3*al*ev(sp.diff(B,SS))/u**4-2*be*ev(sp.diff(B,TT))/u**5
# use hypotheses: V(H)=B  => Vx p_x + Va p_a + 2 S VS + 3 T^2 VT = B ; differentiate w.r.t. S,T: gives B_S, B_T ; div V = 0.
# Implement by solving: express derivatives of B via the normal identity. Do it by substituting B := Vx p_x + Va p_a + 2 S VS + 3 T^2 VT (as function) and div condition Vx_x + Va_a + VS_S + VT_T = 0 -> replace VT_T.
Bdef=Vx*sp.diff(pp,X)+Va*sp.diff(pp,A)+2*SS*VS+3*TT**2*VT
rhs2=ev(Bdef)/u**7-3*al*ev(sp.diff(Bdef,SS))/u**4-2*be*ev(sp.diff(Bdef,TT))/u**5
expr=sp.simplify(sp.expand(lhs-rhs2))
# impose div V = 0: solve for Derivative(VT, T)
dVT=sp.Derivative(VT,TT)
sol=sp.solve(sp.Eq(sp.diff(Vx,X)+sp.diff(Va,A)+sp.diff(VS,SS)+dVT,0), dVT)[0]
# after substitution S->u^3 alpha the derivative objects appear as Subs; handle by doing the substitution AFTER differentiation:
expr_raw=sp.expand(lhs-rhs2)
expr_raw=expr_raw.doit()
# replace Subs(Derivative(VT,T), T, u^2 beta) etc. : easier approach: differentiate first symbolically with generic jets.
# Simplest robust check: use polynomial test functions for V and p, B := V(H), and check identity exactly when div V = 0.
import random
random.seed(1)
def randpoly(vars_, deg, k=3):
    return sum(random.randint(-5,5)*sp.Mul(*[random.choice(vars_) for _ in range(random.randint(0,deg))]) for _ in range(k))
ok_all=True; residuals=[]
for trial in range(4):
    pt=randpoly([X,A],6,5)
    # divergence-free V: V = curl-type: take generating functions
    g1=randpoly([X,A,SS,TT],3,4); g2=randpoly([X,A,SS,TT],3,4); g3=randpoly([X,A,SS,TT],3,4)
    Vxp= sp.diff(g1,A)-sp.diff(g2,SS); Vap=-sp.diff(g1,X)+sp.diff(g3,TT); VSp=sp.diff(g2,X); VTp=-sp.diff(g3,A)
    assert sp.expand(sp.diff(Vxp,X)+sp.diff(Vap,A)+sp.diff(VSp,SS)+sp.diff(VTp,TT))==0
    Bp=sp.expand(Vxp*sp.diff(pt,X)+Vap*sp.diff(pt,A)+2*SS*VSp+3*TT**2*VTp)
    subs2={SS:u**3*al,TT:u**2*be}
    def Dx2(e): return sp.diff(e,X)+sp.diff(pt,X)/(6*u**5)*sp.diff(e,u)
    def Da2(e): return sp.diff(e,A)+sp.diff(pt,A)/(6*u**5)*sp.diff(e,u)
    fv=6*al*VTp.subs(subs2)/u**3-4*be*VSp.subs(subs2)/u**4
    Cxv=-6*Vxp.subs(subs2)/u; Cav=-6*Vap.subs(subs2)/u
    L=dl(fv)+Dx2(Cxv)+Da2(Cav)
    R=Bp.subs(subs2)/u**7-3*al*sp.diff(Bp,SS).subs(subs2)/u**4-2*be*sp.diff(Bp,TT).subs(subs2)/u**5
    d=sp.expand(sp.numer(sp.together(L-R)))
    # reduce modulo u^6 - p and alpha^2 + beta^3 + 1
    d=sp.rem(sp.Poly(d,u), sp.Poly(u**6-pt,u)).as_expr()
    d=sp.rem(sp.Poly(sp.expand(d),al), sp.Poly(al**2+be**3+1,al)).as_expr()
    d=sp.expand(d); residuals.append(d==0); ok_all=ok_all and d==0
print('(3,2) scalar conversion identity on random divergence-free V, random p:', residuals)
log(what='scalar-conversion-(3,2)-reverified', result={'random_trials_identity_holds':residuals})
