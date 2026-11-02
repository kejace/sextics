# Exact finite-field seed data + high-precision 5-adic Hensel lift for the profile-A seed on W_{23,f}.
import sys, json, time
import sympy as sp
sys.path.insert(0,'/tmp/hodge/results/polarized-k3-remainder')
LOG='/tmp/hodge/results/polarized-k3-remainder/wave4/log.jsonl'
def log(**kw):
    kw['ts']=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()); open(LOG,'a').write(json.dumps(kw, default=str)+'\n')
t=sp.symbols('t'); c=sp.symbols('c0:13')
k0,k1,k2,k3,k4 = 2,3,1,1,1
x=4*t; a=c[0]+c[1]*t+c[2]*t**2; y=k0+c[3]*t+c[4]*t**2
S=k1+c[5]*t+c[6]*t**2+c[7]*t**3+c[8]*t**4+c[9]*t**5+k3*t**6
T=k2+c[10]*t+c[11]*t**2+c[12]*t**3+k4*t**4
f7 = 7*x**6 + 5*x**5*a + 5*x**2*a**3*y + 15*x*y**5 + 7*a**6 + 10*a**5*y + 7*y**6
H = sp.expand(f7 + 7*S**2 + 7*T**3); P=sp.Poly(H,t)
eqs=[sp.expand(P.coeff_monomial(t**i)) for i in range(13)]
Jsym=sp.Matrix([[sp.diff(e,ci) for ci in c] for e in eqs])
seed=[1,1,2, 0,3, 4,0,0,0,4, 0,2,0]
sub=dict(zip(c,seed))
E0=[int(e.subs(sub)) for e in eqs]
assert all(v%5==0 for v in E0), E0
J0=Jsym.subs(sub); J0m=J0.applyfunc(lambda v: int(v)%5)
det5=int(J0.det())%5
print('seed residuals mod 5 all zero; det J mod 5 =', det5)
Jinv5=J0.inv_mod(5)
assert (J0*Jinv5).applyfunc(lambda v: int(v)%5)==sp.eye(13)
# freeness matrices mod 5: columns l*H_x (deg<=1), l*H_a, l*H_y, sigma*H_S (deg<=5), theta*H_T (deg<=3) -> coefficients t^0..t^11
grad={'x':sp.diff(f7, t) if False else None}
X_,A_,Y_,S_,T_=sp.symbols('X A Y S T')
Hw = 7*X_**6 + 5*X_**5*A_ + 5*X_**2*A_**3*Y_ + 15*X_*Y_**5 + 7*A_**6 + 10*A_**5*Y_ + 7*Y_**6 + 7*S_**2 + 7*T_**3
phi={X_:x,A_:a,Y_:y,S_:S,T_:T}
g={v: sp.Poly(sp.expand(sp.diff(Hw,v).subs(phi).subs(sub)), t) for v in (X_,A_,Y_,S_,T_)}
def coeffs(poly,n):
    return [int(poly.coeff_monomial(t**i)) for i in range(n)]
gc={str(v):coeffs(g[v],13) for v in g}
def freeness(with_x=True):
    cols=[]
    for v,w in ((X_,1),(A_,1),(Y_,1),(S_,3),(T_,2)):
        if v==X_ and not with_x: continue
        gv=coeffs(g[v],13)
        for k in range(2*w):
            cols.append([ (gv[i-k] if 0<=i-k<13 else 0) for i in range(12)])
    return sp.Matrix(cols).T
M16=freeness(True); M14=freeness(False)
r16=M16.applyfunc(lambda v:v%5).rank(iszerofunc=lambda v: v%5==0)
# ranks mod 5 properly: use nullspace over GF(5) via sympy's Matrix.rref with modulus not available; do manual gaussian elim
def rank_mod(M,p):
    M=[[int(v)%p for v in row] for row in M.tolist()]; r=0; rows=len(M); cols=len(M[0])
    for cc in range(cols):
        piv=next((i for i in range(r,rows) if M[i][cc]%p),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][cc],p-2,p); M[r]=[(v*inv)%p for v in M[r]]
        for i in range(rows):
            if i!=r and M[i][cc]: M[i]=[(vi-M[i][cc]*vr)%p for vi,vr in zip(M[i],M[r])]
        r+=1
    return r
def right_inverse_mod(M,p):
    # find R with M R = I (M is m x n, rank m) by solving per column over GF(p)
    m,n=M.shape; Mi=[[int(v)%p for v in row] for row in M.tolist()]
    import itertools
    # gaussian elimination on augmented [M | I]
    aug=[Mi[i]+[1 if j==i else 0 for j in range(m)] for i in range(m)]
    r=0; pivcols=[]
    for cc in range(n):
        piv=next((i for i in range(r,m) if aug[i][cc]%p),None)
        if piv is None: continue
        aug[r],aug[piv]=aug[piv],aug[r]; inv=pow(aug[r][cc],p-2,p); aug[r]=[(v*inv)%p for v in aug[r]]
        for i in range(m):
            if i!=r and aug[i][cc]: aug[i]=[(vi-aug[i][cc]*vr)%p for vi,vr in zip(aug[i],aug[r])]
        pivcols.append(cc); r+=1
        if r==m: break
    assert r==m
    R=[[0]*m for _ in range(n)]
    for i,cc in enumerate(pivcols):
        for j in range(m): R[cc][j]=aug[i][n+j]
    Rm=sp.Matrix(R); assert (M*Rm).applyfunc(lambda v:int(v)%p)==sp.eye(m)
    return Rm
rk16=rank_mod(M16,5); rk14=rank_mod(M14,5)
R16=right_inverse_mod(M16,5); R14=right_inverse_mod(M14,5)
print('freeness ranks mod 5:',rk16,rk14)
# endpoint gradient at t=0 and dominance: d phi / d(k0,k1,k2,t) at t=0 -> point (a(0),y(0),S(0),T(0)) with x=0
grad0={k:v[0]%5 for k,v in gc.items()}
print('endpoint gradient (H_x,H_a,H_y,H_S,H_T) at t=0 mod 5:', grad0)
# Hensel lift to 5^64
prec=64; mod=5**prec
cur=[int(v) for v in seed]
import sympy
for k in range(1,8):
    m=5**(2**k) if 2**k<=prec else mod
    E=[int(e.subs(dict(zip(c,cur))))%m for e in eqs]
    Jc=Jsym.subs(dict(zip(c,cur))).applyfunc(lambda v:int(v)%m)
    Ji=Jc.inv_mod(m)
    delta=(Ji*sp.Matrix(E)).applyfunc(lambda v:int(v)%m)
    cur=[(cur[i]-int(delta[i]))%m for i in range(13)]
    res=max(int(e.subs(dict(zip(c,cur))))%m for e in eqs)
    print('precision 5^%d: max residual %d'%(2**k if 2**k<=prec else prec, res), flush=True)
    if m==mod: break
res64=[int(e.subs(dict(zip(c,cur))))%mod for e in eqs]
assert max(res64)==0
json.dump({'params':{'k0':k0,'k1':k1,'k2':k2,'S6':k3,'T4':k4},'seed_mod5':seed,'detJ_mod5':det5,
           'Jacobian_mod5':[[int(v) for v in row] for row in J0m.tolist()],'Jacobian_inverse_mod5':[[int(v) for v in row] for row in Jinv5.tolist()],
           'freeness16_mod5':[[int(v)%5 for v in row] for row in M16.tolist()],'freeness16_right_inverse_mod5':[[int(v) for v in row] for row in R16.tolist()],
           'freeness14_nox_mod5':[[int(v)%5 for v in row] for row in M14.tolist()],'freeness14_right_inverse_mod5':[[int(v) for v in row] for row in R14.tolist()],
           'endpoint_gradient_mod5':grad0,'gradient_coeffs_on_curve':gc,
           'lift_mod_5^64':[str(v) for v in cur],'lift_residual_max':max(res64)},
          open('/tmp/hodge/results/polarized-k3-remainder/wave4/w23_seed_exact_data.json','w'), indent=0)
log(what='w23-seed-exact-data-and-hensel-lift', prime=5, params=[k0,k1,k2,k3,k4], result={'detJ_mod5':det5,'freeness_rank16':rk16,'freeness_rank14':rk14,'right_inverses_verified':True,'lift_precision':'5^64','lift_residual_max':max(res64)}, file='wave4/w23_seed_exact_data.json')
print('done')
