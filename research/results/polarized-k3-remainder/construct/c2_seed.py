"""C2: finite-field incidence seed for a polynomial curve family on the second elliptic companion
   W_{23,f} = {7 f(z) + 7 s^2 + 7 t^3 = 0} in P(1,1,1,3,2)  (book Ch. 3A.5), affine cone coordinates (z0,z1,z2,s,t).
Pattern (crossed-ternary style): z0 = kappa*tau exactly; the other four coordinates are polynomials in tau whose
constant terms are base parameters (k0 for z2, k1 for s, k2 for t) except z1 whose constant c0 is an unknown.
Equations: all coefficients of tau^0..tau^{6e} in H(phi(tau)).  Unknown count is chosen = 6e+1 (square).
Checks: (i) H(phi)=0 exactly over F_p, (ii) incidence Jacobian invertible, (iii) polynomial normal psi with
grad H(phi).psi = 1 from the same matrix, (iv) 5x5 determinant J(tau) of [phi', dphi/dk0, dphi/dk1, dphi/dk2, psi]
nonzero (dominance of the family + transversal normal), (v) coordinates have no common root (cone vertex avoided).
Usage: python3 c2_seed.py p e pattern seed_count
"""
import sys, itertools, random, json, time
p = int(sys.argv[1]); e = int(sys.argv[2]); pattern = sys.argv[3]; nseeds = int(sys.argv[4]) if len(sys.argv)>4 else 1
# 7f integral: 7z0^6+5z0^5z1+5z0^2z1^3z2+15z0z2^5+7z1^6+10z1^5z2+7z2^6 ; H = 7f + 7s^2 + 7t^3
# polynomials in tau as lists mod p
def padd(a,b):
    n=max(len(a),len(b)); return [((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0))%p for i in range(n)]
def pmul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%p
    return r
def psc(c,a): return [(c*x)%p for x in a]
def ppow(a,n):
    r=[1]
    for _ in range(n): r=pmul(r,a)
    return r
def H(z0,z1,z2,s,t):
    terms=[psc(7,ppow(z0,6)), psc(5,pmul(ppow(z0,5),z1)), psc(5,pmul(pmul(ppow(z0,2),ppow(z1,3)),z2)),
           psc(15,pmul(z0,ppow(z2,5))), psc(7,ppow(z1,6)), psc(10,pmul(ppow(z1,5),z2)), psc(7,ppow(z2,6)),
           psc(7,ppow(s,2)), psc(7,ppow(t,3))]
    r=[0]
    for u in terms: r=padd(r,u)
    return r
def gradH(z0,z1,z2,s,t):
    g0=padd(padd(psc(42,ppow(z0,5)),psc(25,pmul(ppow(z0,4),z1))),padd(psc(10,pmul(pmul(z0,ppow(z1,3)),z2)),psc(15,ppow(z2,5))))
    g1=padd(padd(psc(5,ppow(z0,5)),psc(15,pmul(pmul(ppow(z0,2),ppow(z1,2)),z2))),padd(psc(42,ppow(z1,5)),psc(50,pmul(ppow(z1,4),z2))))
    g2=padd(padd(psc(5,pmul(ppow(z0,2),ppow(z1,3))),psc(75,pmul(z0,ppow(z2,4)))),padd(psc(10,ppow(z1,5)),psc(42,ppow(z2,5))))
    gs=psc(14,s); gt=psc(21,ppow(t,2))
    return [g0,g1,g2,gs,gt]
# degree patterns: dict coord -> max degree ; coords 0..4 = z0,z1,z2,s,t
patterns = {
 # e=1 candidates (square: 7 unknowns)
 "e1a": {1:1, 2:1, 3:3, 4:1},   # z1 c0+c1 tau ; z2 k0+c2 tau ; s k1+c3..c5 ; t k2+c6 tau  -> unknowns 2+1+3+1=7
 "e1b": {1:1, 2:1, 3:2, 4:2},   # s quadratic, t quadratic -> 2+1+2+2 = 7
 # e=2 candidates (square: 13 unknowns)
 "e2a": {1:2, 2:2, 3:5, 4:3},   # 3+2+5+3 = 13
 "e2b": {1:2, 2:2, 3:4, 4:4},   # 3+2+4+4 = 13
 "e2c": {1:2, 2:2, 3:6, 4:2},   # 3+2+6+2 = 13
}
deg = patterns[pattern]
kappa = 1
def build(c0, coef, k):
    """coef: dict (coord,degree)->value for degree>=1; k=(k0,k1,k2)."""
    z0=[0,kappa]; z1=[c0]+[coef[(1,j)] for j in range(1,deg[1]+1)]
    z2=[k[0]]+[coef[(2,j)] for j in range(1,deg[2]+1)]
    s=[k[1]]+[coef[(3,j)] for j in range(1,deg[3]+1)]
    t=[k[2]]+[coef[(4,j)] for j in range(1,deg[4]+1)]
    return z0,z1,z2,s,t
unk = [(0,0)] + [(c,j) for j in range(1,max(deg.values())+1) for c in (1,2,3,4) if j<=deg[c]]  # (0,0) denotes c0
nunk=len(unk); neq=6*e+1
print(f"p={p} e={e} pattern={pattern} deg={deg} unknowns={nunk} equations={neq}")
assert nunk==neq, "not square"
def eqs(c0,coef,k):
    h=H(*build(c0,coef,k)); h=h+[0]*(neq-len(h)); return h[:neq], h[neq:]
def coef_from_vec(v):
    c0=v[0]; coef={unk[i]:v[i] for i in range(1,nunk)}; return c0,coef
def residual(v,k):
    c0,coef=coef_from_vec(v); h,tail=eqs(c0,coef,k)
    assert all(x==0 for x in tail)
    return h
def jacobian(v,k):
    """columns: d(eqs)/d(unknown) by exact finite differences? No: use linearity: eqs are polynomial; use
       derivative via gradient: d/dc of coeff vector = coefficients of gradH(phi).(basis perturbation)."""
    c0,coef=coef_from_vec(v); phi=build(c0,coef,k); g=gradH(*phi)
    J=[]
    for (c,j) in unk:
        coord = 1 if (c,j)==(0,0) else c; dgr = 0 if (c,j)==(0,0) else j
        col=pmul(g[coord],[0]*dgr+[1]); col=col+[0]*(neq-len(col)); J.append(col[:neq])
    # J[i] is column i ; return matrix rows x cols
    return [[J[i][r] for i in range(nunk)] for r in range(neq)]
def solve_mod(M,b):
    """Gauss-Jordan mod p; returns (solution or None, rank)."""
    n=len(M); m=len(M[0]); A=[row[:]+[bb] for row,bb in zip(M,b)]
    r=0; piv=[]
    for c in range(m):
        pr=None
        for i in range(r,n):
            if A[i][c]%p: pr=i;break
        if pr is None: continue
        A[r],A[pr]=A[pr],A[r]; inv=pow(A[r][c],p-2,p); A[r]=[(x*inv)%p for x in A[r]]
        for i in range(n):
            if i!=r and A[i][c]%p:
                f=A[i][c]; A[i]=[(x-f*y)%p for x,y in zip(A[i],A[r])]
        piv.append(c); r+=1
    for i in range(r,n):
        if A[i][m]%p: return None,r
    x=[0]*m
    for i,c in enumerate(piv): x[c]=A[i][m]
    return x,r
def det_poly(Mp):
    """determinant of a matrix of polynomials (lists) mod p via fraction-free? Use evaluation/interpolation at many points of an extension?  Simpler: expand Leibniz for 5x5 (120 terms)."""
    n=len(Mp); import itertools as it
    det=[0]
    for perm in it.permutations(range(n)):
        sgn=1
        for i in range(n):
            for j in range(i+1,n):
                if perm[i]>perm[j]: sgn=-sgn
        term=[1]
        for i in range(n): term=pmul(term,Mp[i][perm[i]])
        det=padd(det,psc(sgn%p,term))
    while len(det)>1 and det[-1]==0: det.pop()
    return det
random.seed(int(sys.argv[5]) if len(sys.argv)>5 else 1)
found=[]; tried=0; t0=time.time()
maxdeg=max(deg.values())
free_per_degree={j:[(c,j) for c in (1,2,3,4) if j<=deg[c]] for j in range(1,maxdeg+1)}
while len(found)<nseeds and time.time()-t0<float(sys.argv[6]) if len(sys.argv)>6 else len(found)<nseeds:
    k=tuple(random.randrange(p) for _ in range(3))
    # endpoint: z0=0,z1=c0,z2=k0,s=k1,t=k2 on H=0
    for c0 in range(p):
        h0=H([0],[c0],[k[0]],[k[1]],[k[2]])[0]
        if h0%p: continue
        g=gradH([0],[c0],[k[0]],[k[1]],[k[2]]); g=[x[0] for x in g]
        if all(x%p==0 for x in g): continue
        # enumerate degree by degree: at degree j, the tau^j equation is affine-linear in degree-j coefficients
        # (coefficient = g . (deg-j coeffs) + lower-order stuff).  Enumerate all assignments satisfying it.
        partial=[{}]
        for j in range(1,maxdeg+1):
            new=[]
            vars_j=free_per_degree[j]
            for coef in partial:
                # compute lower-order contribution: evaluate tau^j coefficient with degree-j coeffs = 0
                base=dict(coef); 
                for v in vars_j: base[v]=0
                # pad higher degrees with 0 for evaluation
                full=dict(base)
                for jj in range(j+1,maxdeg+1):
                    for v in free_per_degree[jj]: full[v]=0
                hj=H(*build(c0,full,k)); hj=hj+[0]*(6*e+1-len(hj)); lo=hj[j]
                lin=[g[c] for (c,_) in vars_j]
                # need lo + sum lin_i x_i = 0
                for vals in itertools.product(range(p),repeat=len(vars_j)):
                    if (lo+sum(l*x for l,x in zip(lin,vals)))%p==0:
                        d=dict(coef); d.update(dict(zip(vars_j,vals))); new.append(d)
            partial=new
        for coef in partial:
            tried+=1
            v=[c0]+[coef[u] for u in unk[1:]]
            if any(residual(v,k)): continue
            J=jacobian(v,k)
            sol,rank=solve_mod(J,[1]+[0]*(neq-1))
            if rank<nunk: continue
            # psi: normal with gradH(phi).psi = 1 : psi has the same degree pattern as unknowns
            psi_coef=sol
            phi=build(*coef_from_vec(v),k)
            # dphi/dk_m: implicit: J dc/dk = -dE/dk ; dE/dk_m = coeffs of g_coord(m)*1
            gfull=gradH(*phi)
            cols=[]
            dphi=[[ [0] for _ in range(5)] for _ in range(3)]
            for m,coord in enumerate((2,3,4)):
                rhs=gfull[coord]+[0]*neq; rhs=[(-x)%p for x in rhs[:neq]]
                dc,_=solve_mod(J,rhs)
                vec=[[0] for _ in range(5)]; vec[coord]=[1]
                for i,(c,j) in enumerate(unk):
                    cc = 1 if (c,j)==(0,0) else c; dd = 0 if (c,j)==(0,0) else j
                    vec[cc]=padd(vec[cc],[0]*dd+[dc[i]])
                dphi[m]=vec
            psi=[[0] for _ in range(5)]
            for i,(c,j) in enumerate(unk):
                cc = 1 if (c,j)==(0,0) else c; dd = 0 if (c,j)==(0,0) else j
                psi[cc]=padd(psi[cc],[0]*dd+[psi_coef[i]])
            dtau=[[ (i*x)%p for i,x in enumerate(P)][1:] or [0] for P in phi]
            M=[[dtau[r],dphi[0][r],dphi[1][r],dphi[2][r],psi[r]] for r in range(5)]
            Jdet=det_poly(M)
            if Jdet==[0]: continue
            # check normal identity
            chk=[0]
            for r in range(5): chk=padd(chk,pmul(gfull[r],psi[r]))
            assert chk[0]==1 and all(x==0 for x in chk[1:]), chk
            # common root check of coordinates over F_p (vertex avoidance, weak)
            found.append({"p":p,"e":e,"pattern":pattern,"k":k,"c0":c0,"coef":{f"{c},{j}":coef[(c,j)] for (c,j) in coef},
                          "phi":{"z0":phi[0],"z1":phi[1],"z2":phi[2],"s":phi[3],"t":phi[4]},"psi":psi,"Jdet":Jdet,"jac_rank":rank,
                          "tried":tried,"seconds":round(time.time()-t0,1)})
            print("SEED", json.dumps(found[-1])); sys.stdout.flush()
            if len(found)>=nseeds: break
        if len(found)>=nseeds: break
print(f"done: found={len(found)} tried={tried} in {time.time()-t0:.1f}s")
json.dump(found, open(f'/tmp/hodge/results/polarized-k3-remainder/construct/c2_seeds_p{p}_{pattern}.json','w'))
