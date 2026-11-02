#!/usr/bin/env python3
"""Referee replay for off-paired (q=1, A=X0^3 Y0^3): INDEPENDENT code (python-flint 0.9 only; no hodge engine, no w2lib).
Exact ranks over Q via fmpz_mat (columns scaled to integers), optional modular ranks at a fresh prime.
usage: python3 ref_replay.py <check> [...]   checks: D0 f6 blocks d10 e1d11 sectors
"""
import sys, json, time, itertools, math
from flint import fmpq_mpoly_ctx, Ordering, fmpq, fmpz_mat, nmod_mat
VARS = ('X0','X1','X2','Y0','Y1','Y2'); NV = 6
ctx = fmpq_mpoly_ctx.get(VARS, ordering=Ordering.lex)
Z = list(ctx.gens()); X0,X1,X2,Y0,Y1,Y2 = Z
h = fmpq
F = X0**6 + X0**3*X1*Y0*Y1*h(1,2) + X0**2*X1*X2*Y0*Y2*h(1,3) + X1**6 + X1**3*X2*Y1*Y2*h(1,2) + X2**6 + Y0**6 + Y1**6 + Y2**6
A = X0**3*Y0**3
Fi = [F.derivative(i) for i in range(NV)]
WG=(0,2,0,1,3,3); WH=(0,0,1,0,0,5)
def char(e): return (sum(a*b for a,b in zip(WG,e))%6, sum(a*b for a,b in zip(WH,e))%6)
def cadd(a,b): return ((a[0]+b[0])%6,(a[1]+b[1])%6)
VARCHAR = [char(tuple(1 if j==i else 0 for j in range(NV))) for i in range(NV)]
def char_of(P):
    cs = {char(e) for e in P.to_dict()}; assert len(cs)==1, cs; return cs.pop()
_mc={}
def monos(deg, c=None):
    if deg < 0: return []
    key=(deg,c)
    if key in _mc: return _mc[key]
    out=[]
    for comb in itertools.combinations(range(deg+NV-1), NV-1):
        e=[]; prev=-1
        for x in comb: e.append(x-prev-1); prev=x
        e.append(deg+NV-2-prev); e=tuple(e)
        if c is None or char(e)==c: out.append(e)
    _mc[key]=out; return out
def mono(e): return ctx.from_dict({e: h(1)})
def degree(P): return max(sum(e) for e in P.to_dict())
# group sanity: g, h preserve F; A has character (3,0)
assert all(char(e)==(0,0) for e in F.to_dict()) and char_of(A)==(3,0)

def flat_matrix(cols):
    """columns = fmpq_mpoly, each scaled to integers (column scaling keeps rank/consistency). Returns (flat, R, C)."""
    rows={}; dcols=[]
    for P in cols:
        d=P.to_dict(); dcols.append(d)
        for e in d: rows.setdefault(e, len(rows))
    R=len(rows); C=len(dcols); flat=[0]*(R*C)
    for j,d in enumerate(dcols):
        L=1
        for v in d.values(): L = L*int(v.q)//math.gcd(L,int(v.q))
        for e,v in d.items(): flat[rows[e]*C+j] = int(v.p)*(L//int(v.q))
    return flat, R, C

def ranks(flat, R, C, p=None):
    """(rank of first C-1 columns, rank of all C columns)."""
    sub=[flat[i*C+j] for i in range(R) for j in range(C-1)]
    if p is None:
        return int(fmpz_mat(R, C-1, sub).rank()), int(fmpz_mat(R, C, flat).rank())
    return int(nmod_mat(R, C-1, [x%p for x in sub], p).rank()), int(nmod_mat(R, C, [x%p for x in flat], p).rank())

def explicit_kernel(d, cD):
    cL=cadd(char_of(A), cD)
    return len(monos(d,cL)) + sum(len(monos(d-5, cadd(cL,VARCHAR[i]))) for i in range(NV)) - len(monos(d-6,cL))

def fixedD(D, label, p=None):
    """Euler-normalized fixed-D system in D's character sector: -D N(F) + F B_D(N) - F^2 Q = -A D^2."""
    t0=time.time(); d=degree(D); cD=char_of(D); cA=char_of(A)
    cols=[]
    for i in range(NV):
        for e in monos(d+1, cadd(cadd(cA,cD),VARCHAR[i])):
            m=mono(e); cols.append(-D*(m*Fi[i]) + F*(D*m.derivative(i) - m*D.derivative(i)))
    nN=len(cols); F2=F*F
    for e in monos(2*d-6, cadd(cA,cadd(cD,cD))): cols.append(-F2*mono(e))
    nQ=len(cols)-nN
    cols.append(-A*D*D)
    flat,R,C = flat_matrix(cols)
    r0, ra = ranks(flat,R,C,p)
    rec={'check':'fixedD','label':label,'d':d,'char':cD,'field':('Q' if p is None else f'F_{p}'),'rows':R,'cols':C-1,'nN':nN,'nQ':nQ,
         'rank_M':r0,'rank_aug':ra,'nullity':C-1-r0,'explicit_kernel':explicit_kernel(d,cD),'inconsistent':ra==r0+1,'seconds':round(time.time()-t0,1)}
    if p is None: rec['verdict'] = 'EXACT: no certificate with this denominator' if ra==r0+1 else 'EXACT: consistent (!!)'
    else: rec['verdict'] = ('MODULAR inconsistent; rank sandwich ' + ('HOLDS (nullity == explicit kernel) => exact exclusion' if C-1-r0==rec['explicit_kernel'] else 'does NOT hold')) if ra==r0+1 else 'MODULAR consistent'
    print(json.dumps(rec, default=int), flush=True); return rec

def e1(D, label, p=None):
    t0=time.time(); d=degree(D); cD=char_of(D); cA=char_of(A)
    cols=[mono(e)*Fi[i] for i in range(NV) for e in monos(d+1, cadd(cadd(cA,cD),VARCHAR[i]))]
    cols.append(A*D)
    flat,R,C=flat_matrix(cols); r0,ra=ranks(flat,R,C,p)
    rec={'check':'E1','label':label,'d':d,'char':cD,'field':('Q' if p is None else f'F_{p}'),'rows':R,'cols':C-1,'rank_J':r0,'rank_aug':ra,'e1_passes':r0==ra,'seconds':round(time.time()-t0,1)}
    rec['verdict']=('A*D in J(F)' if r0==ra else 'A*D NOT in J(F)') + (' (exact over Q)' if p is None else ' (modular only)')
    print(json.dumps(rec, default=int), flush=True); return rec

def colon_sector(d, c):
    """Exact rational kernel of N(F) - A D = 0 with D in S_d[c]; returns (list of D polys spanning the colon sector mod gauge, colon_dim, syz_dim, gauge_dim)."""
    cA=char_of(A); Dm=monos(d,c)
    cols=[-A*mono(e) for e in Dm]
    for i in range(NV):
        for e in monos(d+1, cadd(cadd(cA,c),VARCHAR[i])): cols.append(mono(e)*Fi[i])
    flat,R,C=flat_matrix(cols)      # D columns have integer coefficients (scale 1), so kernel D-parts are exact
    M=fmpz_mat(R,C,flat); X,nul=M.nullspace(); nul=int(nul)
    nD=len(Dm)
    Dvecs=[[int(X[j,k]) for j in range(nD)] for k in range(nul)]
    colon_dim = int(fmpz_mat(nul,nD,[v for r in Dvecs for v in r]).rank()) if nul else 0
    gauge=[F*mono(e) for e in monos(d-6,c)]
    gidx={e:i for i,e in enumerate(Dm)}
    gvecs=[]
    for g in gauge:
        v=[0]*nD
        for e,x in g.to_dict().items(): v[gidx[e]]=int(x.p)*1 if int(x.q)==1 else None
        gvecs.append(v)
    gauge_dim=int(fmpz_mat(len(gvecs),nD,[v for r in gvecs for v in r]).rank()) if gvecs else 0
    # choose D's independent modulo gauge
    chosen=[]; cur=gauge_dim
    for v in Dvecs:
        if not any(v): continue
        rows_=gvecs+chosen+[v]
        rk=int(fmpz_mat(len(rows_),nD,[x for r in rows_ for x in r]).rank())
        if rk>cur+len(chosen): chosen.append(v)
    Ds=[]
    for v in chosen:
        g=0
        for x in v: g=math.gcd(g,x)
        Ds.append(ctx.from_dict({Dm[j]: h(x//g) for j,x in enumerate(v) if x}))
    return Ds, colon_dim, nul-colon_dim, gauge_dim

M_=X2**6-Y2**6; S_=X2**6+Y2**6; H_=X0**6+X1**6+Y0**6+Y1**6+X0**3*X1*Y0*Y1*h(1,2); C_=X1**3*Y1*h(1,2)+X0**2*X1*Y0*h(1,3)
D0=X0**4+X0*X1*Y0*Y1*h(1,4)+X1*X2*Y0*Y2*h(1,9); f6=X0**6+X1**6+X2**6; g6=Y0**6+Y1**6+Y2**6
FX0,FX1,FX2,FY0,FY1,FY2=Fi
PNEW=1000000009   # fresh prime (lane used 1000003, 2147483647, 1073741827, 2147483629, 2147483587; aggregator 1000000007, 998244353)

if __name__=='__main__':
    for chk in sys.argv[1:]:
        if chk=='D0':
            assert 6*X0*D0 == FX0, 'D0 != F_X0/(6 X0)'; print('D0 == F_X0/(6*X0): exact identity holds')
            fixedD(D0,'D0'); e1(D0,'D0')
        elif chk=='f6':
            e1(f6,'f6'); e1(g6,'g6'); e1(f6+g6,'f6+g6')
        elif chk=='blocks':
            fixedD(X1*X2*FX1,'X1*X2*FX1'); fixedD(FX0*X1*Y1,'FX0*X1*Y1')
        elif chk=='d10':
            fixedD(X1*D0*FX2,'X1*D0*FX2', p=PNEW)
        elif chk=='e1d11':
            e1(X0*H_*C_,'X0*H*C', p=PNEW); e1(X0*H_*C_,'X0*H*C')
        elif chk=='sectors':
            for c in [(1,2),(2,2),(1,4),(2,4),(0,3)]:
                Ds,cd,sd,gd=colon_sector(7,c)
                print(json.dumps(default=int, obj={'check':'colon-sector','d':7,'char':c,'colon_dim':cd,'syz_dim':sd,'gauge_dim':gd,'nt':len(Ds)}), flush=True)
                if len(Ds)==1:
                    print('  unique D (up to scalar):', str(Ds[0])[:160])
                    fixedD(Ds[0], f'sector d7 {c} unique D')
        else:
            print('unknown check', chk)
