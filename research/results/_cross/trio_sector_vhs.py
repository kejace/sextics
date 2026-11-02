#!/usr/bin/env python3
"""Sector-variation test (generalisation of ctr construct c5_sector_vhs.py) for a lane with a diagonal symmetry group
and optional exchange symmetries.  For the symmetric family F + sum u_k m_k (m_k = invariant sextics under the full
visible group Gamma) and the target residue class xi = [A Omega/F^2] in its Gamma-character sector S_A:
  (a) dim Def.xi in R_12^chi      (= h^{2,2}(S_A) if the sector is 'generated' by xi at first order)
  (b) dim Def^2.xi in R_18^chi    (= h^{1,3}(S_A) ...)
  (c) dim K = {v in R_6^chi (epsilon part) : Def.v = 0 in R_12}
Ranks mod p (python-flint); rank_p <= rank_Q, dim ker_p >= dim ker_Q are the rigorous directions.
Usage: python3 trio_sector_vhs.py <lane> <out.json>
"""
import os, sys, json, time, itertools, re
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[v]='1'
from fractions import Fraction
from flint import nmod_mat
import sympy as sp

lane, out_path = sys.argv[1], sys.argv[2]
X = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
def read_rec(lane):
    txt = open(f'/home/kejace/devel/sextics/data/{lane}/instance.rec').read()
    def field(name):
        m = re.search(r'^'+name+r': (.*(?:\n  .*)*)', txt, re.M)
        return m.group(1).replace('\n', ' ')
    return field('F'), field('A')
Fs, As = read_rec(lane)
Fp = sp.Poly(sp.sympify(Fs), *X); Ap = sp.Poly(sp.sympify(As), *X)
den = sp.ilcm(*[sp.fraction(sp.Rational(c))[1] for c in Fp.coeffs()])
Fz = {tuple(e): int(c*den) for e, c in zip(Fp.monoms(), Fp.coeffs())}
A = {tuple(e): int(c) for e, c in zip(Ap.monoms(), Ap.coeffs())}
# lane-specific visible symmetry: diagonal weights mod 6 and exchange permutations (as coordinate permutations, with
# the sign the target residue class has under them: +1/-1).  Only the invariant family and the chi-sector are used.
CFG = {
 'cyclic-triple-residual': dict(W=[(0,0,2,0,0,0),(3,3,3,0,0,0),(2,4,0,4,4,0)], EX=[]),
 'reflection-coupled-residual': dict(W=[(0,0,0,0,0,3),(0,0,0,3,3,0)], EX=[]),
 'quadric-contact-residual': dict(W=[(3,0,0,3,0,0),(0,0,3,0,0,3),(0,1,0,0,5,0),(0,0,0,3,3,3)],
                                  EX=[((3,1,2,0,4,5),+1), ((0,4,2,3,1,5),+1)]),   # s0: X0<->Y0, r1: X1<->Y1; POLYNOMIAL sign of A (+1; residue sign -1 comes from det)
}[lane]
W, EX = CFG['W'], CFG['EX']
def chi(e): return tuple(sum(e[i]*w[i] for i in range(6)) % 6 for w in W)
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
            f=list(e); f[i]-=1; out[tuple(f)] = out.get(tuple(f),0) + v*e[i]
    return out
def padd(P, Q, s=1):
    out=dict(P)
    for e,v in Q.items(): out[e]=out.get(e,0)+s*v
    return {e:v for e,v in out.items() if v}
def perm_poly(P, perm):   # x_i -> x_{perm[i]}
    return {tuple(e[perm.index(j)] for j in range(6)): v for e,v in P.items()}
def mulpoly(P, Q):
    out={}
    for e,v in P.items():
        for f,w in Q.items():
            g=tuple(a+b for a,b in zip(e,f)); out[g]=out.get(g,0)+v*w
    return {e:v for e,v in out.items() if v}
dF = [pderiv(Fz,i) for i in range(6)]
chiF = chi(next(iter(Fz))); assert all(chi(e)==chiF for e in Fz) and chiF == (0,)*len(W), 'F must be invariant'
chiA = chi(next(iter(A)))
# check exchanges are symmetries of F and record the polynomial sign of A
for perm, sgn in EX:
    assert perm_poly(Fz, perm) == Fz, ('exchange not a symmetry', perm)
E6 = [tuple(1 if j==i else 0 for j in range(6)) for i in range(6)]
def csub(a,b): return tuple((x-y)%6 for x,y in zip(a,b))
def cadd(a,b): return tuple((x+y)%6 for x,y in zip(a,b))
# invariant family: chi = 0 monomials, symmetrised over the exchange group (orbit sums), linearly independent
def group_elems():
    els=[tuple(range(6))]
    frontier=list(els)
    while frontier:
        new=[]
        for g in frontier:
            for perm,_ in EX:
                h=tuple(perm[g[i]] for i in range(6))
                if h not in els: els.append(h); new.append(h)
        frontier=new
    return els
GE = group_elems()
def orbit_sum(m, sign_fn=None):
    P={}
    for g in GE:
        P = padd(P, perm_poly({m:1}, g))
    return P
seen=set(); Def=[]
for m in mons(6):
    if chi(m)!=(0,)*len(W): continue
    P = orbit_sum(m); key=tuple(sorted(P.items()))
    if key in seen: continue
    seen.add(key); Def.append(P)
print(lane, 'invariant sextics (orbit sums):', len(Def), 'exchange group order', len(GE), flush=True)
def sector_data(k, chik, p):
    S = [e for e in mons(k) if chi(e)==chik]
    idx = {e:i for i,e in enumerate(S)}
    gens = []
    for i in range(6):
        ci = csub(chiF, chi(E6[i]))
        for m in mons(k-5):
            if cadd(chi(m), ci) == chik:
                gens.append({tuple(a+b for a,b in zip(e,m)): v for e,v in dF[i].items()})
    M = nmod_mat(len(S), len(gens), p)
    for j,g in enumerate(gens):
        for e,v in g.items(): M[idx[e], j] = v % p
    return S, idx, M, M.rank()
def quotient_rank(S, idx, M, rJ, vecs, p):
    N = nmod_mat(len(S), M.ncols()+len(vecs), p)
    for i in range(len(S)):
        for j in range(M.ncols()):
            if M[i,j] != 0: N[i,j] = M[i,j]
    for j,vv in enumerate(vecs):
        for e,v in vv.items(): N[idx[e], M.ncols()+j] = v % p
    return N.rank() - rJ
out = {'lane': lane, 'n_Def': len(Def), 'chiA': chiA}
PRIMES = tuple(int(x) for x in os.environ.get('TRIO_PRIMES','1000003,998244353').split(','))
for p in PRIMES:
    t0=time.time()
    S6, idx6, M6, r6 = sector_data(6, chiA, p)
    S12, idx12, M12, r12 = sector_data(12, chiA, p)
    dims = {'R6_diag_sector': len(S6)-r6, 'R12_diag_sector': len(S12)-r12}
    Am = [mulpoly(A, m) for m in Def]
    a = quotient_rank(S12, idx12, M12, r12, Am, p)
    print({'p':p, 'dims':dims, 'dim_Def_xi_in_R12': a, 't': round(time.time()-t0,1)}, flush=True)
    if os.environ.get('TRIO_SKIP_R18'):
        b = None
    else:
        S18, idx18, M18, r18 = sector_data(18, chiA, p)
        dims['R18_diag_sector'] = len(S18)-r18
        Amm = [mulpoly(Am[i], Def[j]) for i in range(len(Def)) for j in range(i, len(Def))]
        b = quotient_rank(S18, idx18, M18, r18, Amm, p)
        print({'p':p, 'dim_Def2_xi_in_R18': b, 't': round(time.time()-t0,1)}, flush=True)
    # (c) kernel: v in epsilon-part of S6^chi with Def.v in J12.  epsilon-part basis: projector with the target signs.
    def eps_project(m):
        P={}
        for g in GE:
            # sign of g on the class: product of signs of generators along a word; for a Klein/(Z2)^2 exchange group the
            # sign is multiplicative in the coordinates moved: compute via parity of each generator's occurrence
            P = padd(P, perm_poly({m:1}, g), sign_of(g))
        return P
    def sign_of(g):
        s=1
        # express g as product of generator perms by BFS (small group)
        from collections import deque
        start=tuple(range(6)); q=deque([(start,1)]); seen={start:1}
        while q:
            h,sg=q.popleft()
            if h==g: return sg
            for perm,sgn in EX:
                h2=tuple(perm[h[i]] for i in range(6))
                if h2 not in seen: seen[h2]=sg*sgn; q.append((h2,sg*sgn))
        return seen[g]
    basis6=[]; seen6=set()
    for m in S6:
        P=eps_project(m); key=tuple(sorted(P.items()))
        if not P or key in seen6: continue
        seen6.add(key); basis6.append(P)
    # normal form modulo J12 via rref of generators^T
    MT = nmod_mat(M12.ncols(), len(S12), p)
    for i in range(len(S12)):
        for j in range(M12.ncols()):
            if M12[i,j] != 0: MT[j,i] = M12[i,j]
    MTr, rkT = MT.rref()
    piv=[]
    for i in range(rkT):
        j=0
        while j < len(S12) and MTr[i,j]==0: j+=1
        piv.append(j)
    pivset=set(piv); nonpiv=[j for j in range(len(S12)) if j not in pivset]
    # projection P: S12 -> S12/J12 in the nonpivot coordinates, as a matrix (len(nonpiv) x len(S12))
    npv = len(nonpiv)
    P = nmod_mat(npv, len(S12), p)
    for a_, j in enumerate(nonpiv): P[a_, j] = 1
    for i, pi in enumerate(piv):
        for a_, j in enumerate(nonpiv):
            val = MTr[i, j]
            if val != 0: P[a_, pi] = (int(P[a_, pi]) - int(val)) % p
    blocks = []
    for m in Def:
        Mult = nmod_mat(len(S12), len(basis6), p)
        for jcol, b6 in enumerate(basis6):
            prod = mulpoly(b6, m)
            for e, v in prod.items(): Mult[idx12[e], jcol] = (int(Mult[idx12[e], jcol]) + v) % p
        blocks.append(P * Mult)
    big = nmod_mat(len(Def)*npv, len(basis6), p)
    for k, blk in enumerate(blocks):
        for i in range(npv):
            for j in range(len(basis6)):
                if blk[i, j] != 0: big[k*npv+i, j] = int(blk[i, j])
    kerdim = len(basis6) - big.rank()
    # dim of J6 inside the epsilon part: rank of J6 generators projected
    J6eps = quotient_rank(S6, idx6, nmod_mat(len(S6),0,p), 0, [], p)  # placeholder 0
    # epsilon-part of R6: rank of basis6 modulo J6
    r6eps = quotient_rank(S6, idx6, M6, r6, basis6, p)
    allprod = [mulpoly(b6, m) for b6 in basis6 for m in Def]
    defR6 = quotient_rank(S12, idx12, M12, r12, allprod, p)
    rec = {'p':p, 'dims':dims, 'dim_eps_S6': len(basis6), 'dim_eps_R6': r6eps, 'dim_Def_xi_in_R12': a,
           'dim_Def2_xi_in_R18': b, 'dim_ker_in_eps_S6': kerdim, 'dim_K_annihilated_by_Def_in_R6eps': kerdim-(len(basis6)-r6eps), 'dim_Def_R6eps_in_R12': defR6,
           'seconds': round(time.time()-t0,1)}
    print(rec, flush=True); out[str(p)]=rec
    json.dump(out, open(out_path,'w'), indent=1)
json.dump(out, open(out_path,'w'), indent=1)
print('done')
