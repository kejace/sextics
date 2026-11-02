# C5: infinitesimal-variation test for the exact G-sector S_A of xi = [X0^3 Y0^3] over the G-symmetric family (49 invariant
# sextic monomials): (a) dim Def_G.xi in R_12^chi (= S_A^{2,2}, expect 119), (b) dim Def_G^2.xi in R_18^chi (= S_A^{1,3}, 26),
# (c) K = {v in R_6^chi : Def_G.v = 0 in R_12} (expect 0).  (a)=119, (b)=26, (c)=0 (mod p, hence over Q) imply that the
# smallest flat sub-VHS of S_A over the family containing xi is all of S_A (see part5 §5.3 / part6).  Ranks mod p with flint.
import os, sys, json, time, itertools
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[v]='1'
sys.path.insert(0,'/tmp/hodge/results/cyclic-triple-residual'); sys.path.insert(0,'/home/kejace/devel/sextics-work')
from fractions import Fraction
from flint import nmod_mat
from w2core import load, core
F, A, q = load()
Fz = {e: int(7*Fraction(v)) for e, v in F.items()}            # 7F integral
assert all(7*Fraction(v) == int(7*Fraction(v)) for v in F.values())
W = [(0,0,2,0,0,0), (3,3,3,0,0,0), (2,4,0,4,4,0)]              # mu3, tau, g3 (weights mod 6)
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
dF = [pderiv(Fz,i) for i in range(6)]
chiF = chi(next(iter(Fz)))
assert all(chi(e)==chiF for e in Fz), 'F not an eigenvector'
chiA = chi(next(iter(A)))
def cadd(a,b): return tuple((x+y)%6 for x,y in zip(a,b))
def csub(a,b): return tuple((x-y)%6 for x,y in zip(a,b))
E6 = [tuple(1 if j==i else 0 for j in range(6)) for i in range(6)]
def mulmono(P, m):
    return {tuple(a+b for a,b in zip(e,m)): v for e,v in P.items()}
def sector_data(k, chik, p):
    """monomials of degree k with character chik, J_k generators in that sector, rref projection to the quotient."""
    S = [e for e in mons(k) if chi(e)==chik]
    idx = {e:i for i,e in enumerate(S)}
    gens = []
    for i in range(6):
        ci = csub(chiF, chi(E6[i]))              # character of dF_i
        for m in mons(k-5):
            if cadd(chi(m), ci) == chik:
                gens.append(mulmono(dF[i], m))
    M = nmod_mat(len(S), len(gens), p)
    for j,g in enumerate(gens):
        for e,v in g.items(): M[idx[e], j] = v % p
    r = M.rank()
    return S, idx, gens, M, r
def quotient_rank(S, idx, M, rJ, vecs, p):
    """rank of the images of vecs in S/J: rank([J | vecs]) - rank(J)."""
    N = nmod_mat(len(S), M.ncols()+len(vecs), p)
    for i in range(len(S)):
        for j in range(M.ncols()):
            if M[i,j] != 0: N[i,j] = M[i,j]
    for j,vv in enumerate(vecs):
        for e,v in vv.items(): N[idx[e], M.ncols()+j] = v % p
    return N.rank() - rJ
Def = [e for e in mons(6) if chi(e)==(0,0,0)]
print('G-invariant sextic monomials:', len(Def))
out = {'n_Def': len(Def)}
for p in (1000003, 998244353):
    t0=time.time()
    S6, idx6, g6, M6, r6 = sector_data(6, chiA, p)
    S12, idx12, g12, M12, r12 = sector_data(12, chiA, p)
    S18, idx18, g18, M18, r18 = sector_data(18, chiA, p)
    dims = {'R6': len(S6)-r6, 'R12': len(S12)-r12, 'R18': len(S18)-r18}
    Am = [mulmono(A, m) for m in Def]
    a = quotient_rank(S12, idx12, M12, r12, [ {e:int(v) for e,v in P.items()} for P in Am], p)
    Amm = [mulmono(mulmono(A, m), m2) for m, m2 in itertools.combinations_with_replacement(Def, 2)]
    b = quotient_rank(S18, idx18, M18, r18, [ {e:int(v) for e,v in P.items()} for P in Amm], p)
    # (c): K = {v in S6^chi : m v in J12 for all m} / J6^chi.  Build projection P12 from rref of J12.
    Mr, rk = M12.rref()
    # pivot columns of M12^T? we need a basis of the quotient S12/J12: use rref of the TRANSPOSE (row space = span of generators)
    MT = nmod_mat(M12.ncols(), len(S12), p)
    for i in range(len(S12)):
        for j in range(M12.ncols()):
            if M12[i,j] != 0: MT[j,i] = M12[i,j]
    MTr, rkT = MT.rref()
    piv = []
    for i in range(rkT):
        j=0
        while j < len(S12) and MTr[i,j]==0: j+=1
        piv.append(j)
    nonpiv = [j for j in range(len(S12)) if j not in set(piv)]   # quotient coordinates
    assert len(nonpiv) == dims['R12']
    # normal form of a vector w in S12: w - sum_i w[piv_i] * row_i  -> read off nonpiv coordinates
    rows = [[int(MTr[i,j]) for j in range(len(S12))] for i in range(rkT)]
    def nf(w):
        w = [x % p for x in w]
        for i,pi in enumerate(piv):
            if w[pi]:
                c = w[pi]
                row = rows[i]
                for j in range(len(S12)):
                    if row[j]: w[j] = (w[j] - c*row[j]) % p
        return [w[j] for j in nonpiv]
    big = nmod_mat(len(Def)*dims['R12'], len(S6), p)
    for jcol, b6 in enumerate(S6):
        for k, m in enumerate(Def):
            w = [0]*len(S12)
            e = tuple(x+y for x,y in zip(b6, m)); w[idx12[e]] = 1
            v = nf(w)
            for i,val in enumerate(v):
                if val: big[k*dims['R12']+i, jcol] = val
    kerdim = len(S6) - big.rank()
    Kdim = kerdim - (len(S6) - dims['R6'])     # subtract dim J6^chi (which lies in the kernel)
    rec = {'p': p, 'dims': dims, 'dim_DefG_xi_in_R12': a, 'dim_DefG2_xi_in_R18': b, 'dim_K_annihilated_by_Def': Kdim,
           'seconds': round(time.time()-t0,1)}
    print(rec, flush=True); out[str(p)] = rec
json.dump(out, open('/tmp/hodge/results/cyclic-triple-residual/wave2-construct/c5_sector_vhs.json','w'), indent=1)
