"""Referee library for one-reflection-skew (independent of the lane's w2x.py / engine).
Polynomials: dict {exponent tuple: int}. Variables X0,X1,X2,Y0,Y1,Y2 (indices 0..5).  F7 = 7F (integer)."""
import itertools, json, time, sys
from flint import fmpz_mat, nmod_mat, fmpz

NV = 6
# 7F
F7 = {(6,0,0,0,0,0):7, (3,1,0,1,1,0):7, (1,1,0,3,1,0):1, (1,1,0,1,3,0):-1, (0,6,0,0,0,0):7,
      (0,3,1,1,1,0):7, (0,0,6,0,0,0):7, (0,0,0,6,0,0):7, (0,0,0,0,6,0):7, (0,0,0,0,0,6):7}
A = {(3,0,0,3,0,0):1}
G = [(3,3,3,0,0,0),(0,3,0,3,0,0),(4,1,4,4,1,0)]   # diagonal symmetry generators (weights mod 6)

def chi(e):
    return tuple(sum(a*b for a,b in zip(e,g)) % 6 for g in G)
def chi_add(*cs):
    return tuple(sum(c[i] for c in cs) % 6 for i in range(3))
def chi_neg(c):
    return tuple((-x) % 6 for x in c)
def mons(n, nv=NV):
    out = []
    for c in itertools.combinations_with_replacement(range(nv), n):
        e = [0]*NV
        for i in c: e[i] += 1
        out.append(tuple(e))
    return out
def secmons(n, ch, nv=NV):
    return [e for e in mons(n, nv) if chi(e) == ch]
def padd(p, q, s=1):
    r = dict(p)
    for e, c in q.items():
        r[e] = r.get(e, 0) + s*c
        if r[e] == 0: del r[e]
    return r
def pmul(p, q):
    r = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = tuple(a+b for a,b in zip(e1,e2)); r[e] = r.get(e,0) + c1*c2
    return {e:c for e,c in r.items() if c}
def pscale(p, s):
    return {e: c*s for e, c in p.items()} if s else {}
def pderiv(p, i):
    r = {}
    for e, c in p.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1; r[tuple(e2)] = r.get(tuple(e2),0) + c*e[i]
    return {e:c for e,c in r.items() if c}
def mono(e, c=1):
    return {tuple(e): c}
def pdeg(p):
    return max(sum(e) for e in p) if p else -1
def pchi(p):
    cs = {chi(e) for e in p}
    return cs.pop() if len(cs) == 1 else None
def restrict(p, nv):
    """drop variables >= nv (must not occur)"""
    for e in p: assert all(x == 0 for x in e[nv:])
    return {e[:nv]+(0,)*(NV-nv): c for e, c in p.items()}

Fi = [pderiv(F7, i) for i in range(NV)]
f7 = padd(F7, {(0,0,0,0,0,6):7}, -1)          # 7f (Y2-free part)
fi = [pderiv(f7, i) for i in range(5)]
CHI_A = pchi(A)

def build(cols, rowidx, p=None):
    """cols: list of polys; rowidx: dict exponent->row.  Returns fmpz_mat (p None) or nmod_mat."""
    nr = len(rowidx); nc = len(cols)
    M = fmpz_mat(nr, nc) if p is None else nmod_mat(nr, nc, p)
    for j, pol in enumerate(cols):
        for e, c in pol.items():
            M[rowidx[e], j] = c if p is None else c % p
    return M

def e1_kernel_dim(d, ch, nv=NV, p=None, exact=True):
    """dim {D in S_d[ch] : A D in J_{d+6}} where J = (partials of F7) (nv=6) or (partials of f7) (nv=5)."""
    Fs = Fi if nv == 6 else fi
    rows = secmons(d+6, chi_add(ch, CHI_A), nv); rowidx = {e:i for i,e in enumerate(rows)}
    Dm = secmons(d, ch, nv)
    jcols = []
    for i in range(nv):
        for m in secmons(d+1, chi_add(ch, CHI_A, chi(tuple(1 if k==i else 0 for k in range(NV)))), nv):
            jcols.append(pmul(mono(m), Fs[i]))
    acols = [pmul(mono(m), A) for m in Dm]
    if not Dm: return 0, 0, len(rows)
    MJ = build(jcols, rowidx, p); rJ = MJ.rank()
    MA = build(acols + jcols, rowidx, p); rA = MA.rank()
    return len(Dm) - (rA - rJ), rJ, len(rows)

def j_dim(d, ch, nv=NV, p=None):
    if d < 5: return 0
    Fs = Fi if nv == 6 else fi
    rows = secmons(d, ch, nv); rowidx = {e:i for i,e in enumerate(rows)}
    cols = []
    for i in range(nv):
        for m in secmons(d-5, chi_add(ch, chi(tuple(1 if k==i else 0 for k in range(NV)))), nv):
            cols.append(pmul(mono(m), Fs[i]))
    if not cols or not rows: return 0
    return build(cols, rowidx, p).rank()

def all_sectors():
    return [(a,b,c) for a in (0,3) for b in (0,3) for c in range(6)]

def log(path, **rec):
    rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    with open(path, 'a') as fh: fh.write(json.dumps(rec, default=str) + '\n')
