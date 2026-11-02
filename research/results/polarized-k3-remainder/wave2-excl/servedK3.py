"""K3-block projection of the served target space of a mu6-homogeneous denominator D ((U,V) frame), mod p.

For D of character b (all monomials have V-degree == b mod 6) the fixed-D certificate system with a VARIABLE target
A' in the K3 block {B'(U) + B''(V) : B', B'' sextics} (character 0) is mu6-graded: N_i has character b (i<3) / b+1 (i>=3),
Q has character 2b, and the identity lives in character 2b.  We solve the homogeneous system in (N, Q, A') and project
its nullspace to the A'-coordinates: this is the space of K3-block targets served by D (mod p).  Its image in the
Jacobian ring (R_f,6 (x) 1  (+)  1 (x) R_f,6, dimension 19 + 19) is the K3-class projection of served(D).
Rigor: served_Q(D) reduces into served_p(D), so a ZERO modular projection proves the rational projection is zero.
usage: servedK3.py <label> "<D expr in X0..X2 (=U), Y0..Y2 (=V)>" [p] [--out=FILE]
"""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
from flint import nmod_mat, fmpz_mat, fmpq_mat
from fractions import Fraction as Fr
import math
from hodge import core
from hodge.core import monomials, mul, add, sub, scale, pderiv, shift, power
from uvring import TernaryJacobianRing, load_f
NV = 6
def vdeg(e): return e[3] + e[4] + e[5]
def run_exact(label, Dexpr, verbose=True):
    """Exact rational version: integer model G = 14*Fhat, D scaled to integers, targets = K3-block monomials.
    served_Q = projection of the exact nullspace; class dimension modulo J(7f)_6(U) + J(7f)_6(V) over Q."""
    t0 = time.time()
    uv = json.load(open('/tmp/hodge/results/polarized-k3-remainder/uv_forms.json'))
    Fh = core.scale(core.parse(uv['Fhat']), 14)                       # integer coefficients
    D = core.parse(Dexpr)
    den = 1
    for v in D.values():
        fv = Fr(v); den = den * fv.denominator // math.gcd(den, fv.denominator)
    D = {e: int(Fr(v) * den) for e, v in D.items()}
    d = core.degree(D)
    chs = {vdeg(e) % 6 for e in D}
    assert len(chs) == 1, 'D must be mu6-homogeneous'
    b = chs.pop()
    q = 1; degN, degQ = d + 1, 2 * d - 6
    Fi = [pderiv(Fh, i) for i in range(NV)]
    Di = [pderiv(D, i) for i in range(NV)]
    FD = mul(Fh, D)
    Gc = [scale(add(scale(mul(D, Fi[i]), q), mul(Fh, Di[i])), -1) for i in range(NV)]
    Fq1 = scale(power(Fh, q + 1), -1)
    D2 = mul(D, D)
    cols, polys = [], []
    for i in range(NV):
        want = b % 6 if i < 3 else (b + 1) % 6
        for m in monomials(degN):
            if vdeg(m) % 6 != want: continue
            col = shift(Gc[i], m, 1)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i]))
            cols.append(('N', i, m)); polys.append(col)
    for m in monomials(degQ):
        if vdeg(m) % 6 == (2 * b) % 6:
            cols.append(('Q', -1, m)); polys.append(shift(Fq1, m, 1))
    tmons = [m for m in monomials(6) if vdeg(m) in (0, 6)]
    for m in tmons:
        cols.append(('A', -1, m)); polys.append(shift(D2, m, 1))
    rows = {}
    for col in polys:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    M = fmpz_mat(R, C)
    for j, col in enumerate(polys):
        for e, v in col.items(): M[rows[e], j] = int(v)
    X, nul = M.nullspace()           # X: C x C matrix whose first nul columns span the kernel (fmpz)
    a_idx = [j for j, c in enumerate(cols) if c[0] == 'A']
    served = []
    for k in range(nul):
        poly = {cols[j][2]: int(X[j, k]) for j in a_idx if int(X[j, k]) != 0}
        if poly: served.append(poly)
    # K3 classes over Q: reduce modulo J(7f)_6 on each block (U-sextics, V-sextics) via exact ranks
    f7 = {(e[0], e[1], e[2]): int(7 * Fr(v)) for e, v in core.parse(uv['f']).items()}
    def tern_gens():
        gens = []
        fi = [{tuple(ee for ee in (e[0]-(i==0), e[1]-(i==1), e[2]-(i==2))): v * e[i] for e, v in f7.items() if e[i]} for i in range(3)]
        for i in range(3):
            for m in [(a, bb, 1 - a - bb) for a in range(2) for bb in range(2 - a)]:
                gens.append({(e[0]+m[0], e[1]+m[1], e[2]+m[2]): v for e, v in fi[i].items()})
        return gens
    gens = tern_gens()                 # 9 generators of J(7f)_6
    mons6 = [(a, bb, 6 - a - bb) for a in range(7) for bb in range(7 - a)]
    idx6 = {m: i for i, m in enumerate(mons6)}
    def block_rank(vecs):
        # rank of [J-gens ; vecs] over Q minus rank of J-gens
        rowsM = [[0]*28 for _ in range(len(gens) + len(vecs))]
        for r, g in enumerate(gens):
            for e, v in g.items(): rowsM[r][idx6[e]] = v
        for r, vec in enumerate(vecs):
            for e, v in vec.items(): rowsM[len(gens) + r][idx6[e]] = v
        rJ = fmpz_mat([[0]*28 for _ in range(0)] or [rowsM[r] for r in range(len(gens))]).rank()
        return fmpz_mat(rowsM).rank() - rJ, rJ
    Uvecs = [{(e[0], e[1], e[2]): v for e, v in s_.items() if vdeg(e) == 0} for s_ in served]
    Vvecs = [{(e[3], e[4], e[5]): v for e, v in s_.items() if vdeg(e) == 6} for s_ in served]
    # joint class space: U-block and V-block classes are independent coordinates (19 + 19)
    cU, rJ = block_rank([u for u in Uvecs if u])
    cV, _ = block_rank([v for v in Vvecs if v])
    # joint rank (U and V coordinates together)
    rowsJ = []
    for g in gens:
        rowU = [0]*56; rowV = [0]*56
        for e, v in g.items(): rowU[idx6[e]] = v; rowV[28 + idx6[e]] = v
        rowsJ.append(rowU); rowsJ.append(rowV)
    rowsS = []
    for u, v in zip(Uvecs, Vvecs):
        row = [0]*56
        for e, c in u.items(): row[idx6[e]] = c
        for e, c in v.items(): row[28 + idx6[e]] = c
        rowsS.append(row)
    rankJ = fmpz_mat(rowsJ).rank()
    rankAll = fmpz_mat(rowsJ + rowsS).rank() if rowsS else rankJ
    res = {'label': label, 'D': Dexpr, 'd': d, 'char': b, 'arith': 'exact Q (fmpz nullspace)', 'rows': R, 'cols': C, 'nullity': nul,
           'served_K3_monomial_vectors': len(served), 'J6_rank_each_block': rJ, 'served_K3_class_dim_Q': rankAll - rankJ,
           'class_dim_U_block': cU, 'class_dim_V_block': cV, 'seconds': round(time.time() - t0, 1),
           'note': 'exact over Q: class dim 0 => the K3-block projection of served(D) is zero over Q (rigorous)'}
    if verbose: print('[servedK3-exact]', json.dumps(res), flush=True)
    return res

def run(label, Dexpr, p=2147483647, verbose=True):
    t0 = time.time()
    uv = json.load(open('/tmp/hodge/results/polarized-k3-remainder/uv_forms.json'))
    Fh = core.parse(uv['Fhat'], p)
    D = core.parse(Dexpr, p)
    d = core.degree(D)
    chs = {vdeg(e) % 6 for e in D}
    assert len(chs) == 1, 'D must be mu6-homogeneous'
    b = chs.pop()
    q = 1
    degN, degQ = d + 1, 2 * d - 6
    Fi = [pderiv(Fh, i, p) for i in range(NV)]
    Di = [pderiv(D, i, p) for i in range(NV)]
    FD = mul(Fh, D, p)
    Gc = [scale(add(scale(mul(D, Fi[i], p), q, p), mul(Fh, Di[i], p), p), p - 1, p) for i in range(NV)]
    Fq1 = scale(power(Fh, q + 1, p), p - 1, p)
    D2 = mul(D, D, p)
    cols, polys = [], []
    for i in range(NV):
        want = b % 6 if i < 3 else (b + 1) % 6
        for m in monomials(degN):
            if vdeg(m) % 6 != want: continue
            col = shift(Gc[i], m, 1, p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i], p), p)
            cols.append(('N', i, m)); polys.append(col)
    for m in monomials(degQ):
        if vdeg(m) % 6 == (2 * b) % 6:
            cols.append(('Q', -1, m)); polys.append(shift(Fq1, m, 1, p))
    tmons = [m for m in monomials(6) if vdeg(m) in (0, 6)]      # K3 block: pure U-sextics and pure V-sextics
    for m in tmons:
        cols.append(('A', -1, m)); polys.append(shift(D2, m, 1, p))
    rows = {}
    for col in polys:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(polys):
        for e, v in col.items(): flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nul = M.nullspace()
    a_idx = [j for j, c in enumerate(cols) if c[0] == 'A']
    if nul:
        proj = nmod_mat(len(a_idx), nul, [int(X[j, k]) for j in a_idx for k in range(nul)], p)
        Tr, rk = proj.transpose().rref()
        served = []
        for i in range(rk):
            poly = {cols[j][2]: int(Tr[i, jj]) for jj, j in enumerate(a_idx) if int(Tr[i, jj])}
            served.append(poly)
    else:
        rk, served = 0, []
    # class dimension modulo J(f)(U) (+) J(f)(V)
    f = load_f(p); Rf = TernaryJacobianRing(f, p=p, maxdeg=6)
    def classvec(poly):
        pu = {(e[0], e[1], e[2]): v for e, v in poly.items() if vdeg(e) == 0}
        pv = {(e[3], e[4], e[5]): v for e, v in poly.items() if vdeg(e) == 6}
        return Rf.nf(pu, 6) + Rf.nf(pv, 6)
    vecs = [classvec(s) for s in served]
    cdim = nmod_mat(len(vecs), 38, [x for v in vecs for x in v], p).rank() if vecs else 0
    res = {'label': label, 'D': Dexpr, 'd': d, 'char': b, 'prime': p, 'rows': R, 'cols': C, 'nullity': nul,
           'served_K3_monomial_dim': rk, 'served_K3_class_dim': cdim, 'seconds': round(time.time() - t0, 1),
           'note': 'class dim 0 => K3-projection of served(D) vanishes over Q (rigorous); classes live in R_f,6(U)+R_f,6(V) (19+19)'}
    if verbose: print('[servedK3]', json.dumps(res), flush=True)
    return res
if __name__ == '__main__':
    label, Dexpr = sys.argv[1], sys.argv[2]
    p = 2147483647; out = None
    for s in sys.argv[3:]:
        if s.startswith('--out='): out = s.split('=', 1)[1]
        elif s.isdigit(): p = int(s)
        elif s == '--exact': pass
    res = run_exact(label, Dexpr) if '--exact' in sys.argv else run(label, Dexpr, p)
    if out:
        with open(out, 'a') as fh: fh.write(json.dumps(res) + '\n')
