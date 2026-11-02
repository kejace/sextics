"""Wave-2 EXCLUSION/RIGOR library for lane klein-coupled-character (q = 1, A = X0^3 Y0^3).

Conventions (engine v2, hodge.core): polys are dicts {exponent 6-tuple: coeff}; variables (X0,X1,X2,Y0,Y1,Y2);
certificate identity  E = A D^2 - q D sum_i N_i F_i + F (D sum_i d_i N_i - sum_i N_i d_i D) = F^(q+1) Q.
All exact systems use F7 = 7F (integer coefficients); consistency for (7F,A,q,D) <=> for (F,A,q,D)
(N -> N/7, Q -> Q/7^(q+1)).

sigma = diag(zeta,zeta^2,zeta^4,1,1,1), zeta^7 = 1, weights W = (1,2,4,0,0,0); F has character 0, A character 3.

Rigor conventions used throughout:
 * rank_p(M) <= rank_Q(M) for an integer matrix M (reduction mod p can only lower the rank);
 * dim_Fp J(F)_k = dim_Q J(F)_k for every k when the partials of F mod p form a regular sequence
   (checked by comparing dim_Fp J_k with the Hilbert series of the complete intersection);
 * Farkas certificate: exact rational y with y^T M = 0 and y^T b != 0 proves the linear system M x = b has no
   solution over any field of characteristic 0.
"""
import sys, os, json, time, hashlib, math, itertools
from fractions import Fraction
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, symmetry, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from flint import nmod_mat, fmpz_mat, fmpq_mat, fmpz, fmpq

LANE = '/tmp/hodge/results/klein-coupled-character'
LOG = f'{LANE}/wave2-excl.log.jsonl'
PACKET = f'{LANE}/wave2-excl'
W = (1, 2, 4, 0, 0, 0)
PI = (0, 0, 0, 1, 1, 1)
# wave-2 primes (wave-1 used 2147483647; the dossier used 32003 and 101)
P1, P2, P3 = 1000000007, 998244353, 2147483629

def log(**kw):
    kw.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    kw.setdefault('engine', 2)
    kw.setdefault('lens', 'excl')
    with open(LOG, 'a') as f:
        f.write(json.dumps(kw, default=str) + '\n')

def load():
    F, A, q, rec = core.load_instance('klein-coupled-character')
    return F, A, q

def wt(e):
    return (e[0] + 2 * e[1] + 4 * e[2]) % 7

def char_of(P):
    ws = {wt(e) for e in P}
    return ws.pop() if len(ws) == 1 else None

def wmonomials(deg, par=None, k=None):
    return [e for e in monomials(deg, par) if k is None or wt(e) == k]

def integral(D):
    den = 1
    for v in D.values():
        fv = Fraction(v)
        den = den * fv.denominator // math.gcd(den, fv.denominator)
    return {e: int(Fraction(v) * den) for e, v in D.items()}, den

def hilbert_ci(k, n=6, a=5):
    """dim_Q (S/J)_k for a complete intersection of n forms of degree a in n variables: coeff of t^k in ((1-t^a)/(1-t))^n."""
    c = [0] * (k + 1)
    c[0] = 1
    for _ in range(n):
        new = [0] * (k + 1)
        for i in range(k + 1):
            if c[i]:
                for j in range(a):
                    if i + j <= k:
                        new[i + j] += c[i]
        c = new
    return c[k]

def dimS(k):
    return math.comb(k + 5, 5)

# ------------------------------------------------------------------------------------------------
# character/parity blocks of S_k, J_k, R_k mod p

def jac_block(F, deg, p, par=None, k=None):
    """columns m*F_i spanning J_deg in the (par, k) block (mod p)."""
    Fp = core.reduce_mod(F, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    keys, polys = [], []
    for i in range(NV):
        pm = (par + PI[i]) % 2 if par is not None else None
        km = (k + W[i]) % 7 if k is not None else None
        for m in wmonomials(deg - 5, pm, km):
            keys.append((i, m)); polys.append(shift(Fi[i], m, 1, p))
    return keys, polys

def dense(colpolys, rows, p, rhs=None):
    R = len(rows); C = len(colpolys) + (1 if rhs is not None else 0)
    flat = [0] * (R * C)
    for j, col in enumerate(colpolys):
        for e, v in col.items():
            flat[rows[e] * C + j] = v % p
    if rhs is not None:
        for e, v in rhs.items():
            flat[rows[e] * C + len(colpolys)] = v % p
    return nmod_mat(R, C, flat, p)

def e1_kernel_block(F, A, d, p, parD=0, k=None, want_basis=True):
    """{D in S_d^(parD,k) : A D in J_{d+6}} modulo J_d^(parD,k)  (mod p).
    Returns dict(dimS, dimJ_d, dimR_d, dimJ_t, dimR_t (target block), kernel_total, kernel_modJ, reps).
    Method: functionals on S_{d+6} block vanishing on J block; multiplication matrix; nullspace; reduce mod J_d."""
    pA = core.parity_of(A); a = char_of(A)
    Dm = wmonomials(d, parD, k)
    out = {'d': d, 'parD': parD, 'k': k, 'p': p, 'dimS': len(Dm)}
    if not Dm:
        out.update(kernel_total=0, kernel_modJ=0, reps=[]); return out
    tpar = (parD + pA) % 2; tk = (k + a) % 7 if k is not None else None
    tgt = wmonomials(d + 6, tpar, tk)
    rows = {m: i for i, m in enumerate(tgt)}
    jkeys, jpolys = jac_block(F, d + 6, p, tpar, tk)
    MJ = dense(jpolys, rows, p)
    L, nl = MJ.transpose().nullspace()       # functionals vanishing on J block: first nl columns
    rJ = len(tgt) - nl
    out.update(dimS_t=len(tgt), dimJ_t=rJ, dimR_t=nl)
    Ap = core.reduce_mod(A, p)
    AD = dense([shift(Ap, m, 1, p) for m in Dm], rows, p)
    Lt = nmod_mat(nl, len(rows), [int(L[r, c]) for c in range(nl) for r in range(len(rows))], p)
    Mult = Lt * AD
    K, nk = Mult.nullspace()
    kervecs = [[int(K[r, c]) for r in range(len(Dm))] for c in range(nk)]
    jdkeys, jdpolys = jac_block(F, d, p, parD, k)
    Drows = {m: i for i, m in enumerate(Dm)}
    jd = dense(jdpolys, Drows, p).rank() if jdpolys else 0
    out.update(dimJ_d=jd, dimR_d=len(Dm) - jd, kernel_total=nk, kernel_modJ=nk - jd)
    reps = []
    if want_basis and nk > jd:
        nD = len(Dm)
        allcols = ([[int(c) for c in col] for col in zip(*[[int(dense([pp], Drows, p)[r, 0]) for r in range(nD)] for pp in jdpolys])] if False else None)
        # simpler: build matrix [J_d cols | kernel vecs]
        Jm = dense(jdpolys, Drows, p) if jdpolys else None
        C = len(jdpolys) + nk
        flat = [0] * (nD * C)
        if Jm is not None:
            for r in range(nD):
                for j in range(len(jdpolys)):
                    v = int(Jm[r, j])
                    if v: flat[r * C + j] = v
        for c, vec in enumerate(kervecs):
            for r, v in enumerate(vec):
                if v: flat[r * C + len(jdpolys) + c] = v
        Mall = nmod_mat(nD, C, flat, p)
        Mr, rk = Mall.rref()
        piv = []
        for i in range(rk):
            j = 0
            while j < C and int(Mr[i, j]) == 0: j += 1
            piv.append(j)
        nJ = len(jdpolys)
        for j in piv:
            if j >= nJ:
                reps.append({Dm[r]: v for r, v in enumerate(kervecs[j - nJ]) if v})
    out['reps'] = reps
    return out

# ------------------------------------------------------------------------------------------------
# fixed-D test (engine v2, character-restricted) + exact Farkas certificate

def exact_columns(F, A, q, D, cols):
    """exact integer column polys for F7 = 7F and integralized D, same order as `cols`; rhs = -A D^2."""
    F7 = scale(F, 7)
    Dz, _ = integral(D)
    Fi = [pderiv(F7, i) for i in range(NV)]
    Di = [pderiv(Dz, i) for i in range(NV)]
    FD = mul(F7, Dz)
    G = [scale(add(scale(mul(Dz, Fi[i]), q), mul(F7, Di[i])), -1) for i in range(NV)]
    Fq1 = scale(power(F7, q + 1), -1)
    colpolys = []
    for (kind, i, m) in cols:
        if kind == 'N':
            col = shift(G[i], m, 1)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i]))
        else:
            col = shift(Fq1, m, 1)
        colpolys.append({e: int(v) for e, v in col.items()})
    rhs = scale(mul(A, mul(Dz, Dz)), -1)
    rhs = {e: int(v) for e, v in rhs.items()}
    return colpolys, rhs, Dz

_G = None
def group(F):
    global _G
    if _G is None:
        _G = symmetry.diagonal_group(F)
    return _G

def solve_modp(F, A, q, D, p, compressed=False):
    S, sub = symmetry.restricted_system(F, A, q, D, p=p, group=group(F))
    r = fixed_d.solve_compressed(S) if compressed else S.solve()
    out = {'consistent': r['consistent'], 'rank_aug': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows,
           'cols': S.ncols, 'nN': S.nN, 'nQ': S.nQ, 'full_cols': S.full_cols, 'stab': len(sub), 'p': p, 'd': S.d,
           'charD': char_of(D), 'parD': core.parity_of(D)}
    out['rank_M'] = r['rank'] - (0 if r['consistent'] else 1)
    return out, S, r

def farkas(F, A, q, D, p, save=None, Dstr=None, verbose=False, max_store=4_000_000):
    """Exact Farkas certificate for the sector-restricted fixed-D system (status: consistent-mod-p | certified | failed)."""
    t0 = time.time()
    out, S, r = solve_modp(F, A, q, D, p)
    res = dict(out); res['D'] = Dstr or (core.to_str(D) if len(core.to_str(D)) < 2000 else core.to_str(D)[:2000] + '...')
    if out['consistent']:
        res['status'] = 'consistent-mod-p'; res['seconds'] = round(time.time() - t0, 2)
        return res
    n = S.ncols; piv = r['pivots']
    assert piv[-1] == n
    C = piv[:-1]; rk = len(C)
    colpolys, rhs, Dz = exact_columns(F, A, q, D, S.cols)
    rows = {}
    for col in colpolys:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    mrows = len(rows)
    # pivot rows of [M|b] mod p via rref of the transpose
    flatT = [0] * ((n + 1) * mrows)
    for j, col in enumerate(colpolys):
        base = j * mrows
        for e, v in col.items(): flatT[base + rows[e]] = v % p
    base = n * mrows
    for e, v in rhs.items(): flatT[base + rows[e]] = v % p
    MT = nmod_mat(n + 1, mrows, flatT, p)
    MTr, rkT = MT.rref()
    assert rkT == rk + 1, (rkT, rk)
    Srows = []
    for i in range(rkT):
        j = 0
        while j < mrows and int(MTr[i, j]) == 0: j += 1
        Srows.append(j)
    rowpos = {ri: k for k, ri in enumerate(Srows)}
    Bt = [[0] * (rk + 1) for _ in range(rk)]
    for jj, j in enumerate(C):
        for e, v in colpolys[j].items():
            ri = rows[e]
            if ri in rowpos: Bt[jj][rowpos[ri]] = v
    Tp = nmod_mat(rk, rk + 1, [v % p for row in Bt for v in row], p)
    X, nullity = Tp.nullspace()
    assert nullity == 1, nullity
    j0 = next(k for k in range(rk + 1) if int(X[k, 0]) != 0)
    others = [k for k in range(rk + 1) if k != j0]
    Tsq = fmpz_mat(rk, rk, [Bt[i][k] for i in range(rk) for k in others])
    rhsv = fmpz_mat(rk, 1, [-Bt[i][j0] for i in range(rk)])
    t1 = time.time()
    ysol = Tsq.solve(rhsv)
    t2 = time.time()
    y = [None] * (rk + 1); y[j0] = Fraction(1)
    for idx, k in enumerate(others):
        y[k] = Fraction(int(ysol[idx, 0].p), int(ysol[idx, 0].q))
    yrow = {Srows[k]: y[k] for k in range(rk + 1)}
    bad = 0
    for j, col in enumerate(colpolys):
        s = Fraction(0)
        for e, v in col.items():
            ri = rows[e]
            if ri in yrow: s += yrow[ri] * v
        if s != 0: bad += 1
    yb = Fraction(0)
    for e, v in rhs.items():
        ri = rows[e]
        if ri in yrow: yb += yrow[ri] * v
    ok = (bad == 0 and yb != 0)
    inv_rows = {v: k for k, v in rows.items()}
    cert_rows = [list(inv_rows[ri]) for ri in Srows]
    ytxt = [str(v) for v in y]
    h = hashlib.sha256(('\n'.join(ytxt) + '\n' + json.dumps(cert_rows)).encode()).hexdigest()
    res.update({'status': 'certified' if ok else 'failed', 'farkas_rows': rk + 1, 'bad_cols': bad, 'yb': str(yb),
                'sha256_y': h, 'bits_max': max(max(v.numerator.bit_length(), v.denominator.bit_length()) for v in y),
                't_solve': round(t2 - t1, 2), 'seconds': round(time.time() - t0, 2), 'F_scaled': '7F'})
    if save:
        size = sum(len(s) for s in ytxt)
        cert = {'instance': 'klein-coupled-character', 'F': '7F', 'q': q, 'D': core.to_str(Dz), 'p_used_for_pivots': p,
                'columns_note': 'character-sector restricted columns of engine v2 symmetry.restricted_system (lossless by averaging)',
                'row_monomials': cert_rows, 'y': ytxt if size <= max_store else None,
                'y_omitted_size': None if size <= max_store else size, 'sha256_y': h, 'y_dot_b': str(yb),
                'verified_exact': ok, 'rank_mod_p': rk, 'cols': n, 'rows': mrows,
                'pivot_cols': [[c[0], c[1], list(c[2])] for c in (S.cols[j] for j in C)] if size <= max_store else None, 'j0': j0}
        with open(save, 'w') as f: json.dump(cert, f)
    if verbose: print(res, file=sys.stderr)
    return res

def cert_from_file(path, F, A, q):
    c = json.load(open(path))
    D = core.parse(c['D'])
    S, sub = symmetry.restricted_system(F, A, q, D, p=c['p_used_for_pivots'], group=group(F))
    colpolys, rhs, Dz = exact_columns(F, A, q, D, S.cols)
    y = {tuple(r): Fraction(v) for r, v in zip(c['row_monomials'], c['y'])}
    bad = sum(1 for col in colpolys if sum((y.get(e, 0) * v for e, v in col.items()), Fraction(0)) != 0)
    yb = sum((y.get(e, 0) * v for e, v in rhs.items()), Fraction(0))
    return bad == 0 and yb != 0, {'bad': bad, 'yb': str(yb)}

# polynomials used repeatedly
def polars(F):
    return [pderiv(F, i) for i in range(NV)]

def jhat():
    x, y, z, t, r, s = [core.parse(v) for v in ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2')]
    F, A, q = load()
    def P(k): return power
    t6, r6, s6 = power(t, 6), power(r, 6), power(s, 6)
    inner1 = add(add(t6, scale(r6, -1)), s6)                         # t^6 - r^6 + s^6
    inner1 = scale(add(inner1, scale(F, -1)), 9)                     # 9(t^6 - r^6 + s^6) - 9F
    part1 = scale(mul(mul(x, mul(z, z)), inner1), 7)                 # 7 x z^2 [...]
    inner2 = add(add(scale(t6, 2), scale(r6, -4)), scale(s6, 8))     # 2t^6 - 4r^6 + 8s^6
    inner2 = add(inner2, scale(F, -2))
    part2 = scale(mul(mul(y, mul(s, t)), inner2), 15)
    return add(part1, part2)
