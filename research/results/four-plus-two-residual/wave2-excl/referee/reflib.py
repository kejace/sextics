"""Independent referee library (no dependency on the hodge engine or w2lib).
Polys: dict {exponent 6-tuple: Fraction}; variables X0,X1,X2,Y0,Y1,Y2 = 0..5; Y-parity = (e3+e4+e5) mod 2."""
import itertools, math, time, json, re
from fractions import Fraction
from flint import fmpz_mat, nmod_mat
import sympy
NV = 6
CORE, BIN = (0, 1, 3, 4), (2, 5)

def monomials(deg, par=None, idx=tuple(range(NV))):
    out = []
    for c in itertools.combinations_with_replacement(idx, deg):
        e = [0] * NV
        for i in c: e[i] += 1
        if par is None or (e[3] + e[4] + e[5]) % 2 == par:
            out.append(tuple(e))
    return out

def padd(P, Q, s=1):
    R = dict(P)
    for e, c in Q.items():
        v = R.get(e, 0) + s * c
        if v: R[e] = v
        else: R.pop(e, None)
    return R

def pscale(P, s):
    return {e: c * s for e, c in P.items()} if s else {}

def pmul(P, Q):
    R = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            v = R.get(e, 0) + c1 * c2
            if v: R[e] = v
            else: R.pop(e, None)
    return R

def pderiv(P, i):
    R = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            R[tuple(e2)] = c * e[i]
    return R

def pshift(P, m):
    return {tuple(a + b for a, b in zip(e, m)): c for e, c in P.items()}

def parity(P):
    ps = {(e[3] + e[4] + e[5]) % 2 for e in P}
    return ps.pop() if len(ps) == 1 else None

def degree(P):
    return next(iter(P)) and sum(next(iter(P)))

def z(i):
    return {tuple(1 if j == i else 0 for j in range(NV)): Fraction(1)}

def load_instance(path=None):
    import os
    path = path or os.environ.get('HODGE_INSTANCE_REC', '/home/kejace/devel/sextics/data/four-plus-two-residual/instance.rec')
    txt = open(path).read()
    def grab(key):
        m = re.search(r'^' + key + r':\s*(.*?)(?=^\w+:|\Z)', txt, re.S | re.M)
        return m.group(1).replace('\n', ' ')
    X0, X1, X2, Y0, Y1, Y2 = sympy.symbols('X0 X1 X2 Y0 Y1 Y2')
    syms = (X0, X1, X2, Y0, Y1, Y2)
    def todict(expr):
        P = sympy.Poly(sympy.sympify(expr, dict(zip('X0 X1 X2 Y0 Y1 Y2'.split(), syms))), *syms)
        return {tuple(int(x) for x in e): Fraction(int(c.p), int(c.q)) for e, c in zip(P.monoms(), P.coeffs())}
    return todict(grab('F')), todict(grab('A'))

# ---------- matrices ----------
def _rows(cols, extra):
    rows = {}
    for col in list(cols) + (list(extra) if extra else []):
        for e in col: rows.setdefault(e, len(rows))
    return rows

def build_exact(cols, extra=None):
    """fmpz_mat = L * [cols | extra] with ONE global scalar L (nullspace preserved exactly)."""
    allc = list(cols) + (list(extra) if extra else [])
    rows = _rows(cols, extra)
    L = 1
    for col in allc:
        for c in col.values():
            L = L * c.denominator // math.gcd(L, c.denominator)
    R, C = len(rows), len(allc)
    flat = [0] * (R * C)
    for j, col in enumerate(allc):
        for e, c in col.items():
            flat[rows[e] * C + j] = int(c * L)
    return fmpz_mat(R, C, flat), rows, L

def build_modp(cols, p, extra=None):
    allc = list(cols) + (list(extra) if extra else [])
    rows = _rows(cols, extra)
    R, C = len(rows), len(allc)
    flat = [0] * (R * C)
    for j, col in enumerate(allc):
        for e, c in col.items():
            flat[rows[e] * C + j] = (c.numerator * pow(c.denominator, -1, p)) % p
    return nmod_mat(R, C, flat, p), rows

def rank_exact(cols, extra=None):
    if not cols and not extra: return 0
    M, _, _ = build_exact(cols, extra); return int(M.rank())

def rank_modp(cols, p, extra=None):
    if not cols and not extra: return 0
    M, _ = build_modp(cols, p, extra); return int(M.rank())

# ---------- Jacobian data on a variable subset ----------
def jac_cols(P, idx, k, par=None):
    """Spanning set of J(P)_k on variables idx (optional Y-parity of the result)."""
    cols = []
    for i in idx:
        Pi = pderiv(P, i)
        if not Pi or k < 5: continue
        pi = parity(Pi)
        for m in monomials(k - 5, None if par is None else (par - pi) % 2, idx):
            cols.append(pshift(Pi, m))
    return cols

def mult_kernel_dims(P, idx, a, k, par, p=None):
    """dim{D in S_k^par(idx): a D in J(P)_{k+deg a}} and its dimension modulo J(P)_k^par (exact if p None)."""
    da = sum(next(iter(a)))
    msD = monomials(k, par, idx)
    tpar = (par + parity(a)) % 2
    Acols = [pmul(a, {m: Fraction(1)}) for m in msD]
    Jc = jac_cols(P, idx, k + da, tpar)
    rk = rank_exact if p is None else (lambda c: rank_modp(c, p))
    rJ = rk(Jc) if Jc else 0
    rAll = rk(Acols + Jc)
    dim_full = len(msD) + rJ - rAll
    Jk = jac_cols(P, idx, k, par)
    dJk = rk(Jk) if Jk else 0
    img_dim = len(monomials(k + da, tpar, idx)) - rJ  # dim R_{k+da}^{tpar}
    return {'k': k, 'par': par, 'n_mons': len(msD), 'dim_full': dim_full, 'dim_J': dJk, 'dim_modJ': dim_full - dJk,
            'dim_R_k': len(msD) - dJk, 'dim_R_target': img_dim, 'surjective': (len(msD) - dJk) - (dim_full - dJk) == img_dim}

# ---------- fixed-D certificate system (q = 1) ----------
def numerator_parities(D, A):
    pD, pA = parity(D), parity(A)
    return None if pD is None else [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(NV)]

GENS = [(0,0,1,0,0,1), (0,0,3,0,0,0), (0,0,0,3,3,0)]   # h, r2, s as weight vectors mod 6 on X0,X1,X2,Y0,Y1,Y2
def weight(e): return tuple(sum(a*b for a, b in zip(e, g)) % 6 for g in GENS)
def poly_weight(P):
    ws = {weight(e) for e in P}
    assert len(ws) == 1, 'not a G-eigenvector'
    return ws.pop()
def sector_filter(F, A, D, q=1):
    """allowed(kind, i, m): monomial m is in the character sector forced for N_i / Q (G = <h, r2, s>, order 24)."""
    cF, cA, cD = poly_weight(F), poly_weight(A), poly_weight(D)
    zw = [weight(tuple(1 if j == i else 0 for j in range(NV))) for i in range(NV)]
    reqN = [tuple((cA[g] + cD[g] - cF[g] + zw[i][g]) % 6 for g in range(3)) for i in range(NV)]
    reqQ = tuple((cA[g] + 2 * cD[g] - (q + 1) * cF[g]) % 6 for g in range(3))
    def allowed(kind, i, m):
        return weight(m) == (reqN[i] if kind == 'N' else reqQ)
    return allowed

def fixed_D_system(F, A, D, parity_free=False, allowed=None):
    """Columns of  -D N(F) + F (D div N - N(D)) - F^2 Q  (unknown coefficients), rhs = -A D^2 (q = 1)."""
    d = sum(next(iter(D)))
    pars = None if parity_free else numerator_parities(D, A)
    Fi = [pderiv(F, i) for i in range(NV)]
    Di = [pderiv(D, i) for i in range(NV)]
    FD = pmul(F, D)
    G = [pscale(padd(pmul(D, Fi[i]), pmul(F, Di[i])), Fraction(-1)) for i in range(NV)]
    F2 = pscale(pmul(F, F), Fraction(-1))
    cols, info = [], []
    for i in range(NV):
        for m in monomials(d + 1, None if parity_free else pars[i]):
            if allowed is not None and not allowed('N', i, m): continue
            col = pshift(G[i], m)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = padd(col, pscale(pshift(FD, tuple(m2)), Fraction(m[i])))
            cols.append(col); info.append(('N', i, m))
    nN = len(cols)
    for m in monomials(2 * d - 6, None if parity_free else parity(A)):
        if allowed is not None and not allowed('Q', -1, m): continue
        cols.append(pshift(F2, m)); info.append(('Q', m))
    rhs = pscale(pmul(A, pmul(D, D)), Fraction(-1))
    return cols, info, rhs, nN, len(cols) - nN

def koszul_nullity_exact(F, A, D, parity_free=False, allowed=None):
    """Exact nullity of the fixed-D matrix for d <= 8: n_m + dim{(k,Q): D div K - K(D) = F Q}, K = sum k_ij (F_j d_i - F_i d_j)."""
    d = sum(next(iter(D)))
    assert d <= 8
    pars = None if parity_free else numerator_parities(D, A)
    z0 = tuple(1 if j == 0 else 0 for j in range(NV))
    n_m = len([m for m in monomials(d, None if parity_free else pars[0]) if allowed is None or allowed('N', 0, tuple(a + b for a, b in zip(m, z0)))])
    Fi = [pderiv(F, i) for i in range(NV)]
    Di = [pderiv(D, i) for i in range(NV)]
    cols = []
    for i in range(NV):
        for j in range(i + 1, NV):
            park = None if parity_free else (pars[i] - (0 if j < 3 else 1)) % 2
            for m in monomials(d - 4, park):
                if allowed is not None and not allowed('N', i, tuple(a + b for a, b in zip(m, next(iter(Fi[j]))))): continue
                mp = {m: Fraction(1)}
                t1 = pmul(D, padd(pmul(pderiv(mp, i), Fi[j]), pmul(pderiv(mp, j), Fi[i]), -1))
                t2 = pmul(mp, padd(pmul(Fi[j], Di[i]), pmul(Fi[i], Di[j]), -1))
                cols.append(padd(t1, t2, -1))
    n_k = len(cols)
    for m in monomials(2 * d - 6, None if parity_free else parity(A)):
        if allowed is not None and not allowed('Q', -1, m): continue
        cols.append(pscale(pshift(F, m), Fraction(-1)))
    r = rank_exact(cols)
    return {'n_m': n_m, 'n_k': n_k, 'n_Q': len(cols) - n_k, 'koszul_rank': r, 'koszul_nullity': len(cols) - r, 'nullity_total': n_m + len(cols) - r}

def fixed_D_test(F, A, D, p, parity_free=False, exact_nullity=True, allowed=None):
    t0 = time.time()
    cols, info, rhs, nN, nQ = fixed_D_system(F, A, D, parity_free, allowed=allowed)
    M, rows = build_modp(cols, p, extra=[rhs]); rMb = int(M.rank())
    M0, _ = build_modp(cols, p); rM = int(M0.rank())
    out = {'rows': len(rows), 'd': sum(next(iter(D))), 'p': p, 'cols': len(cols), 'nN': nN, 'nQ': nQ, 'rank_M_p': rM, 'rank_Mb_p': rMb, 'nullity_p': len(cols) - rM,
           'consistent_p': rM == rMb}
    if exact_nullity and not out['consistent_p'] and out['d'] <= 8:
        kz = koszul_nullity_exact(F, A, D, parity_free, allowed=allowed)
        out['nullity_exact_Q'] = kz['nullity_total']; out['koszul'] = kz
        out['verdict'] = 'inconsistent over C [RANK]' if kz['nullity_total'] == out['nullity_p'] else 'inconsistent mod p only (nullity mismatch)'
    elif out['consistent_p']:
        out['verdict'] = 'CONSISTENT mod p'
    else:
        out['verdict'] = 'inconsistent mod p [MODP]'
    out['seconds'] = round(time.time() - t0, 1)
    return out

def log(path, ent):
    ent = dict(ent); ent.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(path, 'a') as f: f.write(json.dumps(ent, default=str) + '\n')
    print(json.dumps(ent, default=str)[:600], flush=True)
