"""wave-2 EXCLUSION lens, lane klein-coupled-character (q = 1, A = X0^3 Y0^3).

Shared helpers (all exact mod p unless stated; python-flint nmod_mat):
  * sector bookkeeping for the Z/7 x Z/2 grading (sigma weights (1,2,4,0,0,0), Y-parity),
  * colon_sector: dim (J:A)_d, dim J_d, kernel representatives per sector (E1-admissible denominators),
  * hilbert_R_sector: combinatorial Hilbert function of S/J (regular sequence => Koszul exact),
  * hess_system: the E2-mod-J ("Hessian") relaxation for D = grad(F).w, gauge-fixed, as quadrics in the
    w-coefficients, written in msolve format,
  * irrelevance check of a leading ideal / Macaulay-matrix rank,
  * exact rational dual certificates of inconsistency for fixed-D systems (engine RestrictedFixedDSystem).
Derivation of the relaxation (q = 1): E1 N(F) = A D with D = sum_i w_i F_i  =>  N = A w + K(p) (Koszul);
E2: B_D(N) = D div N - N(D) in (F) subset J.  B_D(K(p)) in J for every D (p antisymmetric, Hess symmetric);
B_D(A w) = A w^T Hess(F) w  mod J.  So E2 implies  w^T Hess(F) w in (J : A)_{2d-6}.  The gauge w -> w + l z
(D -> D + 6 l F) preserves solvability (book ch.1 denominator gauge), so w may be gauge-fixed.
"""
import sys, json, time, os, hashlib, itertools, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, jacobian, symmetry
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from flint import nmod_mat, fmpq_mat, fmpz_mat, fmpq, fmpz
from fractions import Fraction

LANE = '/tmp/hodge/results/klein-coupled-character'
PACKET = LANE + '/wave2-excl'
LOG = LANE + '/wave2-excl.log.jsonl'
W = (1, 2, 4, 0, 0, 0)
PI = (0, 0, 0, 1, 1, 1)
PRIMES_NEW = [1000003, 1073741827, 998244353, 1000000007, 2147483647]

def wt(e):
    return (e[0] + 2 * e[1] + 4 * e[2]) % 7

def char_of(P):
    ws = {wt(e) for e in P}
    return ws.pop() if len(ws) == 1 else None

def par_of(e):
    return (e[3] + e[4] + e[5]) % 2

def wmonomials(deg, par=None, k=None):
    return [e for e in monomials(deg, par) if k is None or wt(e) == k]

def log(entry, path=LOG):
    entry = dict(entry)
    entry.setdefault('time', time.strftime('%Y-%m-%dT%H:%M:%S'))
    entry.setdefault('lens', 'excl')
    with open(path, 'a') as f:
        f.write(json.dumps(entry, default=str) + '\n')

def load():
    F, A, q, rec = core.load_instance('klein-coupled-character')
    return F, A, q

def sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]

# ---------------------------------------------------------------- dense helpers
def dense_cols(colpolys, rowidx, p, rhs=None):
    R, C = len(rowidx), len(colpolys)
    Wd = C + (1 if rhs is not None else 0)
    flat = [0] * (R * Wd)
    for j, col in enumerate(colpolys):
        for e, v in col.items():
            flat[rowidx[e] * Wd + j] = v % p
    if rhs is not None:
        for e, v in rhs.items():
            flat[rowidx[e] * Wd + C] = v % p
    return nmod_mat(R, Wd, flat, p)

def pivots_of(Mr, rank, ncols):
    piv = []
    for i in range(rank):
        j = 0
        while j < ncols and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
    return piv

def rows_of_nmod(M):
    R, C = M.nrows(), M.ncols()
    return [[int(M[i, j]) for j in range(C)] for i in range(R)]

# ---------------------------------------------------------------- Hilbert function of R = S/J (combinatorial)
def hilbert_R_sector(d, par, k):
    """dim (S/J)_d in the (parity, character) sector, from the Koszul complex of the regular sequence
    F_0..F_5 (deg 5, parity PI[i], character -W[i])."""
    tot = 0
    for r in range(7):
        for Ssub in itertools.combinations(range(NV), r):
            dd = d - 5 * r
            if dd < 0:
                continue
            sp = sum(PI[i] for i in Ssub) % 2
            sk = sum(W[i] for i in Ssub) % 7
            cnt = 0
            for m in monomials(dd):
                if (par_of(m) + sp) % 2 == par and (wt(m) - sk) % 7 == k:
                    cnt += 1
            tot += (-1) ** r * cnt
    return tot

# ---------------------------------------------------------------- Jacobian / colon ideal per sector
class Ctx:
    def __init__(self, F, A, p):
        self.p = p
        self.F = core.reduce_mod(F, p); self.A = core.reduce_mod(A, p)
        self.Fi = [pderiv(self.F, i, p) for i in range(NV)]
        self.H = [[pderiv(self.Fi[i], j, p) for j in range(NV)] for i in range(NV)]
        self.pA = core.parity_of(A); self.kA = char_of(A)

def jac_cols(ctx, deg, par, k):
    """generators m*F_i of J_deg in the (par, k) sector. char(m F_i) = char(m) - W[i]; parity(m F_i) = par(m) + PI[i]."""
    keys, polys = [], []
    for i in range(NV):
        for m in wmonomials(deg - 5, None if par is None else (par + PI[i]) % 2, None if k is None else (k + W[i]) % 7):
            keys.append((i, m)); polys.append(shift(ctx.Fi[i], m, 1, ctx.p))
    return keys, polys

def colon_sector(ctx, d, par, k, want_reps=True):
    """(J : A)_d in the (par, k) sector, mod p.
    Returns dict: dimS, dimJ, dimColon, kernel (= dimColon - dimJ), reps (kernel representatives reduced mod J_d),
    colon_basis (rref rows over Dm), Dm (monomial list), jkeys/jpolys (J_d generators)."""
    p = ctx.p
    Dm = wmonomials(d, par, k)
    out = {'d': d, 'par': par, 'k': k, 'p': p, 'dimS': len(Dm)}
    if not Dm:
        out.update({'dimJ': 0, 'dimColon': 0, 'kernel': 0, 'reps': [], 'colon_basis': [], 'Dm': Dm, 'jkeys': [], 'jpolys': []})
        return out
    tpar = (par + ctx.pA) % 2; tk = (k + ctx.kA) % 7
    T = wmonomials(d + 6, tpar, tk)
    ridx = {m: i for i, m in enumerate(T)}
    jk, jp = jac_cols(ctx, d + 6, tpar, tk)
    cols = [shift(ctx.A, m, 1, p) for m in Dm] + jp
    M = dense_cols(cols, ridx, p)
    X, nul = M.nullspace()
    nD = len(Dm)
    if nul == 0:
        colon_basis = []
    else:
        proj = nmod_mat(nul, nD, [int(X[j, c]) for c in range(nul) for j in range(nD)], p)
        Pr, rk = proj.rref()
        colon_basis = [[int(Pr[i, j]) for j in range(nD)] for i in range(rk)]
    out['dimColon'] = len(colon_basis)
    jdk, jdp = jac_cols(ctx, d, par, k)
    didx = {m: i for i, m in enumerate(Dm)}
    if jdp:
        JM = dense_cols(jdp, didx, p)
        dimJ = JM.rank()
    else:
        dimJ = 0
    out['dimJ'] = dimJ
    out['kernel'] = out['dimColon'] - dimJ
    out['Dm'] = Dm; out['jkeys'] = jdk; out['jpolys'] = jdp; out['colon_basis'] = colon_basis
    reps = []
    if want_reps and out['kernel'] > 0:
        # reduce colon basis vectors modulo J_d: stack [J_d gens | colon vecs] as columns, rref, keep colon pivots
        allcols = []
        for P in jdp:
            allcols.append([P.get(m, 0) % p for m in Dm])
        allcols += colon_basis
        C = len(allcols)
        flat = [0] * (nD * C)
        for j, col in enumerate(allcols):
            for r, v in enumerate(col):
                if v: flat[r * C + j] = v
        Mall = nmod_mat(nD, C, flat, p)
        Mr, rk = Mall.rref()
        piv = pivots_of(Mr, rk, C)
        nJ = len(jdp)
        for j in piv:
            if j >= nJ:
                v = colon_basis[j - nJ]
                reps.append({Dm[r]: v[r] for r in range(nD) if v[r]})
        assert len(reps) == out['kernel'], (len(reps), out['kernel'])
    out['reps'] = reps
    return out

# ---------------------------------------------------------------- normal form modulo a subspace (rref rows over a monomial list)
class SubspaceReducer:
    """Given rref rows (over monomial list mons) spanning a subspace U, reduce vectors modulo U; coordinates on
    the standard (non-pivot) monomials give S/U."""
    def __init__(self, rows, mons, p):
        self.p = p; self.mons = mons; self.idx = {m: i for i, m in enumerate(mons)}
        n = len(mons)
        if rows:
            B = nmod_mat(len(rows), n, [v for row in rows for v in row], p)
            Br, rk = B.rref()
            self.piv = pivots_of(Br, rk, n)
            self.B = nmod_mat(rk, n, [int(Br[i, j]) for i in range(rk) for j in range(n)], p) if rk else None
        else:
            self.piv = []; self.B = None
        pset = set(self.piv)
        self.std = [j for j in range(n) if j not in pset]
        self.rank = len(self.piv)

    def reduce_matrix(self, V):
        """V: nmod_mat (nvec x n). returns nmod_mat (nvec x len(std)) of normal-form coordinates."""
        p = self.p; n = len(self.mons)
        if self.B is not None and self.piv:
            # V - V[:, piv] * B
            nv = V.nrows()
            Vp = nmod_mat(nv, len(self.piv), [int(V[i, j]) for i in range(nv) for j in self.piv], p)
            V = V - Vp * self.B
        nv = V.nrows()
        return nmod_mat(nv, len(self.std), [int(V[i, j]) for i in range(nv) for j in self.std], p)

def colon_reducer(ctx, deg, par, k, use_colon=True):
    """Reducer modulo (J:A)_deg (or modulo J_deg if use_colon=False) in the (par,k) sector."""
    cs = colon_sector(ctx, deg, par, k, want_reps=False)
    Dm = cs['Dm']
    if use_colon:
        rows = cs['colon_basis']
    else:
        rows = [[P.get(m, 0) % ctx.p for m in Dm] for P in cs['jpolys']]
    red = SubspaceReducer(rows, Dm, ctx.p)
    red.info = {'dimS': cs['dimS'], 'dimJ': cs['dimJ'], 'dimColon': cs['dimColon'], 'kernel': cs['kernel']}
    return red

# ---------------------------------------------------------------- Hessian relaxation system
def hess_system(ctx, d, Dpars, Dchars, gauge=True, use_colon=True, verbose=False):
    """Quadratic system in the w-coefficients for D = sum w_i F_i with D parity in Dpars, character in Dchars.
    Unknown (i, m): m of degree d-5 with (par(m)+PI[i]) in Dpars and (wt(m)-W[i]) in Dchars.
    Gauge: remove (0, m'*z0) for m' of degree d-6 with parity in Dpars, char in Dchars.
    Equations: coordinates of  w^T Hess(F) w  modulo (J:A)_{2d-6} (use_colon) or modulo J_{2d-6}.
    Returns dict with 'vars' (list of (i,m)), 'eqs' (list of dict {(a,b): coeff} with a<=b), sector bookkeeping."""
    p = ctx.p
    Dpars = set(Dpars); Dchars = set(Dchars)
    unknowns = []
    for i in range(NV):
        for m in monomials(d - 5):
            if (par_of(m) + PI[i]) % 2 in Dpars and (wt(m) - W[i]) % 7 in Dchars:
                unknowns.append((i, m))
    removed = set()
    if gauge and d >= 6:
        for mp in monomials(d - 6):
            if par_of(mp) in Dpars and wt(mp) in Dchars:
                m0 = list(mp); m0[0] += 1
                removed.add((0, tuple(m0)))
    vars_ = [u for u in unknowns if u not in removed]
    assert len(removed) == len(unknowns) - len(vars_)
    n = len(vars_)
    deg2 = 2 * d - 6
    # target sectors: parity p1+p2, char k1+k2 where (p_a, k_a) is the D-sector of unknown a
    def dsec(u):
        i, m = u
        return ((par_of(m) + PI[i]) % 2, (wt(m) - W[i]) % 7)
    secs = sorted({((dsec(a)[0] + dsec(b)[0]) % 2, (dsec(a)[1] + dsec(b)[1]) % 7) for a in vars_ for b in vars_})
    reducers = {}
    for (tp, tk) in secs:
        reducers[(tp, tk)] = colon_reducer(ctx, deg2, tp, tk, use_colon=use_colon)
        if verbose:
            print(f'   target sector par={tp} k={tk}: {reducers[(tp,tk)].info} -> {len(reducers[(tp,tk)].std)} equations', file=sys.stderr, flush=True)
    # products grouped by target sector
    groups = {}
    for a in range(n):
        ia, ma = vars_[a]
        for b in range(a, n):
            ib, mb = vars_[b]
            sa, sb = dsec(vars_[a]), dsec(vars_[b])
            key = ((sa[0] + sb[0]) % 2, (sa[1] + sb[1]) % 7)
            groups.setdefault(key, []).append((a, b))
    eqs = []
    nprod = 0
    for key, pairs in groups.items():
        red = reducers[key]
        if not red.std:
            continue
        idx = red.idx; nm = len(red.mons)
        flat = [0] * (len(pairs) * nm)
        for r, (a, b) in enumerate(pairs):
            ia, ma = vars_[a]; ib, mb = vars_[b]
            P = ctx.H[ia][ib]
            if not P:
                continue
            mm = tuple(x + y for x, y in zip(ma, mb))
            c = 1 if a == b else 2
            base = r * nm
            for e, v in P.items():
                ee = tuple(x + y for x, y in zip(e, mm))
                flat[base + idx[ee]] = (flat[base + idx[ee]] + c * v) % p
        V = nmod_mat(len(pairs), nm, flat, p)
        NF = red.reduce_matrix(V)
        nprod += len(pairs)
        rows = rows_of_nmod(NF.transpose())  # each row = one equation over the pairs
        for row in rows:
            eq = {pairs[r]: v for r, v in enumerate(row) if v}
            if eq:
                eqs.append(eq)
    out = {'d': d, 'p': p, 'Dpars': sorted(Dpars), 'Dchars': sorted(Dchars), 'gauge': gauge, 'use_colon': use_colon,
           'vars': vars_, 'n': n, 'n_unknowns_before_gauge': len(unknowns), 'n_gauge': len(removed),
           'eqs': eqs, 'n_eqs': len(eqs), 'n_products': nprod,
           'target_info': {f'{tp},{tk}': reducers[(tp, tk)].info | {'n_eq': len(reducers[(tp, tk)].std)} for (tp, tk) in secs}}
    return out

def hess_check_gauge(ctx, d, Dpars, Dchars):
    """Sanity: the Euler/gauge directions w = l z must satisfy the un-gauged relaxation (w^T H w = 30 l^2 F in J)."""
    sysg = hess_system(ctx, d, Dpars, Dchars, gauge=False)
    vidx = {u: a for a, u in enumerate(sysg['vars'])}
    bad = 0; tested = 0
    for mp in monomials(d - 6):
        if par_of(mp) in set(Dpars) and wt(mp) in set(Dchars):
            w = {}
            for i in range(NV):
                m = list(mp); m[i] += 1; m = tuple(m)
                if (i, m) in vidx: w[vidx[(i, m)]] = 1
            for eq in sysg['eqs']:
                s = sum(c * w.get(a, 0) * w.get(b, 0) for (a, b), c in eq.items()) % ctx.p
                if s: bad += 1
            tested += 1
    return tested, bad

def write_msolve(sysd, path):
    p = sysd['p']; n = sysd['n']
    names = [f'u{a}' for a in range(n)]
    lines = [','.join(names), str(p)]
    polys = []
    for eq in sysd['eqs']:
        terms = []
        for (a, b), c in sorted(eq.items()):
            terms.append(f'{c}*u{a}*u{b}' if a != b else f'{c}*u{a}^2')
        polys.append('+'.join(terms))
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n' + ',\n'.join(polys) + '\n')
    return path

def macaulay_rank(sysd, deg=3):
    """rank mod p of the degree-`deg` Macaulay matrix of the quadrics (rows: eq * monomial of degree deg-2).
    Full column rank (= #monomials of degree deg in n vars) proves the ideal contains all forms of that degree,
    i.e. the projective variety is empty (mod p, hence over Qbar by the rank inequality)."""
    p = sysd['p']; n = sysd['n']
    mons = [c for c in itertools.combinations_with_replacement(range(n), deg)]
    midx = {m: i for i, m in enumerate(mons)}
    mults = list(itertools.combinations_with_replacement(range(n), deg - 2))
    rows = []
    for eq in sysd['eqs']:
        for mu in mults:
            row = {}
            for (a, b), c in eq.items():
                key = tuple(sorted((a, b) + mu))
                row[midx[key]] = (row.get(midx[key], 0) + c) % p
            rows.append(row)
    R, C = len(rows), len(mons)
    flat = [0] * (R * C)
    for r, row in enumerate(rows):
        for j, v in row.items():
            flat[r * C + j] = v
    M = nmod_mat(R, C, flat, p)
    rk = M.rank()
    return {'deg': deg, 'rows': R, 'cols': C, 'rank': rk, 'full_column_rank': rk == C}

def parse_msolve_leading(path):
    """parse msolve -g 1 output: list of leading monomials as exponent dicts {var: e}."""
    txt = open(path).read()
    body = txt[txt.index('[') + 1: txt.rindex(']')]
    mons = []
    for tok in body.replace('\n', '').split(','):
        tok = tok.strip()
        if not tok:
            continue
        if tok == '1':
            mons.append({})
            continue
        m = {}
        for fac in tok.split('*'):
            if '^' in fac:
                v, e = fac.split('^'); m[v] = int(e)
            else:
                m[fac] = 1
        mons.append(m)
    return mons

def leading_ideal_irrelevant(mons, nvars):
    """True iff every variable u0..u{n-1} has a pure power among the leading monomials (=> quotient finite-dim
    => projective variety empty), or 1 is in the ideal."""
    if any(len(m) == 0 for m in mons):
        return True, 'unit ideal'
    pure = {}
    for m in mons:
        if len(m) == 1:
            v, e = next(iter(m.items()))
            pure[v] = min(pure.get(v, 10**9), e)
    missing = [f'u{a}' for a in range(nvars) if f'u{a}' not in pure]
    return (not missing), {'pure_powers': pure, 'missing': missing}

def standard_monomial_count(mons, nvars, maxdeg=60):
    """number of standard monomials (finite if irrelevant) by brute force on bounded exponents."""
    pure = {}
    for m in mons:
        if len(m) == 1:
            v, e = next(iter(m.items())); pure[v] = min(pure.get(v, 10**9), e)
    bounds = [pure.get(f'u{a}', None) for a in range(nvars)]
    if any(b is None for b in bounds):
        return None
    gens = [tuple(m.get(f'u{a}', 0) for a in range(nvars)) for m in mons]
    cnt = 0
    for e in itertools.product(*[range(b) for b in bounds]):
        if not any(all(e[a] >= g[a] for a in range(nvars)) for g in gens):
            cnt += 1
    return cnt

# ---------------------------------------------------------------- fixed-D systems (engine, character-restricted) + exact duals
def restricted_system(F, A, q, D, p, group=None):
    group = group or symmetry.diagonal_group(F)
    S, sub = symmetry.restricted_system(F, A, q, D, p, group)
    return S, sub, group

def exact_restricted_columns(F, A, q, D, S):
    """Exact (Fraction) column polynomials and rhs for the same unknown set as the mod-p restricted system S."""
    d = core.degree(D)
    Fi = [pderiv(F, i) for i in range(NV)]
    Di = [pderiv(D, i) for i in range(NV)]
    FD = mul(F, D)
    G = [scale(add(scale(mul(D, Fi[i]), q), mul(F, Di[i])), -1) for i in range(NV)]
    Fq1 = scale(power(F, q + 1), -1)
    cols = []
    for (kind, i, m) in S.cols:
        if kind == 'N':
            col = shift(G[i], m, 1)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i]))
        else:
            col = shift(Fq1, m, 1)
        cols.append(col)
    rhs = scale(mul(A, mul(D, D)), -1)
    return cols, rhs

def exact_dual(F, A, q, D, p, verbose=False, group=None):
    """For an inconsistent restricted system, produce an exact rational vector y (on a small row set) with
    y^T M = 0 (all columns) and y^T b = 1, verified exactly over Q.  Returns dict (ok, support, hash, ...)."""
    t0 = time.time()
    S, sub, group = restricted_system(F, A, q, D, p, group)
    Aug = S.dense_augmented()
    R, C = S.nrows, S.ncols
    Mr, rk = Aug.rref()
    piv = pivots_of(Mr, rk, C + 1)
    if not piv or piv[-1] != C:
        return {'ok': False, 'reason': 'system consistent mod p', 'rank': rk, 'rows': R, 'cols': C, 'p': p}
    Jcols = piv  # r pivot columns of M + the b column
    r1 = len(Jcols)
    # choose rows: pivots of the transpose of Aug[:, Jcols]
    sub_flat = [int(Aug[i, j]) for j in Jcols for i in range(R)]
    T = nmod_mat(r1, R, sub_flat, p)
    Tr, rk2 = T.rref()
    Irows = pivots_of(Tr, rk2, R)
    assert rk2 == r1, (rk2, r1)
    # exact square system: B = Aug[Irows, Jcols] (r1 x r1); solve B^T y = e_last
    cols, rhs = exact_restricted_columns(F, A, q, D, S)
    rowmon = [None] * R
    for e, i in S.rows.items():
        rowmon[i] = e
    def entry(i, j):
        e = rowmon[i]
        if j == C:
            return rhs.get(e, Fraction(0))
        return cols[j].get(e, Fraction(0))
    Bq = fmpq_mat(r1, r1)
    for a, i in enumerate(Irows):
        for b, j in enumerate(Jcols):
            v = Fraction(entry(i, j))
            if v:
                Bq[b, a] = fmpq(v.numerator, v.denominator)   # transpose: B^T[b,a] = B[a,b]
    e = fmpq_mat(r1, 1)
    e[r1 - 1, 0] = fmpq(1)
    t1 = time.time()
    y = Bq.solve(e)
    t2 = time.time()
    yv = [Fraction(int(y[a, 0].p), int(y[a, 0].q)) for a in range(r1)]
    # exact verification over ALL columns
    bad = 0; maxbits = 0
    for j in range(C + 1):
        s = Fraction(0)
        for a, i in enumerate(Irows):
            v = entry(i, j)
            if v:
                s += yv[a] * v
        if (j < C and s != 0) or (j == C and s != 1):
            bad += 1
    for v in yv:
        maxbits = max(maxbits, v.numerator.bit_length(), v.denominator.bit_length())
    t3 = time.time()
    yrep = {str(rowmon[i]): f'{yv[a].numerator}/{yv[a].denominator}' for a, i in enumerate(Irows)}
    h = hashlib.sha256(json.dumps(yrep, sort_keys=True).encode()).hexdigest()
    return {'ok': bad == 0, 'bad_columns': bad, 'support': r1, 'rows': R, 'cols': C, 'rank_mod_p': rk - 1, 'p': p,
            'max_bits': maxbits, 'subgroup': len(sub), 'y': yrep, 'sha256': h,
            't_modp': round(t1 - t0, 1), 't_solve': round(t2 - t1, 1), 't_verify': round(t3 - t2, 1)}

def test_D_restricted(F, A, q, D, p, name, family, note='', compressed=None, group=None, want_log=True, extra=None):
    """E1 (A D in J mod p) then the full fixed-D system restricted to the (parity, character) sector. Logs."""
    t0 = time.time()
    d = core.degree(D); k = char_of(D); par = core.parity_of(D)
    e1, info = fixed_d.e1_solvable(F, A, q, D, p)
    entry = {'exp': 'fixedD', 'D': name if len(name) < 300 else name[:300] + '...', 'd': d, 'par': par, 'k': k, 'prime': p,
             'family': family, 'e1': e1, 'engine': 'v2-restricted', 'result': None, 'rank': None, 'nullity': None, 'note': note}
    if extra: entry.update(extra)
    if e1:
        S, sub, group = restricted_system(F, A, q, D, p, group)
        if compressed is None:
            compressed = S.nrows * (S.ncols + 1) > 6e7
        r = fixed_d.solve_compressed(S) if compressed else S.solve()
        entry.update({'result': 'consistent' if r['consistent'] else 'inconsistent', 'rank': r['rank'], 'nullity': r['nullity'],
                      'rows': S.nrows, 'cols': S.ncols, 'compressed': bool(compressed), 'subgroup': len(sub)})
        entry['_S'] = S; entry['_r'] = r
    else:
        entry['result'] = 'E1-fail'
        entry.update({'e1_rank': info['rank'], 'e1_rank_aug': info['rank_aug']})
    entry['seconds'] = round(time.time() - t0, 2)
    if want_log:
        log({kk: v for kk, v in entry.items() if not kk.startswith('_')})
    return entry

def rand_combo(polys, p, rng):
    out = {}
    for P in polys:
        c = rng.randrange(1, p)
        out = add(out, scale(P, c, p), p)
    return out

def write_singular(sysd, path, extra=''):
    """Singular script: ideal I of the relaxation quadrics over ZZ/p, prints dim/vdim of std(I) (affine cone)."""
    p = sysd['p']; n = sysd['n']
    names = ','.join(f'u{a}' for a in range(n))
    polys = []
    for eq in sysd['eqs']:
        terms = [(f'{c}*u{a}*u{b}' if a != b else f'{c}*u{a}^2') for (a, b), c in sorted(eq.items())]
        polys.append('+'.join(terms))
    body = f'ring R = {p}, ({names}), dp;\nideal I = ' + ',\n'.join(polys) + ';\n'
    body += 'ideal SI = std(I);\n"n_gens", size(I), "gb_size", size(SI), "dim", dim(SI), "vdim", vdim(SI);\n' + extra + 'quit;\n'
    open(path, 'w').write(body)
    return path

# ---------------------------------------------------------------- exact (rational) relaxation system for small d
class ExactReducer:
    def __init__(self, rows, mons):
        self.mons = mons; self.idx = {m: i for i, m in enumerate(mons)}
        n = len(mons)
        if rows:
            B = fmpq_mat(len(rows), n, [fmpq(Fraction(v).numerator, Fraction(v).denominator) for row in rows for v in row])
            Br, rk = B.rref()
            piv = []
            for i in range(rk):
                j = 0
                while j < n and Br[i, j] == 0: j += 1
                piv.append(j)
            self.piv = piv; self.rank = rk
            self.B = [[Fraction(int(Br[i, j].p), int(Br[i, j].q)) for j in range(n)] for i in range(rk)]
        else:
            self.piv = []; self.rank = 0; self.B = []
        ps = set(self.piv); self.std = [j for j in range(n) if j not in ps]
    def reduce(self, vec):
        v = list(vec)
        for row, j in zip(self.B, self.piv):
            c = v[j]
            if c:
                v = [a - c * b for a, b in zip(v, row)]
        return [v[j] for j in self.std]

def colon_reducer_exact(F, A, deg, par, k, use_colon=False):
    """exact reducer modulo J_deg (use_colon=False) in the (par,k) sector; use_colon=True not implemented (needs exact colon)."""
    assert not use_colon
    Fi = [pderiv(F, i) for i in range(NV)]
    Dm = wmonomials(deg, par, k)
    rows = []
    for i in range(NV):
        for m in wmonomials(deg - 5, (par + PI[i]) % 2, (k + W[i]) % 7):
            P_ = shift(Fi[i], m, 1); rows.append([Fraction(P_.get(mm, 0)) for mm in Dm])
    return ExactReducer(rows, Dm)

def hess_system_exact(F, A, d, Dpars, Dchars, gauge=True):
    """Exact rational relaxation system modulo J_{2d-6} (valid when (J:A)_{2d-6} = J_{2d-6}, i.e. 2d-6 <= 8)."""
    assert 2 * d - 6 <= 8
    Fi = [pderiv(F, i) for i in range(NV)]
    H = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
    Dpars = set(Dpars); Dchars = set(Dchars)
    unknowns = [(i, m) for i in range(NV) for m in monomials(d - 5) if (par_of(m) + PI[i]) % 2 in Dpars and (wt(m) - W[i]) % 7 in Dchars]
    removed = set()
    if gauge and d >= 6:
        for mp in monomials(d - 6):
            if par_of(mp) in Dpars and wt(mp) in Dchars:
                m0 = list(mp); m0[0] += 1; removed.add((0, tuple(m0)))
    vars_ = [u for u in unknowns if u not in removed]
    n = len(vars_); deg2 = 2 * d - 6
    def dsec(u):
        i, m = u; return ((par_of(m) + PI[i]) % 2, (wt(m) - W[i]) % 7)
    groups = {}
    for a in range(n):
        for b in range(a, n):
            sa, sb = dsec(vars_[a]), dsec(vars_[b])
            groups.setdefault(((sa[0] + sb[0]) % 2, (sa[1] + sb[1]) % 7), []).append((a, b))
    eqs = []
    for key, pairs in groups.items():
        red = colon_reducer_exact(F, A, deg2, key[0], key[1])
        if not red.std: continue
        cols = {}
        for (a, b) in pairs:
            ia, ma = vars_[a]; ib, mb = vars_[b]
            P_ = H[ia][ib]
            if not P_: continue
            mm = tuple(x + y for x, y in zip(ma, mb)); c = 1 if a == b else 2
            vec = [Fraction(0)] * len(red.mons)
            for e, v in P_.items():
                ee = tuple(x + y for x, y in zip(e, mm)); vec[red.idx[ee]] += c * Fraction(v)
            nf = red.reduce(vec)
            for r, v in enumerate(nf):
                if v: cols.setdefault(r, {})[(a, b)] = v
        for r, eq in cols.items():
            eqs.append(eq)
    return {'d': d, 'p': 0, 'Dpars': sorted(Dpars), 'Dchars': sorted(Dchars), 'gauge': gauge, 'vars': vars_, 'n': n,
            'n_unknowns_before_gauge': len(unknowns), 'n_gauge': len(removed), 'eqs': eqs, 'n_eqs': len(eqs)}

def write_singular_exact(sysd, path, extra=''):
    n = sysd['n']; names = ','.join(f'u{a}' for a in range(n))
    polys = []
    import math
    for eq in sysd['eqs']:
        den = 1
        for v in eq.values(): den = den * v.denominator // math.gcd(den, v.denominator)
        terms = []
        for (a, b), c in sorted(eq.items()):
            ci = int(c * den)
            terms.append(f'({ci})*u{a}*u{b}' if a != b else f'({ci})*u{a}^2')
        polys.append('+'.join(terms))
    body = f'ring R = 0, ({names}), dp;\nideal I = ' + ',\n'.join(polys) + ';\n'
    body += 'ideal SI = std(I);\n"n_gens", size(I), "gb_size", size(SI), "dim", dim(SI), "vdim", vdim(SI);\n' + extra + 'quit;\n'
    open(path, 'w').write(body)
    return path

# ---------------------------------------------------------------- Koszul-form fixed-D system for D = grad(F).w (E1 lift N0 = A w)
def koszul_system(ctx, w, Dchar=None, Dpar=None, restrict=True):
    """Unknowns: p_ij m (i<j, deg m = d-4) and Q' monomials (deg 2d-6);  equation
         sum_{i<j} p_ij L_ij(m)  -  F Q'  =  -B_D(A w),   L_ij(m) = D[(d_i m)F_j - (d_j m)F_i] - m(F_j D_i - F_i D_j),
         B_D(N) = D div N - N(D).
    w: list of 6 polys (mod p) of degree d-5.  If D is a sigma-eigenvector of pure parity and restrict=True the unknowns are
    restricted to the sector (lossless by averaging).  Returns (cols, colpolys, rhs, info)."""
    p = ctx.p
    D = {}
    for i in range(NV):
        if w[i]: D = add(D, mul(w[i], ctx.Fi[i], p), p)
    d = core.degree(D)
    kD = char_of(D) if Dchar is None else Dchar
    pD = core.parity_of(D) if Dpar is None else Dpar
    restrict = restrict and kD is not None and pD is not None
    Di = [pderiv(D, i, p) for i in range(NV)]
    N0 = [mul(ctx.A, w[i], p) if w[i] else {} for i in range(NV)]
    B = {}
    for i in range(NV):
        if N0[i]:
            B = add(B, mul(D, pderiv(N0[i], i, p), p), p)
            B = sub(B, mul(N0[i], Di[i], p), p)
    cols, polys = [], []
    for i in range(NV):
        for j in range(i + 1, NV):
            G = sub(mul(ctx.Fi[j], Di[i], p), mul(ctx.Fi[i], Di[j], p), p)
            DFj = mul(D, ctx.Fi[j], p); DFi = mul(D, ctx.Fi[i], p)
            par = (ctx.pA + pD + PI[i] + PI[j]) % 2 if restrict else None
            kp = (ctx.kA + kD + W[i] + W[j]) % 7 if restrict else None
            for m in wmonomials(d - 4, par, kp):
                col = scale(shift(G, m, 1, p), p - 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1; col = add(col, shift(DFj, tuple(m2), m[i], p), p)
                if m[j]:
                    m2 = list(m); m2[j] -= 1; col = sub(col, shift(DFi, tuple(m2), m[j], p), p)
                cols.append(('p', i, j, m)); polys.append(col)
    nP = len(cols)
    if 2 * d - 6 >= 0:
        par = ctx.pA if restrict else None
        kq = (ctx.kA + 2 * kD) % 7 if restrict else None
        for m in wmonomials(2 * d - 6, par, kq):
            cols.append(('Q', m)); polys.append(scale(shift(ctx.F, m, 1, p), p - 1, p))
    rhs = scale(B, p - 1, p)
    return cols, polys, rhs, {'d': d, 'kD': kD, 'pD': pD, 'nP': nP, 'nQ': len(cols) - nP, 'D': D}

def solve_cols(polys, rhs, p):
    rows = {}
    for col in polys:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    M = dense_cols(polys, rows, p, rhs)
    C = len(polys)
    Mr, rank = M.rref()
    piv = pivots_of(Mr, rank, C + 1)
    consistent = not (piv and piv[-1] == C)
    rank0 = rank if consistent else rank - 1
    x = None
    if consistent:
        x = [0] * C
        for i, j in enumerate(piv): x[j] = int(Mr[i, C])
    return {'consistent': consistent, 'rank': rank0, 'nullity': C - rank0, 'rows': len(rows), 'cols': C, 'x': x, 'pivots': piv, 'rowidx': rows}

def koszul_test(ctx, w, restrict=True):
    cols, polys, rhs, info = koszul_system(ctx, w, restrict=restrict)
    r = solve_cols(polys, rhs, ctx.p)
    r.update(info); r.pop('x', None); r.pop('pivots', None); r.pop('rowidx', None)
    return r

def w_from_vars(vars_, u, d):
    """w (6 polys mod p) from the variable list of hess_system and a coordinate vector u."""
    w = [{} for _ in range(NV)]
    for (i, m), v in zip(vars_, u):
        if v:
            w[i][m] = w[i].get(m, 0) + v
    return w

# ---------------------------------------------------------------- lifted family L = S_1 . C (gauge-fixed coordinates) and radical membership
def lifted_family_generators(vars_, lift_pars, C_quads, nlam_index=(1, 2, 3, 4, 5)):
    """Generators (as Singular polynomial strings in u-variables) of the ideal of the lifted family
    L = { w_i = lambda_i Z_i l : i in nlam_index, l linear of parity lift_pars (0: span(X), 1: span(Y)), lambda in C },
    inside the coordinate space given by vars_ (list of (i, m) with m of degree 2).  C_quads: list of dicts
    {(i,i'): Fraction} (exact quadrics in lambda_1..lambda_5).  Returns (linear_gens, minor_gens, cquad_gens, Mindex)."""
    lin_vars = [j for j in range(3)] if lift_pars == 0 else [3, 4, 5]
    vidx = {u: a for a, u in enumerate(vars_)}
    Mindex = {}   # (i, j) -> variable index of coefficient of Z_i * Z_j in w_i
    for i in nlam_index:
        for j in lin_vars:
            m = [0] * NV; m[i] += 1; m[j] += 1; m = tuple(m)
            if (i, m) in vidx: Mindex[(i, j)] = vidx[(i, m)]
    used = set(Mindex.values())
    linear = [f'u{a}' for a in range(len(vars_)) if a not in used]
    minors = []
    items = sorted(Mindex.items())
    for (i1, j1), a1 in items:
        for (i2, j2), a2 in items:
            if i1 < i2 and j1 < j2:
                b1 = Mindex.get((i1, j2)); b2 = Mindex.get((i2, j1))
                if b1 is not None and b2 is not None:
                    minors.append(f'u{a1}*u{a2}-u{b1}*u{b2}')
    cq = []
    import math
    for qd in C_quads:
        den = 1
        for v in qd.values(): den = den * v.denominator // math.gcd(den, v.denominator)
        for j in lin_vars:
            for k in lin_vars:
                if j > k: continue
                terms = []
                for (i, ip), v in qd.items():
                    a = Mindex.get((i, j)); b = Mindex.get((ip, k))
                    if a is None or b is None: continue
                    terms.append(f'({int(v*den)})*u{a}*u{b}')
                if terms: cq.append('+'.join(terms))
    return linear, minors, cq, Mindex

def singular_radical_membership(gb_file_lead_or_ms, sysd, gens_by_group, out_path, p, maxk=6, from_msolve_gb=None):
    """Singular script: compute std(I) (or load an msolve reduced GB), then for each generator g test g^k in I for k<=maxk.
    gens_by_group: dict name -> list of polynomial strings."""
    n = sysd['n']; names = ','.join(f'u{a}' for a in range(n))
    polys = []
    for eq in sysd['eqs']:
        polys.append('+'.join((f'{c}*u{a}*u{b}' if a != b else f'{c}*u{a}^2') for (a, b), c in sorted(eq.items())))
    body = f'ring R = {p}, ({names}), dp;\nideal I = ' + ',\n'.join(polys) + ';\n'
    if from_msolve_gb:
        body += f'ideal SI = ' + from_msolve_gb + ';\nattrib(SI, "isSB", 1);\n'
    else:
        body += 'ideal SI = std(I);\n'
    body += '"gb_size", size(SI), "dim", dim(SI);\nint k; int bad = 0;\n'
    for name, gens in gens_by_group.items():
        body += f'"group {name}: {len(gens)} generators";\n'
        for g in gens:
            body += f'k = 1; while (k <= {maxk}) {{ if (reduce(({g})^k, SI) == 0) {{ break; }} k = k + 1; }} if (k > {maxk}) {{ bad = bad + 1; "FAIL", "{g[:60]}"; }}\n'
    body += '"total failures", bad;\nquit;\n'
    open(out_path, 'w').write(body)
    return out_path

def msolve_gb_to_singular(path):
    """read msolve -g 2 output; return the GB as a Singular ideal-list string (polynomials comma-separated)."""
    txt = open(path).read()
    body = txt[txt.index('[') + 1: txt.rindex(']')]
    polys = [t.strip() for t in body.replace('\n', '').split(',') if t.strip()]
    return ','.join(polys)

def lifted_family_generators_general(vars_, mult_monomials, C_quads, nlam_index=(1, 2, 3, 4, 5)):
    """Ideal generators of L = { w_i = lambda_i Z_i m : m in span(mult_monomials), lambda in C } in the coordinates vars_
    (list of (i, m) with deg m = 1 + deg multiplier).  Returns (linear, minors, cquads, Mindex) as Singular strings."""
    import math
    vidx = {u: a for a, u in enumerate(vars_)}
    Mindex = {}
    for i in nlam_index:
        for jm, mm in enumerate(mult_monomials):
            m = list(mm); m[i] += 1; m = tuple(m)
            if (i, m) in vidx: Mindex[(i, jm)] = vidx[(i, m)]
    used = set(Mindex.values())
    linear = [f'u{a}' for a in range(len(vars_)) if a not in used]
    minors = []
    items = sorted(Mindex.items())
    for (i1, j1), a1 in items:
        for (i2, j2), a2 in items:
            if i1 < i2 and j1 < j2:
                b1 = Mindex.get((i1, j2)); b2 = Mindex.get((i2, j1))
                if b1 is not None and b2 is not None:
                    minors.append(f'u{a1}*u{a2}-u{b1}*u{b2}')
    cq = []
    js = sorted({j for (i, j) in Mindex})
    for qd in C_quads:
        den = 1
        for v in qd.values(): den = den * v.denominator // math.gcd(den, v.denominator)
        for j in js:
            for k in js:
                if j > k: continue
                terms = []
                for (i, ip), v in qd.items():
                    a = Mindex.get((i, j)); b = Mindex.get((ip, k))
                    if a is None or b is None: continue
                    terms.append(f'({int(v*den)})*u{a}*u{b}')
                if terms: cq.append('+'.join(terms))
    return linear, minors, cq, Mindex

def e1_lift_modp(ctx, D, restrict=False):
    """solve N(F) = A D mod p (unrestricted); returns N or None."""
    p = ctx.p; D = core.reduce_mod(D, p); d = core.degree(D); pD = core.parity_of(D)
    cols, keys = [], []
    for i in range(NV):
        par = (ctx.pA - pD - PI[i]) % 2 if pD is not None else None
        for m in monomials(d + 1, par):
            cols.append(shift(ctx.Fi[i], m, 1, p)); keys.append((i, m))
    rhs = mul(ctx.A, D, p)
    r = solve_cols(cols, rhs, p)
    if not r['consistent']: return None
    N = [{} for _ in range(NV)]
    for (i, m), v in zip(keys, r['x']):
        if v: N[i][m] = v
    return N

def koszul_test_general(ctx, D, restrict=False):
    """Koszul-form fixed-D test for arbitrary D (E1 lift solved mod p)."""
    p = ctx.p; D = core.reduce_mod(D, p); d = core.degree(D)
    N0 = e1_lift_modp(ctx, D)
    if N0 is None: return {'consistent': False, 'e1': False}
    Di = [pderiv(D, i, p) for i in range(NV)]
    B = {}
    for i in range(NV):
        if N0[i]:
            B = add(B, mul(D, pderiv(N0[i], i, p), p), p); B = sub(B, mul(N0[i], Di[i], p), p)
    cols, polys = [], []
    pD = core.parity_of(D)
    for i in range(NV):
        for j in range(i + 1, NV):
            G = sub(mul(ctx.Fi[j], Di[i], p), mul(ctx.Fi[i], Di[j], p), p)
            DFj = mul(D, ctx.Fi[j], p); DFi = mul(D, ctx.Fi[i], p)
            par = (ctx.pA + pD + PI[i] + PI[j]) % 2 if pD is not None else None
            for m in monomials(d - 4, par):
                col = scale(shift(G, m, 1, p), p - 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1; col = add(col, shift(DFj, tuple(m2), m[i], p), p)
                if m[j]:
                    m2 = list(m); m2[j] -= 1; col = sub(col, shift(DFi, tuple(m2), m[j], p), p)
                polys.append(col)
    if 2 * d - 6 >= 0:
        for m in monomials(2 * d - 6, ctx.pA):
            polys.append(scale(shift(ctx.F, m, 1, p), p - 1, p))
    r = solve_cols(polys, scale(B, p - 1, p), p)
    r['e1'] = True; r.pop('x', None); r.pop('pivots', None); r.pop('rowidx', None)
    return r
