"""wave-2 EXCLUSION library for one-reflection-skew (q = 1), six homogeneous variables X0,X1,X2,Y0,Y1,Y2.

Everything exact unless a function name ends in `_modp`.  Polynomials are python-flint fmpz_mpoly
(integer coefficients).  We work with F7 = 7*F (integer coefficients); scaling F does not change
the existence of a certificate for fixed D (replace (N,Q) by (N/c, Q/c^2)).

Certificate (normalized, q=1, deg D = d):   E1: N(F) = A D  (degree d+6),   E2: B_D(N) = F Q (degree 2d),
B_D(N) = D div N - N(D).  E1+E2 for some (N,Q) <=> a certificate (D,N,Q) with this D exists.

Gauges used (all valid at fixed degree):  D ~ D + F L (L in S_{d-6}),  N ~ N + K (Koszul syzygy), N ~ N + F u.

Diagonal symmetry group: 24 elements, generators g1=(3,3,3,0,0,0), g2=(0,3,0,3,0,0), g3=(4,1,4,4,1,0) (weights mod 6,
projective normalization w5 = 0).  Character of a monomial e: (e.g1, e.g2, e.g3) mod 6.
"""
import sys, time, json, itertools, math, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from hodge import core
from flint import fmpz_mpoly_ctx, fmpz_mpoly, Ordering, fmpz_mat, fmpq_mat, nmod_mat, fmpz, fmpq

VARS = ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2')
NV = 6
ctx = fmpz_mpoly_ctx.get(VARS, Ordering.degrevlex)
Z = ctx.gens()
GENS = [(3, 3, 3, 0, 0, 0), (0, 3, 0, 3, 0, 0), (4, 1, 4, 4, 1, 0)]
LANE = '/tmp/hodge/results/one-reflection-skew'
LOG = os.environ.get('W2X_LOG', LANE + '/wave2-excl.log.jsonl')

def _jsdefault(o):
    try:
        return int(o)
    except Exception:
        return str(o)

def logexp(**rec):
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(LOG, 'a') as fh:
        fh.write(json.dumps(rec, default=_jsdefault) + '\n')

def from_dict(d):
    return ctx.from_dict({tuple(e): int(c) for e, c in d.items() if int(c)})

def load_instance(iid='one-reflection-skew'):
    Fd, Ad, q, rec = core.load_instance(iid)
    den = 1
    for c in Fd.values():
        den = den * Fraction(c).denominator // math.gcd(den, Fraction(c).denominator)
    Fint = from_dict({e: Fraction(c) * den for e, c in Fd.items()})
    Aint = from_dict(Ad)
    return Fint, Aint, q, den

F7, A, q, FDEN = load_instance()
Fi = [F7.derivative(i) for i in range(NV)]
assert q == 1

# ---------------------------------------------------------------- monomials and characters
_mons = {}
def mons(n):
    if n < 0: return []
    if n not in _mons:
        _mons[n] = core.monomials(n)
    return _mons[n]

def chr_of(e):
    return tuple(sum(a * b for a, b in zip(e, g)) % 6 for g in GENS)

def chr_add(*cs):
    return tuple(sum(x) % 6 for x in zip(*cs))

def chr_neg(c):
    return tuple((-x) % 6 for x in c)

W = [tuple(g[i] for g in GENS) for i in range(NV)]   # character of Z_i
CHI_F = chr_of(next(iter(F7.to_dict())))
CHI_A = chr_of(next(iter(A.to_dict())))
assert all(chr_of(e) == CHI_F for e in F7.to_dict()), 'F not an eigenvector'
ALLCHR = sorted({chr_of(e) for n in (6, 7) for e in mons(n)} | {chr_of(e) for e in mons(8)})

_sec = {}
def secmons(n, chi):
    key = (n, chi)
    if key not in _sec:
        _sec[key] = [e for e in mons(n) if chr_of(e) == chi]
    return _sec[key]

_secidx = {}
def secidx(n, chi):
    key = (n, chi)
    if key not in _secidx:
        _secidx[key] = {e: k for k, e in enumerate(secmons(n, chi))}
    return _secidx[key]

def poly_chr(P):
    cs = {chr_of(e) for e in P.to_dict()}
    return cs.pop() if len(cs) == 1 else None

def poly_deg(P):
    ds = {sum(e) for e in P.to_dict()}
    assert len(ds) <= 1, 'not homogeneous'
    return ds.pop() if ds else -1

def vec_of(P, n, chi):
    """sparse {index: int} of P inside sector (n, chi); asserts membership."""
    idx = secidx(n, chi)
    out = {}
    for e, c in P.to_dict().items():
        out[idx[e]] = int(c)
    return out

def mat_from_cols(cols, nrows):
    """cols: list of sparse {row: int}; returns fmpz_mat nrows x len(cols)."""
    M = fmpz_mat(nrows, len(cols))
    for j, col in enumerate(cols):
        for r, v in col.items():
            M[r, j] = v
    return M

def nmod_from_cols(cols, nrows, p):
    M = nmod_mat(nrows, len(cols), p)
    for j, col in enumerate(cols):
        for r, v in col.items():
            M[r, j] = v % p
    return M

def monom(e, c=1):
    return ctx.from_dict({tuple(e): int(c)})

# ---------------------------------------------------------------- Jacobian ideal in a sector
def jgens(n, chi):
    """generators m*F_i of J_n in sector chi: list of (i, m, poly)."""
    out = []
    if n < 5: return out
    for i in range(NV):
        if Fi[i].is_zero(): continue
        ci = poly_chr(Fi[i])
        need = chr_add(chi, chr_neg(ci))
        for m in secmons(n - 5, need):
            out.append((i, m, monom(m) * Fi[i]))
    return out

_jdual = {}
def jac_dual(n, chi):
    """Exact integer basis Y (fmpz_mat, nrows_sector x k) of the left kernel of the J_n[chi] generator
    matrix, i.e. of the functionals on S_n[chi] vanishing on J_n[chi]; k = dim (S/J)_n[chi].
    Also returns rank of J_n[chi]."""
    key = (n, chi)
    if key in _jdual: return _jdual[key]
    rows = secmons(n, chi)
    gens = jgens(n, chi)
    cols = [vec_of(g[2], n, chi) for g in gens]
    M = mat_from_cols(cols, len(rows))          # nrows x ngens
    MT = M.transpose()                          # ngens x nrows ; left kernel of M = right kernel of MT
    if MT.nrows() == 0:
        Y = fmpz_mat(len(rows), len(rows))
        for i in range(len(rows)): Y[i, i] = 1
        _jdual[key] = (Y, 0); return _jdual[key]
    Yfull, k = MT.nullspace()                   # nrows x nrows, first k columns span the kernel
    Y = fmpz_mat(len(rows), k)
    for i in range(len(rows)):
        for j in range(k):
            Y[i, j] = Yfull[i, j]
    rank = len(rows) - k
    _jdual[key] = (Y, rank)
    return _jdual[key]

def in_jacobian_exact(P):
    """Exact test P in J (P homogeneous eigenvector): all J-dual functionals vanish."""
    n = poly_deg(P); chi = poly_chr(P)
    if n < 0: return True
    Y, _ = jac_dual(n, chi)
    v = vec_of(P, n, chi)
    for j in range(Y.ncols()):
        s = 0
        for r, c in v.items():
            s += int(Y[r, j]) * c
        if s: return False
    return True

# ---------------------------------------------------------------- exact linear solving helpers
def exact_solve(M, b, p=2147483647):
    """Solve M x = b exactly over Q (M fmpz_mat r x c, b fmpz_mat r x 1) or return None if inconsistent.
    Uses a mod-p rref to find a nonsingular pivot block, Dixon-solves it exactly, verifies on all rows."""
    r, c = M.nrows(), M.ncols()
    Mp = nmod_mat(r, c + 1, p)
    for i in range(r):
        for j in range(c):
            v = int(M[i, j])
            if v: Mp[i, j] = v % p
        v = int(b[i, 0])
        if v: Mp[i, c] = v % p
    R, rank = Mp.rref()
    piv = []
    for i in range(rank):
        j = 0
        while j <= c and int(R[i, j]) == 0: j += 1
        piv.append(j)
    if piv and piv[-1] == c:
        return None                                 # inconsistent mod p => (modular evidence only)
    # pivot rows: rows of Mp whose restriction to pivot columns has rank `rank`: find via rref of the transpose
    pc = piv
    MpT = nmod_mat(len(pc), r, p)
    for a, j in enumerate(pc):
        for i in range(r):
            v = int(M[i, j])
            if v: MpT[a, i] = v % p
    RT, rk2 = MpT.rref()
    assert rk2 == rank
    prow = []
    for i in range(rank):
        j = 0
        while j < r and int(RT[i, j]) == 0: j += 1
        prow.append(j)
    S = fmpz_mat(rank, rank)
    for a, i in enumerate(prow):
        for bb, j in enumerate(pc):
            S[a, bb] = M[i, j]
    rhs = fmpz_mat(rank, 1)
    for a, i in enumerate(prow): rhs[a, 0] = b[i, 0]
    xs = S.solve(rhs)                               # fmpq_mat
    x = [fmpq(0)] * c
    for bb, j in enumerate(pc): x[j] = xs[bb, 0]
    # exact verification on all rows
    for i in range(r):
        s = fmpq(0)
        for j in range(c):
            v = int(M[i, j])
            if v and x[j] != 0: s += v * x[j]
        if s != b[i, 0]:
            return None
    return x

def exact_left_kernel(M):
    """Exact integer basis of {y : y^T M = 0} as list of lists (fmpz_mat nullspace of M^T)."""
    MT = M.transpose()
    Y, k = MT.nullspace()
    out = []
    for j in range(k):
        out.append([int(Y[i, j]) for i in range(M.nrows())])
    return out

# ---------------------------------------------------------------- E1 kernel in a sector, exact
def gauge_monomials(d):
    """monomials X0^6 * m, m in S_{d-6}: coefficient constraints fixing the gauge D ~ D + F L."""
    out = []
    for m in mons(d - 6):
        e = list(m); e[0] += 6
        out.append(tuple(e))
    return out

def e1_kernel(d, chi, gauge=True, verbose=False):
    """Exact basis of K_d[chi] = {D in S_d[chi] : A D in J} (with gauge coefficients zero if gauge).
    Returns dict with 'basis' (list of fmpz_mpoly), 'dimK', 'dimJ' (dim J_d[chi]), 'extra' (dimK - dimJ)."""
    t0 = time.time()
    Dm = secmons(d, chi)
    tchi = chr_add(chi, CHI_A)
    Y, rankJ = jac_dual(d + 6, tchi)                # functionals on S_{d+6}[tchi] killing J
    k = Y.ncols()
    # C = Y^T [A m_D]  (k x nD)
    cols = [vec_of(A * monom(m), d + 6, tchi) for m in Dm]
    nrows = len(secmons(d + 6, tchi))
    Mc = mat_from_cols(cols, nrows)
    C = Y.transpose() * Mc                          # k x nD
    if gauge:
        gm = [m for m in gauge_monomials(d) if chr_of(m) == chi]
        if gm:
            C2 = fmpz_mat(C.nrows() + len(gm), C.ncols())
            for i in range(C.nrows()):
                for j in range(C.ncols()): C2[i, j] = C[i, j]
            idx = secidx(d, chi)
            for a, m in enumerate(gm):
                C2[C.nrows() + a, idx[m]] = 1
            C = C2
    if C.nrows() == 0:
        basis = [monom(m) for m in Dm]
    else:
        K, kd = C.nullspace()
        basis = []
        for j in range(kd):
            P = ctx.from_dict({Dm[i]: int(K[i, j]) for i in range(len(Dm)) if int(K[i, j])})
            basis.append(P)
    # dim J_d[chi] (after gauge: J_d[chi] / F S_{d-6}[chi])
    if d >= 5:
        _, rJ = jac_dual(d, chi)
    else:
        rJ = 0
    ng = len([m for m in gauge_monomials(d) if chr_of(m) == chi]) if gauge else 0
    out = {'d': d, 'chi': chi, 'nD': len(Dm), 'dimK': len(basis), 'dimJ': rJ - ng, 'dimJ_raw': rJ,
           'gauge_removed': ng, 'extra': len(basis) - (rJ - ng), 'dimSJ_target': k, 'seconds': round(time.time() - t0, 2)}
    out['basis'] = basis
    if verbose: print({kk: v for kk, v in out.items() if kk != 'basis'}, file=sys.stderr, flush=True)
    return out

def e1_lift(D, verbose=False):
    """Exact N (6 fmpz_mpoly) and integer den with sum_i N_i F_i = den * A * D, D eigenvector.  None if impossible."""
    d = poly_deg(D); chi = poly_chr(D)
    tchi = chr_add(chi, CHI_A)
    n = d + 6
    gens = jgens(n, tchi)
    rows = secmons(n, tchi)
    cols = [vec_of(g[2], n, tchi) for g in gens]
    M = mat_from_cols(cols, len(rows))
    b = fmpz_mat(len(rows), 1)
    for r, v in vec_of(A * D, n, tchi).items(): b[r, 0] = v
    x = exact_solve(M, b)
    if x is None: return None
    den = 1
    for v in x:
        dd = int(fmpq(v).q)
        den = den * dd // math.gcd(den, dd)
    N = [ctx.from_dict({}) for _ in range(NV)]
    for (i, m, _), v in zip(gens, x):
        if v != 0:
            N[i] += monom(m, int(v * den))
    assert sum((N[i] * Fi[i] for i in range(NV)), ctx.from_dict({})) == A * D * den
    return N, den

def BD(D, N):
    """B_D(N) = D div N - N(D)."""
    divN = sum((N[i].derivative(i) for i in range(NV)), ctx.from_dict({}))
    ND = sum((N[i] * D.derivative(i) for i in range(NV)), ctx.from_dict({}))
    return D * divN - ND

def koszul_cols(d, chiD):
    """Koszul correction basis for deg N = d+1: (i, j, c monomial) with K_i = c F_j, K_j = -c F_i,
    char(c) = chi_A + chi_D + w_i + w_j.  Returns list of (i, j, m)."""
    out = []
    if d - 4 < 0: return out
    for i in range(NV):
        for j in range(i + 1, NV):
            need = chr_add(CHI_A, chiD, W[i], W[j])
            for m in secmons(d - 4, need):
                out.append((i, j, m))
    return out

def BD_koszul(D, i, j, m):
    """B_D(K) for K = c (F_j e_i - F_i e_j), c = monomial m: F_j (D c_i - c D_i) - F_i (D c_j - c D_j)."""
    c = monom(m)
    ci, cj = c.derivative(i), c.derivative(j)
    return Fi[j] * (D * ci - c * D.derivative(i)) - Fi[i] * (D * cj - c * D.derivative(j))

# ---------------------------------------------------------------- fixed-D exact test (sector restricted)
def fixed_D_system(D, p=None, use_e1=True, verbose=False):
    """Normalized system for fixed eigenvector D: unknowns N_i (deg d+1, sector-filtered), Q (deg 2d-6);
    equations E1: N(F) - A D = 0 (deg d+6), E2: B_D(N) - F Q = 0 (deg 2d).  Returns (M, b, meta) exact fmpz."""
    d = poly_deg(D); chi = poly_chr(D)
    assert chi is not None, 'D must be a diagonal eigenvector'
    Ncols = []      # (i, m)
    for i in range(NV):
        need = chr_add(CHI_A, chi, chr_neg(CHI_F), W[i])
        for m in secmons(d + 1, need):
            Ncols.append((i, m))
    Qchi = chr_add(CHI_A, chi, chi, chr_neg(CHI_F), chr_neg(CHI_F))
    Qcols = secmons(2 * d - 6, Qchi) if 2 * d - 6 >= 0 else []
    e1chi = chr_add(CHI_A, chi); e2chi = chr_add(CHI_A, chi, chi)
    r1 = secmons(d + 6, e1chi); r2 = secmons(2 * d, e2chi)
    n1 = len(r1); n2 = len(r2)
    cols = []
    Di = [D.derivative(i) for i in range(NV)]
    for (i, m) in Ncols:
        mm = monom(m)
        col = {}
        if use_e1:
            for r, v in vec_of(mm * Fi[i], d + 6, e1chi).items(): col[r] = v
        # E2 part: D * d_i(m) - m * D_i
        P2 = D * mm.derivative(i) - mm * Di[i]
        for r, v in vec_of(P2, 2 * d, e2chi).items(): col[n1 + r] = v
        cols.append(col)
    for m in Qcols:
        P2 = -(F7 * monom(m))
        cols.append({n1 + r: v for r, v in vec_of(P2, 2 * d, e2chi).items()})
    M = mat_from_cols(cols, n1 + n2)
    b = fmpz_mat(n1 + n2, 1)
    if use_e1:
        for r, v in vec_of(A * D, d + 6, e1chi).items(): b[r, 0] = v
    meta = {'d': d, 'chi': chi, 'nN': len(Ncols), 'nQ': len(Qcols), 'rows': n1 + n2, 'rowsE1': n1, 'rowsE2': n2}
    return M, b, meta, (Ncols, Qcols)

def exact_dual_certificate(M, b, p=2147483647):
    """If M x = b is inconsistent mod p, produce an exact rational y with y^T M = 0 and y^T b != 0 (proof of
    inconsistency over Q and C).  Returns (y as list of fmpq, y^T b) or None if consistent mod p."""
    r, c = M.nrows(), M.ncols()
    # left kernel mod p containing a vector with y.b != 0:  work with [M | b]^T
    MT = nmod_mat(c + 1, r, p)
    for i in range(r):
        for j in range(c):
            v = int(M[i, j])
            if v: MT[j, i] = v % p
        v = int(b[i, 0])
        if v: MT[c, i] = v % p
    # rank of M^T vs [M|b]^T
    rM = nmod_mat(c, r, p)
    for i in range(r):
        for j in range(c):
            v = int(M[i, j])
            if v: rM[j, i] = v % p
    rank_M = rM.rank(); rank_Mb = MT.rank()
    if rank_Mb == rank_M:
        return None
    # find y: y^T M = 0, y^T b = 1.  Solve [M | b]^T y = e_last  over Q: pick pivot structure mod p.
    R, rk = MT.rref()
    # Use exact_solve on the (c+1) x r system MT y = e_{c}
    MTz = fmpz_mat(c + 1, r)
    for i in range(r):
        for j in range(c):
            MTz[j, i] = M[i, j]
        MTz[c, i] = b[i, 0]
    e = fmpz_mat(c + 1, 1); e[c, 0] = 1
    y = exact_solve(MTz, e, p)
    if y is None:
        return None
    # verify
    for j in range(c):
        s = fmpq(0)
        for i in range(r):
            v = int(M[i, j])
            if v and y[i] != 0: s += v * y[i]
        assert s == 0
    s = fmpq(0)
    for i in range(r):
        if y[i] != 0: s += int(b[i, 0]) * y[i]
    assert s == 1
    return y, s

def test_fixed_D(D, name, family, p=2147483647, want_exact_dual=True, verbose=False):
    """Full E1+E2 fixed-D test in D's sector; mod-p rank decision + exact dual certificate when inconsistent;
    exact solve attempt when consistent mod p.  Logs and returns a dict."""
    t0 = time.time()
    d = poly_deg(D); chi = poly_chr(D)
    M, b, meta, (Ncols, Qcols) = fixed_D_system(D, verbose=verbose)
    r, c = M.nrows(), M.ncols()
    Mp = nmod_mat(r, c + 1, p)
    for i in range(r):
        for j in range(c):
            v = int(M[i, j])
            if v: Mp[i, j] = v % p
        v = int(b[i, 0])
        if v: Mp[i, c] = v % p
    Mp0 = nmod_mat(r, c, p)
    for i in range(r):
        for j in range(c):
            v = int(M[i, j])
            if v: Mp0[i, j] = v % p
    rank = int(Mp0.rank()); rank_aug = int(Mp.rank())
    out = dict(meta); out.update({'name': name, 'family': family, 'p': p, 'rank': rank, 'rank_aug': rank_aug,
                                  'nullity': c - rank, 'consistent_modp': rank_aug == rank})
    if rank_aug == rank:
        x = exact_solve(M, b, p)
        if x is not None:
            den = 1
            for v in x:
                dd = int(v.q); den = den * dd // math.gcd(den, dd)
            N = [ctx.from_dict({}) for _ in range(NV)]
            for (i, m), v in zip(Ncols, x[:len(Ncols)]):
                if v != 0: N[i] += monom(m, int(v * den))
            Qp = ctx.from_dict({})
            for m, v in zip(Qcols, x[len(Ncols):]):
                if v != 0: Qp += monom(m, int(v * den))
            # verify cleared identity:  A D^2 den - D N(F) + F B_D(N) = F^2 Q   (all scaled by den)
            lhs = A * D * D * den - D * sum((N[i] * Fi[i] for i in range(NV)), ctx.from_dict({})) + F7 * BD(D, N)
            ok = (lhs == F7 * F7 * Qp)
            out['status'] = 'WITNESS' if ok else 'consistent-modp-lift-failed'
            out['witness'] = {'den': den, 'D': str(D), 'N': [str(n) for n in N], 'Q': str(Qp), 'verified': ok}
        else:
            out['status'] = 'consistent-modp-lift-failed'
    else:
        out['status'] = 'inconsistent-modp'
        if want_exact_dual:
            cert = exact_dual_certificate(M, b, p)
            if cert is not None:
                y, s = cert
                out['status'] = 'EXCLUDED-exact'
                out['dual_height_bits'] = max(int(abs(v.p)).bit_length() + int(v.q).bit_length() for v in y if v != 0)
    out['seconds'] = round(time.time() - t0, 2)
    logexp(exp='fixed-D', D=name, family=family, d=d, chi=list(chi), prime=p, result=out['status'], rank=rank,
           rank_aug=rank_aug, cols=c, rows=r, seconds=out['seconds'], note=f"sector-restricted E1+E2; nN={meta['nN']} nQ={meta['nQ']}")
    if verbose:
        print(name, {k: v for k, v in out.items() if k not in ('witness',)}, file=sys.stderr, flush=True)
    return out

# ---------------------------------------------------------------- Hessian-type relaxations (exact)
class Family:
    """A linear family of E1-admissible denominators: basis D_a (eigenvectors) with exact lifts (N_a, den_a)."""
    def __init__(self, d, items):
        # items: list of (name, D, N, den, chi, kind)
        self.d = d
        self.items = items
        self.n = len(items)
        self.chi = [it[4] for it in items]

    @staticmethod
    def from_e1(d, gauge=True, sectors=None, verbose=False):
        items = []
        for chi in (sectors or ALLCHR):
            r = e1_kernel(d, chi, gauge=gauge)
            # Jacobian part first: identify which basis elements are in J (exactly) -- lift generically
            for k, D in enumerate(r['basis']):
                lift = e1_lift(D)
                assert lift is not None
                N, den = lift
                items.append((f'd{d}:{chi}:{k}', D, N, den, chi, 'J' if in_jacobian_exact(D) else 'extra'))
            if verbose: print('sector', chi, 'dimK', r['dimK'], 'extra', r['extra'], file=sys.stderr, flush=True)
        return Family(d, items)

    def pairs(self):
        for a in range(self.n):
            for b in range(a, self.n):
                yield a, b

def pair_B(fam, a, b):
    """Symmetrized B-term of the quadratic map, scaled to integers: for D = sum lambda_a D_a, N = sum lambda_a N_a/den_a,
    B_D(N) = sum_{a<=b} lambda_a lambda_b * (B_{D_a}(N_b)/den_b + B_{D_b}(N_a)/den_a) [a<b], B_{D_a}(N_a)/den_a [a=b].
    We return the integer polynomial den_a den_b * (that) and the scale (den_a*den_b)."""
    _, Da, Na, dena, _, _ = fam.items[a]
    _, Db, Nb, denb, _, _ = fam.items[b]
    if a == b:
        return BD(Da, Na), dena
    return BD(Da, Nb) * dena + BD(Db, Na) * denb, dena * denb

def quadrics_from_duals(fam, duals_by_sector, verbose=False, which_pairs=None):
    """duals_by_sector: {psi: fmpz_mat Y (nrows x k)} exact functionals on S_{2d}[psi] (e.g. J-duals or relaxed duals).
    Returns list of quadrics, each as dict {(a,b): int} with a<=b, from y^T b_{ab}; coefficients exact integers.
    Rational scales 1/(den_a den_b) are cleared per quadric (ideal unchanged)."""
    n2 = 2 * fam.d
    # group pairs by sector
    bysec = {}
    for a, b in (which_pairs or fam.pairs()):
        psi = chr_add(fam.chi[a], fam.chi[b], CHI_A)
        bysec.setdefault(psi, []).append((a, b))
    quads = []
    for psi, prs in bysec.items():
        if psi not in duals_by_sector: continue
        Y = duals_by_sector[psi]
        k = Y.ncols()
        if k == 0: continue
        nrows = len(secmons(n2, psi))
        cols = []; scales = []
        for a, b in prs:
            P, sc = pair_B(fam, a, b)
            cols.append(vec_of(P, n2, psi) if not P.is_zero() else {})
            scales.append(sc)
        Mb = mat_from_cols(cols, nrows)
        Qm = Y.transpose() * Mb                     # k x npairs, entry = y . b_ab * (den_a den_b)
        L = 1
        for sc in scales: L = L * sc // math.gcd(L, sc)
        for i in range(k):
            qd = {}
            for j, (a, b) in enumerate(prs):
                v = int(Qm[i, j])
                if v: qd[(a, b)] = v * (L // scales[j])
            if qd:
                g = 0
                for v in qd.values(): g = math.gcd(g, abs(v))
                quads.append({'psi': psi, 'coef': {kk: v // g for kk, v in qd.items()}})
        if verbose: print('sector', psi, 'pairs', len(prs), 'duals', k, file=sys.stderr, flush=True)
    return quads

def jac_duals_for(fam):
    n2 = 2 * fam.d
    out = {}
    for psi in ALLCHR:
        Y, _ = jac_dual(n2, psi)
        out[psi] = Y
    return out

def quad_matrix(quads, n):
    """Integer matrix rows=quadrics, cols=pairs (a<=b) indexed; returns (fmpz_mat, pair list)."""
    prs = sorted({k for qd in quads for k in qd['coef']})
    idx = {k: j for j, k in enumerate(prs)}
    M = fmpz_mat(len(quads), len(prs))
    for i, qd in enumerate(quads):
        for k, v in qd['coef'].items():
            M[i, idx[k]] = v
    return M, prs

def span_analysis(quads, n, p=2147483647, exact=True, verbose=False):
    """Which pair-monomials lambda_a lambda_b lie in the Q-span of the quadrics (as vectors in Sym^2)?
    Returns dict: rank, squares_in_span (list a), monomials_in_span (list (a,b)).  Exact ranks via fmpz_mat."""
    M, prs = quad_matrix(quads, n)
    if M.nrows() == 0:
        return {'rank': 0, 'squares': [], 'monos': [], 'npairs': len(prs)}
    # exact rref gives membership test: e_j in rowspace iff rank([M; e_j]) == rank(M)
    R, den, rank = M.rref()
    # rowspace membership for unit vectors: e_j in rowspace iff the rref has a row equal to e_j up to scale
    # (true iff column j is a pivot column and that pivot row has no other nonzero entries)
    squares = []; monos = []
    piv_rows = {}
    for i in range(rank):
        j = 0
        while j < R.ncols() and int(R[i, j]) == 0: j += 1
        piv_rows[j] = i
    for j, (a, b) in enumerate(prs):
        if j in piv_rows:
            i = piv_rows[j]
            if all(int(R[i, jj]) == 0 for jj in range(R.ncols()) if jj != j):
                (squares if a == b else monos).append((a, b))
    return {'rank': rank, 'squares': [a for a, b in squares], 'monos': monos, 'npairs': len(prs), 'nquads': len(quads)}

# ---------------------------------------------------------------- relaxed mod-F system with Koszul corrections (exact duals)
def relaxed_duals(fam, subset=None, verbose=False, include_jduals=False):
    """For each row sector psi of S_{2d}: exact integer left-kernel basis of the column space
    [F*S_{2d-6}] + [B_{D_a}(K) : a in subset, K Koszul correction in D_a's sector].
    Returns {psi: fmpz_mat Y}.  (These duals contain the J-duals; the resulting quadrics y^T b(lambda) are
    necessary conditions for E2 after relaxing the products lambda_a * c_K to independent unknowns.)"""
    d = fam.d; n2 = 2 * d
    subset = list(range(fam.n)) if subset is None else list(subset)
    cols_by = {}
    for psi in ALLCHR:
        cols_by[psi] = []
    # F * S_{2d-6}
    for m in mons(n2 - 6) if n2 - 6 >= 0 else []:
        P = F7 * monom(m)
        psi = chr_add(chr_of(m), CHI_F)
        cols_by[psi].append(vec_of(P, n2, psi))
    ncorr = 0
    # row sectors that carry at least one alive pair (others give no quadrics)
    need = {chr_add(fam.chi[a], fam.chi[b], CHI_A) for a in subset for b in subset if a <= b}
    for a in subset:
        _, Da, _, _, chia, _ = fam.items[a]
        for psi in need:
            # Koszul K with N-character psi - chi_a, i.e. koszul_cols 'chiD' argument = psi - chi_a - chi_A
            chiK = chr_add(psi, chr_neg(chia), chr_neg(CHI_A))
            for (i, j, m) in koszul_cols(d, chiK):
                P = BD_koszul(Da, i, j, m)
                if P.is_zero(): continue
                assert poly_chr(P) == psi
                cols_by[psi].append(vec_of(P, n2, psi)); ncorr += 1
    out = {}
    for psi in ALLCHR:
        rows = secmons(n2, psi)
        if not rows: continue
        cols = cols_by[psi]
        M = mat_from_cols(cols, len(rows))
        if M.ncols() == 0:
            Y = fmpz_mat(len(rows), len(rows))
            for i in range(len(rows)): Y[i, i] = 1
            out[psi] = Y; continue
        Yfull, k = M.transpose().nullspace()
        Y = fmpz_mat(len(rows), k)
        for i in range(len(rows)):
            for j in range(k): Y[i, j] = Yfull[i, j]
        out[psi] = Y
        if verbose: print('relaxed sector', psi, 'rows', len(rows), 'cols', len(cols), 'duals', k, file=sys.stderr, flush=True)
    if verbose: print('relaxed: correction columns', ncorr, file=sys.stderr, flush=True)
    return out

def restrict_quads(quads, alive):
    """Set lambda_a = 0 for a not in alive: drop pairs touching dead indices."""
    alive = set(alive)
    out = []
    for qd in quads:
        c = {k: v for k, v in qd['coef'].items() if k[0] in alive and k[1] in alive}
        if c:
            out.append({'psi': qd['psi'], 'coef': c})
    return out

def kill_squares(quads, alive, verbose=False):
    """Iterate: squares in the exact span => kill those coordinates.  Returns (alive, history)."""
    alive = set(alive); hist = []
    while True:
        qs = restrict_quads(quads, alive)
        an = span_analysis(qs, None)
        hist.append({'alive': len(alive), 'rank': an['rank'], 'npairs': an['npairs'], 'squares': len(an['squares']), 'monos': len(an['monos'])})
        if verbose: print(hist[-1], file=sys.stderr, flush=True)
        if not an['squares']:
            return alive, hist, an
        alive -= set(an['squares'])
        if not alive:
            return alive, hist, an

def quads_to_msolve(quads, alive, path, p=1073741827, extra_linear=()):
    """Write the restricted quadrics (and optional linear forms) as an msolve input file over F_p (or char 0 if p=0)."""
    alive = sorted(alive)
    names = [f'l{a}' for a in alive]
    pos = {a: k for k, a in enumerate(alive)}
    lines = []
    for qd in restrict_quads(quads, alive):
        terms = []
        for (a, b), v in qd['coef'].items():
            vv = v % p if p else v
            if vv == 0: continue
            terms.append(f'{vv}*l{a}*l{b}' if a != b else f'{vv}*l{a}^2')
        if terms: lines.append('+'.join(terms))
    for lf in extra_linear:
        lines.append(lf)
    with open(path, 'w') as fh:
        fh.write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(lines) + '\n')
    return len(lines), len(names)

# ---------------------------------------------------------------- mod-p span analysis + exact certification
def span_analysis_modp(quads, p=2147483647):
    prs = sorted({k for qd in quads for k in qd['coef']})
    idx = {k: j for j, k in enumerate(prs)}
    if not quads or not prs:
        return {'rank': 0, 'squares': [], 'monos': [], 'npairs': len(prs), 'nquads': len(quads), 'prs': prs}
    M = nmod_mat(len(quads), len(prs), p)
    for i, qd in enumerate(quads):
        for k, v in qd['coef'].items():
            M[i, idx[k]] = v % p
    R, rank = M.rref()
    squares = []; monos = []
    for i in range(rank):
        nz = [j for j in range(len(prs)) if int(R[i, j]) != 0]
        if len(nz) == 1:
            a, b = prs[nz[0]]
            (squares if a == b else monos).append((a, b))
    return {'rank': rank, 'squares': [a for a, b in squares], 'monos': monos, 'npairs': len(prs), 'nquads': len(quads), 'prs': prs}

def exact_solve_multi(M, B, p=2147483647):
    """Solve M X = B exactly over Q for a matrix RHS B (fmpz_mat r x t).  Returns list of booleans (consistent per column)
    and the fmpq_mat X (columns for inconsistent RHS are meaningless).  Pivot block from a mod-p rref of M, Dixon solve,
    verification by exact matrix product."""
    r, c = M.nrows(), M.ncols(); t = B.ncols()
    Mp = nmod_mat(r, c, p)
    for i in range(r):
        for j in range(c):
            v = int(M[i, j])
            if v: Mp[i, j] = v % p
    R, rank = Mp.rref()
    pc = []
    for i in range(rank):
        j = 0
        while j < c and int(R[i, j]) == 0: j += 1
        pc.append(j)
    MpT = nmod_mat(len(pc), r, p)
    for a, j in enumerate(pc):
        for i in range(r):
            v = int(M[i, j])
            if v: MpT[a, i] = v % p
    RT, rk2 = MpT.rref(); assert rk2 == rank
    prow = []
    for i in range(rank):
        j = 0
        while j < r and int(RT[i, j]) == 0: j += 1
        prow.append(j)
    S = fmpz_mat(rank, rank)
    for a, i in enumerate(prow):
        for bb, j in enumerate(pc): S[a, bb] = M[i, j]
    rhs = fmpz_mat(rank, t)
    for a, i in enumerate(prow):
        for k in range(t): rhs[a, k] = B[i, k]
    xs = S.solve(rhs)                                    # fmpq_mat rank x t
    X = fmpq_mat(c, t)
    for bb, j in enumerate(pc):
        for k in range(t): X[j, k] = xs[bb, k]
    prod = fmpq_mat(M) * X
    ok = []
    for k in range(t):
        good = all(prod[i, k] == B[i, k] for i in range(r))
        ok.append(good)
    return ok, X

def certify_in_span(quads, targets, p=2147483647):
    """Exact certificate that each target pair-monomial (a,b) lies in the Q-span of the quadrics.
    Returns list of certified targets (combination verified exactly by a rational matrix product)."""
    if not targets: return []
    M, prs = quad_matrix(quads, None)
    idx = {k: j for j, k in enumerate(prs)}
    MT = M.transpose()                                   # npairs x nquads
    B = fmpz_mat(len(prs), len(targets))
    for k, t in enumerate(targets): B[idx[t], k] = 1
    ok, X = exact_solve_multi(MT, B, p)
    return [t for t, g in zip(targets, ok) if g]

def kill_squares(quads, alive, verbose=False, p=2147483647, certify=True):
    """Iterate: squares in the span (mod p, then certified exactly) => kill those coordinates."""
    alive = set(alive); hist = []
    while True:
        qs = restrict_quads(quads, alive)
        an = span_analysis_modp(qs, p)
        sq = [(a, a) for a in an['squares']]
        if certify and sq:
            cert = certify_in_span(qs, sq, p)
            an['squares_certified'] = [a for a, _ in cert]
            killed = set(an['squares_certified'])
        else:
            killed = set(an['squares'])
        hist.append({'alive': len(alive), 'rank': an['rank'], 'npairs': an['npairs'], 'nquads': an['nquads'],
                     'squares_modp': len(an['squares']), 'squares_certified': len(killed), 'monos': len(an['monos'])})
        if verbose: print(hist[-1], file=sys.stderr, flush=True)
        if not killed:
            return alive, hist, an
        alive -= killed
        if not alive:
            return alive, hist, an

def macaulay_full_rank(quads, alive, k, p=2147483647):
    """Macaulay matrix of the (exact integer) quadrics in the alive variables in degree k, rank mod p.
    Full column rank (== number of degree-k monomials) proves V_C(quadrics) ∩ (alive coords) = {0}."""
    alive = sorted(alive); s = len(alive); pos = {a: i for i, a in enumerate(alive)}
    qs = restrict_quads(quads, alive)
    # monomials of degree k in s variables as exponent tuples
    def monos_deg(kk):
        out = []
        for c in itertools.combinations(range(kk + s - 1), s - 1):
            e = []; prev = -1
            for x in c:
                e.append(x - prev - 1); prev = x
            e.append(kk + s - 2 - prev)
            out.append(tuple(e))
        return out
    cols = monos_deg(k); cidx = {m: j for j, m in enumerate(cols)}
    mult = monos_deg(k - 2)
    rows = []
    for qd in qs:
        for mm in mult:
            row = {}
            for (a, b), v in qd['coef'].items():
                e = list(mm); e[pos[a]] += 1; e[pos[b]] += 1
                row[cidx[tuple(e)]] = (row.get(cidx[tuple(e)], 0) + v) % p
            rows.append(row)
    M = nmod_mat(len(rows), len(cols), p)
    for i, row in enumerate(rows):
        for j, v in row.items():
            if v: M[i, j] = v
    r = M.rank()
    return r, len(cols), len(rows)

# ---------------------------------------------------------------- Singular leaf analysis (dim / mult / Hilbert) mod p
def singular_leaf(quads, alive, path, p=1073741827, timeout=3600):
    """Run Singular on the restricted quadrics mod p: returns dict(dim, mult, hilb_numerator, leading_pure_powers)."""
    import subprocess, re
    al = sorted(alive)
    polys = []
    for qd in restrict_quads(quads, al):
        terms = []
        for (a, b), v in qd['coef'].items():
            vv = v % p
            if vv: terms.append(f'{vv}*l{a}*l{b}')
        if terms: polys.append('+'.join(terms))
    with open(path, 'w') as fh:
        fh.write(f'ring r = {p}, (' + ','.join(f'l{a}' for a in al) + '), dp;\n')
        fh.write('ideal I = ' + ',\n'.join(polys) + ';\n')
        fh.write('option(redSB); ideal G = std(I);\n')
        fh.write('print("DIM " + string(dim(G)));\nprint("MULT " + string(mult(G)));\n')
        fh.write('intvec h = hilb(G, 1); print("HILB1 " + string(h));\nintvec h2 = hilb(G, 2); print("HILB2 " + string(h2));\n')
        fh.write('ideal L = lead(G); print("NLEAD " + string(size(L)));\n')
        fh.write('int i; string s = ""; for (i = 1; i <= nvars(r); i++) { if (reduce(var(i)^50, std(L)) == 0) { s = s + "1"; } else { s = s + "0"; } } print("PURE " + s);\n')
        fh.write('quit;\n')
    cp = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=timeout)
    out = {'raw': cp.stdout[-3000:], 'stderr': cp.stderr[-500:]}
    for key in ('DIM', 'MULT', 'HILB1', 'HILB2', 'NLEAD', 'PURE'):
        m = re.search(rf'^{key} (.*)$', cp.stdout, re.M)
        if m: out[key] = m.group(1).strip()
    if 'PURE' in out:
        out['no_pure_power'] = [al[i] for i, ch in enumerate(out['PURE']) if ch == '0']
    return out

def quads_vanish_at_axis(quads, a):
    """Exact: every quadric has zero coefficient on lambda_a^2  <=>  the axis through e_a lies in V_Q(I)."""
    return all(qd['coef'].get((a, a), 0) == 0 for qd in quads)
