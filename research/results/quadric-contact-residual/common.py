"""Shared helpers for the quadric-contact-residual attack.

Variables X0,X1,X2,Y0,Y1,Y2 = indices 0..5.  All polynomials are dicts {exp tuple: coeff}.
"""
import sys, json, time, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from hodge import core, fixed_d
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, power, shift
import flint
from flint import nmod_mat

IID = 'quadric-contact-residual'
RES = '/tmp/hodge/results/quadric-contact-residual'
LOG = RES + '/log.jsonl'
P0 = core.PRIMES[0]

F, A, q, rec = core.load_instance(IID)
assert q == 1

def chi(e):
    """Character of a monomial under fqlx's affine diagonal group G_aff (order 48)."""
    return ((e[0] + e[3]) % 2, (e[2] + e[5]) % 2, (e[1] - e[4]) % 6, (e[3] + e[4] + e[5]) % 2)

def chi_add(a, b):
    return ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2, (a[2] + b[2]) % 6, (a[3] + b[3]) % 2)

def chi_neg(a):
    return ((-a[0]) % 2, (-a[1]) % 2, (-a[2]) % 6, (-a[3]) % 2)

def chi_of(P):
    cs = {chi(e) for e in P}
    return cs.pop() if len(cs) == 1 else None

CHI_A = chi_of(A)
CHI_Z = [chi(tuple(1 if j == i else 0 for j in range(NV))) for i in range(NV)]

def log(D_expr, d, family, e1, consistent, rank, nullity, seconds, note='', **extra):
    rec = {"D": D_expr, "d": d, "family": family, "e1": e1, "consistent": consistent,
           "rank": rank, "nullity": nullity, "seconds": round(seconds, 2), "note": note}
    rec['engine'] = 2
    rec.update(extra)
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec

def short(P, n=90):
    s = core.to_str(P)
    return s if len(s) <= n else s[:n] + '...'

# ----------------------------------------------------------------------------
# sector-restricted fixed-D solver (q=1), mod p
# ----------------------------------------------------------------------------

class SectorSystem:
    """Full fixed-D linear system restricted to the character sectors forced by D (D must be a
    G_aff eigenvector).  Unknowns: coefficients of N_i (degree d+1, character chi(A)+chi(D)+chi(Z_i))
    and Q (degree 2d-6, character chi(A)+2chi(D)).  Equation: A D^2 - D N(F) + F(D div N - N(D)) = F^2 Q.
    """
    def __init__(self, D, p=P0, verbose=False, F_=None, A_=None, chi_fn=None):
        t0 = time.time()
        self.p = p
        global chi
        chi_save = chi
        if chi_fn is not None:
            chi = chi_fn
        F = F_ if F_ is not None else globals()['F']
        A = A_ if A_ is not None else globals()['A']
        CHI_A = chi_of(A)
        CHI_Z = [chi(tuple(1 if j == i else 0 for j in range(NV))) for i in range(NV)]
        cD = chi_of(D)
        if cD is None:
            raise ValueError('D is not a G_aff eigenvector')
        self.cD = cD
        d = core.degree(D)
        self.d = d
        Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
        self.F, self.A, self.D = Fp, Ap, Dp
        Fi = [pderiv(Fp, i, p) for i in range(NV)]
        Di = [pderiv(Dp, i, p) for i in range(NV)]
        FD = mul(Fp, Dp, p)
        # column for N_i = m:  -D F_i m + F (D d_i m - D_i m)
        G = [scale(add(mul(Dp, Fi[i], p), mul(Fp, Di[i], p), p), p - 1, p) for i in range(NV)]
        # careful: -D F_i m - F D_i m + F D d_i m ; G_i = -(D F_i + F D_i)
        F2 = scale(power(Fp, 2, p), p - 1, p)
        cols, colpolys = [], []
        degN = d + 1
        degQ = 2 * d - 6
        cA2D = chi_add(CHI_A, chi_add(cD, cD))
        for i in range(NV):
            ci = chi_add(CHI_A, chi_add(cD, CHI_Z[i]))
            for m in monomials(degN):
                if chi(m) != ci:
                    continue
                col = shift(G[i], m, 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(FD, tuple(m2), m[i], p), p)
                cols.append(('N', i, m)); colpolys.append(col)
        if degQ >= 0:
            for m in monomials(degQ):
                if chi(m) != cA2D:
                    continue
                cols.append(('Q', -1, m)); colpolys.append(shift(F2, m, 1, p))
        rhs = scale(mul(Ap, mul(Dp, Dp, p), p), p - 1, p)
        rows = {}
        for col in colpolys:
            for e in col:
                rows.setdefault(e, len(rows))
        for e in rhs:
            rows.setdefault(e, len(rows))
        # sanity: all rows should have character cA2D
        bad = [e for e in rows if chi(e) != cA2D]
        assert not bad, f'row character mismatch: {bad[:3]}'
        self.cols, self.colpolys, self.rows, self.rhs = cols, colpolys, rows, rhs
        self.nrows, self.ncols = len(rows), len(cols)
        self.nN = sum(1 for c in cols if c[0] == 'N'); self.nQ = self.ncols - self.nN
        self.t_build = time.time() - t0
        chi = chi_save
        if verbose:
            print(f'[sector d={d} chi(D)={cD}] rows={self.nrows} cols={self.ncols} (N {self.nN}, Q {self.nQ}) '
                  f'built {self.t_build:.1f}s', file=sys.stderr)

    def solve(self, verbose=False):
        t0 = time.time()
        R, C, p = self.nrows, self.ncols, self.p
        W = C + 1
        flat = [0] * (R * W)
        for j, col in enumerate(self.colpolys):
            for e, v in col.items():
                flat[self.rows[e] * W + j] = v
        for e, v in self.rhs.items():
            flat[self.rows[e] * W + C] = v
        M = nmod_mat(R, W, flat, p)
        Mr, rank = M.rref()
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0:
                j += 1
            pivots.append(j)
        consistent = not (pivots and pivots[-1] == C)
        x = None
        if consistent:
            x = [0] * C
            for i, j in enumerate(pivots):
                x[j] = int(Mr[i, C])
        # rank of the coefficient matrix alone
        rank0 = rank - (0 if consistent else 1)
        res = {'consistent': consistent, 'rank': rank0, 'nullity': C - rank0, 'x': x, 'pivots': pivots,
               'rows': R, 'cols': C, 't_rref': round(time.time() - t0, 2)}
        if verbose:
            print(f'[sector solve p={p}] consistent={consistent} rank={rank0} nullity={C-rank0} '
                  f'({time.time()-t0:.1f}s)', file=sys.stderr)
        return res

    def unpack(self, x):
        N = [{} for _ in range(NV)]; Q = {}
        for (kind, i, m), v in zip(self.cols, x):
            if not v: continue
            if kind == 'N': N[i][m] = v
            else: Q[m] = v
        return N, Q


def sector_e1(D, p=P0):
    """E1 test (A D in J) restricted to the sector of D.  Returns (ok, info)."""
    cD = chi_of(D)
    d = core.degree(D)
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    cols = []
    for i in range(NV):
        ci = chi_add(CHI_A, chi_add(cD, CHI_Z[i]))
        for m in monomials(d + 1):
            if chi(m) == ci:
                cols.append(shift(Fi[i], m, 1, p))
    rhs = mul(Ap, Dp, p)
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    ra = M.rank()
    M0 = nmod_mat(R, C, [flat[r * (C + 1) + j] for r in range(R) for j in range(C)], p)
    r0 = M0.rank()
    return ra == r0, {'rank': r0, 'rank_aug': ra, 'rows': R, 'cols': C}


def test_D(D, family, note='', verbose=True, p=P0, lift=False):
    """Run E1 (sector) then the full sector system; log; return result dict."""
    t0 = time.time()
    d = core.degree(D)
    Dstr = core.to_str(D)
    cD = chi_of(D)
    if cD is None:
        if core.parity_of(D) is None:
            r = log(Dstr, d, family, None, None, -1, -1, 0.0, note + ' [mixed Y-parity: skipped; use D*tau(D)]')
            if verbose: print('mixed parity, skipped', short(D), file=sys.stderr)
            return r
        # fall back to the generic engine (mixed characters)
        ok, info = fixed_d.e1_solvable(F, A, q, D, p=p)
        if not ok:
            r = log(Dstr, d, family, False, None, info['rank'], info['cols'] - info['rank'], time.time() - t0, note + ' [generic E1 fail]')
            if verbose: print('E1 FAIL (generic)', short(D), info, file=sys.stderr)
            return r
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose, want_exact=lift)
        r = log(Dstr, d, family, True, out['consistent'], out['rank'], out['nullity'], time.time() - t0, note + ' [generic full]', rows=out['rows'], cols=out['cols'])
        if out['consistent']:
            print('!!! CONSISTENT (generic):', short(D), file=sys.stderr)
            r['out'] = out
        return r
    ok, info = sector_e1(D, p=p)
    if not ok:
        r = log(Dstr, d, family, False, None, info['rank'], info['cols'] - info['rank'], time.time() - t0, note + f' [E1 fail chi={cD}]')
        if verbose: print(f'E1 FAIL d={d} chi={cD}', short(D), info, file=sys.stderr)
        return r
    S = SectorSystem(D, p=p, verbose=verbose)
    res = S.solve(verbose=verbose)
    r = log(Dstr, d, family, True, res['consistent'], res['rank'], res['nullity'], time.time() - t0,
            note + f' [sector chi={cD}]', rows=res['rows'], cols=res['cols'])
    if res['consistent']:
        print('!!! CONSISTENT:', short(D), file=sys.stderr)
        r['system'] = S; r['res'] = res
    return r


# polynomial helpers in the (z,w) description -------------------------------------------------
X0 = {(1,0,0,0,0,0): Fraction(1)}; X1 = {(0,1,0,0,0,0): Fraction(1)}; X2 = {(0,0,1,0,0,0): Fraction(1)}
Y0 = {(0,0,0,1,0,0): Fraction(1)}; Y1 = {(0,0,0,0,1,0): Fraction(1)}; Y2 = {(0,0,0,0,0,1): Fraction(1)}
ONE = {(0,)*6: Fraction(1)}
Z = [X0, X1, X2, Y0, Y1, Y2]

def lin(c):
    """linear form sum c_i Z_i"""
    out = {}
    for i, ci in enumerate(c):
        if ci: out = add(out, scale(Z[i], Fraction(ci)))
    return out

def dir_deriv(Pp, c):
    """directional derivative sum c_i d_i P"""
    out = {}
    for i, ci in enumerate(c):
        if ci: out = add(out, scale(pderiv(Pp, i), Fraction(ci)))
    return out

def exact_div(Pp, Qq):
    """exact division of multivariate polynomials (dicts), assumes divisible; returns quotient."""
    import sympy as sp
    xs = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
    def to_sp(P):
        return sum(sp.Rational(c.numerator, c.denominator) * sp.prod([x**k for x, k in zip(xs, e)]) for e, c in P.items())
    qq, rr = sp.div(sp.Poly(to_sp(Pp), *xs), sp.Poly(to_sp(Qq), *xs))
    assert rr.is_zero, 'not divisible'
    out = {}
    for mon, c in qq.terms():
        out[tuple(int(k) for k in mon)] = Fraction(int(c.p), int(c.q))
    return out
