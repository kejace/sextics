"""Diagonal symmetries of F and character sectors of the fixed-D system.

Projective diagonal group
-------------------------
A diagonal element acts by z_i -> zeta^{w_i} z_i.  F is an eigenvector iff all
monomials of F have the same weight m.w (mod the order).  Normalizing w_5 = 0
(projective quotient by scalars) the group is

    G = { w in (Q/Z)^5 : (m - m0).w in Z for all m in supp F }

which is dual to Z^5 / L, L = lattice spanned by the exponent differences.
We compute L by Hermite normal form (python-flint) and enumerate G explicitly
as integer weight vectors modulo n = det(L).  For a sextic with all sixth
powers this reproduces the (Z/6)^5 group of scripts/audit-frontier.py
diagonal_screen; for the Klein block it gives the order-126 lift.

Characters (all modulo n = group.order):
    chi_P(g)      = m.w for any monomial m of P         (None if P is not an eigenvector)
    residue char  = chi_A + sum(w) - (q+1) chi_F         (character of A Omega / F^{q+1})

Sector restriction.  If F, A, D are eigenvectors of g then every term of the
cleared identity transforms by chi_A + 2 chi_D once
    chi(N_i) = chi_A + chi_D - chi_F + w_i,      chi(Q) = chi_A + 2 chi_D - (q+1) chi_F,
and averaging a solution over the stabilizer subgroup keeps it a solution, so
restricting the unknowns to that sector is lossless.  The restricted system is
smaller by roughly the order of the stabilizer.

Public API
----------
diagonal_group(F)                       -> DiagonalGroup (order n, elements, generators)
character(P, w, n)                      -> int or None
residue_character(F, A, q, w, n)        -> int or None
stabilizer(group, polys)                -> list of weight vectors fixing all polys' lines
sector_filter(F, A, q, D, group)        -> (allowed(kind, i, m) predicate, subgroup)
RestrictedFixedDSystem(F, A, q, D, p, allowed) -- FixedDSystem with filtered unknowns
solve_fixed_D_symmetric(F, A, q, D, ...) -> (status, data) like fixed_d.solve_fixed_D
"""
from __future__ import annotations
import itertools, math, sys, time
from fractions import Fraction
from . import core, fixed_d
from .core import NV, monomials, shift, mul, add, sub, scale, pderiv, power, reduce_mod

try:
    from flint import fmpz_mat, fmpq_mat, nmod_mat
except ImportError:  # pragma: no cover
    fmpz_mat = None


class DiagonalGroup:
    def __init__(self, order, elements, generators, hnf):
        self.order = order            # weights are integers modulo this
        self.elements = elements      # list of 6-tuples (last entry 0)
        self.generators = generators  # list of 6-tuples
        self.hnf = hnf

    def __len__(self):
        return len(self.elements)

    def __iter__(self):
        return iter(self.elements)

    def __repr__(self):
        return f'DiagonalGroup(size={len(self.elements)}, order={self.order})'


def diagonal_group(F, max_size=10**7):
    """Projective diagonal symmetry group of the line of F (see module docstring)."""
    supp = sorted(F)
    m0 = supp[0]
    rows = [[m[i] - m0[i] for i in range(5)] for m in supp[1:]]
    if not rows:
        raise ValueError('F must have at least two monomials')
    B = fmpz_mat(rows)
    H = B.hnf()
    hrows = [[int(H[i, j]) for j in range(5)] for i in range(H.nrows())]
    hrows = [r for r in hrows if any(r)]
    if len(hrows) < 5:
        raise ValueError('positive-dimensional diagonal stabilizer (F degenerate for this method)')
    H5 = fmpz_mat(hrows[:5])
    det = abs(int(H5.det()))
    if det == 0:
        raise ValueError('singular HNF block')
    if det > max_size:
        raise ValueError(f'group too large to enumerate ({det})')
    # adj(H5) = det * H5^{-1}, integer matrix
    Hinv = fmpq_mat(H5).inv()
    adj = [[int(Hinv[i, j] * det) for j in range(5)] for i in range(5)]
    diag = [hrows[i][i] for i in range(5)]
    gens = []
    for j in range(5):
        w = tuple((adj[i][j]) % det for i in range(5)) + (0,)
        if any(w):
            gens.append(w)
    elements = []
    for k in itertools.product(*[range(abs(h)) for h in diag]):
        w = tuple(sum(adj[i][j] * k[j] for j in range(5)) % det for i in range(5)) + (0,)
        elements.append(w)
    # reduce the modulus by the common gcd of all weights (Fermat: 7776 -> 6)
    g = det
    for w in elements:
        for v in w:
            g = math.gcd(g, v)
            if g == 1:
                break
        if g == 1:
            break
    if g > 1:
        det //= g
        elements = [tuple(v // g for v in w) for w in elements]
        gens = [tuple(v // g for v in w) for w in gens]
    # sanity: every element fixes the line of F
    for w in elements[:50]:
        assert character(F, w, det) is not None
    return DiagonalGroup(det, elements, gens, hrows[:5])


def character(P, w, n):
    """Weight m.w mod n common to all monomials of P, or None if P is not an eigenvector."""
    vals = {sum(a * b for a, b in zip(m, w)) % n for m in P}
    return vals.pop() if len(vals) == 1 else None


def residue_character(F, A, q, w, n):
    """Character of the residue form A Omega / F^(q+1) under the diagonal element w."""
    cA = character(A, w, n)
    cF = character(F, w, n)
    if cA is None or cF is None:
        return None
    return (cA + sum(w) - (q + 1) * cF) % n


def residue_invariant_elements(F, A, q, group):
    """Elements of the group under which the residue class of A is invariant (character 0)."""
    return [w for w in group.elements if residue_character(F, A, q, w, group.order) == 0]


def stabilizer(group, polys):
    """Subgroup of elements for which every given polynomial is an eigenvector."""
    n = group.order
    return [w for w in group.elements if all(character(P, w, n) is not None for P in polys)]


def _weight_table(mons, ws, n):
    """numpy array (len(mons) x len(ws)) of m.w mod n."""
    import numpy as np
    E = np.array(mons, dtype=np.int64)
    W = np.array(ws, dtype=np.int64).T  # 6 x G
    return (E @ W) % n


def sector_filter(F, A, q, D, group):
    """Predicate allowed(kind, i, m) for fixed-D unknowns in the character sector forced by D.

    kind in {'N','Q'}; i the numerator index for kind 'N'.  Uses the subgroup of
    `group` stabilizing the lines of A and D.  Returns (allowed, subgroup)."""
    n = group.order
    sub = stabilizer(group, [A, D])
    if not sub:
        return (lambda kind, i, m: True), sub
    cF = [character(F, w, n) for w in sub]
    cA = [character(A, w, n) for w in sub]
    cD = [character(D, w, n) for w in sub]
    reqN = [[(cA[g] + cD[g] - cF[g] + sub[g][i]) % n for g in range(len(sub))] for i in range(NV)]
    reqQ = [(cA[g] + 2 * cD[g] - (q + 1) * cF[g]) % n for g in range(len(sub))]

    def allowed(kind, i, m):
        req = reqN[i] if kind == 'N' else reqQ
        for g, w in enumerate(sub):
            if sum(a * b for a, b in zip(m, w)) % n != req[g]:
                return False
        return True
    allowed.reqN, allowed.reqQ, allowed.subgroup, allowed.n = reqN, reqQ, sub, n
    return allowed, sub


def filter_monomials(mons, kind, i, allowed):
    """Vectorized filter of a monomial list by a sector_filter predicate."""
    sub = getattr(allowed, 'subgroup', None)
    if not sub or not mons:
        return [m for m in mons if allowed(kind, i, m)]
    import numpy as np
    req = np.array(allowed.reqN[i] if kind == 'N' else allowed.reqQ, dtype=np.int64)
    E = np.array(mons, dtype=np.int64)
    W = np.array(sub, dtype=np.int64).T
    T = (E @ W) % allowed.n
    ok = np.all(T == req[None, :], axis=1)
    return [m for m, k in zip(mons, ok) if k]


class RestrictedFixedDSystem(fixed_d.FixedDSystem):
    """FixedDSystem whose unknown monomials are filtered by allowed(kind, i, m)."""

    def __init__(self, F, A, q, D, p, allowed, verbose=False):
        t0 = time.time()
        self.p = p; self.q = q
        Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); Dp = reduce_mod(D, p)
        self.F, self.A, self.D = Fp, Ap, Dp
        d = core.degree(D); self.d = d
        pars = fixed_d.numerator_parities(D, A); self.pars = pars
        degN = d + 6 * q - 5
        degQ = 2 * d - 6
        Fi = [pderiv(Fp, i, p) for i in range(NV)]
        Di = [pderiv(Dp, i, p) for i in range(NV)]
        FD = mul(Fp, Dp, p)
        G = [add(scale(mul(Dp, Fi[i], p), q, p), mul(Fp, Di[i], p), p) for i in range(NV)]
        G = [scale(g, p - 1, p) for g in G]
        Fq1 = scale(power(Fp, q + 1, p), p - 1, p)
        cols, colpolys = [], []
        self.full_cols = 0
        for i in range(NV):
            ms = monomials(degN, pars[i])
            self.full_cols += len(ms)
            for m in filter_monomials(ms, 'N', i, allowed):
                col = shift(G[i], m, 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(FD, tuple(m2), m[i], p), p)
                cols.append(('N', i, m)); colpolys.append(col)
        if degQ >= 0:
            ms = monomials(degQ, core.parity_of(A))
            self.full_cols += len(ms)
            for m in filter_monomials(ms, 'Q', -1, allowed):
                cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1, p))
        self.nN = sum(1 for c in cols if c[0] == 'N')
        self.nQ = sum(1 for c in cols if c[0] == 'Q')
        D2 = mul(Dp, Dp, p)
        rhs = scale(mul(Ap, D2, p), p - 1, p)
        self.target_cols = []
        rows = {}
        for col in colpolys:
            for e in col:
                if e not in rows:
                    rows[e] = len(rows)
        for e in rhs:
            if e not in rows:
                rows[e] = len(rows)
        self.cols, self.rows, self.colpolys, self.rhs = cols, rows, colpolys, rhs
        self.nrows, self.ncols = len(rows), len(cols)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[restricted system] d={d} q={q} rows={self.nrows} cols={self.ncols} '
                  f'(of {self.full_cols} unrestricted; N {self.nN}, Q {self.nQ}) '
                  f'nnz={sum(len(c) for c in colpolys)} built in {self.build_time:.1f}s', file=sys.stderr)


def restricted_system(F, A, q, D, p=None, group=None, verbose=False):
    """Build the sector-restricted fixed-D system.  Returns (system, subgroup)."""
    p = p or core.PRIMES[0]
    group = group or diagonal_group(F)
    allowed, sub = sector_filter(F, A, q, D, group)
    S = RestrictedFixedDSystem(F, A, q, D, p, allowed, verbose=verbose)
    return S, sub


def solve_fixed_D_symmetric(F, A, q, D, primes=None, verbose=False, max_primes=8, group=None,
                            compressed=False):
    """Like fixed_d.solve_fixed_D but on the character-restricted system.

    Returns (status, data) with status in {'inconsistent', 'solved', 'lift-failed'}.
    data always has 'subgroup_size', 'cols', 'full_cols', 'rows', 'rank', 'nullity'."""
    primes = primes or core.PRIMES
    group = group or diagonal_group(F)
    acc = None; pivots0 = None
    d = core.degree(D)
    last = None
    for k, p in enumerate(primes[:max_primes]):
        S, sub = restricted_system(F, A, q, D, p, group, verbose=verbose)
        r = fixed_d.solve_compressed(S, verbose=verbose) if compressed else S.solve(verbose=verbose)
        base = {'prime': p, 'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows,
                'cols': S.ncols, 'full_cols': S.full_cols, 'subgroup_size': len(sub),
                'group_size': len(group), 'd': d}
        if not r['consistent']:
            return 'inconsistent', base
        last = base
        if pivots0 is None:
            pivots0 = r['pivots']
        elif pivots0 != r['pivots']:
            acc = None; pivots0 = r['pivots']
        x = r['x']
        if acc is None:
            acc = (x, p)
        else:
            xs, m = acc
            xs = [core.crt_pair(a, m, b, p)[0] for a, b in zip(xs, x)]
            acc = (xs, m * p)
        xs, m = acc
        rat = [core.ratrecon(a, m) for a in xs]
        if all(v is not None for v in rat):
            N = [{} for _ in range(NV)]; Q = {}
            for (kind, i, mon), v in zip(S.cols, rat):
                if v == 0: continue
                if kind == 'N': N[i][mon] = v
                else: Q[mon] = v
            ok, info = core.verify_witness(F, A, q, D, N, Q)
            if ok:
                base.update({'N': N, 'Q': Q, 'primes_used': k + 1})
                return 'solved', base
            elif verbose:
                print(f'[lift] reconstruction failed verification after {k+1} primes', file=sys.stderr)
    last = last or {}
    last['primes_used'] = max_primes
    return 'lift-failed', last
