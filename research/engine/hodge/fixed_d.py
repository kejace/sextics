"""Fixed-denominator certificate solver (exact, modular + rational lifting).

Given F, A, q and a homogeneous denominator D, decide whether numerators N and a
quotient Q exist with  A D^2 - q D N(F) + F (D div N - N(D)) = F^(q+1) Q,
with deg N_i = d + 6q - 5 and deg Q = 2d - 6.

The unknown coefficients of N and Q enter linearly, so this is finite linear
algebra.  We solve mod large primes with python-flint and lift to Q by CRT +
rational reconstruction, then verify exactly.
"""
from __future__ import annotations
import os as _os
import sys, time
from fractions import Fraction
from . import core
from .core import NV, monomials, shift, mul, add, sub, scale, pderiv, power, yparity


ENGINE_VERSION = 3
import sys as _sys
if not _os.environ.get('HODGE_QUIET'):
    _sys.stderr.write('[hodge engine v3] NOTE: engine v2 (2026-11-01 17:25) fixed a sign bug in the F*N(D) term; results from before '
                      'that time must be recomputed.  Engine v3 (2026-11-02 00:55 UTC) reports rank/nullity of the COEFFICIENT matrix '
                      'for inconsistent systems (v2 reported the augmented rank, i.e. rank+1 / nullity-1; consistency verdicts and all '
                      'consistent-case numbers were already correct).  Record "engine": 3 in your logs from now on.\n')

try:
    import flint
    from flint import nmod_mat, fmpq_mat, fmpz_mat
except ImportError:  # pragma: no cover
    flint = None


def set_threads(n=None):
    """Set the number of FLINT threads used by dense nmod_mat rref/nullspace (4x at 8 threads for n=4000).

    Default: FLINT_NUM_THREADS from the environment (rjob submit --cpus N sets it), else 1.  Returns the value."""
    if n is None:
        try:
            n = int(_os.environ.get('FLINT_NUM_THREADS', '1'))
        except ValueError:
            n = 1
    n = max(1, int(n))
    if flint is not None:
        flint.ctx.threads = n
    return n


THREADS = set_threads()


def numerator_parities(D, A, F=None):
    """Parity each N_i must have so that every term of E is Y-odd given D's parity."""
    pD = core.parity_of(D)
    pA = core.parity_of(A)
    if pD is None:
        raise ValueError('D has mixed Y-parity; choose a pure parity denominator')
    # term D N_i F_i must have parity of A D^2 == pA.  parity(F_i) = 0 for i<3, 1 for i>=3.
    return [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(NV)]


class FixedDSystem:
    """Builds the sparse linear system for a fixed D (mod p)."""

    def __init__(self, F, A, q, D, p, with_target_space=False, verbose=False):
        t0 = time.time()
        self.p = p
        self.q = q
        Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
        self.F, self.A, self.D = Fp, Ap, Dp
        d = core.degree(D)
        self.d = d
        pars = numerator_parities(D, A)
        self.pars = pars
        degN = d + 6 * q - 5
        degQ = 2 * d - 6
        Fi = [pderiv(Fp, i, p) for i in range(NV)]
        Di = [pderiv(Dp, i, p) for i in range(NV)]
        FD = mul(Fp, Dp, p)
        G = [add(scale(mul(Dp, Fi[i], p), q, p), mul(Fp, Di[i], p), p) for i in range(NV)]  # q D F_i + F D_i
        # column for N_i = m:  -q D F_i m + F (D d_i m - D_i m) = -(G_i) m + FD * d_i m
        G = [scale(g, p - 1, p) for g in G]
        Fq1 = scale(power(Fp, q + 1, p), p - 1, p)
        cols = []      # list of (kind, i, exponent)
        colpolys = []  # list of dict polys
        for i in range(NV):
            for m in monomials(degN, pars[i]):
                col = shift(G[i], m, 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(FD, tuple(m2), m[i], p), p)
                cols.append(('N', i, m)); colpolys.append(col)
        for m in monomials(degQ, (core.parity_of(A) + 0) % 2 if degQ >= 0 else None):
            cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1, p))
        self.nN = sum(1 for c in cols if c[0] == 'N')
        self.nQ = sum(1 for c in cols if c[0] == 'Q')
        D2 = mul(Dp, Dp, p)
        rhs = scale(mul(Ap, D2, p), p - 1, p)   # move A D^2 to the right side: sum cols = -A D^2
        self.target_cols = []
        if with_target_space:
            # columns for a variable target: + m * D^2 for Y-odd m of degree 6q
            for m in monomials(6 * q, 1):
                cols.append(('A', -1, m)); colpolys.append(shift(D2, m, 1, p))
        rows = {}
        for col in colpolys:
            for e in col:
                if e not in rows:
                    rows[e] = len(rows)
        for e in rhs:
            if e not in rows:
                rows[e] = len(rows)
        self.cols = cols
        self.rows = rows
        self.nrows = len(rows)
        self.ncols = len(cols)
        self.colpolys = colpolys
        self.rhs = rhs
        self.build_time = time.time() - t0
        if verbose:
            print(f'[system] d={d} q={q} rows={self.nrows} cols={self.ncols} (N {self.nN}, Q {self.nQ}) '
                  f'nnz={sum(len(c) for c in colpolys)} built in {self.build_time:.1f}s', file=sys.stderr)

    def dense_augmented(self):
        """Dense nmod_mat [M | b]."""
        R, C, p = self.nrows, self.ncols, self.p
        flat = [0] * (R * (C + 1))
        W = C + 1
        rows = self.rows
        for j, col in enumerate(self.colpolys):
            for e, v in col.items():
                flat[rows[e] * W + j] = v
        for e, v in self.rhs.items():
            flat[rows[e] * W + C] = v
        return nmod_mat(R, C + 1, flat, p)

    def solve(self, verbose=False):
        """Returns dict with keys: consistent, rank, nullity, x (list mod p, free vars 0)."""
        t0 = time.time()
        M = self.dense_augmented()
        t1 = time.time()
        Mr, rank = M.rref()
        t2 = time.time()
        C = self.ncols
        pivots = []
        r = 0
        # identify pivot columns
        for i in range(rank):
            # find first nonzero in row i
            j = 0
            while j <= C and int(Mr[i, j]) == 0:
                j += 1
            pivots.append(j)
        consistent = not (pivots and pivots[-1] == C)
        rank_aug = rank
        if not consistent:
            rank -= 1            # rank of the coefficient matrix M (the last pivot is the rhs column)
            pivots = pivots[:-1]
        x = None
        if consistent:
            x = [0] * C
            for i, j in enumerate(pivots):
                x[j] = int(Mr[i, C])
        res = {'consistent': consistent, 'rank': rank, 'rank_aug': rank_aug, 'nullity': C - rank, 'x': x,
               'pivots': pivots, 'rows': self.nrows, 'cols': C,
               't_build': round(t1 - t0, 2), 't_rref': round(t2 - t1, 2)}
        if verbose:
            print(f'[solve p={self.p}] consistent={consistent} rank={rank} nullity={C-rank} '
                  f'(dense build {t1-t0:.1f}s, rref {t2-t1:.1f}s)', file=sys.stderr)
        return res

    def unpack(self, x):
        """Split a solution vector into (N list, Q) as polys mod p."""
        N = [{} for _ in range(NV)]
        Q = {}
        for (kind, i, m), v in zip(self.cols, x):
            if not v:
                continue
            if kind == 'N':
                N[i][m] = v
            elif kind == 'Q':
                Q[m] = v
        return N, Q


def solve_fixed_D(F, A, q, D, primes=None, verbose=False, max_primes=8):
    """Try to find exact rational (N, Q) for the given D.

    Returns (status, data): status in {'inconsistent', 'solved', 'lift-failed'}.
    """
    primes = primes or core.PRIMES
    d = core.degree(D)
    acc = None  # (residues list, modulus)
    pivots0 = None
    for k, p in enumerate(primes[:max_primes]):
        S = FixedDSystem(F, A, q, D, p, verbose=verbose)
        r = S.solve(verbose=verbose)
        if not r['consistent']:
            return 'inconsistent', {'prime': p, 'rank': r['rank'], 'nullity': r['nullity'],
                                    'rows': S.nrows, 'cols': S.ncols}
        if pivots0 is None:
            pivots0 = r['pivots']
        elif pivots0 != r['pivots']:
            if verbose:
                print('[lift] pivot structure differs between primes; restarting accumulation', file=sys.stderr)
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
            # unpack and verify exactly
            N = [{} for _ in range(NV)]; Q = {}
            for (kind, i, mon), v in zip(S.cols, rat):
                if v == 0: continue
                if kind == 'N': N[i][mon] = v
                elif kind == 'Q': Q[mon] = v
            ok, info = core.verify_witness(F, A, q, D, N, Q)
            if ok:
                return 'solved', {'N': N, 'Q': Q, 'primes_used': k + 1, 'rank': r['rank'],
                                  'nullity': r['nullity'], 'd': d}
            elif verbose:
                print(f'[lift] rational reconstruction verified false after {k+1} primes; continuing', file=sys.stderr)
    return 'lift-failed', {'primes_used': max_primes}


def e1_solvable(F, A, q, D, p=None):
    """Necessary condition: q N(F) = A D solvable, i.e. A*D in the Jacobian ideal (mod p)."""
    p = p or core.PRIMES[0]
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
    d = core.degree(D)
    pars = numerator_parities(D, A)
    degN = d + 6 * q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    cols = []
    for i in range(NV):
        for m in monomials(degN, pars[i]):
            cols.append(shift(Fi[i], m, q, p))
    rhs = mul(Ap, Dp, p)
    rows = {}
    for col in cols:
        for e in col:
            rows.setdefault(e, len(rows))
    for e in rhs:
        rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items():
        flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    rank_aug = M.rank()
    M0 = nmod_mat(R, C, [flat[r * (C + 1) + j] for r in range(R) for j in range(C)], p)
    rank0 = M0.rank()
    return rank_aug == rank0, {'rank': rank0, 'rank_aug': rank_aug, 'rows': R, 'cols': C}


def target_image(F, q, D, p=None, verbose=False):
    """Compute the space of Y-odd targets A of degree 6q served by denominator D (mod p).

    Returns (dim, basis) where basis is a list of polys (mod p) spanning the image, expressed
    in the monomial basis of degree-6q Y-odd polynomials.
    """
    p = p or core.PRIMES[0]
    A0 = {}
    # use a dummy odd target to set parities: pick the first odd monomial
    A0 = {monomials(6 * q, 1)[0]: 1}
    S = FixedDSystem(F, A0, q, D, p, with_target_space=True, verbose=verbose)
    # homogeneous system: columns N, Q, A ; we want nullspace projected to A-coords
    R, C = S.nrows, S.ncols
    flat = [0] * (R * C)
    for j, col in enumerate(S.colpolys):
        for e, v in col.items():
            flat[S.rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    a_idx = [j for j, c in enumerate(S.cols) if c[0] == 'A']
    if nullity == 0:
        return 0, [], S
    # projection of nullspace vectors to A-coordinates
    proj = nmod_mat(len(a_idx), nullity, [int(X[j, k]) for j in a_idx for k in range(nullity)], p)
    rank = proj.rank()
    # extract a basis of the column space: rref of the transpose
    T = proj.transpose()
    Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, j in enumerate(a_idx):
            v = int(Tr[i, jj])
            if v:
                poly[S.cols[j][2]] = v
        basis.append(poly)
    return rank, basis, S


# ----------------------------------------------------------------------------
# compressed (sparse -> random row combinations) solver
# ----------------------------------------------------------------------------
import numpy as np
import scipy.sparse as sp
from fractions import Fraction as _Fr


def _sparse_matrix(S, with_rhs=True):
    """scipy CSR (rows x (cols[+1])) of residues mod S.p as int64, plus exact entries."""
    rows = S.rows
    ri, ci, vals = [], [], []
    for j, col in enumerate(S.colpolys):
        for e, v in col.items():
            ri.append(rows[e]); ci.append(j); vals.append(v)
    C = S.ncols
    if with_rhs:
        for e, v in S.rhs.items():
            ri.append(rows[e]); ci.append(C); vals.append(v)
    M = sp.csr_matrix((np.array(vals, dtype=np.int64), (np.array(ri), np.array(ci))),
                      shape=(S.nrows, C + (1 if with_rhs else 0)))
    return M


def _compress(M, p, k, seed=1, chunk=4096):
    """Return (k x ncols) int64 array = (R @ M) mod p for random R with entries in [0,p).

    R is generated in chunks of `chunk` rows of M so memory stays ~ k*chunk*8 bytes (the old
    one-shot version needed k*nrows*8 bytes: 54 GB at d=12).  Exact: 16-bit split, needs column
    weights < 2^16 (checked)."""
    from .sge import compress_rows, _check_colweights
    _check_colweights(M)
    return compress_rows(M, p, k, seed=seed, chunk=chunk)


def solve_compressed(S, k_extra=40, seed=1, verbose=False):
    """Consistency + rank + particular solution mod p via compressed rows."""
    t0 = time.time()
    M = _sparse_matrix(S)
    C = S.ncols
    k = min(S.nrows, C + 1 + k_extra)
    Cm = _compress(M, S.p, k, seed)
    t1 = time.time()
    Mn = nmod_mat(k, C + 1, Cm.ravel().tolist(), S.p)
    Mr, rank = Mn.rref()
    t2 = time.time()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0:
            j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == C)
    rank_aug = rank
    if not consistent:
        rank -= 1                # rank of the coefficient matrix M (the last pivot is the rhs column)
        pivots = pivots[:-1]
    x = None
    if consistent:
        x = [0] * C
        for i, j in enumerate(pivots):
            x[j] = int(Mr[i, C])
    if verbose:
        print(f'[compressed p={S.p}] k={k} consistent={consistent} rank={rank} nullity={C-rank} '
              f'(compress {t1-t0:.1f}s, rref {t2-t1:.1f}s)', file=sys.stderr)
    return {'consistent': consistent, 'rank': rank, 'rank_aug': rank_aug, 'nullity': C - rank, 'x': x, 'pivots': pivots,
            'rows': S.nrows, 'cols': C, 'k': k, 't_compress': round(t1 - t0, 2), 't_rref': round(t2 - t1, 2)}


def solve_sparse(S, stop_cols=400, verbose=False, want_nullspace=False, seed=1, **kw):
    """Consistency + rank + particular solution mod p via sparse elimination (hodge.sge) with a dense finish.

    Same result schema as solve_compressed (consistent, rank, nullity, x, pivots, rows, cols, timings) plus
    the sparse-phase statistics.  `x` is verified by an exact sparse matvec (key 'verified'); on failure
    consistent=None.  With want_nullspace the full nullspace basis (uint32 array nullity x cols) is returned.
    Extra keyword arguments go to sge.solve_sparse_system (need_rank, max_nnz, cands, extra, ...)."""
    from . import sge
    t0 = time.time()
    M = _sparse_matrix(S)
    r = sge.solve_sparse_system(M, S.p, with_rhs=True, stop_cols=stop_cols, seed=seed, verbose=verbose,
                                want_nullspace=want_nullspace, **kw)
    r.update({'rows': S.nrows, 'cols': S.ncols, 't_matrix': round(time.time() - t0 - r.get('t_total', 0), 2)})
    if verbose:
        print(f'[sparse p={S.p}] consistent={r["consistent"]} rank={r["rank"]} nullity={r["nullity"]} '
              f'(sparse {r["t_sparse"]}s, dense {r.get("t_dense")}s, total {r["t_total"]}s)', file=sys.stderr)
    return r


def exact_solution_from_pivots(F, A, q, D, pivots, verbose=False, seed=7):
    """Exact rational particular solution (free variables zero) using a random square
    compression over Q and flint's fmpq_mat.solve.  Returns (N, Q) or None."""
    # Build exact system (Fractions)
    d = core.degree(D)
    pars = numerator_parities(D, A)
    degN = d + 6 * q - 5
    degQ = 2 * d - 6
    Fi = [pderiv(F, i) for i in range(NV)]
    Di = [pderiv(D, i) for i in range(NV)]
    FD = mul(F, D)
    G = [scale(add(scale(mul(D, Fi[i]), q), mul(F, Di[i])), -1) for i in range(NV)]  # -(q D F_i + F D_i)
    Fq1 = scale(power(F, q + 1), -1)
    cols, colpolys = [], []
    for i in range(NV):
        for m in monomials(degN, pars[i]):
            col = shift(G[i], m, 1)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i]))
            cols.append(('N', i, m)); colpolys.append(col)
    for m in monomials(degQ, core.parity_of(A)):
        cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1))
    rhs = scale(mul(A, mul(D, D)), -1)
    rows = {}
    for col in colpolys:
        for e in col:
            rows.setdefault(e, len(rows))
    for e in rhs:
        rows.setdefault(e, len(rows))
    C = len(cols); R = len(rows)
    piv = pivots
    r = len(piv)
    # random 0/1 compression to r rows, restricted to pivot columns
    rng = np.random.default_rng(seed)
    Rm = rng.integers(0, 2, size=(r, R))
    # accumulate exact rational entries: use common denominator per column via Fractions
    # To keep it fast, scale each column polynomial to integers.
    import math as _m
    def to_int_poly(P):
        den = 1
        for v in P.values():
            den = den * _Fr(v).denominator // _m.gcd(den, _Fr(v).denominator)
        return {e: int(_Fr(v) * den) for e, v in P.items()}, den
    colint = []
    for j in piv:
        P, den = to_int_poly(colpolys[j])
        colint.append((P, den))
    rhsint, rden = to_int_poly(rhs)
    # Build integer matrix entries using python ints (exact)
    from flint import fmpq_mat, fmpq
    Mq = fmpq_mat(r, r)
    bq = fmpq_mat(r, 1)
    # compute compressed columns: for each pivot column, vector = Rm @ colvec
    for jj, (P, den) in enumerate(colint):
        idx = np.array([rows[e] for e in P], dtype=np.int64)
        vals = [P[e] for e in P]
        # Rm[:, idx] @ vals  -> python-int exact
        sub_ = Rm[:, idx]
        acc = [0] * r
        for t, v in enumerate(vals):
            colmask = sub_[:, t]
            nz = np.nonzero(colmask)[0]
            for i in nz:
                acc[i] += v
        for i in range(r):
            if acc[i]:
                Mq[i, jj] = fmpq(acc[i], den)
    idx = np.array([rows[e] for e in rhsint], dtype=np.int64)
    vals = [rhsint[e] for e in rhsint]
    sub_ = Rm[:, idx]
    acc = [0] * r
    for t, v in enumerate(vals):
        nz = np.nonzero(sub_[:, t])[0]
        for i in nz:
            acc[i] += v
    for i in range(r):
        if acc[i]:
            bq[i, 0] = fmpq(acc[i], rden)
    t0 = time.time()
    try:
        xq = Mq.solve(bq)
    except Exception as ex:  # singular compression; caller may retry with another seed
        if verbose:
            print(f'[exact] solve failed: {ex}', file=sys.stderr)
        return None
    if verbose:
        print(f'[exact] fmpq solve {r}x{r} in {time.time()-t0:.1f}s', file=sys.stderr)
    N = [{} for _ in range(NV)]; Qp = {}
    for jj, j in enumerate(piv):
        v = xq[jj, 0]
        v = _Fr(int(v.p), int(v.q))
        if v == 0: continue
        kind, i, m = cols[j]
        if kind == 'N': N[i][m] = v
        else: Qp[m] = v
    ok, info = core.verify_witness(F, A, q, D, N, Qp)
    if not ok:
        if verbose:
            print('[exact] verification FAILED', info, file=sys.stderr)
        return None
    return N, Qp


def solve_fixed_D_fast(F, A, q, D, p=None, verbose=False, want_exact=True, method='auto', sparse_kw=None):
    """Modular consistency test, then exact solve if consistent.

    method: 'auto' (sparse elimination with a fill limit of 25x the initial nnz, the dense finish absorbs the
    rest; compressed rows only if the sparse path raises), 'sparse', or 'compressed' (the pre-v3 behaviour).
    The result carries 'solver' (which path produced the modular answer) and, for sparse, its statistics."""
    p = p or core.PRIMES[0]
    S = FixedDSystem(F, A, q, D, p, verbose=verbose)
    r = None; solver = None
    if method in ('auto', 'sparse'):
        try:
            kw = dict(sparse_kw or {})
            if method == 'auto':
                kw.setdefault('max_nnz', 25 * max(1, sum(len(c) for c in S.colpolys)))
            r = solve_sparse(S, verbose=verbose, **kw); solver = 'sparse'
            if r['consistent'] is None:
                raise RuntimeError(r.get('note', 'sparse solve unverified'))
        except Exception as ex:
            if method == 'sparse':
                raise
            if verbose:
                print(f'[fast] sparse path failed ({ex.__class__.__name__}: {ex}); using compressed rows', file=sys.stderr)
            r = None
    if r is None:
        r = solve_compressed(S, verbose=verbose); solver = 'compressed'
    out = {'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'],
           'rows': S.nrows, 'cols': S.ncols, 'd': S.d, 'solver': solver}
    if solver == 'sparse':
        out['sparse_stats'] = {k: r.get(k) for k in ('status_sparse', 'npiv_sparse', 'active_rows', 'active_cols',
                                                     'active_nnz', 'fill_events', 't_sparse', 't_dense', 't_total')}
    if not r['consistent'] or not want_exact:
        return out
    for seed in range(3):
        sol = exact_solution_from_pivots(F, A, q, D, r['pivots'], verbose=verbose, seed=seed + 7)
        if sol is not None:
            out['N'], out['Q'] = sol
            out['exact'] = True
            return out
    out['exact'] = False
    return out
