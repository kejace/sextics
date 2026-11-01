"""Jacobian-ideal linear algebra mod p for the hodge.computer certificate problem.

J(F) = (F_0, ..., F_5) with F_i = dF/dz_i.  Everything here is finite linear
algebra in a fixed degree over GF(p) (python-flint nmod_mat, dense).

Caveat: membership / kernel dimensions computed mod p agree with the
characteristic-zero answer for all but finitely many primes; with p ~ 2^31 a
discrepancy is extremely unlikely but a modular result is not a proof.

Public API
----------
jacobian_generators(F, deg, p, parity=None)  -> list of (i, m, poly)   # m*F_i, deg m = deg-5
jacobian_matrix(F, deg, p, parity=None)      -> (nmod_mat J (rows=monomials), monomial list)
in_jacobian(F, P, p=None)                    -> (bool, info)
divisible_by(F, P, p=None)                   -> (bool, G or None)        # P = F*G ?
jacobian_ring_basis(F, deg, p=None, parity=None) -> (standard monomials, reducer)
normal_form(F, P, p=None)                    -> poly (representative modulo J(F)_deg)
mult_kernel(F, A, d, p=None, parity=0)       -> dict(dim, basis, ring_dim, jacobian_dim, kernel_dim, p)
"""
from __future__ import annotations
import sys, time
from . import core
from .core import NV, monomials, shift, mul, add, sub, scale, pderiv, reduce_mod, yparity

try:
    from flint import nmod_mat
except ImportError:  # pragma: no cover
    nmod_mat = None


def _p(p):
    return p or core.PRIMES[0]


def jacobian_generators(F, deg, p, parity=None):
    """Spanning set of J(F)_deg: all m * F_i with deg m = deg - 5 (mod p).

    parity: restrict to products of the given Y-parity (F is Y-even, so F_i is
    Y-even for i < 3 and Y-odd for i >= 3)."""
    Fp = reduce_mod(F, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    out = []
    if deg < 5:
        return out
    for i in range(NV):
        if not Fi[i]:
            continue
        par_i = 0 if i < 3 else 1
        mpar = None if parity is None else (parity - par_i) % 2
        for m in monomials(deg - 5, mpar):
            out.append((i, m, shift(Fi[i], m, 1, p)))
    return out


def _vectorize(polys, mons=None):
    """Return (row index dict, list) for the union of supports (or the given monomial list)."""
    if mons is None:
        seen = {}
        for P in polys:
            for e in P:
                if e not in seen:
                    seen[e] = len(seen)
        return seen
    return {e: i for i, e in enumerate(mons)}


def jacobian_matrix(F, deg, p, parity=None, mons=None):
    """Matrix of J(F)_deg generators as ROWS over the monomial basis (columns).

    Returns (nmod_mat, monomial list).  If mons is None, all monomials of the
    degree (and parity, if given) are used."""
    gens = jacobian_generators(F, deg, p, parity)
    if mons is None:
        mons = monomials(deg, parity)
    idx = {e: i for i, e in enumerate(mons)}
    R, C = len(gens), len(mons)
    flat = [0] * (R * C)
    for r, (_, _, P) in enumerate(gens):
        base = r * C
        for e, v in P.items():
            flat[base + idx[e]] = v
    return nmod_mat(R, C, flat, p), mons


def in_jacobian(F, P, p=None):
    """Is the homogeneous polynomial P in J(F) (mod p)?  Returns (bool, info)."""
    p = _p(p)
    Pp = reduce_mod(P, p)
    if not Pp:
        return True, {'note': 'zero'}
    if not core.is_homogeneous(P):
        raise ValueError('P must be homogeneous')
    deg = core.degree(P)
    par = core.parity_of(Pp)
    gens = jacobian_generators(F, deg, p, par)
    rows = _vectorize([g[2] for g in gens] + [Pp])
    R, C = len(rows), len(gens)
    flat = [0] * (R * (C + 1))
    for j, (_, _, G) in enumerate(gens):
        for e, v in G.items():
            flat[rows[e] * (C + 1) + j] = v
    for e, v in Pp.items():
        flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    rank_aug = M.rank()
    M0 = nmod_mat(R, C, [flat[r * (C + 1) + j] for r in range(R) for j in range(C)], p)
    rank0 = M0.rank()
    return rank_aug == rank0, {'rank': rank0, 'rank_aug': rank_aug, 'rows': R, 'cols': C, 'p': p}


def divisible_by(F, P, p=None):
    """Decide whether P = F * G for some G (mod p).  Returns (bool, G mod p or None)."""
    p = _p(p)
    Fp = reduce_mod(F, p); Pp = reduce_mod(P, p)
    if not Pp:
        return True, {}
    dG = core.degree(P) - core.degree(F)
    if dG < 0:
        return False, None
    parG = None
    pP, pF = core.parity_of(Pp), core.parity_of(Fp)
    if pP is not None and pF is not None:
        parG = (pP - pF) % 2
    gm = monomials(dG, parG)
    cols = [shift(Fp, m, 1, p) for m in gm]
    rows = _vectorize(cols + [Pp])
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * (C + 1) + j] = v
    for e, v in Pp.items():
        flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    Mr, rank = M.rref()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0:
            j += 1
        pivots.append(j)
    if pivots and pivots[-1] == C:
        return False, None
    G = {}
    for i, j in enumerate(pivots):
        v = int(Mr[i, C])
        if v:
            G[gm[j]] = v
    return True, G


class _Reducer:
    """Row-reduced J(F)_deg: normal forms modulo the Jacobian ideal in one degree."""

    def __init__(self, F, deg, p, parity):
        self.p, self.deg, self.parity = p, deg, parity
        mons = monomials(deg, parity)
        self.mons = mons
        self.idx = {e: i for i, e in enumerate(mons)}
        J, _ = jacobian_matrix(F, deg, p, parity, mons)
        if J.nrows() == 0:
            self.rank = 0; self.pivots = []; self.rows = []
        else:
            Jr, rank = J.rref()
            self.rank = rank
            pivots = []
            rowsl = []
            C = len(mons)
            for i in range(rank):
                j = 0
                while j < C and int(Jr[i, j]) == 0:
                    j += 1
                pivots.append(j)
                rowsl.append({mons[k]: int(Jr[i, k]) for k in range(j, C) if int(Jr[i, k])})
            self.pivots = pivots
            self.rows = rowsl
        self.pivot_set = {mons[j] for j in self.pivots}
        self.standard = [m for j, m in enumerate(mons) if j not in set(self.pivots)]
        self.jacobian_dim = self.rank
        self.ring_dim = len(mons) - self.rank

    def normal_form(self, P):
        """Reduce P (mod p, same degree/parity) to a combination of standard monomials."""
        P = dict(reduce_mod(P, self.p))
        p = self.p
        for piv, row in zip(self.pivots, self.rows):
            m = self.mons[piv]
            c = P.get(m, 0)
            if c:
                P = sub(P, scale(row, c, p), p)
        return P


def jacobian_ring_basis(F, deg, p=None, parity=None):
    """Standard monomials of R_deg / J(F)_deg (mod p) and the reducer object."""
    p = _p(p)
    red = _Reducer(F, deg, p, parity)
    return red.standard, red


def normal_form(F, P, p=None):
    """Normal form of a homogeneous P modulo J(F) in its degree (mod p)."""
    p = _p(p)
    Pp = reduce_mod(P, p)
    if not Pp:
        return {}
    red = _Reducer(F, core.degree(P), p, core.parity_of(Pp))
    return red.normal_form(Pp)


def mult_kernel(F, A, d, p=None, parity=0, verbose=False, return_reducer=False):
    """Kernel of multiplication by [A] : (R/J)_d -> (R/J)_{d+deg A}  on Y-parity `parity`.

    Returns dict with
      dim          -- dimension of the kernel in the Jacobian ring
      basis        -- list of polys (mod p) that are representatives in standard monomials
      kernel_dim   -- dimension of {D in R_d : A D in J} (before quotienting by J_d)
      jacobian_dim -- dim J(F)_d (this parity)
      ring_dim     -- dim (R/J)_d (this parity)
    Note: deg D = d, so A*D must be in J(F)_{d + deg A}; the certificate necessary
    condition e1 (q N(F) = A D) is exactly A*D in J.
    """
    p = _p(p)
    t0 = time.time()
    Ap = reduce_mod(A, p)
    degA = core.degree(A)
    pA = core.parity_of(Ap)
    if pA is None:
        raise ValueError('A has mixed Y-parity')
    Dm = monomials(d, parity)
    K = len(Dm)
    target_par = (pA + parity) % 2
    gens = jacobian_generators(F, d + degA, p, target_par)
    cols = [shift(Ap, m, 1, p) for m in Dm] + [g[2] for g in gens]
    rows = _vectorize(cols)
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    t1 = time.time()
    X, nullity = M.nullspace()
    t2 = time.time()
    # project nullspace to the D coordinates and take a basis of the projection
    if nullity == 0:
        vecs = []
    else:
        proj = nmod_mat(nullity, K, [int(X[j, k]) for k in range(nullity) for j in range(K)], p)
        Pr, rk = proj.rref()
        vecs = []
        for i in range(rk):
            v = {Dm[j]: int(Pr[i, j]) for j in range(K) if int(Pr[i, j])}
            vecs.append(v)
    kernel_dim = len(vecs)
    # quotient by J(F)_d
    red = _Reducer(F, d, p, parity)
    reduced = [red.normal_form(v) for v in vecs]
    reduced = [v for v in reduced if v]
    basis = []
    if reduced:
        sm = red.standard
        sidx = {m: i for i, m in enumerate(sm)}
        Q = nmod_mat(len(reduced), len(sm), [0] * (len(reduced) * len(sm)), p)
        flat = [0] * (len(reduced) * len(sm))
        for r, v in enumerate(reduced):
            for e, c in v.items():
                flat[r * len(sm) + sidx[e]] = c
        Q = nmod_mat(len(reduced), len(sm), flat, p)
        Qr, rk = Q.rref()
        for i in range(rk):
            basis.append({sm[j]: int(Qr[i, j]) for j in range(len(sm)) if int(Qr[i, j])})
    out = {'dim': len(basis), 'basis': basis, 'kernel_dim': kernel_dim,
           'jacobian_dim': red.jacobian_dim, 'ring_dim': red.ring_dim, 'p': p,
           'rows': R, 'cols': C, 't_build': round(t1 - t0, 2), 't_null': round(t2 - t1, 2),
           't_total': round(time.time() - t0, 2)}
    if return_reducer:
        out['reducer'] = red
    if verbose:
        print(f'[mult_kernel] d={d} par={parity} rows={R} cols={C} kernel={kernel_dim} '
              f'J_d={red.jacobian_dim} (R/J)_d={red.ring_dim} ker[A]={len(basis)} '
              f'({out["t_total"]}s)', file=sys.stderr)
    return out


def symmetric(P, p):
    """Symmetric representatives in (-p/2, p/2] for printing."""
    return {e: (v if v <= p // 2 else v - p) for e, v in P.items()}
