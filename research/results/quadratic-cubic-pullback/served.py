"""Fast served-target space for a fixed denominator, selectable target parity, mod p, compressed rows.

served_space(F, q, D, parity, p) -> (dim, basis polys (mod p), S)
 Nullspace of the compressed homogeneous system [N | Q | A] projected onto the A-coordinates.
 Compression with k = ncols + margin random row combinations preserves the nullspace w.h.p.
"""
import sys, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from hodge.core import mul, add, scale, shift, monomials
import numpy as np
from flint import nmod_mat


def served_space(F, q, D, parity=1, p=None, margin=60, seed=3, verbose=False, mod_J=True):
    p = p or core.PRIMES[0]
    A0 = {monomials(6 * q, parity)[0]: 1}
    t0 = time.time()
    S = fixed_d.FixedDSystem(F, A0, q, D, p, with_target_space=False, verbose=verbose)
    Dp = core.reduce_mod(D, p)
    D2 = mul(Dp, Dp, p)
    amons = monomials(6 * q, parity)
    for m in amons:
        S.cols.append(('A', -1, m)); S.colpolys.append(shift(D2, m, 1, p))
    for col in S.colpolys[S.ncols:]:
        for e in col:
            if e not in S.rows:
                S.rows[e] = len(S.rows)
    S.nrows = len(S.rows); S.ncols = len(S.cols)
    M = fixed_d._sparse_matrix(S, with_rhs=False)
    C = S.ncols
    k = min(S.nrows, C + margin)
    Cm = fixed_d._compress(M, p, k, seed)
    t1 = time.time()
    Mn = nmod_mat(k, C, Cm.ravel().tolist(), p)
    del Cm
    X, nullity = Mn.nullspace()
    t2 = time.time()
    a_idx = list(range(C - len(amons), C))
    if nullity == 0:
        return 0, [], S
    proj = nmod_mat(len(a_idx), nullity, [int(X[j, kk]) for j in a_idx for kk in range(nullity)], p)
    T = proj.transpose()
    Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, m in enumerate(amons):
            v = int(Tr[i, jj])
            if v:
                poly[m] = v
        basis.append(poly)
    if verbose:
        print(f'[served parity={parity}] rows={S.nrows} cols={C} k={k} nullity={nullity} served dim={rk} '
              f'(build+compress {t1-t0:.1f}s, nullspace {t2-t1:.1f}s)', file=sys.stderr)
    return rk, basis, S


def jacobian_part(F, parity, p=None):
    """Basis (mod p) of J(F)_6 of given Y-parity: Z_j * dF/dZ_i."""
    p = p or core.PRIMES[0]
    Fp = core.reduce_mod(F, p)
    out = []
    for i in range(6):
        Fi = core.pderiv(Fp, i, p)
        for j in range(6):
            e = [0]*6; e[j] = 1
            g = shift(Fi, tuple(e), 1, p)
            if core.parity_of(g) == parity:
                out.append(g)
    return out


def rank_of(polys, mons, p):
    col = {e: j for j, e in enumerate(mons)}
    if not polys:
        return 0
    flat = [0] * (len(polys) * len(mons))
    for r, P in enumerate(polys):
        for e, v in P.items():
            flat[r * len(mons) + col[e]] = v % p
    return nmod_mat(len(polys), len(mons), flat, p).rank()


def analyse(F, q, D, parity, extra=None, p=None, verbose=True, margin=60):
    """served dim, dim mod J, and whether each poly in `extra` (dict name->poly) is served."""
    p = p or core.PRIMES[0]
    dim, basis, S = served_space(F, q, D, parity, p, verbose=verbose, margin=margin)
    mons = monomials(6 * q, parity)
    J = jacobian_part(F, parity, p)
    rJ = rank_of(J, mons, p)
    rU = rank_of(J + basis, mons, p)
    res = {'served': dim, 'J': rJ, 'served_mod_J': rU - rJ}
    for nm, P in (extra or {}).items():
        Pp = core.reduce_mod(P, p)
        res['serves_' + nm] = (rank_of(basis + [Pp], mons, p) == dim)
        res['serves_' + nm + '_modJ'] = (rank_of(J + basis + [Pp], mons, p) == rU)
    return res, basis
