"""Numerical joint search for a denominator in a linear family (variable projection).

Given a family D(a) = sum_k a_k M_k of homogeneous Y-even polynomials of degree d
and the instance (F, A, q), the fixed-D certificate system is

    M(a) x = b(a),   M(a) = sum_k a_k M_k  (N-columns, linear in a) | M_Q (Q-columns, constant),
                     b(a) = -A D(a)^2 = sum_{k,l} a_k a_l v_kl,  v_kl = -A M_k M_l,

with x the coefficient vector of (N_0..N_5, Q).  The reduced functional
    f(a) = min_x || M(a) x - b(a) ||^2
is evaluated by dense least squares (float64) and its gradient by the standard
variable-projection formula  df/da_k = 2 r . (M_k x_N - db/da_k),  r = M x - b,
which is exact at the least-squares optimum (Golub-Pereyra; the x-derivative term
vanishes).  Since f(t a) = t^4 f(a), we minimize over the unit sphere, a = u/|u|.

Spurious minima: if F | D then at q = 1 the identity holds trivially (N = 0,
Q = A G^2).  With exclude_F_multiples=True the parameter space is restricted to
the orthogonal complement of {a : D(a) in F * R_{d-6}} and is_spurious() tests
numerically for divisibility by F.

After a small residual is found, reconstruct() proposes rational vectors
(continued fractions / limit_denominator and integer scaling) and confirm()
re-solves the fixed-D system exactly (fixed_d.solve_fixed_D_fast with fallback
to solve_fixed_D) -- only the exact solve is evidence.

Public API
----------
VarProFamily(F, A, q, basis, exclude_F_multiples=True, verbose=False)
    .K, .d, .rows, .ncols, .nN, .nQ
    .value(a)                     -> (f, x, r)           at a/|a|
    .value_and_grad(u)            -> (f, grad_u)         for the optimizer
    .search(n_restarts=10, seed=0, maxiter=300, ...) -> list of dict(value, rel, a, x, spurious)
    .is_spurious(a)               -> (bool, rel_residual_of_division)
    .reconstruct(a, max_den=200, tol=1e-6) -> list of Fraction vectors
    .D_exact(frac_vector)         -> poly with Fraction coefficients
    .confirm(frac_vector, verbose=False) -> (status, data)   # exact re-solve
search_and_confirm(F, A, q, basis, ...) -> dict summary (best residuals, confirmed witnesses)
monomial_family(d, parity=0, pred=None) -> list of monomial polys
"""
from __future__ import annotations
import math, sys, time
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from . import core, fixed_d
from .core import NV, monomials, shift, mul, add, sub, scale, pderiv, power


def monomial_family(d, parity=0, pred=None):
    """All monomials of degree d and Y-parity `parity` as single-term polys."""
    return [{m: Fraction(1)} for m in monomials(d, parity) if pred is None or pred(m)]


def family_from_exprs(exprs):
    return [core.parse(e) for e in exprs]


def _to_float(P):
    return {e: float(v) for e, v in P.items()}


class VarProFamily:
    def __init__(self, F, A, q, basis, exclude_F_multiples=True, verbose=False, cond=1e-9):
        t0 = time.time()
        # relative singular-value cutoff for the least-squares solve: the N-map has a
        # structural kernel (rank-deficient M(a)), and a fixed cutoff keeps f(a) smooth.
        self.cond = cond
        self.F, self.A, self.q = F, A, q
        self.basis = [dict(B) for B in basis]
        K = len(basis)
        self.K = K
        degs = {core.degree(B) for B in basis}
        if len(degs) != 1:
            raise ValueError('family members must share one degree')
        d = degs.pop(); self.d = d
        pars = {core.parity_of(B) for B in basis}
        if len(pars) != 1 or None in pars:
            raise ValueError('family members must share one pure Y-parity')
        self.pD = pars.pop()
        D0 = self.basis[0]
        self.pars = fixed_d.numerator_parities(D0, A)
        degN = d + 6 * q - 5; degQ = 2 * d - 6
        self.degN, self.degQ = degN, degQ
        Ncols = [(i, m) for i in range(NV) for m in monomials(degN, self.pars[i])]
        Qcols = [m for m in monomials(degQ, core.parity_of(A))] if degQ >= 0 else []
        self.Ncols, self.Qcols = Ncols, Qcols
        self.nN, self.nQ = len(Ncols), len(Qcols)
        Fi = [pderiv(F, i) for i in range(NV)]
        Fq1 = scale(power(F, q + 1), -1)
        rows = {}

        def ridx(e):
            if e not in rows:
                rows[e] = len(rows)
            return rows[e]
        # per-basis N blocks
        blocks = []
        for Bk in self.basis:
            FBk = mul(F, Bk)
            Bki = [pderiv(Bk, i) for i in range(NV)]
            G = [scale(add(scale(mul(Bk, Fi[i]), q), mul(F, Bki[i])), -1) for i in range(NV)]
            ri, ci, vals = [], [], []
            for j, (i, m) in enumerate(Ncols):
                col = shift(G[i], m, 1)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(FBk, tuple(m2), m[i]))
                for e, v in col.items():
                    ri.append(ridx(e)); ci.append(j); vals.append(float(v))
            blocks.append((ri, ci, vals))
        # constant Q block
        qri, qci, qvals = [], [], []
        for j, m in enumerate(Qcols):
            col = shift(Fq1, m, 1)
            for e, v in col.items():
                qri.append(ridx(e)); qci.append(j); qvals.append(float(v))
        # rhs pair vectors v_kl = -A M_k M_l
        AB = [mul(A, Bk) for Bk in self.basis]
        pair = {}
        for k in range(K):
            for l in range(k, K):
                v = scale(mul(AB[k], self.basis[l]), -1)
                pair[(k, l)] = v
                for e in v:
                    ridx(e)
        R = len(rows)
        self.rows = R
        self.row_index = rows
        self.Mk = [sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                 shape=(R, self.nN)) for ri, ci, vals in blocks]
        self.MQ = sp.csr_matrix((np.array(qvals), (np.array(qri, dtype=np.int64), np.array(qci, dtype=np.int64))),
                                shape=(R, self.nQ)).toarray() if self.nQ else np.zeros((R, 0))
        # Vk: rows x K, column l = v_{kl}
        self.Vk = []
        for k in range(K):
            ri, ci, vals = [], [], []
            for l in range(K):
                v = pair[(min(k, l), max(k, l))]
                for e, c in v.items():
                    ri.append(rows[e]); ci.append(l); vals.append(float(c))
            self.Vk.append(sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                         shape=(R, K)))
        self.ncols = self.nN + self.nQ
        # parameter restriction (exclude a with F | D(a))
        self.P = np.eye(K)
        self.excluded_dim = 0
        if exclude_F_multiples and d >= 6:
            self.P, self.excluded_dim = self._complement_of_F_multiples()
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[varpro] d={d} q={q} K={K} rows={R} cols={self.ncols} (N {self.nN}, Q {self.nQ}) '
                  f'excluded F-multiples dim={self.excluded_dim} built in {self.build_time:.1f}s', file=sys.stderr)

    # ------------------------------------------------------------------
    def _family_matrix(self):
        mons = {}
        for B in self.basis:
            for e in B:
                mons.setdefault(e, len(mons))
        return mons

    def _complement_of_F_multiples(self):
        """Orthonormal basis of {a}^perp where a ranges over the subspace with D(a) in F*R_{d-6}."""
        K = self.K
        dG = self.d - 6
        Gm = monomials(dG, (self.pD - core.parity_of(self.F)) % 2)
        Fcols = [shift(self.F, m, 1) for m in Gm]
        mons = {}
        for B in self.basis + Fcols:
            for e in B:
                mons.setdefault(e, len(mons))
        Mb = np.zeros((len(mons), K + len(Fcols)))
        for k, B in enumerate(self.basis):
            for e, v in B.items():
                Mb[mons[e], k] = float(v)
        for j, C in enumerate(Fcols):
            for e, v in C.items():
                Mb[mons[e], K + j] = -float(v)
        ns = sla.null_space(Mb)
        S = ns[:K, :] if ns.size else np.zeros((K, 0))
        if S.shape[1] == 0:
            return np.eye(K), 0
        # subspace of a's spanned by S (columns); complement via SVD
        U, s, _ = np.linalg.svd(S, full_matrices=True)
        r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
        return U[:, r:], r

    # ------------------------------------------------------------------
    def assemble(self, a):
        M = self.Mk[0] * a[0]
        for k in range(1, self.K):
            M = M + self.Mk[k] * a[k]
        M = M.toarray()
        if self.nQ:
            M = np.hstack([M, self.MQ])
        b = np.zeros(self.rows)
        for k in range(self.K):
            b += a[k] * (self.Vk[k] @ a)
        return M, b

    def value(self, a):
        """f(a/|a|), x, r at the normalized point."""
        a = np.asarray(a, dtype=float)
        na = np.linalg.norm(a)
        if na == 0:
            raise ValueError('a = 0')
        a = a / na
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        self._nfev += 1
        return float(r @ r), x, r

    def value_and_grad(self, u):
        """Objective and gradient in reduced coordinates u (a = P u / |P u|)."""
        u = np.asarray(u, dtype=float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        f = float(r @ r)
        xN = x[:self.nN]
        g = np.zeros(self.K)
        for k in range(self.K):
            dM = self.Mk[k] @ xN
            db = 2.0 * (self.Vk[k] @ a)
            g[k] = 2.0 * float(r @ (dM - db))
        # chain rule through normalization a = a_full/|a_full|
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def bnorm(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        _, b = self.assemble(a)
        return float(np.linalg.norm(b))

    # ------------------------------------------------------------------
    def is_spurious(self, a, tol=1e-8):
        """Numerically: is D(a) divisible by F?  Returns (bool, relative residual)."""
        if self.d < 6:
            return False, 1.0
        a = np.asarray(a, float)
        D = {}
        for k, B in enumerate(self.basis):
            for e, v in B.items():
                D[e] = D.get(e, 0.0) + a[k] * float(v)
        dG = self.d - 6
        Gm = monomials(dG, (self.pD - core.parity_of(self.F)) % 2)
        Fcols = [shift(self.F, m, 1) for m in Gm]
        mons = {}
        for C in Fcols + [D]:
            for e in C:
                mons.setdefault(e, len(mons))
        Mb = np.zeros((len(mons), len(Fcols)))
        for j, C in enumerate(Fcols):
            for e, v in C.items():
                Mb[mons[e], j] = float(v)
        bv = np.zeros(len(mons))
        for e, v in D.items():
            bv[mons[e]] = v
        x, *_ = sla.lstsq(Mb, bv, check_finite=False)
        rel = float(np.linalg.norm(Mb @ x - bv) / max(np.linalg.norm(bv), 1e-300))
        return rel < tol, rel

    def search(self, n_restarts=10, seed=0, maxiter=300, gtol=1e-10, verbose=False, init=None,
               stop_rel=1e-12):
        """Random-restart L-BFGS minimization on the sphere.  Returns candidates sorted by residual."""
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-16})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            f, x, r = self.value(a)
            bn = self.bnorm(a)
            rel = f / max(bn * bn, 1e-300)
            spur, srel = self.is_spurious(a)
            cand = {'value': f, 'rel': rel, 'a': a, 'x': x, 'spurious': spur, 'div_rel': srel,
                    'iters': int(res.nit), 'restart': it}
            out.append(cand)
            if verbose:
                print(f'[varpro] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} spurious={spur} '
                      f'({time.time()-t0:.1f}s)', file=sys.stderr)
            if rel < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    # ------------------------------------------------------------------
    def reconstruct(self, a, max_den=200, tol=1e-6):
        """Rational candidates for a (normalized so that the largest |a_k| = 1)."""
        a = np.asarray(a, float)
        j = int(np.argmax(np.abs(a)))
        a = a / a[j]
        cands = []
        # continued fractions per coordinate
        fr = [Fraction(float(v)).limit_denominator(max_den) for v in a]
        if max(abs(float(f) - v) for f, v in zip(fr, a)) < tol:
            cands.append(tuple(fr))
        # common denominator scaling
        for den in range(1, max_den + 1):
            v = a * den
            if np.max(np.abs(v - np.round(v))) < tol * den:
                cands.append(tuple(Fraction(int(round(x)), den) for x in v))
                break
        # dedupe
        uniq = []
        for c in cands:
            if c not in uniq:
                uniq.append(c)
        return uniq

    def D_exact(self, frac):
        D = {}
        for k, B in enumerate(self.basis):
            if frac[k] == 0:
                continue
            D = add(D, scale(B, Fraction(frac[k])))
        return D

    def confirm(self, frac, verbose=False):
        """Exact fixed-D solve for the rational candidate.  Returns (status, data)."""
        D = self.D_exact(frac)
        if not D:
            return 'inconsistent', {'note': 'D = 0'}
        ok, info = fixed_d.e1_solvable(self.F, self.A, self.q, D)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails', 'e1': info, 'D': D}
        try:
            out = fixed_d.solve_fixed_D_fast(self.F, self.A, self.q, D, verbose=verbose)
            if not out['consistent']:
                return 'inconsistent', {'rank': out['rank'], 'nullity': out['nullity'], 'D': D}
            if out.get('exact'):
                return 'solved', {'N': out['N'], 'Q': out['Q'], 'rank': out['rank'],
                                  'nullity': out['nullity'], 'D': D}
        except Exception as ex:  # fall back to the dense solver
            if verbose:
                print(f'[varpro] fast solver failed ({ex}); dense fallback', file=sys.stderr)
        status, data = fixed_d.solve_fixed_D(self.F, self.A, self.q, D, verbose=verbose)
        data['D'] = D
        return status, data


def search_and_confirm(F, A, q, basis, n_restarts=10, seed=0, rel_tol=1e-8, max_den=200,
                       verbose=False, exclude_F_multiples=True, maxiter=300):
    """Convenience driver: numerical search, rational reconstruction, exact confirmation."""
    fam = VarProFamily(F, A, q, basis, exclude_F_multiples=exclude_F_multiples, verbose=verbose)
    cands = fam.search(n_restarts=n_restarts, seed=seed, verbose=verbose, maxiter=maxiter)
    confirmed = []
    tried = []
    for c in cands:
        if c['spurious'] or c['rel'] > rel_tol:
            continue
        for fr in fam.reconstruct(c['a'], max_den=max_den):
            Dstr = core.to_str(fam.D_exact(fr))
            if Dstr in tried:
                continue
            tried.append(Dstr)
            status, data = fam.confirm(fr, verbose=verbose)
            if verbose:
                print(f'[varpro] exact confirm of {Dstr[:80]} -> {status}', file=sys.stderr)
            if status == 'solved':
                confirmed.append({'D': data['D'], 'N': data['N'], 'Q': data['Q'], 'a': fr,
                                  'rel': c['rel']})
                break
        if confirmed:
            break
    return {'family': fam, 'candidates': cands, 'confirmed': confirmed, 'tried': tried,
            'best_rel': cands[0]['rel'] if cands else None}
