#!/usr/bin/env python3
"""w2s.py -- wave-2 SEARCH library for quadratic-cubic-pullback (q = 1).

Reduced certificate system
--------------------------
For a denominator D (degree d, Y-even) that passes E1, i.e. A*D = sum_i u_i F_i with deg u_i = d+1
(the "E1 lift"; F_i = dF/dZ_i), every numerator vector of a certificate has the form
    N_i = u_i + sum_j c_ij F_j + G Z_i / 6,     c_ij = -c_ji of degree d-4,  G of degree d,
(syzygies of the regular sequence F_0..F_5 are Koszul; F = sum Z_i F_i / 6).  Writing the cleared identity
    E = A D^2 - D sum N_i F_i + F (D div N - N.grad D) = F^2 Q
with sum N_i F_i = A D + G F gives  E = F (D div N - N.grad D - D G), and the G-terms cancel
(D div(G Z/6) - (G Z/6).grad D = D G).  Hence a certificate exists iff there is an antisymmetric c with
    P(D, c) := D div(u + K(c)) - (u + K(c)).grad D  ==  0  (mod F),      K(c)_i = sum_j c_ij F_j,
a polynomial identity of degree 2d (instead of 2d+6) with 15 blocks of c-unknowns (instead of 6 N blocks + Q).

Point form (numerical search).  P in (F)  <=>  P vanishes on X = V(F).  F has no real zeros (min on the real
sphere ~ 0.004), so we sample complex points z_s of X, |z_s| = 1, and minimise the variable-projection functional
    f(a) = min_c sum_s |P(D(a), c)(z_s)|^2     over the sphere |a| = 1,
for a linear family D(a) = sum_k a_k M_k with E1 lifts u(a) = sum_k a_k u^(k).  For fixed a this is a complex
linear least-squares problem in the real unknowns c (2S real equations); the gradient in a is the Golub-Pereyra
formula (exact at the least-squares optimum).  With S >= (n_c + K)/2 + margin the sketched problem has the same
zero set as the full one for generic points (confirmed anyway by the EXACT solver below).

Spurious zeros: D in F * S_{d-6} (then N = 0, Q = A G^2) -- excluded by projecting the a-space onto the
orthogonal complement of the F-multiples (as in hodge.varpro).

Exact confirmation: `exact_reduced_solve` builds the same reduced system over GF(p) with Q-unknowns kept
(columns -F m) and runs nmod_mat.rref; `e1_lift` computes an exact lift u mod p.  Rational lifting of a
consistent candidate is left to hodge.fixed_d (CRT) or to the caller.

Families
--------
jacobian_family(d)        : D = sum_i w_i F_i, basis m*F_i, exact lifts u = A m e_i  (K = dim J_d^even for d <= 9)
numerical_kernel_family   : float basis of K_d = {D : A D in J} from an SVD nullspace of [A*m | J-generators]
                            (used at d = 10 where K_10/J_10 has dimension 294)

All polynomials are hodge.core dicts {exponent tuple: Fraction}.  Variables X0,X1,X2,Y0,Y1,Y2.
"""
from __future__ import annotations
import sys, os, json, time, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from scipy.optimize import minimize
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv

LANE = 'quadratic-cubic-pullback'
WD = '/tmp/hodge/results/quadratic-cubic-pullback'
BATCHLOG = os.path.join(WD, 'wave2-search.log.jsonl')

F, A, q, REC = core.load_instance(LANE)
assert q == 1
Fi = [pderiv(F, i) for i in range(NV)]
PAR_A = core.parity_of(A)          # 1 (Y-odd)


def configure(F_, A_):
    """Switch the instance polynomials (used for the Fermat control)."""
    global F, A, Fi, PAR_A, _FE, _FC, _FiE, _FiC
    F, A = F_, A_
    Fi = [pderiv(F, i) for i in range(NV)]
    PAR_A = core.parity_of(A)
    _FE = np.array(list(F.keys()), dtype=np.int64)
    _FC = np.array([float(c) for c in F.values()])
    _FiE = [np.array(list(Fi[i].keys()), dtype=np.int64) for i in range(NV)]
    _FiC = [np.array([float(c) for c in Fi[i].values()]) for i in range(NV)]


def log_batch(rec):
    rec = dict(rec)
    rec.setdefault('lane', LANE)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('engine', 'w2s-reduced-v1')
    with open(BATCHLOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec


# ----------------------------------------------------------------------------------------------------------
# monomial bookkeeping
# ----------------------------------------------------------------------------------------------------------
_MON = {}


def mons(deg, par=None):
    """(list of exponent tuples, index dict) for degree deg and Y-parity par (None = all)."""
    key = (deg, par)
    if key not in _MON:
        if deg < 0:
            lst = []
        else:
            lst = list(monomials(deg, par))
        _MON[key] = (lst, {e: i for i, e in enumerate(lst)})
    return _MON[key]


def vec(P, deg, par):
    """Dense float vector of P over monomials(deg, par)."""
    lst, idx = mons(deg, par)
    v = np.zeros(len(lst))
    for e, c in P.items():
        v[idx[e]] = float(c)
    return v


def diff_matrix(deg, par, i):
    """Sparse matrix of d/dZ_i : S_{deg,par} -> S_{deg-1, par'} (par' = par ^ (i>=3))."""
    src, _ = mons(deg, par)
    par2 = None if par is None else (par ^ (1 if i >= 3 else 0))
    dst, didx = mons(deg - 1, par2)
    r, c, v = [], [], []
    for j, e in enumerate(src):
        if e[i]:
            e2 = list(e); e2[i] -= 1
            r.append(didx[tuple(e2)]); c.append(j); v.append(float(e[i]))
    return sp.csr_matrix((v, (r, c)), shape=(len(dst), len(src)))


def numerator_parities():
    """Parities of N_i (hence of u_i) for Y-even D and Y-odd A."""
    return [(PAR_A - 0 - (0 if i < 3 else 1)) % 2 for i in range(NV)]


def c_parities():
    """Parity of c_ij for i<j: par(N_i) - par(F_j)."""
    pN = numerator_parities()
    out = {}
    for i in range(NV):
        for j in range(i + 1, NV):
            out[(i, j)] = (pN[i] - (0 if j < 3 else 1)) % 2
    return out


# ----------------------------------------------------------------------------------------------------------
# points on X(C)
# ----------------------------------------------------------------------------------------------------------
_FE = np.array(list(F.keys()), dtype=np.int64)
_FC = np.array([float(c) for c in F.values()])
_FiE = [np.array(list(Fi[i].keys()), dtype=np.int64) for i in range(NV)]
_FiC = [np.array([float(c) for c in Fi[i].values()]) for i in range(NV)]


def F_at(Z):
    """F at points Z (S x 6 complex)."""
    return (np.prod(Z[:, None, :] ** _FE[None, :, :], axis=2) * _FC[None, :]).sum(axis=1)


def gradF_at(Z):
    return np.stack([(np.prod(Z[:, None, :] ** _FiE[i][None, :, :], axis=2) * _FiC[i][None, :]).sum(axis=1)
                     for i in range(NV)], axis=1)


def sample_points(S, seed=0, newton=6):
    """S complex points of X with |z| = 1: roots of F along random complex lines, Newton-refined."""
    rng = np.random.default_rng(seed)
    pts = []
    while len(pts) < S:
        z0 = rng.standard_normal(6) + 1j * rng.standard_normal(6)
        v = rng.standard_normal(6) + 1j * rng.standard_normal(6)
        ts = np.exp(2j * np.pi * np.arange(7) / 7) * 1.5
        vals = F_at(z0[None, :] + ts[:, None] * v[None, :])
        coef = np.polyfit(ts, vals, 6)
        roots = np.roots(coef)
        Z = z0[None, :] + roots[:, None] * v[None, :]
        for _ in range(newton):
            f = F_at(Z); g = gradF_at(Z)
            gv = (g * v[None, :]).sum(axis=1)
            Z = Z - (f / gv)[:, None] * v[None, :]
        Z = Z / np.linalg.norm(Z, axis=1)[:, None]
        res = np.abs(F_at(Z))
        for k in range(Z.shape[0]):
            if res[k] < 1e-13 and len(pts) < S:
                pts.append(Z[k])
    Z = np.array(pts)
    return Z


def mono_table(Z, deg, par=None):
    """S x n complex table of monomial values z^e for monomials(deg, par)."""
    lst, _ = mons(deg, par)
    if not lst:
        return np.zeros((Z.shape[0], 0), dtype=complex)
    E = np.array(lst, dtype=np.int64)
    maxe = int(E.max()) if E.size else 0
    pw = np.ones((Z.shape[0], NV, maxe + 1), dtype=complex)
    for k in range(1, maxe + 1):
        pw[:, :, k] = pw[:, :, k - 1] * Z
    T = np.ones((Z.shape[0], len(lst)), dtype=complex)
    for i in range(NV):
        T *= pw[:, i, E[:, i]]
    return T


# ----------------------------------------------------------------------------------------------------------
# families
# ----------------------------------------------------------------------------------------------------------
class Family:
    """Linear family of Y-even denominators of degree d with E1 lifts.

    Mmat : (n_{d,0} x K) float matrix of basis polynomials over monomials(d, 0)
    Umat : list of 6 float matrices (n_{d+1, par_i} x K) of lift coefficients, A*M_k = sum_i U_i[:,k] F_i
    labels: names of the basis elements
    exact: optional list of exact dict polys (for rational reconstruction) or None
    """

    def __init__(self, d, Mmat, Umat, labels=None, exact=None, exact_u=None, name=''):
        self.d = d
        self.K = Mmat.shape[1]
        self.Mmat = Mmat
        self.Umat = Umat
        self.labels = labels or [f'b{k}' for k in range(self.K)]
        self.exact = exact
        self.exact_u = exact_u
        self.name = name

    def D_float(self, a):
        return self.Mmat @ a

    def D_poly_float(self, a):
        lst, _ = mons(self.d, 0)
        v = self.Mmat @ a
        return {lst[j]: v[j] for j in range(len(lst)) if abs(v[j]) > 0}


def jacobian_family(d):
    """D = sum_i w_i F_i (Y-even): basis m*F_i with m in S_{d-5} of parity (0 if i<3 else 1); lifts u = A m e_i."""
    lst, idx = mons(d, 0)
    pN = numerator_parities()
    cols, labels, exact, exact_u = [], [], [], []
    ucols = [[] for _ in range(NV)]
    for i in range(NV):
        mpar = 0 if i < 3 else 1
        for m in monomials(d - 5, mpar):
            P = shift(Fi[i], m, 1)
            cols.append(vec(P, d, 0))
            labels.append((i, m))
            exact.append(P)
            Am = shift(A, m, 1)
            u = [{} for _ in range(NV)]
            u[i] = Am
            exact_u.append(u)
            for jj in range(NV):
                ucols[jj].append(vec(u[jj], d + 1, pN[jj]))
    Mmat = np.array(cols).T
    Umat = [np.array(ucols[jj]).T for jj in range(NV)]
    return Family(d, Mmat, Umat, labels, exact, exact_u, name=f'J{d}')


def family_from_polys(d, polys, lifts, labels=None, name=''):
    """Family from exact polys (dicts) and exact lifts (list of 6 dicts each)."""
    pN = numerator_parities()
    Mmat = np.array([vec(P, d, 0) for P in polys]).T
    Umat = [np.array([vec(u[jj], d + 1, pN[jj]) for u in lifts]).T for jj in range(NV)]
    return Family(d, Mmat, Umat, labels, polys, lifts, name=name)


def stacked_matrix_exact_cols(d):
    """Exact columns (dict polys over Fractions) of [A*m | m' F_i] in degree d+6, odd parity, and their labels."""
    Dm, _ = mons(d, 0)
    tpar = (PAR_A + 0) % 2
    cols, labels = [], []
    for m in Dm:
        cols.append(shift(A, m, 1)); labels.append(('D', m))
    for i in range(NV):
        par_i = 0 if i < 3 else 1
        mpar = (tpar - par_i) % 2
        for m in monomials(d + 1, mpar):
            cols.append(shift(Fi[i], m, 1)); labels.append(('u', i, m))
    return cols, labels


# ----------------------------------------------------------------------------------------------------------
# point-evaluation variable projection
# ----------------------------------------------------------------------------------------------------------
class PointVarPro:
    def __init__(self, fam: Family, S, seed=0, exclude_F_multiples=True, ridge=1e-12, verbose=False,
                 Z=None):
        t0 = time.time()
        self.fam = fam
        d = fam.d
        self.d = d
        self.K = fam.K
        self.Z = sample_points(S, seed=seed) if Z is None else Z
        Z = self.Z
        self.S = Z.shape[0]
        pN = numerator_parities()
        self.pN = pN
        # tables of the family at the points: TM (K x S), TdM[i] (K x S), Tu[i] (K x S), Tdivu (K x S)
        Td = mono_table(Z, d, 0)
        self.TM = (Td @ fam.Mmat).T
        self.TdM = []
        for i in range(NV):
            Di = diff_matrix(d, 0, i)
            Tdm1 = mono_table(Z, d - 1, 0 ^ (1 if i >= 3 else 0))
            self.TdM.append((Tdm1 @ (Di @ fam.Mmat)).T)
        self.Tu = []
        Tdivu = np.zeros((self.K, self.S), dtype=complex)
        for i in range(NV):
            Tu_i = mono_table(Z, d + 1, pN[i])
            self.Tu.append((Tu_i @ fam.Umat[i]).T)
            Di = diff_matrix(d + 1, pN[i], i)
            Tdiv = mono_table(Z, d, pN[i] ^ (1 if i >= 3 else 0))
            Tdivu += (Tdiv @ (Di @ fam.Umat[i])).T
        self.Tdivu = Tdivu
        # F_j at the points
        self.TF = gradF_at(Z).T  # 6 x S
        # c blocks
        self.cpar = c_parities()
        self.pairs = [(i, j) for i in range(NV) for j in range(i + 1, NV)]
        dc = d - 4
        self.dc = dc
        self.Tc = {}; self.Tdc = {}
        for par in (0, 1):
            lst, _ = mons(dc, par)
            self.Tc[par] = mono_table(Z, dc, par)
            self.Tdc[par] = []
            for i in range(NV):
                Di = diff_matrix(dc, par, i)
                T1 = mono_table(Z, dc - 1, par ^ (1 if i >= 3 else 0))
                self.Tdc[par].append(T1 @ Di.toarray())  # S x n_c : values of d_i m
        self.cblocks = [(ij, self.cpar[ij], len(mons(dc, self.cpar[ij])[0])) for ij in self.pairs]
        self.nc = sum(n for _, _, n in self.cblocks)
        self.coff = np.cumsum([0] + [n for _, _, n in self.cblocks])
        # structural kernel of c -> K(c): c_ij = sum_k e_ijk F_k (deg e = d-9), antisymmetric e.  Remove it.
        self.cP = None
        if d - 9 >= 0:
            self.cP = self._c_complement()
        # F-multiples exclusion
        self.P = np.eye(self.K)
        self.excluded_dim = 0
        if exclude_F_multiples and d >= 6:
            self.P, self.excluded_dim = self._complement_of_F_multiples()
        self.dimP = self.P.shape[1]
        self.ridge = ridge
        self.nfev = 0
        self.build_time = time.time() - t0
        if verbose:
            print(f'[w2s] family {fam.name} d={d} K={self.K} S={self.S} nc={self.nc} (eff {self.nc_eff}) '
                  f'excludedF={self.excluded_dim} built in {self.build_time:.1f}s', file=sys.stderr)

    @property
    def nc_eff(self):
        return self.nc if self.cP is None else self.cP.shape[1]

    def _c_vec_index(self):
        out = {}
        for b, (ij, par, n) in enumerate(self.cblocks):
            lst, idx = mons(self.dc, par)
            out[ij] = (self.coff[b], idx)
        return out

    def _c_complement(self):
        """Orthonormal basis of the complement of the Koszul kernel {c_ij = sum_k e_ijk F_k}."""
        de = self.d - 9
        cidx = self._c_vec_index()
        vecs = []
        triples = [(i, j, k) for i in range(NV) for j in range(i + 1, NV) for k in range(j + 1, NV)]
        for (i, j, k) in triples:
            # parity of e_ijk: c_ij = e F_k  -> par(e) = par(c_ij) - par(F_k)
            pe = (self.cpar[(i, j)] - (0 if k < 3 else 1)) % 2
            for m in monomials(de, pe):
                v = np.zeros(self.nc)
                # c_ij += e F_k, c_jk += e F_i, c_ik -= e F_j  (antisymmetric e_ijk)
                for (a, b, cc, sgn) in (((i, j), k, None, 1), ((j, k), i, None, 1), ((i, k), j, None, -1)):
                    off, idx = cidx[a]
                    P = shift(Fi[b], m, 1)
                    for e, c in P.items():
                        v[off + idx[e]] += sgn * float(c)
                vecs.append(v)
        if not vecs:
            return None
        Sm = np.array(vecs).T
        U, s, _ = np.linalg.svd(Sm, full_matrices=True)
        r = int(np.sum(s > 1e-9 * s[0]))
        return U[:, r:]

    def _complement_of_F_multiples(self):
        dG = self.d - 6
        Gm, _ = mons(dG, 0)
        lst, idx = mons(self.d, 0)
        Fcols = np.array([vec(shift(F, m, 1), self.d, 0) for m in Gm]).T  # n_d x nG
        # solve Mmat a = F m  for each m (least squares); keep those that are exact
        sol, res, rk, sv = sla.lstsq(self.fam.Mmat, Fcols, check_finite=False)
        good = []
        for j in range(Fcols.shape[1]):
            r = np.linalg.norm(self.fam.Mmat @ sol[:, j] - Fcols[:, j]) / max(np.linalg.norm(Fcols[:, j]), 1e-300)
            if r < 1e-8:
                good.append(sol[:, j])
        if not good:
            return np.eye(self.K), 0
        Sm = np.array(good).T
        U, s, _ = np.linalg.svd(Sm, full_matrices=True)
        r = int(np.sum(s > 1e-9 * s[0]))
        return U[:, r:], r

    # --------------------------------------------------------------------------------------------------
    def _state(self, a):
        D = a @ self.TM                       # S
        dD = [a @ self.TdM[i] for i in range(NV)]
        u = [a @ self.Tu[i] for i in range(NV)]
        divu = a @ self.Tdivu
        return D, dD, u, divu

    def _Ematrix(self, D, dD):
        """S x nc complex matrix of the c-columns."""
        E = np.empty((self.S, self.nc), dtype=complex)
        TF = self.TF
        for b, ((i, j), par, n) in enumerate(self.cblocks):
            o = self.coff[b]
            DFj = (D * TF[j])[:, None]; DFi = (D * TF[i])[:, None]
            w = (TF[j] * dD[i] - TF[i] * dD[j])[:, None]
            E[:, o:o + n] = DFj * self.Tdc[par][i] - DFi * self.Tdc[par][j] - w * self.Tc[par]
        return E

    def solve_c(self, a):
        """Least-squares c for unit a.  Returns (f, c, Pres, state)."""
        D, dD, u, divu = self._state(a)
        p0 = D * divu - sum(u[i] * dD[i] for i in range(NV))
        E = self._Ematrix(D, dD)
        if self.cP is not None:
            E = E @ self.cP
        Er = np.vstack([E.real, E.imag])
        pr = np.concatenate([p0.real, p0.imag])
        G = Er.T @ Er
        dg = np.mean(np.diag(G)) if G.size else 1.0
        G[np.diag_indices_from(G)] += self.ridge * max(dg, 1e-300)
        rhs = -(Er.T @ pr)
        try:
            cf = sla.solve(G, rhs, assume_a='pos', check_finite=False)
        except Exception:
            cf = sla.lstsq(Er, -pr, check_finite=False)[0]
        Pres = E @ cf + p0
        c = cf if self.cP is None else self.cP @ cf
        f = float(np.sum(np.abs(Pres) ** 2))
        self.nfev += 1
        return f, c, Pres, (D, dD, u, divu, p0)

    def value_and_grad(self, uu):
        a_full = self.P @ uu
        na = np.linalg.norm(a_full)
        a = a_full / na
        f, c, Pres, (D, dD, u, divu, p0) = self.solve_c(a)
        # K(c) and h at the points
        g = [np.zeros(self.S, dtype=complex) for _ in range(NV)]
        h = np.zeros(self.S, dtype=complex)
        TF = self.TF
        for b, ((i, j), par, n) in enumerate(self.cblocks):
            o = self.coff[b]
            cv = c[o:o + n]
            cij = self.Tc[par] @ cv
            dci = self.Tdc[par][i] @ cv
            dcj = self.Tdc[par][j] @ cv
            g[i] += cij * TF[j]; g[j] -= cij * TF[i]
            h += TF[j] * dci - TF[i] * dcj
        Pc = np.conj(Pres)
        grad = (self.TM @ (Pc * (divu + h)) + self.Tdivu @ (Pc * D)
                - sum(self.Tu[i] @ (Pc * dD[i]) for i in range(NV))
                - sum(self.TdM[i] @ (Pc * (u[i] + g[i])) for i in range(NV)))
        grad = 2.0 * grad.real
        g_full = (grad - float(grad @ a) * a) / na
        return f, self.P.T @ g_full

    def rel(self, a):
        """Relative residual f / (sum of squared magnitudes of the terms of p0) -- nonzero scale whenever D,u != 0."""
        f, c, Pres, (D, dD, u, divu, p0) = self.solve_c(a)
        scale_ = float(np.sum(np.abs(D * divu) ** 2) + sum(np.sum(np.abs(u[i] * dD[i]) ** 2) for i in range(NV)))
        return f / max(scale_, 1e-300), f, c

    def is_spurious(self, a, tol=1e-8):
        """Numerically: D(a) divisible by F?"""
        if self.d < 6:
            return False, 1.0
        Dv = self.fam.Mmat @ a
        Gm, _ = mons(self.d - 6, 0)
        Fcols = np.array([vec(shift(F, m, 1), self.d, 0) for m in Gm]).T
        x = sla.lstsq(Fcols, Dv, check_finite=False)[0]
        r = float(np.linalg.norm(Fcols @ x - Dv) / max(np.linalg.norm(Dv), 1e-300))
        return r < tol, r

    def search(self, n_restarts=10, seed=0, maxiter=500, gtol=1e-12, verbose=False, stop_rel=1e-13,
               init=None, callback=None):
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
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 30})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            rel, f, c = self.rel(a)
            spur, srel = self.is_spurious(a)
            cand = {'value': f, 'rel': rel, 'a': a, 'spurious': spur, 'div_rel': srel, 'iters': int(res.nit),
                    'restart': it, 'nfev': int(res.nfev), 'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if verbose:
                print(f'[w2s] {self.fam.name} restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} '
                      f'spurious={spur} ({time.time() - t0:.1f}s)', file=sys.stderr, flush=True)
            if callback:
                callback(cand)
            if rel < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out


# ----------------------------------------------------------------------------------------------------------
# rational reconstruction of candidates
# ----------------------------------------------------------------------------------------------------------
def reconstruct_vector(a, max_dens=(1, 2, 3, 4, 6, 12, 24, 60, 120, 1000, 10000), tol=1e-7):
    """Rational candidates for a (normalised so the largest |a_k| = 1): per-coordinate limit_denominator and
    common-denominator scaling."""
    a = np.asarray(a, float)
    j = int(np.argmax(np.abs(a)))
    a = a / a[j]
    cands = []
    for md in max_dens:
        fr = [Fraction(float(v)).limit_denominator(md) for v in a]
        err = max(abs(float(f) - v) for f, v in zip(fr, a))
        if err < tol:
            t = tuple(fr)
            if t not in cands:
                cands.append(t)
    return cands


def D_exact_from_family(fam: Family, frac):
    D = {}
    for k, B in enumerate(fam.exact):
        if frac[k] == 0:
            continue
        D = add(D, scale(B, Fraction(frac[k])))
    return D


# ----------------------------------------------------------------------------------------------------------
# exact (mod p) reduced solver
# ----------------------------------------------------------------------------------------------------------
def e1_lift(D, p=None, verbose=False):
    """Exact lift mod p: u with A D = sum u_i F_i (deg u_i = d+1), or None if A D not in J.
    Returns (u list of 6 dicts mod p, info)."""
    from flint import nmod_mat
    p = p or core.PRIMES[0]
    t0 = time.time()
    Dp = core.reduce_mod(D, p); Ap = core.reduce_mod(A, p)
    d = core.degree(D)
    rhs = mul(Ap, Dp, p)
    tpar = core.parity_of(rhs)
    gens = jacobian.jacobian_generators(F, d + 6, p, tpar)
    rows = {}
    for g in gens:
        for e in g[2]:
            rows.setdefault(e, len(rows))
    for e in rhs:
        rows.setdefault(e, len(rows))
    R, C = len(rows), len(gens)
    flat = [0] * (R * (C + 1))
    for j, g in enumerate(gens):
        for e, v in g[2].items():
            flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items():
        flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    del flat
    Mr, rank = M.rref()
    del M
    ent = Mr.entries()
    del Mr
    W = C + 1
    u = [{} for _ in range(NV)]
    pivots = []
    for i in range(rank):
        base = i * W
        j = 0
        while j <= C and ent[base + j] == 0:
            j += 1
        if j == C:
            return None, {'note': 'A D not in J', 'rank': rank, 'seconds': round(time.time() - t0, 1)}
        v = ent[base + C]
        if v:
            gi, m, _ = gens[j]
            u[gi][m] = int(v)
    return u, {'rank': rank, 'rows': R, 'cols': C, 'seconds': round(time.time() - t0, 1), 'p': p}


def reduced_columns_exact(D, u, p, with_Q=True):
    """Columns (dict polys mod p) of the reduced system P(D,c) - F Q = -p0, plus rhs."""
    d = core.degree(D)
    Dp = core.reduce_mod(D, p)
    Fp = core.reduce_mod(F, p)
    Fip = [pderiv(Fp, i, p) for i in range(NV)]
    Di = [pderiv(Dp, i, p) for i in range(NV)]
    cpar = c_parities()
    dc = d - 4
    cols, labels = [], []
    DFj = [mul(Dp, Fip[j], p) for j in range(NV)]
    W = {}
    for i in range(NV):
        for j in range(i + 1, NV):
            W[(i, j)] = sub(mul(Fip[j], Di[i], p), mul(Fip[i], Di[j], p), p)
    for (i, j), par in cpar.items():
        for m in monomials(dc, par):
            # D F_j d_i m - D F_i d_j m - m W_ij
            col = {}
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(DFj[j], tuple(m2), m[i], p), p)
            if m[j]:
                m2 = list(m); m2[j] -= 1
                col = sub(col, shift(DFj[i], tuple(m2), m[j], p), p)
            col = sub(col, shift(W[(i, j)], m, 1, p), p)
            cols.append(col); labels.append(('c', i, j, m))
    if with_Q:
        for m in monomials(2 * d - 6, PAR_A):
            cols.append(scale(shift(Fp, m, 1, p), p - 1, p)); labels.append(('Q', m))
    # rhs: -p0 = -(D div u - u . grad D)
    divu = {}
    for i in range(NV):
        divu = add(divu, pderiv(u[i], i, p), p)
    p0 = mul(Dp, divu, p)
    for i in range(NV):
        p0 = sub(p0, mul(u[i], Di[i], p), p)
    rhs = scale(p0, p - 1, p)
    return cols, labels, rhs


def exact_reduced_solve(D, p=None, u=None, verbose=False, return_solution=False):
    """Exact mod-p consistency of the reduced system for denominator D.  Returns dict."""
    from flint import nmod_mat
    p = p or core.PRIMES[0]
    t0 = time.time()
    if u is None:
        u, info = e1_lift(D, p)
        if u is None:
            return {'e1': False, 'consistent': False, 'info': info, 'seconds': round(time.time() - t0, 1)}
    cols, labels, rhs = reduced_columns_exact(D, u, p)
    rows = {}
    for col in cols:
        for e in col:
            rows.setdefault(e, len(rows))
    for e in rhs:
        rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    t1 = time.time()
    flat = [0] * (R * (C + 1))
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items():
        flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    del flat
    if verbose:
        print(f'[w2s exact] d={core.degree(D)} rows={R} cols={C} build {t1 - t0:.1f}s', file=sys.stderr, flush=True)
    Mr, rank_aug = M.rref()
    t2 = time.time()
    ent = Mr.entries()
    W = C + 1
    consistent = True
    pivots = []
    sol = {}
    for i in range(rank_aug):
        base = i * W
        j = 0
        while j <= C and ent[base + j] == 0:
            j += 1
        pivots.append(j)
        if j == C:
            consistent = False
            break
        if return_solution:
            v = ent[base + C]
            if v:
                sol[labels[j]] = int(v)
    rank = rank_aug - (0 if consistent else 1)
    out = {'e1': True, 'consistent': consistent, 'rank': rank, 'rank_aug': rank_aug, 'rows': R, 'cols': C,
           'nullity': C - rank, 'p': p, 'seconds': round(time.time() - t0, 1), 't_rref': round(t2 - t1, 1)}
    if return_solution and consistent:
        out['solution'] = sol
        out['u'] = u
    return out


def numerators_from_reduced(D, u, sol, p):
    """Assemble N_i = u_i + sum_j c_ij F_j (mod p) and Q from a reduced solution."""
    Fp = core.reduce_mod(F, p)
    Fip = [pderiv(Fp, i, p) for i in range(NV)]
    N = [dict(u[i]) for i in range(NV)]
    Q = {}
    for lab, v in sol.items():
        if lab[0] == 'c':
            _, i, j, m = lab
            N[i] = add(N[i], shift(Fip[j], m, v, p), p)
            N[j] = sub(N[j], shift(Fip[i], m, v, p), p)
        else:
            Q[lab[1]] = v
    return N, Q


# ----------------------------------------------------------------------------------------------------------
# exact check of the reduction (used by the self-test)
# ----------------------------------------------------------------------------------------------------------
def reduced_identity_check(D, u, c, G=None):
    """Exact: with N = u + K(c) + G Z/6, E(D,N) == F * (D div N - N.grad D - D G)."""
    N = [dict(u[i]) for i in range(NV)]
    for (i, j), cij in c.items():
        N[i] = add(N[i], mul(cij, Fi[j]))
        N[j] = sub(N[j], mul(cij, Fi[i]))
    if G:
        for i in range(NV):
            e = [0] * NV; e[i] = 1
            N[i] = add(N[i], scale(shift(G, tuple(e), 1), Fraction(1, 6)))
    E = core.defect(F, A, D, 1, N)
    divN = {}
    for i in range(NV):
        divN = add(divN, pderiv(N[i], i))
    P = mul(D, divN)
    for i in range(NV):
        P = sub(P, mul(N[i], pderiv(D, i)))
    if G:
        P = sub(P, mul(D, G))
    return sub(E, mul(F, P)) == {}
