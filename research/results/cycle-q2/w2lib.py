"""wave-2 library for lane cycle-q2: sector-restricted, Koszul-reduced, complex variable-projection
joint search for a denominator D in a linear family, plus exact confirmation.

Formulation (docs/05 §5.2 "full correction map").  For D in the E1 kernel K_kappa = {D in (S_d)_kappa : A D in J(F)}
pick an E1 lift N0(D) (linear in D, q N0(F) = A D).  Every E1 lift is N0 + K(c) with
    K(c) = sum_{i<j} c_ij (F_j d_i - F_i d_j),   deg c_ij = d + 6q - 10,
and E2 is the affine membership problem
    -B_D(N0) in F^q S_{2d-6} + im T_D,   B_D(N) = D div N - N(D),
    T_D(c) = sum_{i<j} D [(d_i c_ij) F_j - (d_j c_ij) F_i] - c_ij (D_i F_j - D_j F_i).
With D = D(a) = sum_k a_k D_k:  M(a) z = b(a),  M(a) = [T_{D(a)} | -F^q .] linear in a,  b(a) = -B_{D(a)}(N0(a)) quadratic.
f(a) = min_z |M(a) z - b(a)|^2 (complex least squares), minimized over the unit sphere by L-BFGS with the
variable-projection gradient.  Everything is restricted to the character sector forced by kappa (lossless by averaging).
Only an exact re-solve (eq.EqSystem / eq.solve_exact + core.verify_witness) counts as evidence.

Conventions: eq.py (lane dir) group = all w in (Z/6)^6 fixing the line of F (order 432 incl. scalars; 72 projective
characters); keys = base-6 digit vectors of weights on a greedy generating set; kF = 0; char(F_j) = -kx[j].
"""
from __future__ import annotations
import sys, os, time, json, math, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
from flint import nmod_mat, fmpz_mat, fmpq_mat, fmpq
from hodge import core
import eq

NV = 6
BASE = eq.BASE
POW = eq.POW
LANE = '/tmp/hodge/results/cycle-q2'
LOG = os.path.join(LANE, 'wave2-search.log.jsonl')
PRIME = eq.PRIME  # 2^31 - 1
ENGINE = 2


def log(**kw):
    kw.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    kw.setdefault('lane', 'cycle-q2')
    kw.setdefault('engine', ENGINE)
    with open(LOG, 'a') as f:
        f.write(json.dumps(kw, default=str) + '\n')


# ----------------------------------------------------------------------------
# polynomials as (codes int64, vals float64) arrays
# ----------------------------------------------------------------------------

def enc(e):
    return int(sum(int(e[i]) * BASE ** i for i in range(NV)))


def dec(c):
    e = []
    for i in range(NV):
        e.append(int(c % BASE)); c //= BASE
    return tuple(e)


def parr(P):
    """exact dict poly -> (codes, float vals)"""
    if not P:
        return np.zeros(0, dtype=np.int64), np.zeros(0)
    codes = np.array([enc(e) for e in P], dtype=np.int64)
    vals = np.array([float(v) for v in P.values()], dtype=float)
    return codes, vals


def carr(P):
    """exact or float dict poly -> (codes, complex vals)"""
    if not P:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=complex)
    codes = np.array([enc(e) for e in P], dtype=np.int64)
    vals = np.array([complex(v) for v in P.values()], dtype=complex)
    return codes, vals


def acc_arr(codes, vals):
    """sum duplicate codes"""
    if len(codes) == 0:
        return codes, vals
    u, inv = np.unique(codes, return_inverse=True)
    out = np.zeros(len(u), dtype=vals.dtype)
    np.add.at(out, inv, vals)
    keep = out != 0
    return u[keep], out[keep]


def mul_arr(c1, v1, c2, v2):
    if len(c1) == 0 or len(c2) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.result_type(v1, v2))
    cc = (c1[:, None] + c2[None, :]).ravel()
    vv = (v1[:, None] * v2[None, :]).ravel()
    return acc_arr(cc, vv)


def deriv_arr(codes, vals, i):
    """d/dx_i of a (codes, vals) poly"""
    ei = (codes // (BASE ** i)) % BASE
    has = ei > 0
    return codes[has] - BASE ** i, vals[has] * ei[has]


def arr_to_dict(codes, vals):
    return {dec(int(c)): v for c, v in zip(codes, vals)}


def dict_from_vec(mons, vec):
    """exponent-tuple list + coefficient vector -> dict (drops zeros)"""
    return {tuple(int(x) for x in m): v for m, v in zip(mons, vec) if v != 0}


# ----------------------------------------------------------------------------
# sectors (character bookkeeping)
# ----------------------------------------------------------------------------

class Sectors:
    def __init__(self, F, A):
        self.F, self.A = F, A
        self.Gw = eq.symmetry_group(F)
        self.gens = eq.generators(self.Gw)
        self.ng = self.gens.shape[0]
        self.order = len(self.Gw)
        self.kA = eq.key_of_poly(A, self.gens)
        self.kx = [int(eq.keys_of(np.eye(NV, dtype=np.int64)[i:i + 1], self.gens)[0]) for i in range(NV)]
        self._cache = {}

    def add(self, k1, k2):
        return eq.key_add(k1, k2, self.ng)

    def scale(self, k, c):
        return eq.key_scale(k, c, self.ng)

    def neg(self, k):
        return eq.key_scale(k, 5, self.ng)

    def key(self, P):
        return eq.key_of_poly(P, self.gens)

    def all_keys(self, deg):
        ex = eq.monomial_exps(deg)
        return sorted(set(eq.keys_of(ex, self.gens).tolist()))

    def mons(self, deg, key):
        """exponent array (n x 6) of the monomials of degree deg with character key"""
        ck = (deg, key)
        if ck not in self._cache:
            ex = eq.monomial_exps(deg)
            ks = eq.keys_of(ex, self.gens)
            self._cache[ck] = ex[ks == key]
        return self._cache[ck]

    def kN(self, kappa):
        return [self.add(self.add(self.kA, kappa), self.kx[i]) for i in range(NV)]

    def kC(self, kappa):
        kN = self.kN(kappa)
        return [[self.add(kN[i], self.kx[j]) for j in range(NV)] for i in range(NV)]

    def kE(self, kappa):
        return self.add(self.kA, self.scale(kappa, 2))

    def parity_of_key(self, deg, key):
        ms = self.mons(deg, key)
        if len(ms) == 0:
            return None
        return int((ms[0, 3] + ms[0, 4] + ms[0, 5]) % 2)


# ----------------------------------------------------------------------------
# E1 kernel in a sector: mod p dimensions and numeric (SVD) basis
# ----------------------------------------------------------------------------

def e1_matrices(S, d, kappa, q=2, dtype=float):
    """Returns (Dm, Amat, Jcols) for the E1 problem  A D = q sum N_i F_i in degree d+deg A, sector kappa:
    Dm: D monomials (n x 6); row codes: monomials of degree d+degA with key kA+kappa (complete sector);
    Amat: rows x nD float matrix of A*m; Jmat: rows x nN matrix of q F_i * n (N monomials by sector)."""
    F, A = S.F, S.A
    degA = core.degree(A)
    Dm = S.mons(d, kappa)
    degN = d + 6 * q - 5
    kN = S.kN(kappa)
    rows_ex = S.mons(d + degA, S.add(S.kA, kappa))
    row_codes = np.sort(rows_ex @ POW)
    ridx = {int(c): i for i, c in enumerate(row_codes)}
    R = len(row_codes)
    Ac, Av = parr(A)
    nD = len(Dm)
    Amat = np.zeros((R, nD))
    if nD:
        mc = Dm @ POW
        cc = (Ac[:, None] + mc[None, :]).ravel()
        vv = np.broadcast_to(Av[:, None], (len(Ac), nD)).ravel()
        jj = np.broadcast_to(np.arange(nD)[None, :], (len(Ac), nD)).ravel()
        rr = np.searchsorted(row_codes, cc)
        assert np.all(row_codes[rr] == cc)
        np.add.at(Amat, (rr, jj), vv)
    Ncols = []
    blocks = []
    ncol = 0
    Fi = [core.pderiv(F, i) for i in range(NV)]
    for i in range(NV):
        sel = S.mons(degN, kN[i])
        if len(sel) == 0:
            continue
        Fc, Fv = parr(core.scale(Fi[i], Fraction(q)))
        mc = sel @ POW
        cc = (Fc[:, None] + mc[None, :]).ravel()
        vv = np.broadcast_to(Fv[:, None], (len(Fc), len(sel))).ravel()
        jj = np.broadcast_to(np.arange(ncol, ncol + len(sel))[None, :], (len(Fc), len(sel))).ravel()
        rr = np.searchsorted(row_codes, cc)
        assert np.all(row_codes[rr] == cc), 'E1 column outside sector rows'
        blocks.append((rr, jj, vv))
        for m in sel:
            Ncols.append((i, tuple(int(x) for x in m)))
        ncol += len(sel)
    Jmat = np.zeros((R, ncol))
    for rr, jj, vv in blocks:
        np.add.at(Jmat, (rr, jj), vv)
    return Dm, Amat, Jmat, Ncols, row_codes


def e1_kernel_modp(S, d, kappa, q=2, p=PRIME):
    """dim K_kappa (E1 kernel), dim (J_d)_kappa, dim (S_d)_kappa, and a basis of K_kappa mod p (list of dict polys).
    Built from the exact polynomials (reduced mod p), so arbitrary rational coefficients are fine."""
    F, A = S.F, S.A
    degA = core.degree(A)
    Dm = S.mons(d, kappa)
    nD = len(Dm)
    if nD == 0:
        return {'dimS': 0, 'dimK': 0, 'dimJ': 0, 'basis': [], 'p': p}
    degN = d + 6 * q - 5
    kN = S.kN(kappa)
    rows_ex = S.mons(d + degA, S.add(S.kA, kappa))
    row_codes = np.sort(rows_ex @ POW)
    R = len(row_codes)
    Ap = core.reduce_mod(A, p)
    Fp = core.reduce_mod(F, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(NV)]
    cols = []
    for m in Dm:
        cols.append(core.shift(Ap, tuple(int(x) for x in m), 1, p))
    for i in range(NV):
        for n in S.mons(degN, kN[i]):
            cols.append(core.shift(Fi[i], tuple(int(x) for x in n), (p - q) % p, p))
    C = len(cols)
    flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items():
            code = enc(e)
            r = int(np.searchsorted(row_codes, code))
            assert r < R and row_codes[r] == code
            flat[r * C + j] = v
    Mn = nmod_mat(R, C, flat, p)
    X, nul = Mn.nullspace()
    vecs = []
    for k in range(nul):
        v = [int(X[j, k]) for j in range(nD)]
        if any(v):
            vecs.append(v)
    basis = []
    dimK = 0
    if vecs:
        P = nmod_mat(len(vecs), nD, [x for v in vecs for x in v], p)
        Pr, rk = P.rref()
        dimK = rk
        for i in range(rk):
            basis.append({tuple(int(x) for x in Dm[j]): int(Pr[i, j]) for j in range(nD) if int(Pr[i, j])})
    dimJ = 0
    if d >= 5:
        jb = jac_basis(S, d, kappa)
        if jb:
            idx = {tuple(int(x) for x in m): j for j, m in enumerate(Dm)}
            flat = []
            for c in jb:
                row = [0] * nD
                for e, v in core.reduce_mod(c, p).items():
                    row[idx[e]] = v
                flat.extend(row)
            Jm = nmod_mat(len(jb), nD, flat, p)
            dimJ = Jm.rank()
    return {'dimS': nD, 'dimK': dimK, 'dimJ': dimJ, 'basis': basis, 'p': p}


def e1_kernel_numeric(S, d, kappa, q=2, tol=1e-9):
    """Orthonormal real basis (columns, nD x K) of K_kappa = {D : A D in J} via SVD. Also returns the
    numeric codimension check and the singular value gap."""
    Dm, Amat, Jmat, Ncols, row_codes = e1_matrices(S, d, kappa, q)
    nD = len(Dm)
    if nD == 0:
        return Dm, np.zeros((0, 0)), {}
    # kernel of [Amat | -Jmat] projected to D coordinates: D-part of null space.
    # Equivalent: D with Amat D in col(Jmat):  (I - P_J) Amat D = 0  where P_J = orthogonal projector on col(Jmat).
    U, s, Vt = np.linalg.svd(Jmat, full_matrices=False)
    r = int(np.sum(s > tol * s[0])) if len(s) else 0
    UJ = U[:, :r]
    Res = Amat - UJ @ (UJ.T @ Amat)          # (I - P_J) A
    Ur, sr, Vrt = np.linalg.svd(Res, full_matrices=True)
    smax = sr[0] if len(sr) else 1.0
    rk = int(np.sum(sr > 1e-9 * max(smax, 1.0)))
    K = Vrt[rk:, :].T   # nD x (nD - rk)
    gap = (sr[rk - 1] / smax if rk > 0 else None, sr[rk] / smax if rk < len(sr) else 0.0)
    return Dm, K, {'dimS': nD, 'dimK': K.shape[1], 'rankJ': r, 'gap': gap}


def f_multiples_basis(S, d, kappa):
    """basis (list of dict polys) of F * (S_{d-6})_kappa (key(F) = 0)."""
    if d < 6:
        return []
    out = []
    for m in S.mons(d - 6, kappa):
        out.append(core.shift(S.F, tuple(int(x) for x in m), 1))
    return out


def jac_basis(S, d, kappa):
    """basis of (J_d)_kappa as dict polys: n * F_i with key(n) = kappa + kx[i]"""
    out = []
    if d < 5:
        return out
    Fi = [core.pderiv(S.F, i) for i in range(NV)]
    for i in range(NV):
        for m in S.mons(d - 5, S.add(kappa, S.kx[i])):
            out.append(core.shift(Fi[i], tuple(int(x) for x in m), 1))
    return out


def vecs_from_dicts(Dm, polys):
    idx = {tuple(int(x) for x in m): j for j, m in enumerate(Dm)}
    V = np.zeros((len(Dm), len(polys)))
    for k, P in enumerate(polys):
        for e, v in P.items():
            V[idx[e], k] = float(v)
    return V


def orth_basis(V, tol=1e-10):
    """orthonormal basis of the column space of V"""
    if V.size == 0:
        return V
    U, s, _ = np.linalg.svd(V, full_matrices=False)
    r = int(np.sum(s > tol * s[0]))
    return U[:, :r]


def intersect_subspaces(K1, K2, tol=1e-9):
    """orthonormal basis of col(K1) ∩ col(K2) (both orthonormal column bases)"""
    if K1.shape[1] == 0 or K2.shape[1] == 0:
        return np.zeros((K1.shape[0], 0))
    # vectors x = K1 a = K2 b  <=> [K1 -K2][a;b] = 0
    M = np.hstack([K1, -K2])
    U, s, Vt = np.linalg.svd(M, full_matrices=True)
    rk = int(np.sum(s > tol))
    N = Vt[rk:, :].T
    X = K1 @ N[:K1.shape[1], :]
    return orth_basis(X) if X.shape[1] else X


def subspace_minus(K, Wsub):
    """orthonormal basis of the orthogonal complement of col(Wsub) inside col(K)"""
    if Wsub.shape[1] == 0:
        return K
    P = K @ (K.T @ Wsub)
    Pb = orth_basis(P)
    # complement within col(K)
    C = K - Pb @ (Pb.T @ K)
    return orth_basis(C, tol=1e-8)


# ----------------------------------------------------------------------------
# Koszul-reduced sector varpro
# ----------------------------------------------------------------------------

class KoszulVarPro:
    """Family D(a) = sum_k a_k D_k with D_k = Dm-monomial combinations given by the columns of Bmat (nD x K)."""

    def __init__(self, S: Sectors, d, kappa, Dm, Bmat, q=2, verbose=False, cond=1e-9, lapack='gelsy'):
        t0 = time.time()
        self.S, self.d, self.kappa, self.q = S, d, kappa, q
        F, A = S.F, S.A
        self.F, self.A = F, A
        self.Dm = Dm
        self.nD = len(Dm)
        self.Bmat = np.asarray(Bmat, dtype=float)
        self.K = self.Bmat.shape[1]
        self.cond = cond
        self.lapack = lapack
        degC = d + 6 * q - 10
        degQ = 2 * d - 6
        degE = 2 * d + 6 * q - 6
        self.degC, self.degQ, self.degE = degC, degQ, degE
        kC = S.kC(kappa); kE = S.kE(kappa)
        self.kE = kE
        rows_ex = S.mons(degE, kE)
        row_codes = np.sort(rows_ex @ POW)
        self.row_codes = row_codes
        R = len(row_codes); self.R = R
        Fi = [core.pderiv(F, i) for i in range(NV)]
        Fa = [parr(Fi[i]) for i in range(NV)]
        # c_ij columns
        cols = []  # (i, j, exps array)
        ncol = 0
        self.cmeta = []
        for i in range(NV):
            for j in range(i + 1, NV):
                sel = S.mons(degC, kC[i][j]) if degC >= 0 else np.zeros((0, NV), dtype=np.int64)
                if len(sel) == 0:
                    continue
                cols.append((i, j, sel, ncol))
                for n in sel:
                    self.cmeta.append(('c', i, j, tuple(int(x) for x in n)))
                ncol += len(sel)
        self.nC = ncol
        # Q columns: -F^q * n
        Qm = S.mons(degQ, kE) if degQ >= 0 else np.zeros((0, NV), dtype=np.int64)
        self.nQ = len(Qm)
        self.ncols = self.nC + self.nQ
        self.Qmeta = [tuple(int(x) for x in n) for n in Qm]
        # ---- monomial T-matrices
        self.Mm = []
        mcodes = Dm @ POW

        def rows_of(cc):
            rr = np.searchsorted(row_codes, cc)
            ok = np.all(row_codes[np.minimum(rr, R - 1)] == cc)
            assert ok, 'code outside sector row set'
            return rr
        for mi_, m in enumerate(Dm):
            cm = int(mcodes[mi_])
            ri, ci, vv = [], [], []
            for (i, j, sel, c0) in cols:
                nc_ = sel @ POW
                ni = sel[:, i]; nj = sel[:, j]
                mi = int(m[i]); mj = int(m[j])
                Fjc, Fjv = Fa[j]; Fic, Fiv = Fa[i]
                # term1: n_i * (m n / x_i) F_j
                has = ni > 0
                if np.any(has):
                    cc = (Fjc[:, None] + (cm + nc_[has] - BASE ** i)[None, :]).ravel()
                    v = (Fjv[:, None] * ni[has][None, :]).ravel()
                    jj = np.broadcast_to((np.arange(len(sel))[has] + c0)[None, :], (len(Fjc), int(has.sum()))).ravel()
                    ri.append(rows_of(cc)); ci.append(jj); vv.append(v)
                # term2: - n_j * (m n / x_j) F_i
                has = nj > 0
                if np.any(has):
                    cc = (Fic[:, None] + (cm + nc_[has] - BASE ** j)[None, :]).ravel()
                    v = -(Fiv[:, None] * nj[has][None, :]).ravel()
                    jj = np.broadcast_to((np.arange(len(sel))[has] + c0)[None, :], (len(Fic), int(has.sum()))).ravel()
                    ri.append(rows_of(cc)); ci.append(jj); vv.append(v)
                # term3: - n * m_i (m/x_i) F_j
                if mi > 0:
                    cc = (Fjc[:, None] + (cm - BASE ** i + nc_)[None, :]).ravel()
                    v = np.broadcast_to((-mi * Fjv)[:, None], (len(Fjc), len(sel))).ravel()
                    jj = np.broadcast_to((np.arange(len(sel)) + c0)[None, :], (len(Fjc), len(sel))).ravel()
                    ri.append(rows_of(cc)); ci.append(jj); vv.append(v)
                # term4: + n * m_j (m/x_j) F_i
                if mj > 0:
                    cc = (Fic[:, None] + (cm - BASE ** j + nc_)[None, :]).ravel()
                    v = np.broadcast_to((mj * Fiv)[:, None], (len(Fic), len(sel))).ravel()
                    jj = np.broadcast_to((np.arange(len(sel)) + c0)[None, :], (len(Fic), len(sel))).ravel()
                    ri.append(rows_of(cc)); ci.append(jj); vv.append(v)
            if ri:
                ri = np.concatenate(ri); ci = np.concatenate(ci); vv = np.concatenate(vv)
                Mk = sp.csr_matrix((vv, (ri, ci)), shape=(R, self.nC))
                Mk.sum_duplicates(); Mk.eliminate_zeros()
            else:
                Mk = sp.csr_matrix((R, self.nC))
            self.Mm.append(Mk)
        # Q block (dense, constant)
        self.MQ = np.zeros((R, self.nQ))
        if self.nQ:
            Fq = core.power(F, q)
            Fqc, Fqv = parr(Fq)
            qc = Qm @ POW
            cc = (Fqc[:, None] + qc[None, :]).ravel()
            v = np.broadcast_to((-Fqv)[:, None], (len(Fqc), self.nQ)).ravel()
            jj = np.broadcast_to(np.arange(self.nQ)[None, :], (len(Fqc), self.nQ)).ravel()
            np.add.at(self.MQ, (rows_of(cc), jj), v)
        # ---- E1 lifts N0^{(k)} for the family basis (min-norm least squares on the sector E1 system)
        Dm2, Amat, Jmat, Ncols, e1rows = e1_matrices(S, d, kappa, q)
        assert len(Dm2) == self.nD
        rhs = Amat @ self.Bmat  # rows x K
        if Jmat.shape[1] == 0:
            raise ValueError('empty N sector')
        N0, res, rk, sv = sla.lstsq(Jmat, rhs, cond=1e-12, lapack_driver='gelsd', check_finite=False)
        e1res = np.linalg.norm(Jmat @ N0 - rhs, axis=0) / np.maximum(np.linalg.norm(rhs, axis=0), 1e-300)
        self.e1_residual = float(e1res.max()) if len(e1res) else 0.0
        self.N0 = N0  # nN x K
        self.Ncols = Ncols
        # ---- W tensor: W[r, m, k] = -B_m(N0^{(k)}) = -( m div N0 - sum_i d_i m * N0_i )
        nN = len(Ncols)
        Ncodes = np.array([enc(n) for (_, n) in Ncols], dtype=np.int64)
        Nvar = np.array([i for (i, _) in Ncols], dtype=np.int64)
        W = np.zeros((R, self.nD, self.K))
        # div N0: for each column (i, n) with n_i > 0 -> code n - BASE^i, coefficient n_i
        ni_all = (Ncodes // (BASE ** np.arange(NV))[Nvar]) % BASE
        for mi_, m in enumerate(Dm):
            cm = int(mcodes[mi_])
            # term a: - m * div N0
            has = ni_all > 0
            cc = cm + Ncodes[has] - (BASE ** Nvar[has])
            coef = -ni_all[has].astype(float)
            rr = rows_of(cc)
            contrib = coef[:, None] * N0[has, :]
            np.add.at(W[:, mi_, :], rr, contrib)
            # term b: + sum_i m_i (m/x_i) N0_i
            for i in range(NV):
                mi = int(m[i])
                if mi == 0:
                    continue
                sel = Nvar == i
                cc = cm - BASE ** i + Ncodes[sel]
                rr = rows_of(cc)
                np.add.at(W[:, mi_, :], rr, mi * N0[sel, :])
        self.W = W
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            nnz = sum(int(Mk.nnz) for Mk in self.Mm)
            print(f'[kvp] d={d} key={kappa} K={self.K} nD={self.nD} rows={R} cols={self.ncols} (c {self.nC}, Q {self.nQ}) '
                  f'nnz(T)={nnz} e1res={self.e1_residual:.1e} built {self.build_time:.1f}s', file=sys.stderr)

    # ------------------------------------------------------------------
    def assemble(self, a):
        """dense complex M(a) (R x ncols) and b(a)"""
        at = self.Bmat @ a                     # nD complex
        M = None
        for mi_, Mk in enumerate(self.Mm):
            if at[mi_] == 0 or Mk.nnz == 0:
                continue
            M = Mk * at[mi_] if M is None else M + Mk * at[mi_]
        Md = np.zeros((self.R, self.ncols), dtype=complex)
        if M is not None:
            Md[:, :self.nC] = M.toarray()
        if self.nQ:
            Md[:, self.nC:] = self.MQ
        b = np.einsum('rmk,m,k->r', self.W, at, a)
        return Md, b, at

    def value(self, a):
        a = np.asarray(a, dtype=complex)
        na = np.linalg.norm(a)
        a = a / na
        M, b, at = self.assemble(a)
        x = sla.lstsq(M, b, cond=self.cond, lapack_driver=self.lapack, check_finite=False)[0]
        r = M @ x - b
        self._nfev += 1
        return float(np.vdot(r, r).real), x, r, b

    def _fg(self, a):
        """f and complex 'gradient' g_k = r^H dr/da_k at unit a (no normalization chain rule)."""
        M, b, at = self.assemble(a)
        x = sla.lstsq(M, b, cond=self.cond, lapack_driver=self.lapack, check_finite=False)[0]
        r = M @ x - b
        f = float(np.vdot(r, r).real)
        xc = x[:self.nC]
        # dM/da_k x = sum_m B[m,k] (M_m xc)
        ym = np.zeros((self.R, self.nD), dtype=complex)
        for mi_, Mk in enumerate(self.Mm):
            if Mk.nnz:
                ym[:, mi_] = Mk @ xc
        rH_ym = np.conj(r) @ ym                     # nD
        gM = self.Bmat.T @ rH_ym                    # K
        # db/da_k = sum_m B[m,k] (W[:,m,:] a) + sum_m at_m W[:,m,k]
        Wa = np.einsum('rmk,k->rm', self.W, a)      # R x nD
        Wat = np.einsum('rmk,m->rk', self.W, at)    # R x K
        rH_Wa = np.conj(r) @ Wa                     # nD
        gb = self.Bmat.T @ rH_Wa + np.conj(r) @ Wat
        g = gM - gb                                 # complex, K
        self._nfev += 1
        return f, g, x, r, b

    def value_and_grad_real(self, z, complex_mode=True):
        """z in R^{2K} (re, im) or R^K; returns f, grad_z (sphere normalization applied)."""
        K = self.K
        if complex_mode:
            a_full = z[:K] + 1j * z[K:]
        else:
            a_full = z.astype(complex)
        na = np.linalg.norm(a_full)
        a = a_full / na
        f, g, x, r, b = self._fg(a)
        # real gradient wrt (re a, im a): d f/d re = 2 Re(g), d f / d im = -2 Im(g)
        G = np.concatenate([2 * g.real, -2 * g.imag])
        zunit = np.concatenate([a.real, a.imag])
        Gt = (G - float(G @ zunit) * zunit) / na
        if complex_mode:
            return f, Gt
        return f, Gt[:K]

    def search(self, n_restarts=10, seed=0, maxiter=200, gtol=1e-12, complex_mode=True, verbose=False,
               init=None, stop_rel=1e-13, time_budget=None, callback=None):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        K = self.K
        dim = 2 * K if complex_mode else K
        for it in range(n_restarts):
            if time_budget and time.time() - t0 > time_budget:
                break
            if init is not None and it == 0:
                z0 = np.asarray(init, float)
            else:
                z0 = rng.standard_normal(dim)
            z0 /= np.linalg.norm(z0)
            res = minimize(lambda z: self.value_and_grad_real(z, complex_mode), z0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            z = res.x
            a = (z[:K] + 1j * z[K:]) if complex_mode else z.astype(complex)
            a /= np.linalg.norm(a)
            f, x, r, b = self.value(a)
            bn = float(np.linalg.norm(b))
            rel = f / max(bn * bn, 1e-300)
            spur, srel = self.is_spurious(a)
            cand = {'f': f, 'rel': rel, 'bnorm': bn, 'a': a, 'iters': int(res.nit), 'restart': it,
                    'spurious': spur, 'div_rel': srel, 'nfev': int(res.nfev)}
            out.append(cand)
            if verbose:
                print(f'[kvp] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} spur={spur} ({time.time()-t0:.1f}s)',
                      file=sys.stderr)
            if callback:
                callback(cand)
            if rel < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    # ------------------------------------------------------------------
    def D_coeffs(self, a):
        """coefficient vector of D(a) over Dm (complex)"""
        return self.Bmat @ np.asarray(a, dtype=complex)

    def is_spurious(self, a, tol=1e-8):
        """is D(a) numerically divisible by F? (compare against F*(S_{d-6})_kappa)"""
        if self.d < 6:
            return False, 1.0
        Fb = f_multiples_basis(self.S, self.d, self.kappa)
        if not Fb:
            return False, 1.0
        V = vecs_from_dicts(self.Dm, Fb)
        dv = self.D_coeffs(a)
        x = sla.lstsq(V.astype(complex), dv, check_finite=False)[0]
        rel = float(np.linalg.norm(V @ x - dv) / max(np.linalg.norm(dv), 1e-300))
        return rel < tol, rel

    def D_dict_float(self, a):
        dv = self.D_coeffs(a)
        return {tuple(int(x) for x in m): complex(v) for m, v in zip(self.Dm, dv) if abs(v) > 1e-14}


# ----------------------------------------------------------------------------
# rational reconstruction and exact confirmation
# ----------------------------------------------------------------------------

def normalize_phase(dv):
    """rotate a complex coefficient vector so the largest entry is real positive; return (vec, imag_rel)."""
    j = int(np.argmax(np.abs(dv)))
    ph = dv[j] / abs(dv[j])
    v = dv / ph
    imag_rel = float(np.linalg.norm(v.imag) / max(np.linalg.norm(v), 1e-300))
    return v, imag_rel


def rational_candidates(v, max_dens=(1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 256, 512, 1024, 4096, 2 ** 14, 2 ** 17, 2 ** 20), tol=1e-7):
    """real vector (max entry 1) -> list of Fraction tuples (distinct), by common-denominator scaling and per-entry CF."""
    v = np.asarray(v, float)
    cands = []
    for den in max_dens:
        w = v * den
        if np.max(np.abs(w - np.round(w))) < tol * den:
            fr = tuple(Fraction(int(round(x)), den) for x in w)
            if fr not in cands:
                cands.append(fr)
    for md in (10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6):
        fr = tuple(Fraction(float(x)).limit_denominator(md) for x in v)
        if max(abs(float(f) - x) for f, x in zip(fr, v)) < tol:
            if fr not in cands:
                cands.append(fr)
    return cands


def lll_integer_relation(v, scale=10 ** 9):
    """integer vector c (small) with c proportional to v, via LLL on [I | scale*v]; returns Fraction tuple or None"""
    v = np.asarray(v, float)
    n = len(v)
    rows = []
    for i in range(n):
        row = [1 if j == i else 0 for j in range(n)] + [int(round(scale * v[i]))]
        rows.append(row)
    # we want integer combos lambda with sum lambda_i v_i small: that is a relation, not proportionality.
    # proportionality: find integer c with c_i/c_j = v_i/v_j: use the lattice of (c_i v_j - c_j v_i) -- instead do
    # common denominator search on v / v_j for the pivot j.
    j = int(np.argmax(np.abs(v)))
    w = v / v[j]
    best = None
    for den in range(1, 20001):
        x = w * den
        if np.max(np.abs(x - np.round(x))) < 1e-7 * den:
            best = tuple(Fraction(int(round(t)), den) for t in x)
            break
    return best


def exact_D_from_fracs(Dm, fr):
    D = {}
    for m, c in zip(Dm, fr):
        if c != 0:
            D[tuple(int(x) for x in m)] = Fraction(c)
    return D


def exact_confirm(F, A, q, D, verbose=False, do_lift=True):
    """E1 check, full sector system consistency mod 2^31-1 (two primes), CRT lift + exact verification, F | D check."""
    from hodge import jacobian
    out = {'d': core.degree(D)}
    r1 = eq.e1_system(F, A, q, D)
    out['e1'] = bool(r1['consistent']); out['e1_rank'] = r1['rank']
    if not r1['consistent']:
        out['status'] = 'e1-fails'
        return out
    S1 = eq.EqSystem(F, A, q, D, p=PRIME)
    res = S1.solve()
    out.update({'consistent_p1': bool(res['consistent']), 'rank': res['rank'], 'nullity': res['nullity'],
                'rows': res['rows'], 'cols': res['cols'], 'order': S1.order})
    if not res['consistent']:
        out['status'] = 'inconsistent'
        return out
    S2 = eq.EqSystem(F, A, q, D, p=core.PRIMES[1])
    res2 = S2.solve()
    out['consistent_p2'] = bool(res2['consistent'])
    divF, _ = jacobian.divisible_by(F, D)
    out['F_divides_D'] = bool(divF)
    if not do_lift:
        out['status'] = 'consistent-modp'
        return out
    status, data = eq.solve_exact(F, A, q, D, max_primes=12, verbose=verbose)
    out['lift'] = status
    if status == 'solved':
        ok, info = core.verify_witness(F, A, q, D, data['N'], data['Q'])
        out['verified'] = bool(ok)
        out['N'] = [core.to_str(n) for n in data['N']]
        out['Q'] = core.to_str(data['Q'])
        out['status'] = 'WITNESS' if (ok and not divF) else ('spurious-F-multiple' if divF else 'verify-failed')
    else:
        out['status'] = 'consistent-but-lift-failed'
    return out


def family_label(kind, **kw):
    return kind + ('' if not kw else ':' + ','.join(f'{k}={v}' for k, v in kw.items()))
