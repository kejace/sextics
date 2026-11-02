#!/usr/bin/env python3
"""w2_jvp.py -- wave 2: Jacobian-type variable-projection search for quadratic-cubic-pullback (q = 1, A = X0^3 Y0^3).

Family of Y-even denominators of degree d:
    D(u) = sum_{l<6} w_l F_l  (+ sum_k a_k B_k),     u = (w, a),
w_l in S_{d-5}^{par_l} (par_l = 0 for l < 3, 1 for l >= 3), B_k = optional real basis of a non-Jacobian part of the E1
kernel with precomputed E1 lifts n_k (A B_k = sum_i n_ki F_i).  E1 lift:  N0(u) = A w + sum_k a_k n_k  (linear in u).

docs/05-computation.md §5.2 (q = 1): every E1 lift is N0 + K with K = sum_{i<j} c_ij (F_j d_i - F_i d_j), and E2 holds iff
    -B_D(N0) in F S_{2d-6} + im T_D,     B_D(N) = D div N - N(D),     T_D(c) = B_D(K)  (every term in J(F)).
Stage 1 ("Hessian"/Jacobian-quotient relaxation, necessary):  r1(u) = Z^T B_D(N0)(u), Z = orthonormal basis of the
    complement of J(F)_{2d}^odd in S_{2d}^odd (dim R_{2d}^odd: 866, 744, 468, 208, 60 for d = 6..10).
Stage 2 (full E2, variable projection):  r2(u) = min_{c,Q} || S (B_D(N0)(u) + T_D(c) - F Q) ||  where S = random row
    compression of the projector onto the complement of F S_{2d-6}^odd (k = nc + extra rows); exact residual is
    recomputed without compression for the best points.
Spurious directions (D = F G, and the kernel of u -> D) are projected out of the parameter space; |u| = 1 is kept by the
extra residual (u.u - 1).  Verified on the Fermat control (F0 = sum Z^6, D = X0^6+X1^6+X2^6 = sum_{i<3} Z_i F_i / 6).
"""
from __future__ import annotations
import sys, os, json, time, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import least_squares, minimize
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, pderiv, mul, add, scale, shift

LANE = 'quadratic-cubic-pullback'


# ----------------------------------------------------------------------------------------------------------------------
# dense monomial bases and product tables
# ----------------------------------------------------------------------------------------------------------------------
def code(e):
    return sum(int(e[i]) << (6 * i) for i in range(6))


class Basis:
    def __init__(self, deg, par):
        self.deg, self.par = deg, par
        mons = monomials(deg, par) if deg >= 0 else []
        codes = np.array([code(m) for m in mons], dtype=np.int64)
        order = np.argsort(codes)
        self.mons = [mons[i] for i in order]
        self.codes = codes[order]
        self.index = {m: i for i, m in enumerate(self.mons)}
        self.n = len(self.mons)

    def vec(self, P):
        v = np.zeros(self.n)
        for e, c in P.items():
            v[self.index[e]] += float(c)
        return v

    def poly(self, v, tol=0.0):
        return {m: float(v[i]) for i, m in enumerate(self.mons) if abs(v[i]) > tol}

    def lookup(self, codes):
        idx = np.searchsorted(self.codes, codes)
        idx = np.minimum(idx, self.n - 1)
        ok = self.codes[idx] == codes
        if not np.all(ok):
            raise ValueError('monomial outside basis')
        return idx.astype(np.int32)


def prod_table(B1, B2, Bout):
    """int32 array T[a,b] = index in Bout of mon_a * mon_b."""
    if B1.n == 0 or B2.n == 0:
        return np.zeros((B1.n, B2.n), dtype=np.int32)
    codes = B1.codes[:, None] + B2.codes[None, :]
    return Bout.lookup(codes.ravel()).reshape(B1.n, B2.n)


def tmul(T, p, q, nout):
    """product of dense polys p (over B1) and q (over B2) via table T -> dense over Bout."""
    return np.bincount(T.ravel(), weights=np.outer(p, q).ravel(), minlength=nout)


def deriv_matrix(Bin, Bout, i):
    """sparse (Bout.n x Bin.n) matrix of d/dz_i."""
    ri, ci, vals = [], [], []
    for j, m in enumerate(Bin.mons):
        if m[i]:
            m2 = list(m); m2[i] -= 1
            ri.append(Bout.index[tuple(m2)]); ci.append(j); vals.append(float(m[i]))
    return sp.csr_matrix((vals, (ri, ci)), shape=(Bout.n, Bin.n))


def mult_matrix(P, Bin, Bout):
    """sparse (Bout.n x Bin.n) matrix of multiplication by the exact poly P."""
    ri, ci, vals = [], [], []
    for j, m in enumerate(Bin.mons):
        for e, c in P.items():
            s = tuple(e[t] + m[t] for t in range(6))
            ri.append(Bout.index[s]); ci.append(j); vals.append(float(c))
    return sp.csr_matrix((vals, (ri, ci)), shape=(Bout.n, Bin.n))


def mult_matrix_dense_poly(v, T, nout, nin):
    """sparse (nout x nin) matrix of multiplication by the dense poly v (over B1) using table T (B1 x B2)."""
    nz = np.nonzero(v)[0]
    if len(nz) == 0:
        return sp.csr_matrix((nout, nin))
    rows = T[nz, :].ravel()
    cols = np.tile(np.arange(nin, dtype=np.int32), len(nz))
    vals = np.repeat(v[nz], nin)
    return sp.csr_matrix((vals, (rows, cols)), shape=(nout, nin))


def complement_basis(G, rdim, verbose=False):
    """Orthonormal basis Z (n x rdim) of the orthogonal complement of the column space of the sparse (n x m) matrix G,
    assuming the complement has dimension exactly rdim (checked by the eigenvalue gap)."""
    t0 = time.time()
    M = (G @ G.T).toarray()
    w, V = sla.eigh(M, driver='evr')
    gap = (w[rdim - 1], w[rdim]) if rdim < len(w) else (w[-1], None)
    Z = V[:, :rdim]
    resid = np.abs(G.T @ Z).max() if G.shape[1] else 0.0
    if verbose:
        print(f'[complement] n={G.shape[0]} m={G.shape[1]} rdim={rdim} eig[rdim-1]={gap[0]:.3e} eig[rdim]={gap[1]} '
              f'max|G^T Z|={resid:.2e} ({time.time()-t0:.1f}s)', file=sys.stderr)
    return Z, gap, resid


# ----------------------------------------------------------------------------------------------------------------------
class JVP:
    def __init__(self, F, A, d, extra=None, stage2=True, compress_extra=200, seed=12345, verbose=True,
                 stage1=True, parD=0):
        """extra: optional dict with 'B' (list of dense vectors over S_d^even), 'n' (list of 6 dense vectors per B over
        the N0 bases) for non-Jacobian kernel directions."""
        t0 = time.time()
        self.F, self.A, self.d = F, A, d
        self.parD = parD
        self.verbose = verbose
        Fi = [pderiv(F, i) for i in range(NV)]
        self.Fi = Fi
        parN = [(1 - parD) if i < 3 else parD for i in range(NV)]   # Y-parity of N0_i (A Y-odd, D parity parD)
        self.parN = parN
        self.BD = Basis(d, parD)
        self.BN = [Basis(d + 1, parN[i]) for i in range(NV)]
        self.BDi = [Basis(d - 1, parD if i < 3 else 1 - parD) for i in range(NV)]
        self.Bdiv = Basis(d, 1 - parD)
        self.Bout = Basis(2 * d, 1)
        self.Bw = [Basis(d - 5, parD if l < 3 else 1 - parD) for l in range(NV)]
        nout = self.Bout.n
        # --- parameter columns
        cols_D, cols_N = [], [[] for _ in range(NV)]
        labels = []
        for l in range(NV):
            for mu in self.Bw[l].mons:
                mono = {mu: Fraction(1)}
                cols_D.append(self.BD.vec(mul(mono, Fi[l])))
                for i in range(NV):
                    cols_N[i].append(self.BN[i].vec(mul(A, mono)) if i == l else np.zeros(self.BN[i].n))
                labels.append(('w', l, mu))
        self.nw = len(cols_D)
        if extra is not None:
            for k, Bk in enumerate(extra['B']):
                cols_D.append(np.asarray(Bk, float))
                for i in range(NV):
                    cols_N[i].append(np.asarray(extra['n'][k][i], float))
                labels.append(('a', k, None))
        self.labels = labels
        self.K = len(cols_D)
        self.U_D = np.array(cols_D).T                        # nD x K
        self.U_N = [np.array(cols_N[i]).T for i in range(NV)]  # nN_i x K
        # derivative maps
        self.Dmat_N = [deriv_matrix(self.BN[i], self.Bdiv, i) for i in range(NV)]       # d_i N0_i -> S_d^odd
        self.Dmat_D = [deriv_matrix(self.BD, self.BDi[i], i) for i in range(NV)]        # d_i D -> S_{d-1}
        self.U_div = sum(self.Dmat_N[i] @ self.U_N[i] for i in range(NV))               # ndiv x K (dense)
        self.U_div = np.asarray(self.U_div)
        self.U_Di = [np.asarray(self.Dmat_D[i] @ self.U_D) for i in range(NV)]          # nDi x K
        # product tables
        self.T1 = prod_table(self.BD, self.Bdiv, self.Bout)
        self.T2 = [prod_table(self.BN[i], self.BDi[i], self.Bout) for i in range(NV)]
        # --- spurious directions: u with D(u) in F * S_{d-6}^even, and ker(U_D)
        self.P = self._reduced_param_space()
        self.dimP = self.P.shape[1]
        # --- stage 1 projector
        self.Z = None
        if stage1:
            gens = []
            for i in range(NV):
                par_m = (1 - (0 if i < 3 else 1)) % 2
                Bm = Basis(2 * d - 5, par_m)
                gens.append(mult_matrix(Fi[i], Bm, self.Bout))
            G = sp.hstack(gens).tocsr()
            self.G_J = G
            rJ = self._hilbert_odd(2 * d)
            self.rJ = rJ
            self.Z, self.gap, self.zres = complement_basis(G, rJ, verbose=verbose)
        # --- stage 2 data
        self.stage2 = stage2
        if stage2:
            self._build_stage2(compress_extra, seed)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[jvp d={d}] K={self.K} (w {self.nw}) dimP={self.dimP} nout={nout} '
                  f'{"rJ=%d" % self.rJ if stage1 else ""} '
                  f'{"nc=%d nQ=%d k=%d" % (self.nc, self.nQ, self.kS) if stage2 else ""} built in {self.build_time:.1f}s',
                  file=sys.stderr)

    # Hilbert function of R = S/J(F) for a regular sequence of six quintics, Y-odd part, via Gorenstein duality
    # with the parity-split monomial counts: dim R_k^odd = |S_k^odd| - dim J_k^odd.  We compute dim J_k^odd exactly by
    # the parity-split Hilbert series (1 - t^5)^3 (1 - u t^5)^3 / ((1 - t)^3 (1 - u t)^3), u^2 = 1.
    @staticmethod
    def _hilbert_odd(k):
        # series in t with coefficients (even, odd) pairs
        N = k + 1
        even = np.zeros(N); odd = np.zeros(N)
        even[0] = 1.0
        def mul_series(e1, o1, e2, o2):
            e = np.zeros(N); o = np.zeros(N)
            for i in range(N):
                if e1[i] == 0 and o1[i] == 0:
                    continue
                e[i:] += e1[i] * e2[:N - i] + o1[i] * o2[:N - i]
                o[i:] += e1[i] * o2[:N - i] + o1[i] * e2[:N - i]
            return e, o
        # 1/(1-t)^3 : coefficients C(n+2,2) all even; 1/(1-ut)^3: coefficient C(n+2,2) u^n
        inv_t = np.array([math.comb(n + 2, 2) for n in range(N)], float)
        e, o = mul_series(even, odd, inv_t, np.zeros(N))
        ue = np.array([math.comb(n + 2, 2) if n % 2 == 0 else 0 for n in range(N)], float)
        uo = np.array([math.comb(n + 2, 2) if n % 2 == 1 else 0 for n in range(N)], float)
        e, o = mul_series(e, o, ue, uo)
        # (1 - t^5)^3 = 1 - 3t^5 + 3t^10 - t^15 (even); (1 - u t^5)^3 = 1 - 3u t^5 + 3 t^10 - u t^15
        p1e = np.zeros(N); p1e[0] = 1
        for pw, cf in ((5, -3), (10, 3), (15, -1)):
            if pw < N: p1e[pw] += cf
        p2e = np.zeros(N); p2o = np.zeros(N); p2e[0] = 1
        for pw, cf, par in ((5, -3, 1), (10, 3, 0), (15, -1, 1)):
            if pw < N:
                (p2o if par else p2e)[pw] += cf
        e, o = mul_series(e, o, p1e, np.zeros(N))
        e, o = mul_series(e, o, p2e, p2o)
        return int(round(o[k]))

    def _reduced_param_space(self):
        """Orthonormal basis of the complement of {u : D(u) in F S_{d-6}} + ker U_D inside R^K."""
        d = self.d
        K = self.K
        blocks = [self.U_D]
        if d >= 6:
            BG = Basis(d - 6, self.parD)
            MF = mult_matrix(self.F, BG, self.BD).toarray()
            blocks.append(-MF)
        Mb = np.hstack(blocks)
        ns = sla.null_space(Mb)
        S = ns[:K, :]
        if S.size == 0 or S.shape[1] == 0:
            return np.eye(K)
        U, s, _ = np.linalg.svd(S, full_matrices=True)
        r = int(np.sum(s > 1e-9 * max(1.0, s[0])))
        self.spurious_dim = r
        return U[:, r:]

    def _build_stage2(self, compress_extra, seed):
        d = self.d
        Fi = self.Fi
        # Koszul coefficient columns c_ij (i<j), degree d-4, parity: N_i parity = parN[i]; c_ij F_j must have parN[i]
        self.cpairs = []
        P1_cols = []      # over S_d^odd : d_i(mu) F_j - d_j(mu) F_i
        P2_cols = [[] for _ in range(NV)]  # over BN[t]? no: mu*F_j lives in degree d+1 with parity parN[i]
        Bc = {}
        for i in range(NV):
            for j in range(i + 1, NV):
                par_c = (self.parN[i] - (0 if j < 3 else 1)) % 2
                Bc[(i, j)] = Basis(d - 4, par_c)
        self.Bc = Bc
        # column polys: T_D(c) = sum_{i<j} [ D * (d_i c F_j - d_j c F_i) - c (D_i F_j - D_j F_i) ]
        # = Mult_D(P1 c) - sum_t Mult_{D_t}(P2_t c)  with P2_t c = sum over pairs containing t of (+/-) c F_other
        self.Bodd_d = Basis(d, 1 - self.parD)
        P1 = sp.lil_matrix((self.Bodd_d.n, 0))
        p1_cols = []
        p2_cols = [[] for _ in range(NV)]   # over BN[t] (degree d+1, parity parN[t])
        clabels = []
        for (i, j), B in Bc.items():
            for mu in B.mons:
                mono = {mu: Fraction(1)}
                p1 = add(mul(pderiv(mono, i), Fi[j]), scale(mul(pderiv(mono, j), Fi[i]), -1))
                p1_cols.append(self.Bodd_d.vec(p1))
                for t in range(NV):
                    if t == i:
                        p2_cols[t].append(self.BN[i].vec(mul(mono, Fi[j])))        # - c F_j D_i  -> coefficient of D_i is c F_j
                    elif t == j:
                        p2_cols[t].append(self.BN[j].vec(scale(mul(mono, Fi[i]), -1)))  # + c F_i D_j
                    else:
                        p2_cols[t].append(np.zeros(self.BN[t].n))
                clabels.append((i, j, mu))
        self.clabels = clabels
        self.nc = len(clabels)
        self.P1 = sp.csr_matrix(np.array(p1_cols).T) if self.nc else sp.csr_matrix((self.Bodd_d.n, 0))
        self.P2 = [sp.csr_matrix(np.array(p2_cols[t]).T) for t in range(NV)]
        # product tables for D (even, deg d) x S_d^odd and D_t (deg d-1) x BN[t] (deg d+1)
        self.T_c1 = prod_table(self.BD, self.Bodd_d, self.Bout)
        self.T_c2 = [prod_table(self.BDi[t], self.BN[t], self.Bout) for t in range(NV)]
        # F Q block and its orthogonal complement compression
        BQ = Basis(2 * d - 6, 1)
        self.BQ = BQ
        self.nQ = BQ.n
        MF = mult_matrix(self.F, BQ, self.Bout)
        self.MF = MF
        QF, _ = sla.qr(MF.toarray(), mode='economic')
        self.QF = QF
        rng = np.random.default_rng(seed)
        k = min(self.Bout.n, self.nc + compress_extra)
        self.kS = k
        R = rng.standard_normal((k, self.Bout.n)) / math.sqrt(self.Bout.n)
        self.S = R - (R @ QF) @ QF.T          # S annihilates F*S_{2d-6}

    # ------------------------------------------------------------------------------------------------------------------
    def polys(self, u):
        D = self.U_D @ u
        N = [self.U_N[i] @ u for i in range(NV)]
        div = self.U_div @ u
        Di = [self.U_Di[i] @ u for i in range(NV)]
        return D, N, div, Di

    def B_vec(self, u):
        """B_D(N0)(u) = D div N0 - sum_i N0_i d_i D  as a dense vector over Bout."""
        D, N, div, Di = self.polys(u)
        nout = self.Bout.n
        B = tmul(self.T1, D, div, nout)
        for i in range(NV):
            B -= tmul(self.T2[i], N[i], Di[i], nout)
        return B, (D, N, div, Di)

    def B_grad_T(self, u, r, data=None):
        """gradient of r . B(u,u) with respect to u (B bilinear: d/du [B(u,u)] v = B(v,u) + B(u,v))."""
        if data is None:
            _, data = self.B_vec(u)
        D, N, div, Di = data
        g = np.zeros(self.K)
        R1 = r[self.T1]                      # nD x ndiv
        g += self.U_D.T @ (R1 @ div)         # dB/dD direction
        g += self.U_div.T @ (R1.T @ D)       # dB/ddiv
        for i in range(NV):
            R2 = r[self.T2[i]]               # nN_i x nDi
            g -= self.U_N[i].T @ (R2 @ Di[i])
            g -= self.U_Di[i].T @ (R2.T @ N[i])
        return g

    # ---- stage 1 ----
    def r1(self, v):
        u = self.P @ v
        B, data = self.B_vec(u)
        r = self.Z.T @ B
        return np.concatenate([r, [u @ u - 1.0]]), B, data

    def r1_jac(self, v):
        u = self.P @ v
        B, data = self.B_vec(u)
        D, N, div, Di = data
        nout = self.Bout.n
        # d/du B(u,u) = Mult_div U_D + Mult_D U_div - sum_i (Mult_{D_i} U_N[i] + Mult_{N_i} U_Di[i])  (nout x K)
        # Mult_q built from product tables: for T1 (BD x Bdiv): mult by div acts on D-vectors; mult by D acts on div-vectors
        M = mult_matrix_dense_poly(div, self.T1.T.copy(), nout, self.BD.n) @ self.U_D
        M = M + mult_matrix_dense_poly(D, self.T1, nout, self.Bdiv.n) @ self.U_div
        for i in range(NV):
            M = M - mult_matrix_dense_poly(Di[i], self.T2[i].T.copy(), nout, self.BN[i].n) @ self.U_N[i]
            M = M - mult_matrix_dense_poly(N[i], self.T2[i], nout, self.BDi[i].n) @ self.U_Di[i]
        M = np.asarray(M)
        cols = self.Z.T @ M
        J = np.vstack([cols, 2.0 * u[None, :]])
        return J @ self.P

    def f1_and_grad(self, v):
        u = self.P @ v
        B, data = self.B_vec(u)
        r = self.Z.T @ B
        f = float(r @ r)
        rfull = self.Z @ r
        g = 2.0 * self.B_grad_T(u, rfull, data)
        nrm = u @ u - 1.0
        f += nrm * nrm
        g += 4.0 * nrm * u
        return f, self.P.T @ g

    def stage1_rel(self, v):
        u = self.P @ v
        B, _ = self.B_vec(u)
        r = self.Z.T @ B
        return float(r @ r) / max(float(B @ B), 1e-300), float(np.sqrt(B @ B))

    # ---- stage 2 ----
    def T_matrix(self, D, Di):
        """dense (nout x nc) matrix of c -> T_D(c) for the given dense D, D_t."""
        nout = self.Bout.n
        M = mult_matrix_dense_poly(D, self.T_c1, nout, self.Bodd_d.n) @ self.P1
        for t in range(NV):
            M = M - mult_matrix_dense_poly(Di[t], self.T_c2[t], nout, self.BN[t].n) @ self.P2[t]
        return M

    def r2_value(self, v, exact=False, cond=1e-10):
        """stage-2 residual: min_c || S (B + T_D c) || (compressed) or exact min over (c, Q) without compression."""
        u = self.P @ v
        B, data = self.B_vec(u)
        D, N, div, Di = data
        T = self.T_matrix(D, Di)
        if exact:
            Bp = B - self.QF @ (self.QF.T @ B)
            Tp = np.asarray(T.todense()) if sp.issparse(T) else T
            Tp = Tp - self.QF @ (self.QF.T @ Tp)
            c, *_ = sla.lstsq(Tp, -Bp, cond=cond, lapack_driver='gelsd', check_finite=False)
            r = Bp + Tp @ c
            return float(r @ r), float(B @ B), c, r
        SB = self.S @ B
        ST = np.asarray(self.S @ T) if sp.issparse(T) else self.S @ T
        c, *_ = sla.lstsq(ST, -SB, cond=cond, lapack_driver='gelsd', check_finite=False)
        r = SB + ST @ c
        return float(r @ r), float(B @ B), c, r

    def f2_and_grad(self, v, cond=1e-10):
        """variable-projection objective ||S(B + T c*)||^2 + (|u|^2-1)^2 and its gradient (Golub-Pereyra)."""
        u = self.P @ v
        B, data = self.B_vec(u)
        D, N, div, Di = data
        T = self.T_matrix(D, Di)
        SB = self.S @ B
        ST = np.asarray(self.S @ T) if sp.issparse(T) else self.S @ T
        c, *_ = sla.lstsq(ST, -SB, cond=cond, lapack_driver='gelsd', check_finite=False)
        r = SB + ST @ c
        f = float(r @ r)
        rfull = self.S.T @ r                 # nout
        # gradient: 2 r . d/du [S (B(u,u) + T_D(u) c)]  (c fixed: VarPro)
        g = 2.0 * self.B_grad_T(u, rfull, data)
        # T_D(u) c = Mult_D (P1 c) - sum_t Mult_{D_t} (P2_t c): linear in (D, D_t)
        p1c = self.P1 @ c
        R1 = rfull[self.T_c1]                # nD x nodd
        g += 2.0 * (self.U_D.T @ (R1 @ p1c))
        for t in range(NV):
            p2c = self.P2[t] @ c
            R2 = rfull[self.T_c2[t]]         # nDi x nN_t
            g -= 2.0 * (self.U_Di[t].T @ (R2 @ p2c))
        nrm = u @ u - 1.0
        f += nrm * nrm
        g += 4.0 * nrm * u
        return f, self.P.T @ g

    # ---- linear restriction of the family: D(u) in W, W = invariants of index permutations of the exponent tuple ----
    def restrict_to_invariants(self, perms, tol=1e-9):
        """Restrict the parameter space to u with D(u) invariant under the given exponent-index permutations."""
        rows = []
        for j, m in enumerate(self.BD.mons):
            for pm in perms:
                y = tuple(m[pm[i]] for i in range(6))
                k = self.BD.index[y]
                if k != j:
                    r = np.zeros(self.BD.n); r[j] = 1.0; r[k] = -1.0; rows.append(r)
        if not rows:
            return self.dimP
        C = np.array(rows)
        M = C @ (self.U_D @ self.P)
        N = sla.null_space(M, rcond=tol)
        self.P = self.P @ N
        self.dimP = self.P.shape[1]
        return self.dimP

    # ---- exact reconstruction ----
    def D_exact_from_u(self, u, max_den=1000):
        """rational D from a parameter vector u via per-coordinate continued fractions (normalized max |u| = 1)."""
        u = np.asarray(u, float)
        j = int(np.argmax(np.abs(u)))
        un = u / u[j]
        fr = [Fraction(float(x)).limit_denominator(max_den) for x in un]
        D = {}
        for k, (lab, f) in enumerate(zip(self.labels, fr)):
            if f == 0:
                continue
            kind, l, mu = lab
            if kind == 'w':
                D = add(D, scale(mul({mu: Fraction(1)}, self.Fi[l]), f))
            else:
                raise ValueError('exact reconstruction of non-J directions needs the exact kernel basis')
        return D, fr


# ----------------------------------------------------------------------------------------------------------------------
# drivers
# ----------------------------------------------------------------------------------------------------------------------
def run_stage1(jvp, n_restarts, seed, out_path, tag, maxiter=200, method='lm', verbose=False, log_every=1,
               time_limit=None):
    rng = np.random.default_rng(seed)
    results = []
    t0 = time.time()
    with open(out_path, 'a') as f:
        for it in range(n_restarts):
            if time_limit and time.time() - t0 > time_limit:
                break
            v0 = rng.standard_normal(jvp.dimP); v0 /= np.linalg.norm(v0)
            t1 = time.time()
            try:
                res = least_squares(jvp.r1_fun, v0, jac=jvp.r1_jacfun, method=method, max_nfev=maxiter,
                                    xtol=1e-15, ftol=1e-15, gtol=1e-15)
                v = res.x
            except Exception as ex:
                print(f'[stage1] restart {it} failed: {ex}', file=sys.stderr)
                continue
            v = v / np.linalg.norm(v)
            rel, bnorm = jvp.stage1_rel(v)
            rec = {'stage': 1, 'd': jvp.d, 'parity': jvp.parD, 'tag': tag, 'dimP': jvp.dimP, 'seed': seed, 'restart': it, 'rel': rel, 'bnorm': bnorm,
                   'nfev': int(res.nfev), 'seconds': round(time.time() - t1, 2), 'v': v.tolist()}
            results.append(rec)
            f.write(json.dumps(rec) + '\n'); f.flush()
            if verbose and it % log_every == 0:
                print(f'[stage1 d={jvp.d}] restart {it}: rel={rel:.3e} |B|={bnorm:.3e} nfev={res.nfev} '
                      f'({time.time()-t0:.0f}s)', file=sys.stderr)
    return results


def _attach_funs(jvp):
    def r1_fun(v):
        r, _, _ = jvp.r1(v)
        return r
    jvp.r1_fun = r1_fun
    jvp.r1_jacfun = jvp.r1_jac


def run_stage2(jvp, starts, seed, out_path, tag, maxiter=300, verbose=False, exact_check=True, time_limit=None):
    results = []
    t0 = time.time()
    with open(out_path, 'a') as f:
        for it, v0 in enumerate(starts):
            if time_limit and time.time() - t0 > time_limit:
                break
            t1 = time.time()
            res = minimize(jvp.f2_and_grad, v0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': 1e-14, 'ftol': 1e-18, 'maxcor': 30})
            v = res.x / np.linalg.norm(res.x)
            f2, b2, c, r = jvp.r2_value(v)
            rel = f2 / max(b2, 1e-300)
            rec = {'stage': 2, 'd': jvp.d, 'parity': jvp.parD, 'tag': tag, 'dimP': jvp.dimP, 'seed': seed, 'restart': it, 'rel_compressed': rel,
                   'f2': f2, 'bnorm2': b2, 'nit': int(res.nit), 'nfev': int(res.nfev),
                   'seconds': round(time.time() - t1, 2), 'v': v.tolist()}
            if jvp.Z is not None:
                rec['rel_stage1'] = jvp.stage1_rel(v)[0]
            if exact_check and rel < 1e-6:
                f2e, b2e, ce, re_ = jvp.r2_value(v, exact=True)
                rec['rel_exact'] = f2e / max(b2e, 1e-300)
            results.append(rec)
            f.write(json.dumps(rec) + '\n'); f.flush()
            if verbose:
                print(f'[stage2 d={jvp.d}] start {it}: rel={rel:.3e} (stage1 {rec.get("rel_stage1", float("nan")):.2e}) '
                      f'nit={res.nit} ({time.time()-t0:.0f}s)', file=sys.stderr)
    return results


_SHARED = {}


def _limit_threads(n=1):
    try:
        import threadpoolctl
        threadpoolctl.threadpool_limits(n)
    except Exception:
        pass


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--instance', default=LANE)
    ap.add_argument('--fermat', action='store_true', help='control: use F0 = sum Z^6 instead of the instance F')
    ap.add_argument('--stage', default='1', help='1, 2, or 12 (stage 2 from stage-1 minima + random)')
    ap.add_argument('--restarts', type=int, default=10)
    ap.add_argument('--restarts2', type=int, default=5, help='random stage-2 restarts (stage 12/2), total over workers')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--tag', default='test')
    ap.add_argument('--out', default=None)
    ap.add_argument('--maxiter', type=int, default=200)
    ap.add_argument('--maxiter2', type=int, default=300)
    ap.add_argument('--thresh1', type=float, default=1e-3, help='stage-1 rel below which stage 2 is run')
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--extra', default=None, help='npz with real non-J kernel basis B (nD x m) and lifts n_i (nN_i x m)')
    ap.add_argument('--quiet', action='store_true')
    ap.add_argument('--time-limit', type=float, default=None, help='seconds per stage per worker')
    ap.add_argument('--seeds-from', default=None, help='jsonl of stage-1 records to seed stage 2 (sorted by rel, below thresh1)')
    ap.add_argument('--max-seeds', type=int, default=1000)
    ap.add_argument('--worker-threads', type=int, default=1)
    ap.add_argument('--refine1', action='store_true', help='with --seeds-from: refine recorded stage-1 minima by long LM runs')
    ap.add_argument('--parity', type=int, default=0, help='Y-parity of D (0 even, 1 odd)')
    ap.add_argument('--restrict', default=None, help='comma list of symmetry generators (rho, sigma, taux) forcing D invariant')
    a = ap.parse_args()
    F, A, q, _ = core.load_instance(a.instance)
    if a.fermat:
        F = {tuple(6 if j == i else 0 for j in range(6)): Fraction(1) for i in range(6)}
    out = a.out or f'w2_jvp_d{a.d}_{"odd" if a.parity else "even"}_{a.tag}.jsonl'
    if a.restrict:
        out = out.replace('.jsonl', f'_{a.restrict.replace(",", "")}.jsonl')
    extra = None
    if a.extra:
        z = np.load(a.extra)
        extra = {'B': list(z['B'].T), 'n': [[z[f'n{i}'][:, k] for i in range(6)] for k in range(z['B'].shape[1])]}
    verbose = not a.quiet
    jvp = JVP(F, A, a.d, extra=extra, stage2=('2' in a.stage), stage1=('1' in a.stage), verbose=verbose, parD=a.parity)
    _attach_funs(jvp)
    if a.restrict:
        gens = {'rho': (1, 2, 0, 4, 5, 3), 'sigma': (3, 4, 5, 0, 1, 2), 'taux': (1, 0, 2, 4, 3, 5)}
        perms = [gens[s.strip()] for s in a.restrict.split(',') if s.strip()]
        dp = jvp.restrict_to_invariants(perms)
        chk = float('nan')
        if dp > 0:
            rng0 = np.random.default_rng(0); u0 = jvp.P @ rng0.standard_normal(dp); D0 = jvp.U_D @ u0; Dp0 = jvp.BD.poly(D0)
            chk = max(abs(Dp0.get(tuple(m[pm[i]] for i in range(6)), 0.0) - c) for m, c in Dp0.items() for pm in perms) / max(np.linalg.norm(D0), 1e-300)
        print(f'[jvp] restricted to {a.restrict}-invariant denominators: parameter dim {dp}; invariance check {chk:.2e}', file=sys.stderr)
        if dp == 0:
            print('[jvp] empty restricted family; nothing to search', file=sys.stderr); return
    _SHARED['jvp'] = jvp
    if a.procs > 1:
        import multiprocessing as mp
        ctx = mp.get_context('fork')
        seeds = [a.seed * 1000 + s for s in range(a.procs)]
        per = max(1, a.restarts // a.procs)
        per2 = max(1, a.restarts2 // a.procs)
        args = [(a, s, per, per2, out.replace('.jsonl', f'_s{s}.jsonl')) for s in seeds]
        with ctx.Pool(a.procs) as pool:
            pool.starmap(_worker, args, chunksize=1)
        return
    _worker(a, a.seed, a.restarts, a.restarts2, out)


def _worker(a, seed, restarts, restarts2, out):
    _limit_threads(a.worker_threads)
    verbose = not a.quiet
    jvp = _SHARED['jvp']
    rng = np.random.default_rng(seed + 777)
    s1 = []
    if '1' in a.stage and a.seeds_from and a.refine1:
        # seeded stage-1 refinement: long LM runs from the best recorded stage-1 points
        recs = []
        for path in a.seeds_from.split(','):
            for line in open(path):
                try:
                    r = json.loads(line)
                except Exception:
                    continue
                if r.get('stage') == 1 and r.get('d') == a.d and r.get('parity', 0) == jvp.parD and r['rel'] < a.thresh1:
                    recs.append(r)
        recs.sort(key=lambda r: r['rel'])
        with open(out, 'a') as f:
            for it, r0 in enumerate(recs[:a.max_seeds]):
                v0 = np.array(r0['v']); t1 = time.time()
                res = least_squares(jvp.r1_fun, v0, jac=jvp.r1_jacfun, method='lm', max_nfev=a.maxiter, xtol=1e-15, ftol=1e-15, gtol=1e-15)
                v = res.x / np.linalg.norm(res.x); rel, bnorm = jvp.stage1_rel(v)
                rec = {'stage': 1, 'd': jvp.d, 'parity': jvp.parD, 'tag': a.tag + '-refine', 'dimP': jvp.dimP, 'seed': seed, 'restart': it,
                       'rel': rel, 'rel_seed': r0['rel'], 'bnorm': bnorm, 'nfev': int(res.nfev), 'status': int(res.status),
                       'seconds': round(time.time() - t1, 2), 'v': v.tolist()}
                s1.append(rec); f.write(json.dumps(rec) + '\n'); f.flush()
                if verbose:
                    print(f'[refine1 d={jvp.d}] seed {it}: rel {r0["rel"]:.3e} -> {rel:.3e} nfev={res.nfev} status={res.status} ({time.time()-t1:.0f}s)', file=sys.stderr)
    elif '1' in a.stage:
        s1 = run_stage1(jvp, restarts, seed, out, a.tag, maxiter=a.maxiter, verbose=verbose, time_limit=a.time_limit)
    if '2' in a.stage:
        starts = []
        if s1:
            good = sorted([r for r in s1 if r['rel'] < a.thresh1], key=lambda r: r['rel'])
            starts += [np.array(r['v']) for r in good]
        if a.seeds_from:
            recs = []
            for path in a.seeds_from.split(','):
                for line in open(path):
                    try:
                        r = json.loads(line)
                    except Exception:
                        continue
                    if r.get('stage') == 1 and r.get('d') == a.d and r['rel'] < a.thresh1:
                        recs.append(r)
            recs.sort(key=lambda r: r['rel'])
            mine = [r for k, r in enumerate(recs[:a.max_seeds]) if k % max(1, a.procs) == (seed % 1000) % max(1, a.procs)]
            starts += [np.array(r['v']) for r in mine]
        for _ in range(restarts2):
            v0 = rng.standard_normal(jvp.dimP); v0 /= np.linalg.norm(v0); starts.append(v0)
        run_stage2(jvp, starts, seed, out, a.tag, maxiter=a.maxiter2, verbose=verbose, time_limit=a.time_limit)


if __name__ == '__main__':
    main()
