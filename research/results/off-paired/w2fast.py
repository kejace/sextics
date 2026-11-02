"""Fast reduced sector-restricted variable projection (vectorized; replaces w2lib.ReducedVarPro).

Bilinear map  B(D, N) = D div N - N . grad D  on coefficient vectors:
   for a D-monomial m_D and an N-monomial (slot i, m_N):  coefficient (e_i(m_N) - e_i(m_D)) at row m_D + m_N - e_i.
Tables tab_row / tab_coef (nD x nN) make B(D, N) a 300k-term scatter-add and T(D) = B(D, .) a sparse matrix.
Residual  rho(a) = B(D(a), S s) - F Q + B(D(a), N(a)),  minimized over (s, Q) by least squares; D(a) = DB^T a,
N(a) = NB^T a (E1 lifts).  The constant Q-block is optionally projected out with a thin QR (qproj=True).
Parameters u live on the unit sphere of an orthonormalized D-space (P = Gram^{-1/2}); rel = f / bref2 with
bref2 = mean |B(D(u), N(u))|^2 over random unit u (fixed reference scale; the raw |b(a)| may vanish at a solution).
"""
import sys, time, math
sys.path.insert(0, '/tmp/hodge/results/off-paired')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
import w2lib as W
from hodge import core, fixed_d, symmetry, varpro
from hodge.core import NV, monomials, mul, add, sub, scale, shift, pderiv


def _keys(E, base):
    E = np.asarray(E, dtype=np.int64)
    k = np.zeros(E.shape[0], dtype=np.int64)
    for i in range(NV):
        k = k * base + E[:, i]
    return k


class FastReducedVarPro(varpro.VarProFamily):
    def __init__(self, F, A, q, basis, nlifts, koszul, c, verbose=False, cond=1e-9, label='', qproj=True, nref=16, seed=0):
        t0 = time.time()
        self.cond = cond; self.F, self.A, self.q = F, A, q
        self.basis = [dict(B) for B in basis]; self.label = label; self.c = c
        K = len(basis); self.K = K
        degs = {core.degree(B) for B in basis}; assert len(degs) == 1, degs
        d = degs.pop(); self.d = d
        self.pD = core.parity_of(basis[0])
        cA = W.poly_chi(A); cQ = W.cadd(cA, c, c)
        self.degN, self.degQ = d + 6 * q - 5, 2 * d - 6
        # D monomials / N monomials / rows
        Dm = W.mons_chi(d, c); nD = len(Dm); Didx = {m: j for j, m in enumerate(Dm)}
        Nm = []; slot = []
        for i in range(NV):
            for m in W.mons_chi(d + 1, W.cadd(cA, c, W.WZ[i])):
                Nm.append((i, m)); slot.append(i)
        nN = len(Nm); Nidx = {im: j for j, im in enumerate(Nm)}
        Rm = W.mons_chi(2 * d, cQ); R = len(Rm)
        base = 2 * d + 2
        rkeys = _keys(Rm, base); order = np.argsort(rkeys); rkeys_sorted = rkeys[order]
        self.rows = R; self.row_mons = Rm
        ED = np.array(Dm, dtype=np.int64); EN = np.array([m for _, m in Nm], dtype=np.int64); slot = np.array(slot)
        # tables
        # e_i(m_N) - e_i(m_D) with i = slot of the N column: shape nD x nN
        coef = (EN[np.arange(nN), slot][None, :] - ED[:, slot]).astype(float)   # ED[:, slot] is nD x nN
        # row monomial m_D + m_N - e_slot
        Esum = ED[:, None, :] + EN[None, :, :]                                   # nD x nN x 6
        Esum[np.arange(nD)[:, None], np.arange(nN)[None, :], slot[None, :]] -= 1
        flat = Esum.reshape(-1, NV)
        valid = np.all(flat >= 0, axis=1)
        keys = _keys(np.where(valid[:, None], flat, 0), base)
        pos = np.searchsorted(rkeys_sorted, keys)
        pos = np.minimum(pos, R - 1)
        found = valid & (rkeys_sorted[pos] == keys)
        rowidx = np.where(found, order[pos], 0)
        coef_flat = coef.reshape(-1) * found
        self.tab_row = rowidx.astype(np.int64); self.tab_coef = coef_flat
        self.nD, self.nN = nD, nN
        self.tab_col = np.tile(np.arange(nN), nD)
        # D basis and N lifts as dense matrices
        DB = np.zeros((K, nD)); NB = np.zeros((K, nN))
        for k, B in enumerate(self.basis):
            for e, v in B.items():
                DB[k, Didx[e]] = float(v)
        for k, N in enumerate(nlifts):
            for i in range(NV):
                for e, v in N[i].items():
                    NB[k, Nidx[(i, e)]] = float(v)
        self.DB, self.NB = DB, NB
        # Koszul generators as sparse nS x nN
        nS = len(koszul); self.nS = nS; self.nN_cols = nS
        ri, ci, vals = [], [], []
        for j, S in enumerate(koszul):
            for i in range(NV):
                for e, v in S[i].items():
                    ri.append(j); ci.append(Nidx[(i, e)]); vals.append(float(v))
        self.SB = sp.csr_matrix((vals, (ri, ci)), shape=(nS, nN)) if nS else None
        self.SBT = self.SB.T.tocsr() if nS else None
        # Q block (dense rows x nQ): -F m
        Qm = W.mons_chi(self.degQ, cQ) if self.degQ >= 0 else []
        self.Qmons = Qm; nQ = len(Qm); self.nQ = nQ
        ridx = {m: j for j, m in enumerate(Rm)}
        MQ = np.zeros((R, nQ))
        for j, m in enumerate(Qm):
            for e, v in shift(F, m, -1).items():
                MQ[ridx[e], j] = float(v)
        self.MQ = MQ
        self.qproj = qproj and nQ > 0
        if self.qproj:
            Qq, Rr = np.linalg.qr(MQ, mode='reduced')
            # drop dependent columns (F m are independent, but be safe)
            self.Qq = Qq
        self.ncols = nS + nQ
        # orthonormalized parameters: P = Gram^{-1/2}
        G = DB @ DB.T
        w, V = np.linalg.eigh(G)
        w = np.maximum(w, 1e-14 * w.max())
        self.P = V @ np.diag(w ** -0.5) @ V.T
        self.dimP = K; self.excluded_dim = 0
        # reference scale
        rng = np.random.default_rng(seed)
        bs = []
        for _ in range(nref):
            u = rng.standard_normal(K); a = self.P @ u; a /= np.linalg.norm(a)
            bs.append(np.sum(self._b(a) ** 2))
        self.bref2 = float(np.mean(bs)) if bs else 1.0
        self.build_time = time.time() - t0; self._nfev = 0
        if verbose:
            print(f'[frvp {label}] d={d} chi={c} K={K} rows={R} nD={nD} nN={nN} S={nS} Q={nQ} tab={nD*nN} '
                  f'bref2={self.bref2:.3e} built in {self.build_time:.1f}s', file=sys.stderr)

    # ---- bilinear machinery
    def _vecs(self, a):
        return self.DB.T @ a, self.NB.T @ a

    def _B(self, Dv, Nv):
        vals = (Dv[:, None] * Nv[None, :]).reshape(-1) * self.tab_coef
        return np.bincount(self.tab_row, weights=vals, minlength=self.rows)

    def _b(self, a):
        Dv, Nv = self._vecs(a)
        return self._B(Dv, Nv)          # rho(a) = M(a) x + b(a)  with b = B(D(a), N(a))

    def _T(self, Dv):
        data = (Dv[:, None] * self.tab_coef.reshape(self.nD, self.nN)).reshape(-1)
        return sp.csr_matrix((data, (self.tab_row, self.tab_col)), shape=(self.rows, self.nN))

    def _G1(self, Nv, r):
        """G1[mD] = sum_mN Nv[mN] coef r[row]  -> gradient wrt D-coefficients."""
        vals = (self.tab_coef * r[self.tab_row]).reshape(self.nD, self.nN)
        return vals @ Nv

    def _G2(self, Dv, r):
        vals = (self.tab_coef * r[self.tab_row]).reshape(self.nD, self.nN)
        return Dv @ vals

    def assemble(self, a):
        Dv, Nv = self._vecs(a)
        b = -self._B(Dv, Nv)
        if self.nS:
            MS = (self._T(Dv) @ self.SBT).toarray()
            M = np.hstack([MS, self.MQ]) if self.nQ else MS
        else:
            M = self.MQ
        return M, b

    def _solve(self, a):
        """Least squares; returns (f, xS, r) with r the full residual M x - b (Q optimal)."""
        Dv, Nv = self._vecs(a)
        b = -self._B(Dv, Nv)
        if self.nS:
            MS = (self._T(Dv) @ self.SBT).toarray()
        else:
            MS = np.zeros((self.rows, 0))
        if self.qproj:
            Qq = self.Qq
            MSp = MS - Qq @ (Qq.T @ MS) if self.nS else MS
            bp = b - Qq @ (Qq.T @ b)
            if self.nS:
                xS, *_ = sla.lstsq(MSp, bp, cond=self.cond, lapack_driver='gelsd', check_finite=False)
                r = MSp @ xS - bp
            else:
                xS = np.zeros(0); r = -bp
        else:
            M = np.hstack([MS, self.MQ]) if self.nQ else MS
            x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
            r = M @ x - b; xS = x[:self.nS]
        self._nfev += 1
        return float(r @ r), xS, r, Dv, Nv

    def value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        f, xS, r, _, _ = self._solve(a)
        return f, xS, r

    def value_and_grad(self, u):
        u = np.asarray(u, float)
        a_full = self.P @ u; na = np.linalg.norm(a_full); a = a_full / na
        f, xS, r, Dv, Nv = self._solve(a)
        Ns = self.SBT @ xS if self.nS else np.zeros(self.nN)
        # d rho / d a_k = B(D_k, Ns + N(a)) + B(D(a), N_k)
        G1 = self._G1(Ns + Nv, r); G2 = self._G2(Dv, r)
        g = 2.0 * (self.DB @ G1 + self.NB @ G2)
        g_full = (g - float(g @ a) * a) / na
        return f, self.P.T @ g_full

    def bnorm(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        return float(np.linalg.norm(self._b(a)))

    def rel(self, f):
        return f / self.bref2

    def search(self, n_restarts=10, seed=0, maxiter=300, gtol=1e-12, verbose=False, init=None, stop_rel=1e-14):
        rng = np.random.default_rng(seed); out = []; t0 = time.time()
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = np.linalg.solve(self.P, np.asarray(init, float))
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            f, x, r = self.value(a)
            rel = self.rel(f)
            spur, srel = self.is_spurious(a)
            out.append({'value': f, 'rel': rel, 'a': a, 'x': x, 'spurious': spur, 'div_rel': srel,
                        'iters': int(res.nit), 'restart': it, 'bnorm': self.bnorm(a)})
            if verbose:
                print(f'[frvp] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} spurious={spur} ({time.time()-t0:.1f}s)', file=sys.stderr)
            if rel < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    def residual_vec(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        return self._solve(a)[2]

    def confirm(self, frac, verbose=False, group=None):
        D = self.D_exact(frac)
        if not D:
            return 'inconsistent', {'note': 'D = 0'}
        ok, info = fixed_d.e1_solvable(self.F, self.A, self.q, D)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails', 'e1': info, 'D': D}
        status, data = symmetry.solve_fixed_D_symmetric(self.F, self.A, self.q, D, verbose=verbose, group=group)
        data['D'] = D
        return status, data

    def info(self):
        return {'K': self.K, 'rows': self.rows, 'nS': self.nS, 'nQ': self.nQ, 'nD': self.nD, 'nN': self.nN,
                'd': self.d, 'chi': list(self.c), 'bref2': self.bref2, 'build_seconds': round(self.build_time, 1),
                'qproj': bool(self.qproj)}


def build_family(F, A, q, d, c, kind, hessian=False, verbose=False, spaces=None, qproj=True):
    sp_ = spaces or W.SectorSpaces(F, A, q, d, c, verbose=verbose)
    basis, lifts = [], []
    if kind in ('jac', 'full'):
        for tag, B in zip(sp_.jac_tags, sp_.jac_basis):
            basis.append(B); lifts.append(W.n_lift_exact_jac(A, q, tag))
    if kind in ('ker', 'full') and sp_.ker_basis:
        L = W.E1Lifter(F, A, q, d, c)
        for B in sp_.ker_basis:
            N, res = L.lift(B)
            if res > 1e-8:
                raise RuntimeError(f'E1 lift failed for kernel vector (res {res:.2e})')
            basis.append(B); lifts.append(N)
    if not basis:
        return None, sp_
    kos = [] if hessian else [N for _, N in W.koszul_generators(F, d, c, W.poly_chi(A))]
    fam = FastReducedVarPro(F, A, q, basis, lifts, kos, c, verbose=verbose, qproj=qproj,
                            label=f'{kind}{"-hess" if hessian else ""}')
    return fam, sp_
