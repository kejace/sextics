#!/usr/bin/env python3
"""Koszul-reduced variable projection (full E2) over an E1-admissible family, lane reflection-coupled-residual.

For D in the family with an exact E1 solution N0 (N0(F) = A D), every E1 solution is N0 + K_a with
K_a,i = sum_j a_ij F_j (a antisymmetric, a_ij in S_{d-4}), and
    B_D(K_a) = sum_{i<j} [ D (F_j d_i a_ij - F_i d_j a_ij) - a_ij (F_j D_i - F_i D_j) ].
The certificate exists iff  B_D(N0) + B_D(K_a) in (F)  for some a.  Hence for the family
D(c) = sum c_k D_k, N0(c) = sum c_k N_k:
    f(c) = min_a || Pi_F ( L(c) a - quad(c) ) ||^2,   quad(c) = B(D(c), N0(c)),  L(c) a = B_{D(c)}(K_a),
with Pi_F the orthogonal projector onto the complement of F * S_{2d-6}.  f is homogeneous of degree 4 in c;
we minimise on the unit sphere of the complement of the spurious subspace {c : D(c) in F*S}.
Index arrays make quad(c), L(c) and the gradient cost a few sparse passes; the dense work per evaluation is the
projection of L (rows x nQ x n_a) and one least-squares solve (rows x n_a).

CLI:  python3 kvp.py --d D --sec ab|e|o --family J|K|JK --restarts N --seed S [--maxiter 400] [--thresh 1e-9]
"""
from __future__ import annotations
import argparse, json, os, sys, time
from fractions import Fraction
import numpy as np
import scipy.linalg as sla
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rcr2
from rcr2 import (instance, chi, chi_var, sadd, sec_str, sectors_of, parse_sec, sector_monomials, log_event,
                  jacobian_pairs, kernel_pairs, reconstruct_poly, exact_confirm, _pack, HERE)
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv


STOP = {'flag': False}


def _graceful(signum, frame):
    STOP['flag'] = True


import signal
signal.signal(signal.SIGTERM, _graceful)
signal.signal(signal.SIGINT, _graceful)


class KoszulVarPro:
    def __init__(self, d, pairs, verbose=True, cond=1e-10, reduce='sparse', extra_rows=60, seed=12345):
        F, A, q, Fi, Hs = instance()
        self.d = d; self.pairs = pairs; self.K = len(pairs); self.cond = cond; self.reduce = reduce
        K = self.K
        t0 = time.time()
        dsecs = sorted({chi(e) for p in pairs for e in p[0]})
        self.dsecs = dsecs
        outsecs = sorted({sadd((0, 1), sadd(sk, sl)) for sk in dsecs for sl in dsecs})
        # --- monomial bookkeeping
        Dmons = sorted({e for p in pairs for e in p[0]})
        didx = {e: i for i, e in enumerate(Dmons)}
        nD = len(Dmons)
        Dmat = np.zeros((nD, K))
        for k, (Dk, Nk, tag) in enumerate(pairs):
            for e, v in Dk.items():
                Dmat[didx[e], k] = float(v)
        Ncols = []
        for i in range(NV):
            for s in sorted({sadd(sadd((0, 1), sd), chi_var(i)) for sd in dsecs}):
                Ncols += [(i, m) for m in sector_monomials(d + 1, s)]
        nidx = {c: j for j, c in enumerate(Ncols)}
        nN = len(Ncols)
        Nmat = np.zeros((nN, K))
        for k, (Dk, Nk, tag) in enumerate(pairs):
            for i in range(NV):
                for m, v in Nk[i].items():
                    Nmat[nidx[(i, m)], k] += float(v)
        outmons = []
        for s in outsecs:
            outmons += sector_monomials(2 * d, s)
        R = len(outmons)
        keys = _pack(outmons)
        order = np.argsort(keys); skeys = keys[order]

        def rows_of(E):
            kk = _pack(E)
            pos = np.clip(np.searchsorted(skeys, kk), 0, len(skeys) - 1)
            if not np.all(skeys[pos] == kk):
                raise RuntimeError('monomial outside the output sectors')
            return order[pos]
        # --- quad index arrays: entry (e, (i,m)) -> row(e + m - e_i), coef (m_i - e_i)
        Ed = np.array(Dmons, dtype=np.int64)
        QE, QN, QR, QC = [], [], [], []
        for i in range(NV):
            cols_i = [j for j, (ii, m) in enumerate(Ncols) if ii == i]
            if not cols_i:
                continue
            Em = np.array([Ncols[j][1] for j in cols_i], dtype=np.int64)
            E = Ed[:, None, :] + Em[None, :, :]
            E[:, :, i] -= 1
            coef = (Em[None, :, i] - Ed[:, None, i]).astype(float)
            valid = (E[:, :, i] >= 0) & (coef != 0)
            QE.append(np.broadcast_to(np.arange(nD)[:, None], valid.shape)[valid])
            QN.append(np.broadcast_to(np.array(cols_i)[None, :], valid.shape)[valid])
            QR.append(rows_of(E[valid]))
            QC.append(coef[valid])
        self.QE, self.QN, self.QR, self.QC = (np.concatenate(QE), np.concatenate(QN), np.concatenate(QR), np.concatenate(QC))
        # --- Koszul columns a_ij (i<j), h in S_{d-4} with chi(h) = (0,1) + s + chi(z_i) + chi(z_j), s in dsecs
        acols = []
        for i in range(NV):
            for j in range(i + 1, NV):
                hsecs = sorted({sadd(sadd((0, 1), s), sadd(chi_var(i), chi_var(j))) for s in dsecs})
                for s in hsecs:
                    acols += [(i, j, h) for h in sector_monomials(d - 4, s)] if d - 4 >= 0 else []
        n_a = len(acols)
        self.acols = acols
        Fif = [{e: float(v) for e, v in f.items()} for f in Fi]
        LE, LC, LR, LV = [], [], [], []
        # part 1: e * (F_j d_i h - F_i d_j h)  for every e in Dmons
        for cidx, (i, j, h) in enumerate(acols):
            hp = {h: 1.0}
            G = sub(mul(Fif[j], pderiv(hp, i)), mul(Fif[i], pderiv(hp, j)))
            if not G:
                continue
            Eg = np.array(list(G.keys()), dtype=np.int64); cg = np.array(list(G.values()))
            E = Ed[:, None, :] + Eg[None, :, :]
            rows = rows_of(E.reshape(-1, 6))
            LE.append(np.repeat(np.arange(nD), len(cg))); LC.append(np.full(nD * len(cg), cidx))
            LR.append(rows); LV.append(np.tile(cg, nD))
        # part 2: - h * (F_j d_i e - F_i d_j e) for every e; group columns by (i,j)
        by_ij = {}
        for cidx, (i, j, h) in enumerate(acols):
            by_ij.setdefault((i, j), []).append((cidx, h))
        for (i, j), lst in by_ij.items():
            Hh = np.array([h for _, h in lst], dtype=np.int64); Hc = np.array([c for c, _ in lst])
            for ie, e in enumerate(Dmons):
                ep = {e: 1.0}
                M = sub(mul(Fif[j], pderiv(ep, i)), mul(Fif[i], pderiv(ep, j)))
                if not M:
                    continue
                Em_ = np.array(list(M.keys()), dtype=np.int64); cm = np.array(list(M.values()))
                E = Em_[:, None, :] + Hh[None, :, :]
                rows = rows_of(E.reshape(-1, 6))
                LE.append(np.full(len(cm) * len(Hc), ie)); LC.append(np.tile(Hc, len(cm)))
                LR.append(rows); LV.append(-np.repeat(cm, len(Hc)))
        self.LE, self.LC, self.LR, self.LV = (np.concatenate(LE), np.concatenate(LC), np.concatenate(LR), np.concatenate(LV))
        self.Lflat = self.LR * n_a + self.LC
        # --- F-multiples projector (Q unknowns)
        Qsecs = sorted({sadd((0, 1), sadd(s, s2)) for s in dsecs for s2 in dsecs})
        Fcols = []
        if 2 * d - 6 >= 0:
            for s in Qsecs:
                Fcols += [shift({e: float(v) for e, v in F.items()}, h, 1) for h in sector_monomials(2 * d - 6, s)]
        Fm = np.zeros((R, len(Fcols)))
        for jq, col in enumerate(Fcols):
            for e, v in col.items():
                Fm[rows_of(np.array([e], dtype=np.int64))[0], jq] = v
        self.nQ = len(Fcols)
        if reduce == 'sparse':
            self.QF = None
            self.C = self._sparse_reduction(F, outmons, rows_of, n_a, extra_rows, seed)   # (k x R) sparse: compression o NF mod F
            self.CT = self.C.T.tocsr()
        else:
            self.QF = np.linalg.qr(Fm, mode='reduced')[0] if len(Fcols) else np.zeros((R, 0))
            self.C = None
        # --- spurious complement
        self.Dmat, self.Nmat, self.nD, self.nN, self.n_a, self.R = Dmat, Nmat, nD, nN, n_a, R
        self.P, self.dimZ = self._complement(F)
        self.dimP = self.P.shape[1]
        self._nfev = 0
        self.build_time = time.time() - t0
        if verbose:
            print(f'[kvp] d={d} secs={[sec_str(s) for s in dsecs]} K={K} rows={R} n_a={n_a} nQ={self.nQ} nN={nN} '
                  f'quad entries={len(self.QC)} L entries={len(self.LV)} spurious={self.dimZ} search dim={self.dimP} '
                  f'reduce={reduce}' + (f' k={self.k_rows} nnz(C)={self.C.nnz} nnz(NF)={self.nnz_NF}' if self.C is not None else '') +
                  f' built in {self.build_time:.1f}s', flush=True)

    def _sparse_reduction(self, F, outmons, rows_of, n_a, extra_rows, seed):
        """C = G o (quotient map S_2d -> S_2d/(F)) as a (k x R) matrix, k = n_a + extra_rows.
        Pivot monomials X0^6 h (h in S_{2d-6}, same sectors) form a unit-triangular block P of the F-multiplication
        matrix (ordering by X0-degree); the quotient coordinates of v are v_std - N_std P^{-1} v_piv.  G is a sparse
        random compression (each standard coordinate feeds 4 compressed rows).  Zeros of the compressed residual
        coincide with zeros of the exact one with probability 1 (false zeros are caught by the exact confirmation)."""
        import scipy.sparse as sp
        from scipy.sparse.linalg import spsolve_triangular
        R = len(outmons)
        Ff = {e: float(v) for e, v in F.items()}
        hs = []
        dsecs = self.dsecs
        for s in sorted({sadd((0, 1), sadd(a, b)) for a in dsecs for b in dsecs}):
            hs += sector_monomials(2 * self.d - 6, s) if 2 * self.d - 6 >= 0 else []
        hs.sort(key=lambda h: h[0])                      # ascending X0-degree
        nQ = len(hs)
        midx = {m: i for i, m in enumerate(outmons)}
        piv_rows = [midx[(h[0] + 6,) + h[1:]] for h in hs]
        piv_set = set(piv_rows)
        std_rows = [i for i in range(R) if i not in piv_set]
        std_pos = {r: j for j, r in enumerate(std_rows)}
        piv_pos = {r: j for j, r in enumerate(piv_rows)}
        pr, pc, pv, sr, sc, sv = [], [], [], [], [], []
        for j, h in enumerate(hs):
            for e, v in shift(Ff, h, 1).items():
                r = midx[e]
                if r in piv_pos:
                    pr.append(piv_pos[r]); pc.append(j); pv.append(v)
                else:
                    sr.append(std_pos[r]); sc.append(j); sv.append(v)
        P = sp.csc_matrix((pv, (pr, pc)), shape=(nQ, nQ))          # unit upper triangular (row h' < col h)
        Nstd = sp.csc_matrix((sv, (sr, sc)), shape=(len(std_rows), nQ))
        Rp = len(std_rows)
        k = n_a + extra_rows
        rng = np.random.default_rng(seed)
        if k >= Rp:
            G = sp.identity(Rp, format='csr'); k = Rp
        else:
            per = 4
            gr = rng.integers(0, k, size=Rp * per)
            gc = np.repeat(np.arange(Rp), per)
            gv = rng.standard_normal(Rp * per)
            G = sp.csr_matrix((gv, (gr, gc)), shape=(k, Rp))
        self.k_rows = k
        # C_piv = -(G Nstd) P^{-1}  (dense k x nQ):  solve  X P = G Nstd  <=>  P^T X^T = (G Nstd)^T
        GN = (G @ Nstd).toarray()                                   # k x nQ
        XT = spsolve_triangular(P.T.tocsr(), GN.T, lower=True, unit_diagonal=True)   # nQ x k
        C_piv = -XT.T
        # assemble C as sparse (k x R): standard columns from G, pivot columns from C_piv
        Gcoo = G.tocoo()
        rows = [Gcoo.row]; cols = [np.array(std_rows)[Gcoo.col]]; vals = [Gcoo.data]
        rr, cc = np.nonzero(np.abs(C_piv) > 1e-14)
        rows.append(rr); cols.append(np.array(piv_rows)[cc]); vals.append(C_piv[rr, cc])
        C = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(k, R))
        self.nnz_NF = Nstd.nnz + P.nnz
        return C

    def _complement(self, F):
        K = self.K
        cols = [p[0] for p in self.pairs]
        Fcols = []
        if self.d >= 6:
            for s in self.dsecs:
                Fcols += [shift(F, h, 1) for h in sector_monomials(self.d - 6, s)]
        mons = {}
        for Cc in cols + Fcols:
            for e in Cc:
                mons.setdefault(e, len(mons))
        Mb = np.zeros((len(mons), K + len(Fcols)))
        for k, Cc in enumerate(cols):
            for e, v in Cc.items():
                Mb[mons[e], k] = float(v)
        for j, Cc in enumerate(Fcols):
            for e, v in Cc.items():
                Mb[mons[e], K + j] = -float(v)
        ns = sla.null_space(Mb, rcond=1e-10)
        S = ns[:K, :] if ns.size else np.zeros((K, 0))
        if S.shape[1] == 0:
            return np.eye(K), 0
        U, s, _ = np.linalg.svd(S, full_matrices=True)
        r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
        return np.ascontiguousarray(U[:, r:]), r

    # ------------------------------------------------------------------
    def assemble(self, c):
        dvec = self.Dmat @ c
        nvec = self.Nmat @ c
        quad = np.bincount(self.QR, weights=self.QC * dvec[self.QE] * nvec[self.QN], minlength=self.R)
        L = np.bincount(self.Lflat, weights=self.LV * dvec[self.LE], minlength=self.R * self.n_a).reshape(self.R, self.n_a)
        return dvec, nvec, quad, L

    def project(self, X):
        if self.C is not None:
            return self.C @ X
        if self.nQ == 0:
            return X
        return X - self.QF @ (self.QF.T @ X)

    def value_full(self, c):
        """Returns f, a, r (pulled back to the output monomial space), dvec, nvec, scale for normalised c."""
        c = np.asarray(c, float); c = c / np.linalg.norm(c)
        dvec, nvec, quad, L = self.assemble(c)
        Lp = self.project(L); qp = self.project(quad)
        a, *_ = sla.lstsq(Lp, qp, cond=self.cond, lapack_driver='gelsy', check_finite=False)
        r = Lp @ a - qp
        if self.C is not None:
            r = self.CT @ r          # adjoint pull-back for the gradient passes (f itself is |r_compressed|^2)
            f = float((Lp @ a - qp) @ (Lp @ a - qp))
        else:
            f = float(r @ r)
        self._nfev += 1
        # scale for the relative residual: |D|^2 |N0|^2, the natural magnitude of quad = B(D, N0) (c normalised)
        return f, a, r, dvec, nvec, float(dvec @ dvec) * float(nvec @ nvec)

    def value_and_grad(self, u):
        c_full = self.P @ u
        na = np.linalg.norm(c_full)
        c = c_full / na
        f, a, r, dvec, nvec, _ = self.value_full(c)
        y = np.bincount(self.LE, weights=r[self.LR] * self.LV * a[self.LC], minlength=self.nD)
        uu = np.bincount(self.QE, weights=r[self.QR] * self.QC * nvec[self.QN], minlength=self.nD)
        w = np.bincount(self.QN, weights=r[self.QR] * self.QC * dvec[self.QE], minlength=self.nN)
        g = 2.0 * (self.Dmat.T @ (y - uu) - self.Nmat.T @ w)
        g_full = (g - float(g @ c) * c) / na
        return f, self.P.T @ g_full

    # ------------------------------------------------------------------
    # joint Levenberg-Marquardt on x = (c, a): residual r = C (L(c) a - quad(c)), plus the sphere residual |c|^2 - 1
    def _joint_residual(self, x):
        K = self.K
        c = x[:K]; a = x[K:]
        dvec = self.Dmat @ c; nvec = self.Nmat @ c
        quad = np.bincount(self.QR, weights=self.QC * dvec[self.QE] * nvec[self.QN], minlength=self.R)
        L = np.bincount(self.Lflat, weights=self.LV * dvec[self.LE], minlength=self.R * self.n_a).reshape(self.R, self.n_a)
        v = L @ a - quad
        r = self.project(v)
        return np.concatenate([r, [float(c @ c) - 1.0]]), (dvec, nvec, L)

    def _joint_jacobian(self, x, cache):
        K = self.K
        c = x[:K]; a = x[K:]
        dvec, nvec, L = cache
        R, nD, nN = self.R, self.nD, self.nN
        # d r / d a = C L
        Ja = self.project(L)
        # d v / d c_j = sum_e Dmat[e,j] (L^(e) a)  -  sum_e Dmat[e,j] B(e, N(c))  -  sum_m Nmat[m,j] B(D(c), m)
        Y = np.bincount(self.LR * nD + self.LE, weights=self.LV * a[self.LC], minlength=R * nD).reshape(R, nD)
        Z1 = np.bincount(self.QR * nD + self.QE, weights=self.QC * nvec[self.QN], minlength=R * nD).reshape(R, nD)
        T = (Y - Z1) @ self.Dmat
        if R * nN <= 80_000_000:
            Z2 = np.bincount(self.QR * nN + self.QN, weights=self.QC * dvec[self.QE], minlength=R * nN).reshape(R, nN)
            T -= Z2 @ self.Nmat
        else:
            for j in range(K):
                T[:, j] -= np.bincount(self.QR, weights=self.QC * dvec[self.QE] * self.Nmat[self.QN, j], minlength=R)
        Jc = self.project(T)
        J = np.zeros((Ja.shape[0] + 1, K + self.n_a))
        J[:-1, :K] = Jc; J[:-1, K:] = Ja
        J[-1, :K] = 2.0 * c
        return J

    def search_lm(self, n_restarts, seed, maxfev=200, verbose=False, stop_rel=1e-13, progress_path=None):
        """Joint LM over (c, a) in the reduced coordinates c = P u (spurious directions removed)."""
        from scipy.optimize import least_squares
        rng = np.random.default_rng(seed)
        K = self.K
        out = []
        t0 = time.time()
        Pm = self.P
        for it in range(n_restarts):
            if STOP['flag']:
                print(f'[kvp-lm] SIGTERM/SIGINT received: stopping after {it} restarts', flush=True)
                break
            u0 = rng.standard_normal(self.dimP); u0 /= np.linalg.norm(u0)
            c0 = Pm @ u0
            # warm start a by least squares at c0
            f0, a0, r0, dv, nv, qq = self.value_full(c0)
            x0 = np.concatenate([u0, a0])
            cache = {}

            def fun(x):
                xf = np.concatenate([Pm @ x[:self.dimP], x[self.dimP:]])
                r, cc = self._joint_residual(xf)
                cache['x'] = x; cache['cc'] = cc
                return r

            def jac(x):
                xf = np.concatenate([Pm @ x[:self.dimP], x[self.dimP:]])
                if cache.get('x') is None or not np.array_equal(cache['x'], x):
                    _, cc = self._joint_residual(xf)
                else:
                    cc = cache['cc']
                J = self._joint_jacobian(xf, cc)
                # chain rule c = P u
                Ju = np.hstack([J[:, :K] @ Pm, J[:, K:]])
                return Ju
            m = self.k_rows + 1 if self.C is not None else self.R + 1
            method = 'lm' if m >= self.dimP + self.n_a else 'trf'
            res = least_squares(fun, x0, jac=jac, method=method, max_nfev=maxfev, xtol=1e-15, ftol=1e-15, gtol=1e-15)
            u = res.x[:self.dimP]; c = Pm @ u; c /= np.linalg.norm(c)
            f, a, r, dvec, nvec, qq = self.value_full(c)
            rel = f / max(qq, 1e-300)
            out.append({'value': f, 'rel': rel, 'c': c, 'iters': int(res.nfev), 'restart': it})
            if verbose and (it < 5 or it % 20 == 0):
                print(f'[kvp-lm] restart {it}: f={f:.3e} rel={rel:.3e} nfev={res.nfev} njev={res.njev} status={res.status} ({time.time()-t0:.1f}s)', flush=True)
            if progress_path:
                with open(progress_path, 'a') as pf:
                    pf.write(json.dumps({'restart': it, 'rel': rel, 'f': f, 'nfev': int(res.nfev), 'opt': 'lm', 't': round(time.time() - t0, 1)}) + '\n')
            if rel < stop_rel:
                break
        out.sort(key=lambda z: z['rel'])
        return out

    def search(self, n_restarts, seed, maxiter=400, verbose=False, stop_rel=1e-13, progress_path=None):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            if STOP['flag']:
                print(f'[kvp] SIGTERM/SIGINT received: stopping after {it} restarts', flush=True)
                break
            u0 = rng.standard_normal(self.dimP); u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': 1e-13, 'ftol': 1e-16})
            c = self.P @ res.x; c /= np.linalg.norm(c)
            f, a, r, dvec, nvec, qq = self.value_full(c)
            rel = f / max(qq, 1e-300)
            out.append({'value': f, 'rel': rel, 'c': c, 'iters': int(res.nit), 'restart': it})
            if verbose and (it < 5 or it % 20 == 0):
                print(f'[kvp] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} nfev={self._nfev} ({time.time()-t0:.1f}s)', flush=True)
            if progress_path:
                with open(progress_path, 'a') as pf:
                    pf.write(json.dumps({'restart': it, 'rel': rel, 'f': f, 'iters': int(res.nit), 't': round(time.time() - t0, 1)}) + '\n')
            if rel < stop_rel:
                break
        out.sort(key=lambda z: z['rel'])
        return out

    def D_of(self, c):
        D = {}
        for k, (Dk, Nk, tag) in enumerate(self.pairs):
            if c[k]:
                D = add(D, scale({e: float(v) for e, v in Dk.items()}, float(c[k])))
        return D


def fixed_D_check(D, N0, kvp_pairs_unused=None):
    """(unused helper placeholder)"""
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--d', type=int, required=True); ap.add_argument('--sec', required=True)
    ap.add_argument('--family', default='JK'); ap.add_argument('--restarts', type=int, default=50)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--maxiter', type=int, default=400)
    ap.add_argument('--thresh', type=float, default=1e-9); ap.add_argument('--keep', type=int, default=20)
    ap.add_argument('--tag', default='')
    ap.add_argument('--exact-metric', action='store_true', help='dense QR projection instead of the sparse NF/compression')
    ap.add_argument('--opt', default='lm', choices=['lm', 'lbfgs'], help='joint LM on (c,a) [default] or L-BFGS variable projection')
    ap.add_argument('--extra-rows', type=int, default=60)
    a = ap.parse_args(argv)
    d = a.d
    secs = sectors_of(parse_sec(a.sec))
    pairs = []
    for s in secs:
        if 'J' in a.family:
            pairs += jacobian_pairs(d, s)
        if 'K' in a.family:
            pairs += kernel_pairs(d, s)
    if not pairs:
        print('[kvp] empty family', flush=True)
        log_event({'event': 'kvp', 'd': d, 'secs': a.sec, 'family': a.family, 'K': 0, 'note': 'empty family'})
        return
    extra = a.extra_rows + (len(pairs) + 1 if a.opt == 'lm' else 0)
    fam = KoszulVarPro(d, pairs, reduce='exact' if a.exact_metric else 'sparse', extra_rows=extra)
    t0 = time.time()
    prog = os.path.join(HERE, f'{rcr2.CACHE_PREFIX}kvp_progress_d{d}_{a.sec}_{a.family}_s{a.seed}.jsonl')
    if a.opt == 'lm':
        res = fam.search_lm(a.restarts, a.seed, maxfev=a.maxiter, verbose=True, progress_path=prog)
    else:
        res = fam.search(a.restarts, a.seed, maxiter=a.maxiter, verbose=True, progress_path=prog)
    if not res:
        log_event({'event': 'kvp', 'd': d, 'secs': a.sec, 'family': a.family, 'K': fam.K, 'restarts': 0, 'note': 'stopped before the first restart finished'})
        return
    rels = np.array([z['rel'] for z in res])
    hist = {f'<1e-{k}': int(np.sum(rels < 10.0 ** (-k))) for k in (2, 4, 6, 8, 10, 12)}
    rec = {'event': 'kvp', 'd': d, 'secs': a.sec, 'family': a.family, 'K': fam.K, 'dimZ': fam.dimZ, 'dimP': fam.dimP,
           'rows': fam.R, 'n_a': fam.n_a, 'nQ': fam.nQ, 'restarts': len(res), 'seed': a.seed, 'best_rel': float(rels.min()),
           'median_rel': float(np.median(rels)), 'hist': hist, 'nfev': fam._nfev, 'seconds': round(time.time() - t0, 1),
           'build_seconds': round(fam.build_time, 1), 'tag': a.tag, 'reduce': fam.reduce, 'k_rows': getattr(fam, 'k_rows', None), 'opt': a.opt}
    print('[kvp] summary', rec, flush=True)
    log_event(rec)
    json.dump([{'rel': z['rel'], 'c': [float(x) for x in z['c']]} for z in res[:a.keep]],
              open(os.path.join(HERE, f'{rcr2.CACHE_PREFIX}kvp_d{d}_{a.sec}_{a.family}_s{a.seed}.json'), 'w'))
    for z in res:
        if z['rel'] > a.thresh:
            break
        Df = fam.D_of(z['c'])
        cands = reconstruct_poly(Df)
        log_event({'event': 'kvp-candidate', 'd': d, 'secs': a.sec, 'family': a.family, 'rel': z['rel'],
                   'n_reconstructions': len(cands), 'top_terms': sorted(((abs(v), core.to_str({e: 1})) for e, v in Df.items()), reverse=True)[:8]})
        for D in cands:
            exact_confirm(D, note=f'kvp d={d} secs={a.sec} fam={a.family} rel={z["rel"]:.2e}')


if __name__ == '__main__':
    main()
