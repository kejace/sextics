"""fplib -- wave-2 search library for the four-plus-two-residual lane (adapted from off-paired/w2lib.py).

Sector-restricted, matrix-free variable projection (VarPro) over numerically computed bases of the
E1-admissible denominator space (J(F):A)_{d,chi} (= non-Jacobian kernel + J(F)_{d,chi}), L-BFGS restarts,
Levenberg-Marquardt polishing, rational reconstruction + exact (CRT-lifted, sector-reduced) re-solve.

Characters: tuple of (m . g mod n) over the generators g of symmetry.diagonal_group(F) (n = 6 here, 3 generators,
24 characters).  Engine: /home/kejace/devel/sextics-work/hodge (read-only).
"""
import sys, os, time, json, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np, scipy.linalg as sla, scipy.sparse as sp
from scipy.optimize import minimize
from fractions import Fraction
from hodge import core, fixed_d, jacobian, symmetry
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from flint import nmod_mat

LANE = '/tmp/hodge/results/four-plus-two-residual'
WLOG = LANE + '/wave2-search.log.jsonl'
P0 = core.PRIMES[0]
BK = 64
POW = np.array([BK ** i for i in range(6)], dtype=np.int64)

_INST = {}


def load(iid='four-plus-two-residual'):
    if iid not in _INST:
        F, A, q, rec = core.load_instance(iid)
        G = symmetry.diagonal_group(F)
        _INST[iid] = (F, A, q, G)
    return _INST[iid]


class Chars:
    """Character bookkeeping for one instance (generators of the projective diagonal group)."""

    def __init__(self, G, elements=None):
        self.n = G.order
        if elements is None:
            self.gens = [tuple(g) for g in G.generators]
            self.size = len(G.elements)
        else:
            self.gens = generating_set([tuple(e) for e in elements], G.order)
            self.size = len(elements)
        self.VARCHAR = [self.char(tuple(1 if j == i else 0 for j in range(6))) for i in range(6)]
        self._mc = {}

    def char(self, e):
        return tuple(sum(a * b for a, b in zip(e, g)) % self.n for g in self.gens)

    def char_of(self, P):
        cs = {self.char(e) for e in P}
        return cs.pop() if len(cs) == 1 else None

    def cadd(self, c1, c2):
        return tuple((a + b) % self.n for a, b in zip(c1, c2))

    def cneg(self, c):
        return tuple((-a) % self.n for a in c)

    def cscale(self, c, k):
        return tuple((a * k) % self.n for a in c)

    def monos_char(self, deg, c):
        k = (deg, c)
        if k not in self._mc:
            self._mc[k] = [m for m in monomials(deg) if self.char(m) == c]
        return self._mc[k]

    def all_chars(self, deg):
        out = {}
        for m in monomials(deg):
            out.setdefault(self.char(m), 0)
            out[self.char(m)] += 1
        return out


def generating_set(elements, n):
    """Greedy small generating set of a finite abelian group of weight vectors mod n."""
    gens = []
    span = {tuple([0] * 6)}
    for w in elements:
        w = tuple(x % n for x in w)
        if w in span:
            continue
        gens.append(w)
        new = set(span)
        for s in span:
            cur = s
            for k in range(1, n):
                cur = tuple((a + b) % n for a, b in zip(cur, w))
                new.add(cur)
        span = new
    return gens


def is_irreducible_sector(ch, d, c):
    """True iff monomials of degree d and character c are even in X2 and in Y2 (otherwise X2|D or Y2|D)."""
    ms = ch.monos_char(d, c)
    return bool(ms) and all(m[2] % 2 == 0 and m[5] % 2 == 0 for m in ms)


def key(e):
    return int(sum(int(a) * int(b) for a, b in zip(e, POW)))


def keys_of(mons):
    return (np.array(mons, dtype=np.int64) @ POW) if len(mons) else np.zeros(0, np.int64)


def poly_arrays(P):
    ks = np.array([key(e) for e in P], dtype=np.int64)
    vs = np.array([float(v) for v in P.values()], dtype=float)
    return ks, vs


def poly_to_vec(P, mons_index, n):
    v = np.zeros(n)
    for e, c in P.items():
        v[mons_index[e]] += float(c)
    return v


# ----------------------------------------------------------------------------------------------------------
# E1-admissible sector spaces (numerical bases + modular dimension checks)
# ----------------------------------------------------------------------------------------------------------
def sector_dims_modp(F, A, q, d, c, ch, p=P0):
    """mod-p dimensions: dim (J:A)_{d,c} (all admissible D incl. J_d), dim J(F)_{d,c}, #monomials."""
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    cA = ch.char_of(Ap)
    Dm = ch.monos_char(d, c)
    nD = len(Dm)
    if nD == 0:
        return {'nD': 0, 'adm': 0, 'jac': 0, 'nonjac': 0}
    degN = d + 6 * q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    Jcols = []
    for i in range(NV):
        ci = ch.cadd(ch.cadd(c, ch.VARCHAR[i]), cA)
        for m in ch.monos_char(degN, ci):
            Jcols.append(shift(Fi[i], m, q, p))
    Acols = [shift(Ap, m, 1, p) for m in Dm]
    rows = {}
    for col in Acols + Jcols:
        for e in col:
            rows.setdefault(e, len(rows))
    R = len(rows)

    def mat(cols):
        C = len(cols)
        flat = [0] * (R * C)
        for j, col in enumerate(cols):
            for e, v in col.items():
                flat[rows[e] * C + j] = v
        return nmod_mat(R, C, flat, p)
    rJ = mat(Jcols).rank() if Jcols else 0
    rT = mat(Acols + Jcols).rank()
    adm = nD - rT + rJ
    # J(F)_{d,c}
    jac = 0
    if d >= 5:
        gens = []
        for i in range(NV):
            ci = ch.cadd(c, ch.VARCHAR[i])
            for m in ch.monos_char(d - 5, ci):
                gens.append(shift(Fi[i], m, 1, p))
        if gens:
            rows2 = {}
            for g in gens:
                for e in g:
                    rows2.setdefault(e, len(rows2))
            C = len(gens)
            flat = [0] * (len(rows2) * C)
            for j, g in enumerate(gens):
                for e, v in g.items():
                    flat[rows2[e] * C + j] = v
            jac = nmod_mat(len(rows2), C, flat, p).rank()
    return {'nD': nD, 'adm': adm, 'jac': jac, 'nonjac': adm - jac, 'rows': R, 'rankJ': rJ, 'rankT': rT}


def sector_basis_numeric(F, A, q, d, c, ch, tol=1e-9, verbose=False):
    """Orthonormal float basis (rows, in the coordinates of ch.monos_char(d,c)) of {D : A D in J(F)} in the sector.
    Returns (B (r x nD), Dm, info)."""
    t0 = time.time()
    cA = ch.char_of(A)
    Dm = ch.monos_char(d, c)
    nD = len(Dm)
    if nD == 0:
        return np.zeros((0, 0)), Dm, {'nD': 0, 'adm': 0}
    degN = d + 6 * q - 5
    Fi = [pderiv(F, i) for i in range(NV)]
    Jcols = []
    for i in range(NV):
        ci = ch.cadd(ch.cadd(c, ch.VARCHAR[i]), cA)
        for m in ch.monos_char(degN, ci):
            Jcols.append(shift(Fi[i], m, q))
    Acols = [shift(A, m, 1) for m in Dm]
    rows = {}
    for col in Acols + Jcols:
        for e in col:
            rows.setdefault(e, len(rows))
    R = len(rows)
    MJ = np.zeros((R, len(Jcols)))
    for j, col in enumerate(Jcols):
        for e, v in col.items():
            MJ[rows[e], j] = float(v)
    MA = np.zeros((R, nD))
    for j, col in enumerate(Acols):
        for e, v in col.items():
            MA[rows[e], j] = float(v)
    # orthonormal basis of col(MJ)
    if MJ.shape[1]:
        # scale columns for conditioning
        cn = np.linalg.norm(MJ, axis=0); cn[cn == 0] = 1
        U, s, Vt = np.linalg.svd(MJ / cn, full_matrices=False)
        rJ = int(np.sum(s > tol * s[0]))
        QJ = U[:, :rJ]
        T = MA - QJ @ (QJ.T @ MA)
    else:
        rJ = 0
        T = MA
    cn = np.linalg.norm(MA, axis=0); cn[cn == 0] = 1
    T = T / cn
    Ut, st, Vtt = np.linalg.svd(T, full_matrices=True)
    rT = int(np.sum(st > max(tol * st[0], 1e-8))) if st.size else 0
    Bt = Vtt[rT:, :]            # nullspace of T in (scaled) D coords
    B = Bt / cn[None, :]        # unscale: D = Bt / cn
    # re-orthonormalize
    if B.shape[0]:
        Ub, sb, Vb = np.linalg.svd(B, full_matrices=False)
        B = Vb
    gap = (st[rT - 1] / max(st[rT], 1e-300) if 0 < rT < len(st) else None)
    info = {'nD': nD, 'adm': int(B.shape[0]), 'rows': R, 'rankJ': rJ, 'rankT': rT,
            'sv_gap': (float(gap) if gap else None), 'smin_kept': float(st[rT - 1]) if rT else None,
            'smax_dropped': float(st[rT]) if rT < len(st) else None, 'time': round(time.time() - t0, 2)}
    if verbose:
        print(f'[sector basis] d={d} c={c} nD={nD} adm={B.shape[0]} rankJ={rJ} rankT={rT} ({info["time"]}s)', file=sys.stderr)
    return B, Dm, info


def jacobian_vecs(F, d, c, ch, Dm):
    """Float vectors (in Dm coordinates) of the generators m*F_i of J(F)_{d,c}."""
    idx = {m: i for i, m in enumerate(Dm)}
    out = []
    if d < 5:
        return np.zeros((0, len(Dm)))
    Fi = [pderiv(F, i) for i in range(NV)]
    for i in range(NV):
        ci = ch.cadd(c, ch.VARCHAR[i])
        for m in ch.monos_char(d - 5, ci):
            out.append(poly_to_vec(shift(Fi[i], m, 1), idx, len(Dm)))
    return np.array(out) if out else np.zeros((0, len(Dm)))


def gauge_vecs(F, d, c, ch, Dm):
    idx = {m: i for i, m in enumerate(Dm)}
    if d < 6:
        return np.zeros((0, len(Dm)))
    out = [poly_to_vec(shift(F, m, 1), idx, len(Dm)) for m in ch.monos_char(d - 6, c)]
    return np.array(out) if out else np.zeros((0, len(Dm)))


def orth_rows(V, tol=1e-10):
    if V.shape[0] == 0:
        return V
    U, s, Vt = np.linalg.svd(V, full_matrices=False)
    r = int(np.sum(s > tol * s[0]))
    return Vt[:r]


def restrict_rows(B, sub, tol=1e-8):
    """Rows of sub (in Dm coords) projected into span(B): returns orthonormal basis of span(B) ∩ span(sub)
    (numerically: the part of span(sub) lying inside span(B))."""
    if sub.shape[0] == 0 or B.shape[0] == 0:
        return np.zeros((0, B.shape[1]))
    Proj = sub - (sub @ B.T) @ B
    Uc, sc, Vct = np.linalg.svd(Proj, full_matrices=True)
    rc = int(np.sum(sc > tol * max(np.linalg.norm(sub, axis=1).max(), 1e-300))) if sc.size else 0
    C = Uc[:, rc:].T
    out = C @ sub
    out = out[np.linalg.norm(out, axis=1) > 1e-12]
    return orth_rows(out)


# ----------------------------------------------------------------------------------------------------------
# sector-restricted matrix-free VarPro
# ----------------------------------------------------------------------------------------------------------
class SectorVarPro:
    """D(a) = sum_k a_k B_k, B rows in the coordinates of Dm = monos_char(d, c); f(a) = min_{N,Q} |M(a)x - b(a)|^2,
    b = -A D^2, a on the unit sphere of the orthogonal complement of the gauge (F * R_{d-6}) inside span(B)."""

    def __init__(self, F, A, q, B, Dm, c, ch, gauge=None, verbose=False, cond=1e-9, tag=''):
        t0 = time.time()
        self.F, self.A, self.q, self.c, self.ch = F, A, q, c, ch
        self.tag = tag
        self.Dm = Dm; nD = len(Dm); self.nD = nD
        d = core.degree({Dm[0]: 1}); self.d = d
        cA = ch.char_of(A); self.cA = cA
        self.cond = cond
        B = orth_rows(np.asarray(B, float))
        self.Borth = B; r = B.shape[0]; self.rankB = r
        G = gauge if gauge is not None else gauge_vecs(F, d, c, ch, Dm)
        ng = 0
        if G.shape[0]:
            Gc = G @ B.T
            resid = np.linalg.norm(G - Gc @ B) / max(np.linalg.norm(G), 1e-300)
            if resid > 1e-6:
                # gauge not fully inside the family: keep only the inside part
                Gin = restrict_rows(B, G)
                Gc = Gin @ B.T
            if Gc.shape[0]:
                Ug, sg, Vgt = np.linalg.svd(Gc, full_matrices=True)
                ng = int(np.sum(sg > 1e-10 * max(sg[0], 1e-300)))
                self.P = Vgt[ng:, :].T
            else:
                self.P = np.eye(r)
        else:
            self.P = np.eye(r)
        self.ng = ng
        self.dimP = self.P.shape[1]
        degN = d + 6 * q - 5; degQ = 2 * d - 6
        self.degN, self.degQ = degN, degQ
        crow = ch.cadd(cA, ch.cadd(c, c))
        rowmons = ch.monos_char(2 * d + 6 * q, crow)
        rowkeys = keys_of(rowmons)
        order = np.argsort(rowkeys)
        self._rk_sorted = rowkeys[order]; self._rk_order = order
        R = len(rowmons); self.R = R

        def row_index(ks):
            pos = np.searchsorted(self._rk_sorted, ks)
            pos = np.minimum(pos, len(self._rk_sorted) - 1)
            if not np.all(self._rk_sorted[pos] == ks):
                raise RuntimeError('monomial outside the row sector')
            return self._rk_order[pos]
        self.row_index = row_index
        Ncols = []; self.Ncols_by_i = []
        for i in range(NV):
            ci = ch.cadd(ch.cadd(c, ch.VARCHAR[i]), cA)
            ms = ch.monos_char(degN, ci)
            self.Ncols_by_i.append((len(Ncols), ms))
            Ncols += [(i, m) for m in ms]
        nN = len(Ncols); self.nN = nN; self.Ncols = Ncols
        Fi = [pderiv(F, i) for i in range(NV)]
        ri_l, ci_l, va_l = [], [], []
        for mu_idx, mu in enumerate(Dm):
            mup = {mu: Fraction(1)}
            Fmu = mul(F, mup)
            fk, fv = poly_arrays(Fmu)
            for i in range(NV):
                start, ms = self.Ncols_by_i[i]
                if not ms:
                    continue
                Gi = scale(add(scale(mul(mup, Fi[i]), q), mul(F, pderiv(mup, i))), -1)
                mk = keys_of(ms)
                cols = np.arange(start, start + len(ms), dtype=np.int64)
                if Gi:
                    gk, gv = poly_arrays(Gi)
                    kk = (gk[:, None] + mk[None, :]).ravel()
                    rws = row_index(kk)
                    ri_l.append(rws * nN + np.broadcast_to(cols[None, :], (len(gk), len(ms))).ravel())
                    ci_l.append(np.full(kk.size, mu_idx, dtype=np.int64))
                    va_l.append(np.broadcast_to(gv[:, None], (len(gk), len(ms))).ravel().copy())
                mi = np.array([m[i] for m in ms], dtype=np.int64)
                sel = np.nonzero(mi)[0]
                if sel.size:
                    mk2 = mk[sel] - BK ** i
                    kk = (fk[:, None] + mk2[None, :]).ravel()
                    rws = row_index(kk)
                    ri_l.append(rws * nN + np.broadcast_to(cols[sel][None, :], (len(fk), sel.size)).ravel())
                    ci_l.append(np.full(kk.size, mu_idx, dtype=np.int64))
                    va_l.append((fv[:, None] * mi[sel][None, :].astype(float)).ravel())
        ri = np.concatenate(ri_l); ci = np.concatenate(ci_l); va = np.concatenate(va_l)
        self.TN = sp.csr_matrix((va, (ri, ci)), shape=(R * nN, nD))
        self.TNT = self.TN.T.tocsr()
        Qm = ch.monos_char(degQ, crow) if degQ >= 0 else []
        self.Qm = Qm; nQ = len(Qm); self.nQ = nQ
        if nQ:
            F2n = scale(power(F, q + 1), -1)
            f2k, f2v = poly_arrays(F2n)
            MQ = np.zeros((R, nQ))
            qk = keys_of(Qm)
            kk = (f2k[:, None] + qk[None, :]).ravel()
            rws = row_index(kk)
            cols = np.broadcast_to(np.arange(nQ)[None, :], (len(f2k), nQ)).ravel()
            np.add.at(MQ, (rws, cols), np.broadcast_to(f2v[:, None], (len(f2k), nQ)).ravel())
            self.MQ = MQ
            Qq, Rq = np.linalg.qr(MQ, mode='reduced')
            self.QQ = Qq
        else:
            self.MQ = np.zeros((R, 0)); self.QQ = np.zeros((R, 0))
        assert len(A) == 1, 'A must be a monomial'
        kA = key(next(iter(A)))
        self.cAval = float(A[next(iter(A))])
        dk = keys_of(Dm)
        mu_i, nu_i = np.meshgrid(np.arange(nD), np.arange(nD), indexing='ij')
        self.pair_mu, self.pair_nu = mu_i.ravel(), nu_i.ravel()
        self.pair_row = row_index(kA + dk[self.pair_mu] + dk[self.pair_nu])
        self.ncols = nN + nQ
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[sector varpro {tag}] d={d} c={c} K={r} gauge={ng} dimP={self.dimP} rows={R} '
                  f'N={nN} Q={nQ} nnz(TN)={self.TN.nnz} built in {self.build_time:.1f}s', file=sys.stderr)

    def D_coef(self, u):
        a = self.P @ np.asarray(u, float)
        na = np.linalg.norm(a)
        a = a / na
        return self.Borth.T @ a, a, na

    def assemble(self, Dc):
        MN = (self.TN @ Dc).reshape(self.R, self.nN)
        b = -self.cAval * np.bincount(self.pair_row, weights=Dc[self.pair_mu] * Dc[self.pair_nu], minlength=self.R)
        return MN, b

    def residual(self, Dc):
        MN, b = self.assemble(Dc)
        if self.nQ:
            MNp = MN - self.QQ @ (self.QQ.T @ MN)
            bp = b - self.QQ @ (self.QQ.T @ b)
        else:
            MNp, bp = MN, b
        xN, *_ = sla.lstsq(MNp, bp, cond=self.cond, lapack_driver='gelsy', check_finite=False)
        r = MNp @ xN - bp
        self._nfev += 1
        return r, xN, b

    def value(self, u):
        Dc, a, na = self.D_coef(u)
        r, xN, b = self.residual(Dc)
        f = float(r @ r); bn = float(b @ b)
        return f, bn, r, xN, Dc

    def value_and_grad(self, u):
        Dc, a, na = self.D_coef(u)
        r, xN, b = self.residual(Dc)
        f = float(r @ r)
        gD = 2.0 * (self.TNT @ np.outer(r, xN).ravel())
        gb = -2.0 * self.cAval * np.bincount(self.pair_mu, weights=r[self.pair_row] * Dc[self.pair_nu], minlength=self.nD)
        gD = gD - 2.0 * gb
        ga = self.Borth @ gD
        g_full = (ga - float(ga @ a) * a) / na
        return f, self.P.T @ g_full

    def relres(self, u):
        f, bn, r, xN, Dc = self.value(u)
        return math.sqrt(f / max(bn, 1e-300))

    def fd_check(self, seed=0, h=1e-6):
        rng = np.random.default_rng(seed)
        u = rng.standard_normal(self.dimP); u /= np.linalg.norm(u)
        f0, g = self.value_and_grad(u)
        v = rng.standard_normal(self.dimP); v -= (v @ u) * u; v /= np.linalg.norm(v)
        fp, _ = self.value_and_grad(u + h * v); fm, _ = self.value_and_grad(u - h * v)
        fd = (fp - fm) / (2 * h)
        return float(g @ v), fd

    def minimize_from(self, seed, maxiter=400, gtol=1e-12, init=None):
        rng = np.random.default_rng(seed)
        u0 = rng.standard_normal(self.dimP) if init is None else np.asarray(init, float)
        u0 /= np.linalg.norm(u0)
        t0 = time.time()
        res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                       options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
        u = res.x / np.linalg.norm(res.x)
        f, bn, r, xN, Dc = self.value(u)
        rel = math.sqrt(f / max(bn, 1e-300))
        return {'seed': int(seed), 'f': f, 'bnorm2': bn, 'relres': rel, 'iters': int(res.nit), 'nfev': int(res.nfev),
                'status': int(res.status), 'time': round(time.time() - t0, 2), 'u': u, 'Dc': Dc}

    # --- Levenberg-Marquardt (Kaufman Jacobian) --------------------------------------------------------------
    def _qr_solve(self, Dc):
        MN, b = self.assemble(Dc)
        if self.nQ:
            MNp = MN - self.QQ @ (self.QQ.T @ MN)
            bp = b - self.QQ @ (self.QQ.T @ b)
        else:
            MNp, bp = MN, b
        Qf, Rf, piv = sla.qr(MNp, mode='economic', pivoting=True, check_finite=False)
        dg = np.abs(np.diag(Rf))
        rank = int(np.sum(dg > self.cond * dg[0])) if dg.size and dg[0] > 0 else 0
        Q1 = Qf[:, :rank]
        y = Q1.T @ bp
        xN = np.zeros(self.nN)
        if rank:
            z = sla.solve_triangular(Rf[:rank, :rank], y, check_finite=False)
            xN[piv[:rank]] = z
        r = -(bp - Q1 @ y)
        self._nfev += 1
        return r, xN, b, Q1

    def resid(self, u):
        Dc, a, na = self.D_coef(u)
        r, xN, b, Q1 = self._qr_solve(Dc)
        return r, b

    def resid_jac(self, u):
        u = np.asarray(u, float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        Dc = self.Borth.T @ a
        r, xN, b, Q1 = self._qr_solve(Dc)
        R, nN, nD = self.R, self.nN, self.nD
        Xm = sp.csr_matrix((np.tile(xN, R), np.arange(R * nN, dtype=np.int64), np.arange(0, R * nN + 1, nN, dtype=np.int64)),
                           shape=(R, R * nN))
        Y = np.asarray((Xm @ self.TN).todense())
        Bp = sp.csr_matrix((-2.0 * self.cAval * Dc[self.pair_nu], (self.pair_row, self.pair_mu)), shape=(R, nD)).toarray()
        JD = Y - Bp
        if self.nQ:
            JD = JD - self.QQ @ (self.QQ.T @ JD)
        JD = JD - Q1 @ (Q1.T @ JD)
        Ja = JD @ self.Borth.T
        Ja = (Ja - np.outer(Ja @ a, a)) / na
        Ju = Ja @ self.P
        return r, Ju, b

    def lm_from(self, seed, maxiter=150, init=None, tol=1e-14, lam0=1e-3):
        rng = np.random.default_rng(seed)
        u = rng.standard_normal(self.dimP) if init is None else np.asarray(init, float)
        u = u / np.linalg.norm(u)
        t0 = time.time()
        r, J, b = self.resid_jac(u)
        f = float(r @ r); bn = float(b @ b)
        lam = lam0
        it = 0; njac = 1; nres = 0
        status = 'maxiter'
        while it < maxiter:
            it += 1
            g = J.T @ r; H = J.T @ J
            dH = np.diag(H).copy(); dH[dH < 1e-300] = 1e-300
            accepted = False
            for _ in range(30):
                try:
                    delta = np.linalg.solve(H + lam * np.diag(dH), -g)
                except np.linalg.LinAlgError:
                    lam *= 10; continue
                u_new = u + delta
                u_new /= np.linalg.norm(u_new)
                r_new, b_new = self.resid(u_new); nres += 1
                f_new = float(r_new @ r_new)
                if f_new < f:
                    accepted = True
                    break
                lam *= 4.0
                if lam > 1e14:
                    break
            if not accepted:
                status = 'stalled'; break
            rel_impr = (f - f_new) / max(f, 1e-300)
            u, f = u_new, f_new
            lam = max(lam / 3.0, 1e-15)
            r, J, b = self.resid_jac(u); njac += 1
            f = float(r @ r); bn = float(b @ b)
            if math.sqrt(f / max(bn, 1e-300)) < tol:
                status = 'converged'; break
            if rel_impr < 1e-12 and lam > 1e6:
                status = 'flat'; break
        Dc, a, na = self.D_coef(u)
        rel = math.sqrt(f / max(bn, 1e-300))
        return {'seed': int(seed), 'f': f, 'bnorm2': bn, 'relres': rel, 'iters': it, 'njac': njac, 'nres': nres,
                'status': status, 'time': round(time.time() - t0, 2), 'u': u, 'Dc': Dc}

    # --- reconstruction ------------------------------------------------------------------------------------
    def D_exact_from_coef(self, Dc, max_den, pivot=None):
        j = int(np.argmax(np.abs(Dc))) if pivot is None else pivot
        v = Dc / Dc[j]
        fr = [Fraction(float(x)).limit_denominator(max_den) for x in v]
        err = max(abs(float(a) - b) for a, b in zip(fr, v))
        D = {self.Dm[i]: fr[i] for i in range(self.nD) if fr[i] != 0}
        return D, err

    def D_float_str(self, Dc, ndig=6):
        j = int(np.argmax(np.abs(Dc)))
        v = Dc / Dc[j]
        terms = []
        for i in np.argsort(-np.abs(v)):
            if abs(v[i]) < 1e-12:
                continue
            terms.append(f'{v[i]:.{ndig}g}*{core.to_str({self.Dm[i]: 1})}')
        return ' + '.join(terms)


# ----------------------------------------------------------------------------------------------------------
# exact confirmation
# ----------------------------------------------------------------------------------------------------------
def confirm_exact(F, A, q, D, G=None, verbose=False, max_primes=12):
    if not D:
        return 'inconsistent', {'note': 'D=0'}
    divF, _ = jacobian.divisible_by(F, D)
    if divF:
        return 'spurious', {'note': 'F divides D'}
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    if not ok:
        return 'inconsistent', {'note': 'e1 fails', 'e1': info}
    st, data = symmetry.solve_fixed_D_symmetric(F, A, q, D, verbose=verbose, max_primes=max_primes, group=G)
    if st == 'solved':
        ok2, info2 = core.verify_witness(F, A, q, D, data['N'], data['Q'])
        data['verified'] = bool(ok2)
    return st, data


def try_reconstruct_and_confirm(fam, Dc, G=None, dens=(1, 2, 3, 6, 7, 12, 14, 42, 60, 84, 120, 420, 1000, 2520, 10000, 100000, 10 ** 6),
                                tol=1e-5, verbose=False, pivots=None):
    out = []
    seen = set()
    pivs = [None] if pivots is None else pivots
    for piv in pivs:
        for md in dens:
            try:
                D, err = fam.D_exact_from_coef(Dc, md, pivot=piv)
            except Exception as ex:
                out.append({'max_den': md, 'error': str(ex)}); continue
            if err > tol or not D:
                out.append({'max_den': md, 'recon_err': err, 'status': 'no-fit'}); continue
            s = core.to_str(D)
            if s in seen:
                continue
            seen.add(s)
            t0 = time.time()
            st, data = confirm_exact(fam.F, fam.A, fam.q, D, G=G, verbose=verbose)
            rec = {'max_den': md, 'recon_err': err, 'status': st, 'D': s, 'time': round(time.time() - t0, 1)}
            for k in ('rank', 'nullity', 'note', 'verified'):
                if k in data:
                    rec[k] = data[k]
            if st == 'solved':
                rec['N'] = [core.to_str(n) for n in data['N']]; rec['Q'] = core.to_str(data['Q'])
                out.append(rec)
                return out
            out.append(rec)
    return out


def save_witness(iid, D, N, Q, path):
    F, A, q, G = load(iid)
    Dp = core.parse(D) if isinstance(D, str) else D
    Np = [core.parse(n) if isinstance(n, str) else n for n in N]
    Qp = core.parse(Q) if isinstance(Q, str) else Q
    ok, info = core.verify_witness(F, A, q, Dp, Np, Qp)
    w = {'instance': iid, 'd': core.degree(Dp), 'D': core.to_str(Dp), 'N': [core.to_str(n) for n in Np],
         'Q': core.to_str(Qp), 'verified': bool(ok), 'info': str(info)}
    with open(path, 'w') as f:
        json.dump(w, f, indent=1)
    return ok


# ----------------------------------------------------------------------------------------------------------
# logging
# ----------------------------------------------------------------------------------------------------------
def append_json(path, rec):
    with open(path, 'a') as f:
        f.write(json.dumps(rec) + '\n')


def wlog(rec):
    rec = dict(rec); rec.setdefault('time', now()); rec.setdefault('engine', 2)
    append_json(WLOG, rec)


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%S')


def cstr(c):
    return ''.join(str(x) for x in c)
