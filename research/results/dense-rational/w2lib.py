"""w2lib -- wave-2 search library for the dense-rational lane (q = 1, A = X0^3 Y0^3).

Normalized-E2 variable projection.  Every E1 solution for a denominator D in a linear family is
    N = N0(w) + K . grad F,   D = D(w) = sum_k w_k D_k,  N0(w) = sum_k w_k N0_k   (N0_k.grad F = A D_k),
K skew (15 entries of degree d-4), and the certificate identity is equivalent (Euler normalization, see
dossier 1.2 / attack-1 4.1) to
    E2:  B_D(N) := D div N - N.grad D  =  F C,     C of degree 2d-6.
For fixed w this is LINEAR in (K, C); the reduced functional f(w) = min_{K,C} ||B_D(N0+K.gradF) - F C||^2
is minimized by Kaufman-type Levenberg-Marquardt (variable projection).  f is homogeneous of degree 4 in w.
Spurious zeros: w with D(w) in F*S_{d-6} (gauge subspace G).  We parametrize w = u + g, u in G^perp with
|u| = 1 (sphere), g in G free; then f(u+g) = 0 iff a genuine certificate exists (D(w) is never an F-multiple),
and no denominator of the family is lost (unlike restricting to G^perp).

Assembly of M_K(w) = [B_D(e_j grad F)]_j uses dense GEMMs:  M_K = S_D @ V - sum_i S_{d_i D} @ W_i  where
S_P is the "shift matrix" (rows x monomials-of-the-other-factor) of a polynomial P and V, W_i are constant.

Engine: /home/kejace/devel/sextics-work/hodge (not modified).  Exact confirmation only through the engine.
"""
import sys, os, time, json, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.linalg as sla
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power

LANE = '/tmp/hodge/results/dense-rational'
LOG = LANE + '/wave2-search.log.jsonl'
IID = 'dense-rational'
BK = 64
POW = np.array([BK ** i for i in range(6)], dtype=np.int64)
ZPAR = [0, 0, 0, 1, 1, 1]


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%S')


def append_json(path, rec):
    with open(path, 'a') as f:
        f.write(json.dumps(rec) + '\n')


def log_batch(rec):
    rec = dict(rec); rec.setdefault('ts', now()); rec.setdefault('engine', fixed_d.ENGINE_VERSION)
    append_json(LOG, rec)


def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q


def key(e):
    return int(sum(int(a) * int(b) for a, b in zip(e, POW)))


def keys_of(mons):
    return (np.array(mons, dtype=np.int64) @ POW) if len(mons) else np.zeros(0, np.int64)


def unkey(k):
    e = []
    for i in range(6):
        e.append(int(k % BK)); k //= BK
    return tuple(e)


def poly_arrays(P):
    ks = np.array([key(e) for e in P], dtype=np.int64)
    vs = np.array([float(v) for v in P.values()], dtype=float)
    return ks, vs


class MonSpace:
    """Monomials of one degree (and optional Y-parity) with fast key -> index lookup."""

    def __init__(self, deg, par=None):
        self.deg, self.par = deg, par
        self.mons = monomials(deg, par) if deg >= 0 else []
        self.keys = keys_of(self.mons)
        self.n = len(self.mons)
        order = np.argsort(self.keys)
        self._sorted = self.keys[order]; self._order = order

    def index(self, ks):
        ks = np.asarray(ks, dtype=np.int64)
        if self.n == 0:
            raise RuntimeError('empty space')
        pos = np.searchsorted(self._sorted, ks)
        pos = np.minimum(pos, self.n - 1)
        if not np.all(self._sorted[pos] == ks):
            bad = ks[self._sorted[pos] != ks][:3]
            raise RuntimeError(f'monomial outside space deg={self.deg} par={self.par}: {[unkey(int(b)) for b in bad]}')
        return self._order[pos]

    def vec(self, P):
        v = np.zeros(self.n)
        if P:
            ks, vs = poly_arrays(P)
            v[self.index(ks)] = vs
        return v

    def poly(self, v, tol=0.0):
        return {self.mons[i]: float(v[i]) for i in range(self.n) if abs(v[i]) > tol}


def deriv_map(src: MonSpace, dst: MonSpace, i):
    """Index arrays for d/dz_i : src -> dst.  Returns (src_idx, dst_idx, mult)."""
    E = np.array(src.mons, dtype=np.int64) if src.n else np.zeros((0, 6), np.int64)
    sel = np.nonzero(E[:, i] > 0)[0]
    if sel.size == 0:
        return sel, sel, np.zeros(0)
    mult = E[sel, i].astype(float)
    dk = src.keys[sel] - BK ** i
    return sel, dst.index(dk), mult


def apply_deriv(dm, vec, n_dst):
    s, t, m = dm
    if s.size == 0:
        return np.zeros(n_dst)
    return np.bincount(t, weights=vec[s] * m, minlength=n_dst)


def shift_table(S1: MonSpace, S2: MonSpace, S3: MonSpace):
    """IDX[mu, t] = S3.index(mu + t) for mu in S1, t in S2 (int64 array nS1 x nS2)."""
    if S1.n == 0 or S2.n == 0:
        return np.zeros((S1.n, S2.n), np.int64)
    kk = (S1.keys[:, None] + S2.keys[None, :]).ravel()
    return S3.index(kk).reshape(S1.n, S2.n)


# ----------------------------------------------------------------------------------------------------------
# families of E1 pairs (D_k, N0_k)
# ----------------------------------------------------------------------------------------------------------
class Family:
    """D_k given as vectors over Dspace (nw x nD); N0_k as six arrays (nw x nT_i) over Nspaces."""

    def __init__(self, d, pD, Dmat, N0mats, names, kind):
        self.d, self.pD, self.Dmat, self.N0mats, self.names, self.kind = d, pD, Dmat, N0mats, names, kind
        self.nw = Dmat.shape[0]


def numerator_parities(pD, pA=1):
    if pD is None:
        return [None] * NV
    return [(pA - pD - ZPAR[i]) % 2 for i in range(NV)]


def spaces(d, pD):
    """All monomial spaces used by the normalized system at degree d, denominator parity pD (None = mixed)."""
    pars = numerator_parities(pD)
    sp = {}
    sp['D'] = MonSpace(d, pD)
    sp['rows'] = MonSpace(2 * d, None if pD is None else 1)
    sp['T'] = MonSpace(d, None if pD is None else (1 - pD) % 2)          # div N lives here
    sp['N'] = [MonSpace(d + 1, pars[i]) for i in range(NV)]
    sp['dD'] = [MonSpace(d - 1, None if pD is None else (pD + ZPAR[i]) % 2) for i in range(NV)]
    sp['C'] = MonSpace(2 * d - 6, None if pD is None else 1)
    sp['K'] = []   # list of (i, j, MonSpace)
    for i in range(NV):
        for j in range(i + 1, NV):
            par = None if pD is None else (pars[i] - ZPAR[j]) % 2
            sp['K'].append((i, j, MonSpace(d - 4, par)))
    sp['pars'] = pars
    return sp


def jacobian_family(F, A, d, pD, sp):
    """Jacobian-type family: D_k = m F_i, N0_k = A m e_i, m in S_{d-5} with parity pD + par(z_i)."""
    Fi = [pderiv(F, i) for i in range(NV)]
    Dsp, Nsp = sp['D'], sp['N']
    rows_D, rows_N, names = [], [], []
    for i in range(NV):
        par = None if pD is None else (pD + ZPAR[i]) % 2
        for m in monomials(d - 5, par):
            Dk = shift(Fi[i], m, 1)
            rows_D.append(Dsp.vec(Dk))
            n0 = [np.zeros(Nsp[j].n) for j in range(NV)]
            n0[i] = Nsp[i].vec(shift(A, m, 1))
            rows_N.append(n0)
            names.append((i, m))
    Dmat = np.array(rows_D)
    N0mats = [np.array([r[j] for r in rows_N]) for j in range(NV)]
    return Family(d, pD, Dmat, N0mats, names, 'jacobian')


# ----------------------------------------------------------------------------------------------------------
# the normalized E2 variable-projection system
# ----------------------------------------------------------------------------------------------------------
class NormE2:
    def __init__(self, F, A, q, fam: Family, sp, cond=1e-10, verbose=False, exclude_gauge=True):
        assert q == 1
        t0 = time.time()
        self.F, self.A, self.q, self.fam, self.sp = F, A, q, fam, sp
        d, pD = fam.d, fam.pD
        self.d, self.pD = d, pD
        self.cond = cond
        rows, Dsp, Tsp, Nsp, dDsp, Csp = sp['rows'], sp['D'], sp['T'], sp['N'], sp['dD'], sp['C']
        self.R, self.nD, self.nT, self.nC = rows.n, Dsp.n, Tsp.n, Csp.n
        self.nTi = [s.n for s in Nsp]; self.ndDi = [s.n for s in dDsp]
        Fi = [pderiv(F, i) for i in range(NV)]
        # --- shift tables
        self.IDX_DT = shift_table(Dsp, Tsp, rows)                       # (nD x nT): row(mu + t)
        self.FLAT_SD = (self.IDX_DT * self.nT + np.arange(self.nT)[None, :]).ravel()      # S_D[row, t]
        self.FLAT_SDp = (self.IDX_DT * self.nD + np.arange(self.nD)[:, None]).ravel()     # S'_P[row, mu]
        self.IDX_dN = [shift_table(dDsp[i], Nsp[i], rows) for i in range(NV)]           # (ndD_i x nT_i)
        self.FLAT_Si = [(self.IDX_dN[i] * self.nTi[i] + np.arange(self.nTi[i])[None, :]).ravel() for i in range(NV)]
        self.FLAT_Sip = [(self.IDX_dN[i] * self.ndDi[i] + np.arange(self.ndDi[i])[:, None]).ravel() for i in range(NV)]
        # --- derivative maps
        self.dD_maps = [deriv_map(Dsp, dDsp[i], i) for i in range(NV)]                  # D -> d_i D
        self.divN_maps = [deriv_map(Nsp[i], Tsp, i) for i in range(NV)]                  # N_i -> d_i N_i in T
        # --- K basis: V (nT x nK) = div(e_j grad F), W_i (nT_i x nK) = (e_j grad F)_i
        Kcols = []
        for (i, j, ksp) in sp['K']:
            for m in ksp.mons:
                Kcols.append((i, j, m))
        self.Kcols = Kcols; nK = len(Kcols); self.nK = nK
        V = np.zeros((self.nT, nK))
        W = [np.zeros((self.nTi[i], nK)) for i in range(NV)]
        for c, (i, j, m) in enumerate(Kcols):
            mm = {m: Fraction(1)}
            Ni = mul(mm, Fi[j]); Nj = scale(mul(mm, Fi[i]), -1)        # K_ij = m, K_ji = -m
            W[i][:, c] = Nsp[i].vec(Ni); W[j][:, c] = Nsp[j].vec(Nj)
            dv = add(pderiv(Ni, i), pderiv(Nj, j))
            V[:, c] = Tsp.vec(dv)
        self.V, self.W = V, W
        # --- C block and its orthonormal range
        MC = np.zeros((self.R, self.nC))
        fk, fv = poly_arrays(F)
        ck = Csp.keys
        kk = (fk[:, None] + ck[None, :]).ravel()
        rr = rows.index(kk)
        MC[rr, np.broadcast_to(np.arange(self.nC)[None, :], (len(fk), self.nC)).ravel()] = \
            np.broadcast_to(fv[:, None], (len(fk), self.nC)).ravel()
        self.MC = MC
        self.QQ, _ = np.linalg.qr(MC)   # R x nC
        # --- family data
        self.Dmat = fam.Dmat                                             # nw x nD
        self.nw = fam.nw
        self.N0mats = fam.N0mats                                         # nw x nT_i
        self.dDmat = [np.array([apply_deriv(self.dD_maps[i], self.Dmat[k], self.ndDi[i]) for k in range(self.nw)])
                      for i in range(NV)]                                # nw x ndD_i
        self.divN0mat = np.zeros((self.nw, self.nT))
        for k in range(self.nw):
            for i in range(NV):
                self.divN0mat[k] += apply_deriv(self.divN_maps[i], self.N0mats[i][k], self.nT)
        # --- gauge: w with D(w) in F * S_{d-6}^{pD}
        self.ng = 0; self.G = np.zeros((self.nw, 0)); self.Pp = np.eye(self.nw)
        if exclude_gauge and d >= 6:
            Gsp = MonSpace(d - 6, pD)
            Fcols = np.array([Dsp.vec(shift(F, m, 1)) for m in Gsp.mons]).T if Gsp.n else np.zeros((self.nD, 0))
            Mb = np.hstack([self.Dmat.T, -Fcols])
            ns = sla.null_space(Mb, rcond=1e-11)
            Wg = ns[:self.nw, :]
            if Wg.shape[1]:
                U, s, _ = np.linalg.svd(Wg, full_matrices=True)
                ng = int(np.sum(s > 1e-9 * s[0]))
                self.G = U[:, :ng]; self.Pp = U[:, ng:]; self.ng = ng
        self.dimU = self.Pp.shape[1]
        # --- buffers
        self.SD = np.zeros((self.R, self.nT))
        self.SDp = np.zeros((self.R, self.nD))
        self.Si = [np.zeros((self.R, self.nTi[i])) for i in range(NV)]
        self.Sip = [np.zeros((self.R, self.ndDi[i])) for i in range(NV)]
        self.build_time = time.time() - t0
        self.nfev = 0
        if verbose:
            print(f'[NormE2] d={d} pD={pD} kind={fam.kind} nw={self.nw} gauge={self.ng} dimU={self.dimU} '
                  f'rows={self.R} nK={nK} nC={self.nC} nT={self.nT} built in {self.build_time:.1f}s', file=sys.stderr)

    # --- w <-> (u, g)
    def w_of(self, u, g):
        u = np.asarray(u, float); nu = np.linalg.norm(u)
        w = self.Pp @ (u / nu)
        if self.ng:
            w = w + self.G @ np.asarray(g, float)
        return w

    # --- shift matrices
    def _fill(self, buf, FLAT, coeffs_rowwise, ncols):
        buf.fill(0.0)
        buf.ravel()[FLAT] = np.repeat(coeffs_rowwise, ncols)

    def _fillT(self, buf, FLAT, coeffs_colwise, nrows):
        buf.fill(0.0)
        buf.ravel()[FLAT] = np.tile(coeffs_colwise, nrows)

    # --- evaluation
    def evaluate(self, w, need_jac=False):
        """Inner solve at w.  Returns state dict with f, bnorm2, r, x, (Jw if need_jac)."""
        w = np.asarray(w, float)
        Dc = self.Dmat.T @ w                                             # nD
        N0 = [self.N0mats[i].T @ w for i in range(NV)]                   # nT_i
        dD = [apply_deriv(self.dD_maps[i], Dc, self.ndDi[i]) for i in range(NV)]
        divN0 = self.divN0mat.T @ w                                      # nT
        self._fill(self.SD, self.FLAT_SD, Dc, self.nT)
        for i in range(NV):
            self._fill(self.Si[i], self.FLAT_Si[i], dD[i], self.nTi[i])
        b = self.SD @ divN0
        MK = self.SD @ self.V
        for i in range(NV):
            b -= self.Si[i] @ N0[i]
            MK -= self.Si[i] @ self.W[i]
        bnorm2 = float(b @ b)
        # project out the C block
        MKp = MK - self.QQ @ (self.QQ.T @ MK)
        bp = b - self.QQ @ (self.QQ.T @ b)
        Qf, Rf, piv = sla.qr(MKp, mode='economic', pivoting=True, check_finite=False)
        dg = np.abs(np.diag(Rf))
        rank = int(np.sum(dg > self.cond * dg[0])) if dg.size and dg[0] > 0 else 0
        Q1 = Qf[:, :rank]
        y = Q1.T @ bp
        xp = np.zeros(self.nK)
        if rank:
            xp[:rank] = sla.solve_triangular(Rf[:rank, :rank], y, check_finite=False)
        x = np.zeros(self.nK); x[piv] = xp
        r = (Q1 @ y) - bp
        f = float(r @ r)
        self.nfev += 1
        st = {'f': f, 'bnorm2': bnorm2, 'r': r, 'x': x, 'rank': rank, 'w': w, 'Dc': Dc, 'Q1': Q1}
        if need_jac:
            # N = N0 + K grad F ; divN
            Ni = [N0[i] + self.W[i] @ x for i in range(NV)]
            divN = divN0 + self.V @ x
            self._fillT(self.SDp, self.FLAT_SDp, divN, self.nD)
            Jw = self.SDp @ self.Dmat.T                                   # [D_k divN]_k
            Jw += self.SD @ self.divN0mat.T                               # [D div N0_k]_k
            for i in range(NV):
                self._fillT(self.Sip[i], self.FLAT_Sip[i], Ni[i], self.ndDi[i])
                Jw -= self.Sip[i] @ self.dDmat[i].T                       # [N_i d_i D_k]_k
                Jw -= self.Si[i] @ self.N0mats[i].T                       # [N0_{k,i} d_i D]_k
            st['Jw'] = Jw
        return st

    def reduced_jac(self, st):
        """Kaufman Jacobian: project Jw onto the complement of range([QQ | Q1])."""
        Jw = st['Jw']
        J = Jw - self.QQ @ (self.QQ.T @ Jw)
        J -= st['Q1'] @ (st['Q1'].T @ J)
        return J

    def relres(self, st):
        return math.sqrt(st['f'] / max(st['bnorm2'], 1e-300))

    # --- Levenberg-Marquardt in (u, g) coordinates
    def lm(self, u0, g0=None, maxiter=60, lam0=1e-3, tol_step=1e-11, tol_f=1e-26, verbose=False, stall=8):
        u = np.asarray(u0, float); u = u / np.linalg.norm(u)
        g = np.zeros(self.ng) if g0 is None else np.asarray(g0, float)
        st = self.evaluate(self.w_of(u, g), need_jac=True)
        lam = lam0; it = 0; nstall = 0; hist = [st['f']]
        t0 = time.time()
        while it < maxiter:
            it += 1
            J = self.reduced_jac(st)
            # chain rule: w = Pp u/|u| + G g ; at |u| = 1 : dw/du = Pp (I - u u^T)
            Ju = (J @ self.Pp); Ju = Ju - np.outer(Ju @ u, u)
            Jg = J @ self.G if self.ng else np.zeros((self.R, 0))
            Jt = np.hstack([Ju, Jg])
            JtJ = Jt.T @ Jt; Jtr = Jt.T @ st['r']
            dscale = np.maximum(np.diag(JtJ), 1e-12)
            improved = False
            for _ in range(12):
                try:
                    delta = np.linalg.solve(JtJ + lam * np.diag(dscale), -Jtr)
                except np.linalg.LinAlgError:
                    lam *= 10; continue
                du = delta[:self.dimU]; dg = delta[self.dimU:]
                un = u + du; un /= np.linalg.norm(un)
                gn = g + dg
                stn = self.evaluate(self.w_of(un, gn), need_jac=True)
                if stn['f'] < st['f']:
                    step = float(np.linalg.norm(delta))
                    u, g, st = un, gn, stn
                    lam = max(lam / 5, 1e-15)
                    improved = True
                    break
                lam *= 5
            hist.append(st['f'])
            if verbose:
                print(f'  it {it} f={st["f"]:.3e} rel={self.relres(st):.3e} lam={lam:.1e} rank={st["rank"]} ({time.time()-t0:.1f}s)', file=sys.stderr)
            if not improved:
                break
            if st['f'] < tol_f or step < tol_step:
                break
            if len(hist) > stall and hist[-stall - 1] > 0 and hist[-1] / hist[-stall - 1] > 0.9:
                nstall += 1
                if nstall >= 2:
                    break
        return {'u': u, 'g': g, 'f': st['f'], 'bnorm2': st['bnorm2'], 'relres': self.relres(st), 'iters': it,
                'rank': st['rank'], 'time': round(time.time() - t0, 2), 'w': st['w'], 'Dc': st['Dc'], 'x': st['x']}

    # --- finite-difference check of the Kaufman gradient
    def fd_check(self, seed=0, h=1e-6):
        rng = np.random.default_rng(seed)
        u = rng.standard_normal(self.dimU); u /= np.linalg.norm(u)
        g = rng.standard_normal(self.ng) * 0.3
        st = self.evaluate(self.w_of(u, g), need_jac=True)
        J = self.reduced_jac(st)
        grad_w = 2 * (J.T @ st['r'])          # d f / d w  (varpro gradient; same as 2 Jw^T r)
        v = rng.standard_normal(self.nw); v /= np.linalg.norm(v)
        w = self.w_of(u, g)
        fp = self.evaluate(w + h * v)['f']; fm = self.evaluate(w - h * v)['f']
        return float(grad_w @ v), (fp - fm) / (2 * h), float(2 * (st['Jw'].T @ st['r']) @ v)

    # --- exact side
    def D_poly_float(self, Dc):
        return self.sp['D'].poly(Dc)

    def reconstruct(self, Dc, dens=(1, 2, 3, 4, 6, 8, 12, 24, 60, 120, 360, 720, 2520, 5040, 10 ** 4, 10 ** 5, 10 ** 6), tol=1e-6):
        """Rational candidates: normalize so the largest |coefficient| = 1, then try common denominators."""
        Dc = np.asarray(Dc, float)
        j = int(np.argmax(np.abs(Dc))); a = Dc / Dc[j]
        cands = []
        for den in dens:
            v = a * den
            if np.max(np.abs(v - np.round(v))) < tol * max(1.0, den ** 0.5):
                fr = [Fraction(int(round(x)), den) for x in v]
                if fr not in cands:
                    cands.append(fr)
        fr = [Fraction(float(x)).limit_denominator(10 ** 4) for x in a]
        if max(abs(float(f) - x) for f, x in zip(fr, a)) < tol and fr not in cands:
            cands.append(fr)
        return cands

    def D_exact(self, fr):
        D = {}
        for i, m in enumerate(self.sp['D'].mons):
            if fr[i] != 0:
                D[m] = Fraction(fr[i])
        return D


def confirm_exact(F, A, q, D, verbose=False):
    """Exact fixed-D solve through the engine.  Returns (status, data)."""
    if not D:
        return 'inconsistent', {'note': 'D = 0'}
    divF, _ = jacobian.divisible_by(F, D)
    if divF:
        return 'spurious', {'note': 'F divides D'}
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    if not ok:
        return 'inconsistent', {'note': 'E1 fails', 'e1': info}
    out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose)
    if not out['consistent']:
        return 'inconsistent', {'rank': out['rank'], 'nullity': out['nullity'], 'rows': out['rows'], 'cols': out['cols']}
    if out.get('exact'):
        okv, infov = core.verify_witness(F, A, q, D, out['N'], out['Q'])
        return ('solved' if okv else 'lift-failed'), {'N': out['N'], 'Q': out['Q'], 'rank': out['rank'], 'verify': infov}
    status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
    return status, data


def save_candidate(d, D, N, Q, path):
    rec = {'instance': IID, 'd': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q),
           'engine': fixed_d.ENGINE_VERSION, 'ts': now()}
    with open(path, 'w') as f:
        json.dump(rec, f)
    return rec
