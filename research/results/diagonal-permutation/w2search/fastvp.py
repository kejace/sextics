"""fastvp: fast joint (variable-projection) search for q=1 certificates with D in a linear family.

Reduced E2 formulation.  For q = 1 and an E1-admissible D with a particular solution N0 of
    N0(F) := sum_i N0_i F_i = A * D                     (N0 linear in D)
every solution of E1 (q N(F) = A D + F H, any H) is  N = N0 + z H/6 + K.gradF  with K antisymmetric
of degree d-4 (Koszul syzygies; grad F is a regular sequence), and the certificate identity
    A D^2 - D N(F) + F (D div N - N(D)) = F^2 Q
becomes (the H-terms cancel identically)
    Phi(D) + Psi(D, K) = F Q,   Phi = D div N0 - N0(D),   Psi = D div(K.gradF) - (K.gradF)(D),
all in degree 2d.  Working modulo F means evaluating at points of X = {F = 0}; at a point pt
    Psi_{(ij),m}(pt) = D [d_i m F_j - d_j m F_i] - m [F_j d_i D - F_i d_j D]     (F_ij cancels)
so for a family D(a) = sum_k a_k B_k the problem is   min_k || Psi(a) k + Phi(a) ||^2  with
Psi linear and Phi quadratic in a (variable projection, Golub-Pereyra gradient).  P >= dim S_{2d}^-/F S_{2d-6}^-
generic points of X make the point conditions equivalent to the polynomial identity.
Only an exact re-solve (fixed_d) of a rationally reconstructed candidate counts as evidence.
"""
import sys, os, time, json, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
from hodge import core, fixed_d, jacobian

NV = 6
SIGMA = [1, 0, 2, 4, 3, 5]          # (01): swap (X0,Y0) <-> (X1,Y1)
def pf(i): return 0 if i < 3 else 1  # Y-parity of F_i (F Y-even)

def monos(deg, par=None):
    return list(core.monomials(deg, par)) if deg >= 0 else []

def poly_to_arrays(P):
    E = np.array(list(P.keys()), dtype=np.int64).reshape(-1, NV)
    c = np.array([float(v) for v in P.values()], dtype=float)
    return E, c

def coeff_matrix(polys, mons):
    """S x K matrix of coefficients of the polys (dicts) in the monomial basis mons."""
    idx = {m: i for i, m in enumerate(mons)}
    C = np.zeros((len(mons), len(polys)))
    for k, P in enumerate(polys):
        for e, v in P.items():
            C[idx[e], k] = float(v)
    return C

def shift_down(C, mons_hi, mons_lo, i):
    """Coefficient matrix of d_i applied to polys given by C over mons_hi, in basis mons_lo."""
    idx = {m: j for j, m in enumerate(mons_lo)}
    out = np.zeros((len(mons_lo), C.shape[1]))
    for r, e in enumerate(mons_hi):
        if e[i]:
            e2 = list(e); e2[i] -= 1
            out[idx[tuple(e2)]] += e[i] * C[r]
    return out

# ----------------------------------------------------------------------------- points on X
def random_points(F, P, seed, sigma_pairs=False, minabs=0.08, maxabs_ratio=1.0):
    """P random complex points on X={F=0}, scaled to max|z_i| = 1.  If sigma_pairs, P must be even and
    points come as consecutive (pt, sigma(pt)) pairs."""
    rng = np.random.default_rng(seed)
    E, c = poly_to_arrays(F)
    pts = []
    want = P // 2 if sigma_pairs else P
    while len(pts) < want:
        j = rng.integers(NV)
        z = (rng.uniform(0.5, 1.5, NV) * np.exp(2j * np.pi * rng.uniform(0, 1, NV))).astype(complex)
        # univariate sextic in z_j
        coefs = np.zeros(7, dtype=complex)
        for e, cv in zip(E, c):
            term = cv
            for i in range(NV):
                if i != j and e[i]:
                    term = term * z[i] ** e[i]
            coefs[e[j]] += term
        roots = np.roots(coefs[::-1])
        if len(roots) == 0:
            continue
        z[j] = roots[rng.integers(len(roots))]
        z = z / np.max(np.abs(z))
        if np.min(np.abs(z)) < minabs:
            continue
        pts.append(z)
    Z = np.array(pts)
    if sigma_pairs:
        Zs = Z[:, SIGMA]
        Z = np.stack([Z, Zs], axis=1).reshape(-1, NV)
    # sanity: |F(pt)| small
    return Z

def eval_monos(mons, Z):
    """P x n complex matrix of monomial values."""
    if not mons:
        return np.zeros((Z.shape[0], 0), dtype=complex)
    E = np.array(mons, dtype=float)
    return np.exp(np.log(Z) @ E.T)

# ----------------------------------------------------------------------------- families
def jacobian_family(F, A, d, parD):
    """Basis (i, m) -> B = m F_i with N0 = A m e_i; returns (labels, Bs, N0s) with exact Fractions."""
    Fi = [core.pderiv(F, i) for i in range(NV)]
    labels, Bs, N0s = [], [], []
    for i in range(NV):
        for m in monos(d - 5, (parD + pf(i)) % 2):
            B = core.shift(Fi[i], m, 1)
            N0 = [dict() for _ in range(NV)]
            N0[i] = core.shift(A, m, 1)
            labels.append(('J', i, m)); Bs.append(B); N0s.append(N0)
    return labels, Bs, N0s

def e1_particular(F, A, Ds, verbose=False):
    """Float particular solutions N0 (list of 6 float-dict polys) of sum_i N0_i F_i = A D for each D in Ds
    (same degree/parity), via one dense least squares in coefficient space.  Returns (N0s, max relative residual)."""
    t_start = time.time()
    d = core.degree(Ds[0]); parD = core.parity_of(Ds[0]); parA = core.parity_of(A)
    Fi = [core.pderiv(F, i) for i in range(NV)]
    tgt_par = (parA + parD) % 2
    rows = monos(d + 6, tgt_par)
    ridx = {m: i for i, m in enumerate(rows)}
    cols = []
    for i in range(NV):
        for m in monos(d + 1, (parA + parD + pf(i)) % 2):
            cols.append((i, m))
    M = np.zeros((len(rows), len(cols)))
    for j, (i, m) in enumerate(cols):
        for e, v in core.shift(Fi[i], m, 1).items():
            M[ridx[e], j] = float(v)
    B = np.zeros((len(rows), len(Ds)))
    for k, D in enumerate(Ds):
        for e, v in core.mul(A, D).items():
            B[ridx[e], k] = float(v)
    t0 = time.time()
    if verbose: print(f'[e1_particular] assembled {M.shape} in {time.time()-t_start:.1f}s', file=sys.stderr)
    X, *_ = sla.lstsq(M, B, cond=1e-12, lapack_driver='gelsd', check_finite=False)
    res = np.linalg.norm(M @ X - B, axis=0) / np.maximum(np.linalg.norm(B, axis=0), 1e-300)
    if verbose:
        print(f'[e1_particular] {M.shape} x {len(Ds)} rhs: max rel residual {res.max():.2e} ({time.time()-t0:.1f}s)', file=sys.stderr)
    N0s = []
    for k in range(len(Ds)):
        N0 = [dict() for _ in range(NV)]
        for j, (i, m) in enumerate(cols):
            if abs(X[j, k]) > 1e-14:
                N0[i][m] = X[j, k]
        N0s.append(N0)
    return N0s, float(res.max())

# ----------------------------------------------------------------------------- the family object
class FastFamily:
    def __init__(self, F, A, Bs, N0s, labels=None, P=None, seed=0, sector='all', verbose=False, cond=1e-9,
                 oversample=1.12, extra_pts=40):
        t0 = time.time()
        self.F, self.A = F, A
        self.Bs, self.N0s = Bs, N0s
        self.labels = labels or [str(k) for k in range(len(Bs))]
        self.K = len(Bs)
        d = core.degree(Bs[0]); self.d = d
        parD = core.parity_of(Bs[0]); parA = core.parity_of(A)
        self.parD, self.parA = parD, parA
        self.cond = cond
        self.sector = sector
        # row count needed: dim S_{2d}^{parA} - dim S_{2d-6}^{parA}
        Rp = len(monos(2 * d, parA)) - len(monos(2 * d - 6, parA))
        self.Rprime = Rp
        if P is None:
            P = int(Rp * oversample) + extra_pts
        sig = (sector == 'inv')
        if sig:
            P += P % 2
        self.P = P
        Z = random_points(F, P, seed, sigma_pairs=sig)
        self.Z = Z
        # F_i at points
        Fi = [core.pderiv(F, i) for i in range(NV)]
        self.Fv = np.stack([self._evalp(Fi[i], Z) for i in range(NV)], axis=1)   # P x 6
        Fval = self._evalp(F, Z)
        self.Fres = float(np.max(np.abs(Fval)))
        # D values and gradients
        md = monos(d, parD)
        Vd = eval_monos(md, Z)
        md1 = {par: monos(d - 1, par) for par in (0, 1)}
        Vd1 = {par: eval_monos(md1[par], Z) for par in (0, 1)}
        CD = coeff_matrix(Bs, md)
        self.Dv = Vd @ CD                                                        # P x K
        self.gradD = np.stack([Vd1[(parD + pf(i)) % 2] @ shift_down(CD, md, md1[(parD + pf(i)) % 2], i)
                               for i in range(NV)], axis=1)                      # P x 6 x K
        # N0 values and divergence
        N0v = np.zeros((P, NV, self.K), dtype=complex)
        mdiv = monos(d, (parA + parD) % 2); Vdiv = eval_monos(mdiv, Z)
        divC = np.zeros((len(mdiv), self.K))
        for i in range(NV):
            par_i = (parA + parD + pf(i)) % 2
            mi = monos(d + 1, par_i)
            Vi = eval_monos(mi, Z)
            Ci = coeff_matrix([N0[i] for N0 in N0s], mi)
            N0v[:, i, :] = Vi @ Ci
            divC += shift_down(Ci, mi, mdiv, i)
        self.N0v = N0v
        self.divN0 = Vdiv @ divC                                                 # P x K
        # Koszul unknown monomials per pair
        self.pairs = [(i, j) for i in range(NV) for j in range(i + 1, NV)]
        self.kmons = {}
        Mv = {}; dMv = {}
        for par in (0, 1):
            mm = monos(d - 4, par)
            self.kmons[par] = mm
            Mv[par] = eval_monos(mm, Z)
            Cm = np.eye(len(mm))
            dMv[par] = []
            for i in range(NV):
                mm1 = monos(d - 5, (par + pf(i)) % 2)
                V1 = eval_monos(mm1, Z)
                dMv[par].append(V1 @ shift_down(Cm, mm, mm1, i))
        self.Mv, self.dMv = Mv, dMv
        self.pair_par = [(parA + parD + pf(i) + pf(j)) % 2 for (i, j) in self.pairs]
        self.col_slices = []
        c0 = 0
        for pi, (i, j) in enumerate(self.pairs):
            n = len(self.kmons[self.pair_par[pi]])
            self.col_slices.append(slice(c0, c0 + n)); c0 += n
        self.nK = c0
        self.nQ = len(monos(2 * d - 6, parA))
        # reduced parameter coordinates: row space of Dv (removes F-multiples and dependent members)
        U, s, Vh = np.linalg.svd(self.Dv, full_matrices=False)
        tol = 1e-9 * s[0]
        r = int(np.sum(s > tol))
        self.Dsv = s
        Vr = Vh[:r].conj().T                                                     # K x r (complex in general)
        Vr = np.real_if_close(Vr, tol=1e6)
        if np.iscomplexobj(Vr):
            # take a real basis of the row space: use real and imaginary parts, orthonormalize
            Mreal = np.hstack([Vr.real, Vr.imag])
            Qr, Rr = np.linalg.qr(Mreal)
            dg = np.abs(np.diag(Rr)); keep = dg > 1e-9 * dg.max()
            Vr = Qr[:, :r] if keep.sum() >= r else Qr[:, keep]
        self.Vr = np.asarray(Vr, dtype=float)
        self.excluded = self.K - r
        if sector == 'inv':
            # (01)-invariant D: D(sigma pt) = D(pt); points are consecutive pairs
            Dp = self.Dv.reshape(P // 2, 2, self.K)
            Mi = (Dp[:, 0, :] - Dp[:, 1, :]) @ self.Vr
            Mi = np.vstack([Mi.real, Mi.imag])
            ns = sla.null_space(Mi, rcond=1e-9)
            self.Vr = self.Vr @ ns
            self.Kcols_basis = self._koszul_sector_basis()
        else:
            self.Kcols_basis = None
        self.dimP = self.Vr.shape[1]
        self.nfev = 0
        rng0 = np.random.default_rng(12345)
        a0 = self.Vr @ rng0.standard_normal(self.dimP); a0 /= np.linalg.norm(a0)
        _, phin0 = self.value(a0)
        self.fscale = 1.0 / max(phin0, 1e-300)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[fastvp] d={d} parD={parD} K={self.K} P={P} (R\'={Rp}) nK={self.nK} nQ(eliminated)={self.nQ} '
                  f'dimP={self.dimP} excluded={self.excluded} sector={sector} |F(pt)|max={self.Fres:.1e} '
                  f'built {self.build_time:.1f}s', file=sys.stderr)

    def _evalp(self, Pp, Z):
        E, c = poly_to_arrays(Pp)
        return np.exp(np.log(Z) @ E.T.astype(float)) @ c

    # ---- (01)-equivariant Koszul sector: signed permutation on K-coordinates, anti-eigenspace
    def _koszul_sector_basis(self):
        n = self.nK
        S = np.zeros((n, n))
        for pi, (i, j) in enumerate(self.pairs):
            par = self.pair_par[pi]
            mm = self.kmons[par]; midx = {m: t for t, m in enumerate(mm)}
            si, sj = SIGMA[i], SIGMA[j]
            sign = 1.0
            if si > sj:
                si, sj = sj, si; sign = -1.0
            pj = self.pairs.index((si, sj))
            assert self.pair_par[pj] == par
            sl_from = self.col_slices[pj]; sl_to = self.col_slices[pi]
            for t, m in enumerate(mm):
                sm = tuple(m[SIGMA[v]] for v in range(NV))
                S[sl_to.start + t, sl_from.start + midx[sm]] = sign
        # anti-equivariant: S k = -k  (checked numerically in self.check_sector)
        self._S = S
        ns_minus = sla.null_space(S + np.eye(n))
        ns_plus = sla.null_space(S - np.eye(n))
        self._sector_bases = {'-': ns_minus, '+': ns_plus}
        return ns_minus

    def check_sector(self, seed=1):
        """Numerically decide which Koszul eigenspace is lossless for invariant D and confirm half-row equivalence."""
        rng = np.random.default_rng(seed)
        u = rng.standard_normal(self.dimP); a = self.Vr @ u; a /= np.linalg.norm(a)
        Psi, Phi = self.assemble(a)
        full = self._lsq(Psi, Phi)[1]
        out = {}
        for key, Bk in self._sector_bases.items():
            if Bk.shape[1] == 0:
                out[key] = None; continue
            out[key] = self._lsq(Psi @ Bk, Phi)[1]
        return full, out

    # ---- assembly
    def assemble(self, a, full_rows=True):
        D = self.Dv @ a; gD = self.gradD @ a; N0 = self.N0v @ a; dN0 = self.divN0 @ a
        Phi = D * dN0 - np.sum(N0 * gD, axis=1)
        Psi = np.empty((self.P, self.nK), dtype=complex)
        for pi, (i, j) in enumerate(self.pairs):
            par = self.pair_par[pi]
            Mv = self.Mv[par]; dMi = self.dMv[par][i]; dMj = self.dMv[par][j]
            Fi = self.Fv[:, i]; Fj = self.Fv[:, j]
            block = D[:, None] * (dMi * Fj[:, None] - dMj * Fi[:, None]) - Mv * (Fj * gD[:, i] - Fi * gD[:, j])[:, None]
            Psi[:, self.col_slices[pi]] = block
        return Psi, Phi

    def _lsq(self, Psi, Phi):
        if self.sector == 'inv':
            # rows: one point per sigma-pair (residual is anti-invariant), columns already restricted by caller
            Psi = Psi[0::2]; Phi = Phi[0::2]
        M = np.vstack([Psi.real, Psi.imag]); b = -np.concatenate([Phi.real, Phi.imag])
        k, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ k - b
        return k, float(r @ r)

    def value(self, a, want_k=False):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        Psi, Phi = self.assemble(a)
        Bk = self.Kcols_basis
        if Bk is not None:
            k, f = self._lsq(Psi @ Bk, Phi); k = Bk @ k
        else:
            k, f = self._lsq(Psi, Phi)
        phin = float(np.vdot(Phi, Phi).real)
        if self.sector == 'inv':
            phin = float(np.vdot(Phi[0::2], Phi[0::2]).real)
        return (f, phin, k) if want_k else (f, phin)

    def value_and_grad(self, u, return_jac=False):
        a_full = self.Vr @ u; na = np.linalg.norm(a_full); a = a_full / na
        D = self.Dv @ a; gD = self.gradD @ a; N0 = self.N0v @ a; dN0 = self.divN0 @ a
        Phi = D * dN0 - np.sum(N0 * gD, axis=1)
        Psi, _ = self.assemble(a)
        Bk = self.Kcols_basis
        if Bk is not None:
            kk, f = self._lsq(Psi @ Bk, Phi); k = Bk @ kk
        else:
            k, f = self._lsq(Psi, Phi)
        r = Psi @ k + Phi                                     # complex residual at all points
        if self.sector == 'inv':
            mask = np.zeros(self.P, bool); mask[0::2] = True
            r = np.where(mask, r, 0)
        # dPhi/da_m
        dPhi = (self.Dv * dN0[:, None] - np.einsum('pi,pim->pm', N0, self.gradD)
                + D[:, None] * self.divN0 - np.einsum('pim,pi->pm', self.N0v, gD))
        # dPsi/da_m k:  U = sum_{i<j} (dM_i F_j - dM_j F_i) k_ij ;  c_i = sum_j F_j w_ij, w antisym, w_ij = Mv k_ij
        U = np.zeros(self.P, dtype=complex); w = {}
        for pi, (i, j) in enumerate(self.pairs):
            par = self.pair_par[pi]; kij = k[self.col_slices[pi]]
            U += (self.dMv[par][i] @ kij) * self.Fv[:, j] - (self.dMv[par][j] @ kij) * self.Fv[:, i]
            w[(i, j)] = self.Mv[par] @ kij
        c = np.zeros((self.P, NV), dtype=complex)
        for (i, j), wij in w.items():
            c[:, i] += self.Fv[:, j] * wij
            c[:, j] -= self.Fv[:, i] * wij
        dPsik = self.Dv * U[:, None] - np.einsum('pim,pi->pm', self.gradD, c)
        Ja = dPhi + dPsik                                     # complex P x K Jacobian of the residual wrt a (Kaufman approx.)
        if return_jac:
            return f, r, Ja, a
        g = 2.0 * np.real(np.conj(r) @ Ja)
        g_full = (g - float(g @ a) * a) / na
        self.nfev += 1
        sc = getattr(self, 'fscale', 1.0)
        return f * sc, (self.Vr.T @ g_full) * sc

    def polish(self, a, iters=12, verbose=False):
        """Kaufman/Gauss-Newton refinement of a near-solution: delta = -(P_perp J_a)^+ r with P_perp the projector off col(Psi)."""
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        best_a, best_f = a, self.value(a)[0]
        for it in range(iters):
            u = self.Vr.T @ a
            f, r, Ja, a_n = self.value_and_grad(u, return_jac=True)
            Psi, _ = self.assemble(a_n)
            Bk = self.Kcols_basis
            if Bk is not None: Psi = Psi @ Bk
            if self.sector == 'inv':
                rows = slice(0, None, 2)
                Psi = Psi[rows]; r = r[rows]; Ja = Ja[rows]
            M = np.vstack([Psi.real, Psi.imag]); Jr = np.vstack([Ja.real, Ja.imag]) @ self.Vr; rr = np.concatenate([r.real, r.imag])
            Qm, _ = np.linalg.qr(M)
            Jp = Jr - Qm @ (Qm.T @ Jr); rp = rr - Qm @ (Qm.T @ rr)
            # remove the radial direction (a on the sphere): restrict to tangent space
            du, *_ = sla.lstsq(Jp, -rp, cond=1e-10, lapack_driver='gelsd', check_finite=False)
            a_new = a_n + self.Vr @ du; a_new /= np.linalg.norm(a_new)
            f_new = self.value(a_new)[0]
            if verbose:
                print(f'[polish] it {it}: f {f:.3e} -> {f_new:.3e} |du|={np.linalg.norm(du):.2e}', file=sys.stderr)
            if f_new < best_f:
                best_a, best_f = a_new, f_new
            if f_new > 0.5 * f or np.linalg.norm(du) < 1e-15:
                break
            a = a_new
        return best_a, best_f

    def fd_check(self, seed=0, h=1e-6):
        rng = np.random.default_rng(seed)
        u = rng.standard_normal(self.dimP); u /= np.linalg.norm(u)
        f0, g = self.value_and_grad(u)
        v = rng.standard_normal(self.dimP); v /= np.linalg.norm(v)
        fp = self.value_and_grad(u + h * v)[0]; fm = self.value_and_grad(u - h * v)[0]
        fd = (fp - fm) / (2 * h)
        return f0, float(g @ v), fd

    # ---- search
    def search(self, n_restarts=10, seed=0, maxiter=400, gtol=1e-12, verbose=False, stop_rel=1e-13, log=None, tag='', polish_rel=1e-5):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            u0 = rng.standard_normal(self.dimP); u0 /= np.linalg.norm(u0)
            t1 = time.time()
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            a = self.Vr @ res.x; a /= np.linalg.norm(a)
            f, phin = self.value(a)
            rel = f / max(phin, 1e-300)
            rel_lbfgs = rel
            if rel < polish_rel:
                a, f = self.polish(a)
                phin = self.value(a)[1]
                rel = f / max(phin, 1e-300)
            Dn = float(np.linalg.norm(self.Dv @ a)) / math.sqrt(self.P)
            cand = {'f': f, 'phi2': phin, 'rel': rel, 'rel_lbfgs': rel_lbfgs, 'iters': int(res.nit), 'restart': it, 'seed': seed,
                    'Dnorm_onX': Dn, 'seconds': round(time.time() - t1, 2), 'a': a.tolist(), 'msg': str(res.message)[:60]}
            out.append(cand)
            if log is not None:
                log({'exp': 'fastvp-restart', 'tag': tag, **{k: v for k, v in cand.items() if k != 'a'}})
            if verbose:
                print(f'[fastvp {tag}] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} ({time.time()-t0:.1f}s)', file=sys.stderr)
            if rel < stop_rel:
                break
        out.sort(key=lambda c: c['rel'])
        return out

    # ---- exact side
    def D_exact(self, frac):
        D = {}
        for k, B in enumerate(self.Bs):
            if frac[k]:
                D = core.add(D, core.scale(B, Fraction(frac[k])))
        return D

    def reconstruct(self, a, max_dens=(12, 100, 1000, 10000, 100000), tol=1e-7):
        a = np.asarray(a, float)
        order = np.argsort(-np.abs(a))
        cands = []
        for j in order[:3]:
            b = a / a[j]
            for md in max_dens:
                fr = [Fraction(float(v)).limit_denominator(md) for v in b]
                err = max(abs(float(f) - v) for f, v in zip(fr, b))
                if err < tol:
                    if tuple(fr) not in cands:
                        cands.append(tuple(fr))
                    break
        return cands

    def confirm(self, frac, verbose=False):
        D = self.D_exact(frac)
        if not D:
            return 'inconsistent', {'note': 'D=0'}
        divF, _ = jacobian.divisible_by(self.F, D)
        if divF:
            return 'spurious', {'note': 'F | D'}
        ok, info = fixed_d.e1_solvable(self.F, self.A, 1, D)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails', 'D': core.to_str(D)}
        out = fixed_d.solve_fixed_D_fast(self.F, self.A, 1, D, verbose=verbose)
        if not out['consistent']:
            return 'inconsistent', {'rank': out['rank'], 'nullity': out['nullity'], 'D': core.to_str(D)}
        if out.get('exact'):
            okv, vinfo = core.verify_witness(self.F, self.A, 1, D, out['N'], out['Q'])
            return ('solved' if okv else 'lift-bad'), {'N': out['N'], 'Q': out['Q'], 'D': D, 'verify': vinfo}
        return 'consistent-modp', {'rank': out['rank'], 'nullity': out['nullity'], 'D': core.to_str(D)}
