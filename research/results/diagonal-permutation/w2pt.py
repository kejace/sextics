"""Pointwise variable-projection search for q=1 certificates (wave 2, lane diagonal-permutation).

Reduced E2 formulation.  For a denominator D in a linear family D(a) = sum_k a_k D_k with E1 particular solutions
N^k (N^k(F) = A D_k exactly), every E1 lift is N = N^0(a) + (H/6) E + K_c (E = Euler field, K_c = Koszul field with
skew coefficients c_ij of degree d-4), the H-term drops out of E2 identically, and E2 reads

    B_D(N^0) + T_D(c) = 0  mod F,      B_D(N) = D div N - N(D),      T_D(c) = B_D(K_c)
    T_D(c) = sum_{i<j} D [(d_i c_ij) F_j - (d_j c_ij) F_i] - c_ij (D_i F_j - D_j F_i).

'mod F' is tested by evaluation at random complex points of V(F) (F irreducible, so vanishing on V(F) <=> divisible
by F).  With r > dim(a) + dim(c) generic points the compressed system has no false zeros (dimension count).  The
reduced functional f(a) = min_c |rho(a,c)|^2 is minimized on the unit sphere of a (orthogonal to the spurious
directions with D(a) in F*S_{d-6} or D(a) = 0) by L-BFGS with the exact variable-projection gradient.
Everything is pointwise: no large matrices are stored; per evaluation one dense least squares (2r x nK).
"""
from __future__ import annotations
import sys, time, json, math
import numpy as np
import scipy.linalg as sla
from fractions import Fraction
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from hodge.core import monomials, NV, pderiv, mul, shift, add, scale

PAR = [0, 0, 0, 1, 1, 1]   # Y-parity of the coordinates


def exps_coefs(P):
    if not P:
        return np.zeros((0, NV), dtype=np.int64), np.zeros(0)
    E = np.array(list(P.keys()), dtype=np.int64)
    C = np.array([float(v) for v in P.values()])
    return E, C


class Points:
    """r random complex points on V(F) (unit-modulus coordinates except the solved one)."""
    def __init__(self, F, r, seed, maxdeg=40):
        rng = np.random.default_rng(seed)
        self.r = r
        Z = np.zeros((r, NV), dtype=complex)
        E, C = exps_coefs(F)
        degF = int(E.sum(axis=1).max())
        for t in range(r):
            while True:
                j = int(rng.integers(NV))           # coordinate to solve for
                z = np.exp(2j * np.pi * rng.random(NV))
                others = [i for i in range(NV) if i != j]
                vals = z[others]
                poly = np.zeros(degF + 1, dtype=complex)
                mon = np.prod(vals[None, :] ** E[:, others], axis=1)
                for k in range(len(C)):
                    poly[E[k, j]] += C[k] * mon[k]
                if abs(poly[degF]) < 1e-12:
                    continue
                roots = np.roots(poly[::-1])
                roots = roots[np.isfinite(roots)]
                if len(roots) == 0:
                    continue
                rt = roots[np.argmin(np.abs(np.abs(roots) - 1.0))]
                if abs(rt) < 0.2 or abs(rt) > 5:
                    continue
                z[j] = rt
                Z[t] = z
                break
        self.Z = Z
        self.pows = np.ones((NV, maxdeg + 1, r), dtype=complex)
        for i in range(NV):
            for e in range(1, maxdeg + 1):
                self.pows[i, e] = self.pows[i, e - 1] * Z[:, i]
        self.F = F
        self.Fres = np.max(np.abs(self.eval_poly(F)))

    def eval_exps(self, E):
        """values of the monomials with exponent rows E: (len(E), r)."""
        E = np.asarray(E, dtype=np.int64)
        out = np.ones((len(E), self.r), dtype=complex)
        for i in range(NV):
            out *= self.pows[i, E[:, i], :]
        return out

    def eval_poly(self, P):
        E, C = exps_coefs(P)
        if len(C) == 0:
            return np.zeros(self.r, dtype=complex)
        return C @ self.eval_exps(E)


# ----------------------------------------------------------------------------------------------------------------
class Family:
    """Linear family of denominators with E1 lifts, evaluated at points.

    members: list of dicts with keys
        'D'    : exact poly (dict) or None
        'Dval' : (r,) complex   D_k(z)
        'Dder' : (6, r)         d_i D_k (z)
        'divN' : (r,)           div N^k (z)
        'Nval' : (6, r)         N^k_i (z)
    spurious: (K, s) array whose columns span the directions a with D(a) in F*S_{d-6} or D(a)=0 (to be removed).
    """
    def __init__(self, d, members, spurious, label):
        self.d, self.label = d, label
        self.K = len(members)
        self.members = members
        self.Dval = np.array([m['Dval'] for m in members])            # K x r
        self.Dder = np.array([m['Dder'] for m in members])            # K x 6 x r
        self.divN = np.array([m['divN'] for m in members])            # K x r
        self.Nval = np.array([m['Nval'] for m in members])            # K x 6 x r
        self.spurious = spurious
        self.exact = all(m.get('D') is not None for m in members)

    def D_exact(self, a_frac):
        D = {}
        for k, m in enumerate(self.members):
            if a_frac[k]:
                D = add(D, scale(m['D'], Fraction(a_frac[k])))
        return D


def jacobian_family(F, A, d, pts, extra=None):
    """D = w(F), w a polynomial field of degree d-5 with Y-even D; N^k = A m d_i for D_k = m F_i.
    extra: list of (D poly, N list-of-6 polys) additional members (e.g. non-Jacobian kernel elements with E1 lifts)."""
    Fi = [pderiv(F, i) for i in range(NV)]
    Ai = [pderiv(A, i) for i in range(NV)]
    Fi_val = [pts.eval_poly(P) for P in Fi]
    A_val = pts.eval_poly(A)
    Ai_val = [pts.eval_poly(P) for P in Ai]
    members = []
    wbasis = []
    for i in range(NV):
        for m in monomials(d - 5, PAR[i]):
            wbasis.append((i, m))
            mp = {m: Fraction(1)}
            Dk = mul(mp, Fi[i])
            mval = pts.eval_poly(mp)
            dm = [pts.eval_poly(pderiv(mp, j)) for j in range(NV)]
            Dder = np.array([pts.eval_poly(pderiv(Dk, j)) for j in range(NV)])
            # N^k = A m d_i : div = m d_i A + A d_i m ; N_i = A m
            divN = mval * Ai_val[i] + A_val * dm[i]
            Nval = np.zeros((NV, pts.r), dtype=complex); Nval[i] = A_val * mval
            members.append({'D': Dk, 'Dval': pts.eval_poly(Dk), 'Dder': Dder, 'divN': divN, 'Nval': Nval, 'w': (i, m)})
    K0 = len(members)
    # spurious directions: w = l * Euler (D = 6 l F), l in S_{d-6}^+ ; Koszul fields (D = 0): c_ij constants etc.
    idx = {wb: k for k, wb in enumerate(wbasis)}
    sp = []
    if d >= 6:
        for l in monomials(d - 6, 0):
            v = np.zeros(K0)
            for i in range(NV):
                e = list(l); e[i] += 1
                v[idx[(i, tuple(e))]] = 1.0
            sp.append(v)
    if d >= 10:
        # Koszul: w_i = c F_j, w_j = -c F_i with deg c = d-10, parity so that w_i has parity PAR[i] (Y-even D)
        for i in range(NV):
            for j in range(i + 1, NV):
                for c in monomials(d - 10, None):
                    cp = {c: Fraction(1)}
                    wi = mul(cp, Fi[j]); wj = scale(mul(cp, Fi[i]), -1)
                    if core.parity_of(wi) != PAR[i] or core.parity_of(wj) != PAR[j]:
                        continue
                    v = np.zeros(K0)
                    for mm, cc in wi.items(): v[idx[(i, mm)]] += float(cc)
                    for mm, cc in wj.items(): v[idx[(j, mm)]] += float(cc)
                    sp.append(v)
    if extra:
        for (Dx, Nx) in extra:
            Dder = np.array([pts.eval_poly(pderiv(Dx, j)) for j in range(NV)])
            divN = np.zeros(pts.r, dtype=complex)
            Nval = np.zeros((NV, pts.r), dtype=complex)
            for i in range(NV):
                divN += pts.eval_poly(pderiv(Nx[i], i))
                Nval[i] = pts.eval_poly(Nx[i])
            members.append({'D': Dx if all(isinstance(v, Fraction) for v in Dx.values()) else None,
                            'Dval': pts.eval_poly(Dx), 'Dder': Dder, 'divN': divN, 'Nval': Nval, 'w': None})
    K = len(members)
    S = np.zeros((K, len(sp)))
    for j, v in enumerate(sp):
        S[:len(v), j] = v
    return Family(d, members, S, f'jac-d{d}' + (f'+{len(extra)}extra' if extra else ''))


# ----------------------------------------------------------------------------------------------------------------
class KoszulBlock:
    """Koszul-correction columns T_D(c) evaluated at the points (everything independent of a precomputed)."""
    def __init__(self, F, A, d, pts):
        self.d = d
        pD = 0; pA = core.parity_of(A)
        pars = [(pA - pD - PAR[i]) % 2 for i in range(NV)]       # parity of N_i
        Fi = [pderiv(F, i) for i in range(NV)]
        self.Fval = [pts.eval_poly(P) for P in Fi]
        self.pairs = []
        self.Cval = []     # per pair: (r x n_c) values c(z)
        self.Amat = []     # per pair: (r x n_c) values (d_i c) F_j - (d_j c) F_i
        self.cols = []     # (i, j, monomial)
        for i in range(NV):
            for j in range(i + 1, NV):
                parc = (pars[i] - PAR[j]) % 2
                cs = monomials(d - 4, parc)
                if not cs:
                    continue
                E = np.array(cs, dtype=np.int64)
                Cv = pts.eval_exps(E).T                                    # r x n_c
                # derivatives: d_i c = e_i * c / z_i
                def dmon(E, i):
                    Ei = E.copy(); mask = Ei[:, i] > 0
                    out = np.zeros((len(E), pts.r), dtype=complex)
                    if mask.any():
                        Em = Ei[mask].copy(); Em[:, i] -= 1
                        out[mask] = (Ei[mask, i][:, None]) * pts.eval_exps(Em)
                    return out.T                                            # r x n_c
                Am = dmon(E, i) * self.Fval[j][:, None] - dmon(E, j) * self.Fval[i][:, None]
                self.pairs.append((i, j)); self.Cval.append(Cv); self.Amat.append(Am)
                self.cols += [(i, j, c) for c in cs]
        self.nK = len(self.cols)
        self.offsets = np.cumsum([0] + [C.shape[1] for C in self.Cval])


class VarPro:
    def __init__(self, fam: Family, kb: KoszulBlock, pts: Points, cond=1e-10, csector=None):
        self.fam, self.kb, self.pts, self.cond = fam, kb, pts, cond
        self.r = pts.r
        self.csector = csector                      # nK x nK' matrix (orbit sums) or None
        self.nKc = kb.nK if csector is None else csector.shape[1]
        K = fam.K
        # orthonormal complement of spurious directions
        S = fam.spurious
        if S.size and S.shape[1]:
            U, s, _ = np.linalg.svd(S, full_matrices=True)
            rk = int(np.sum(s > 1e-10 * max(1.0, s[0])))
            self.P = U[:, rk:]
            self.nspur = rk
        else:
            self.P = np.eye(K); self.nspur = 0
        self.dimP = self.P.shape[1]
        self.nfev = 0

    # --- assembly -------------------------------------------------------------------------------------------
    def assemble(self, a):
        fam, kb = self.fam, self.kb
        D = a @ fam.Dval                                   # r
        Di = np.einsum('k,kir->ir', a, fam.Dder)           # 6 x r
        divN = a @ fam.divN
        Ni = np.einsum('k,kir->ir', a, fam.Nval)
        B0 = D * divN - np.sum(Ni * Di, axis=0)            # r
        blocks = []
        G = {}
        for (i, j), Cv, Am in zip(kb.pairs, kb.Cval, kb.Amat):
            g = Di[i] * kb.Fval[j] - Di[j] * kb.Fval[i]
            G[(i, j)] = g
            blocks.append(D[:, None] * Am - g[:, None] * Cv)
        M = np.hstack(blocks)                              # r x nK complex
        return M, B0, D, Di, divN, Ni, G

    def solve(self, a):
        M, B0, D, Di, divN, Ni, G = self.assemble(a)
        Ms = M if self.csector is None else M @ self.csector
        Mr = np.vstack([Ms.real, Ms.imag]); br = np.concatenate([(-B0).real, (-B0).imag])
        c, *_ = sla.lstsq(Mr, br, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        if self.csector is not None:
            c = self.csector @ c
        rho = M @ c + B0
        self.nfev += 1
        return c, rho, (M, B0, D, Di, divN, Ni, G)

    def value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        c, rho, aux = self.solve(a)
        B0 = aux[1]
        f = float(np.vdot(rho, rho).real)
        bn = float(np.vdot(B0, B0).real)
        return f, f / max(bn, 1e-300), c, rho

    def value_and_grad(self, u):
        fam, kb = self.fam, self.kb
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        c, rho, (M, B0, D, Di, divN, Ni, G) = self.solve(a)
        f = float(np.vdot(rho, rho).real)
        # d rho / d a_k = Dval_k * uvec - sum_ij (Dder_k,i F_j - Dder_k,j F_i) * v_ij  +  dB0/da_k
        uvec = np.zeros(self.r, dtype=complex)
        drho = np.zeros((fam.K, self.r), dtype=complex)
        for t, ((i, j), Cv, Am) in enumerate(zip(kb.pairs, kb.Cval, kb.Amat)):
            cc = c[kb.offsets[t]:kb.offsets[t + 1]]
            uvec += Am @ cc
            v = Cv @ cc
            drho -= (fam.Dder[:, i, :] * kb.Fval[j][None, :] - fam.Dder[:, j, :] * kb.Fval[i][None, :]) * v[None, :]
        drho += fam.Dval * uvec[None, :]
        # dB0/da_k = Dval_k divN + D divN_k - sum_i (Nval_k,i Di + Ni Dder_k,i)
        drho += fam.Dval * divN[None, :] + D[None, :] * fam.divN
        drho -= np.einsum('kir,ir->kr', fam.Nval, Di) + np.einsum('ir,kir->kr', Ni, fam.Dder)
        g = 2.0 * (drho @ np.conj(rho)).real          # K
        g_full = (g - float(g @ a) * a) / na
        return f, self.P.T @ g_full

    def bnorm2(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        _, B0, *_ = self.assemble(a)
        return float(np.vdot(B0, B0).real)

    # --- search --------------------------------------------------------------------------------------------
    def search(self, n_restarts, seed, maxiter=400, gtol=1e-12, callback=None, init=None):
        from scipy.optimize import minimize
        rng = np.random.default_rng(seed)
        out = []
        for it in range(n_restarts):
            t0 = time.time()
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 30})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            f, rel, c, rho = self.value(a)
            cand = {'f': f, 'rel': rel, 'a': a, 'c': c, 'iters': int(res.nit), 'restart': it, 'seconds': round(time.time() - t0, 2)}
            out.append(cand)
            if callback:
                callback(cand)
        out.sort(key=lambda z: z['rel'])
        return out

    def refine(self, a, maxiter=2000):
        from scipy.optimize import minimize
        u0 = self.P.T @ np.asarray(a, float)
        res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                       options={'maxiter': maxiter, 'gtol': 1e-16, 'ftol': 1e-20, 'maxcor': 50})
        a = self.P @ res.x; a /= np.linalg.norm(a)
        return a


# ----------------------------------------------------------------------------------------------------------------
def rational_candidates(a, max_den=10 ** 6, tol=1e-7):
    """Rational vectors near a (normalized by the largest entry)."""
    a = np.asarray(a, float)
    j = int(np.argmax(np.abs(a))); a = a / a[j]
    out = []
    for md in (10, 100, 1000, 10 ** 4, 10 ** 5, max_den):
        fr = [Fraction(float(v)).limit_denominator(md) for v in a]
        err = max(abs(float(f) - v) for f, v in zip(fr, a))
        if err < tol:
            t = tuple(fr)
            if t not in out:
                out.append(t)
    return out


def confirm_exact(F, A, q, D, verbose=False):
    """Exact fixed-D solve (engine).  Returns (status, data)."""
    if not D:
        return 'inconsistent', {'note': 'D=0'}
    from hodge import jacobian
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    if not ok:
        return 'inconsistent', {'note': 'e1 fails', 'e1': info}
    divF, _ = jacobian.divisible_by(F, D) if core.degree(D) >= 6 else (False, None)
    if divF:
        return 'spurious', {'note': 'F divides D'}
    out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose, want_exact=True)
    if not out['consistent']:
        return 'inconsistent', {'rank': out['rank'], 'nullity': out['nullity']}
    if out.get('exact'):
        return 'solved', out
    status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
    return status, data
