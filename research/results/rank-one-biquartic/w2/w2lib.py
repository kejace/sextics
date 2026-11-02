"""Reduced-E2 variable-projection search library (wave 2, lane rank-one-biquartic).

Normalized certificate system at q = 1 (book 5.2):
  E1:  N(F) = sum_i N_i F_i = A D                (N_i of degree d+1)
  E2:  B_D(N) := D div N - N(D)  in  F * S_{2d-6}   (degree 2d)
Every E1 lift is N0 + K(c), K(c) = sum_{i<j} c_ij (F_j e_i - F_i e_j), deg c_ij = d-4.
For D in a linear family D(a) = sum_k a_k D_k with linear E1 lifts N0(a) = sum_k a_k N0_k,
E2 is linear in (c, C):  T_{D(a)}(c) - F C = -B_{D(a)}(N0(a)).
We reduce modulo F (monomial reduction by X0^6) so that C disappears:
  R_F T_{D(a)}(c) = -R_F B_{D(a)}(N0(a)),
and minimize f(a) = min_c || R_F ( T_{D(a)} c + B_{D(a)}(N0(a)) ) ||^2 over the unit sphere.
'Hessian relaxation' = same with c = 0 (N = N0 only).

Everything is float64 numpy/scipy; exact checks (mod p / rational) are separate functions.
"""
from __future__ import annotations
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power

IID = 'rank-one-biquartic'
WORK = f'/tmp/hodge/results/{IID}'
W2 = f'{WORK}/w2'
BATCHLOG = f'{WORK}/wave2-search.log.jsonl'
PARF = [0, 0, 0, 1, 1, 1]   # Y-parity of F_i = dF/dz_i


def log_batch(entry, path=BATCHLOG):
    entry.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    entry.setdefault('lane', IID)
    entry.setdefault('engine', 2)
    with open(path, 'a') as fh:
        fh.write(json.dumps(entry, default=str) + '\n')


def to_float(P):
    return {e: float(v) for e, v in P.items()}


def unit(e, i):
    return tuple(e[j] + (1 if j == i else 0) for j in range(NV))


def addexp(e, f):
    return tuple(a + b for a, b in zip(e, f))


class Index:
    def __init__(self, mons):
        self.mons = list(mons)
        self.idx = {m: i for i, m in enumerate(self.mons)}
        self.n = len(self.mons)

    def vec(self, P):
        v = np.zeros(self.n)
        for e, c in P.items():
            v[self.idx[e]] += float(c)
        return v

    def poly(self, v, tol=0.0):
        return {m: float(v[i]) for i, m in enumerate(self.mons) if abs(v[i]) > tol}


# ----------------------------------------------------------------------------
# E1 (annihilator) spaces
# ----------------------------------------------------------------------------
def e1_matrix_float(F, A, d, pD):
    """Dense float matrix [A*m_j | -(m F_i)] in degree d+6 (rows = monomials of parity pA+pD)."""
    pA = core.parity_of(A)
    Dm = Index(monomials(d, pD))
    parsN = [(pA - pD - PARF[i]) % 2 for i in range(NV)]
    Ncols = [(i, m) for i in range(NV) for m in monomials(d + 1, parsN[i])]
    rows = Index(monomials(d + 6, (pA + pD) % 2))
    Fi = [to_float(pderiv(F, i)) for i in range(NV)]
    Af = to_float(A)
    M = np.zeros((rows.n, Dm.n + len(Ncols)))
    for j, m in enumerate(Dm.mons):
        for e, c in Af.items():
            M[rows.idx[addexp(e, m)], j] += c
    for j, (i, m) in enumerate(Ncols):
        for e, c in Fi[i].items():
            M[rows.idx[addexp(e, m)], Dm.n + j] -= c
    return M, Dm, Ncols, rows


def nullspace_float(M, nullity=None, gap_report=True, driver='gesdd'):
    """Right nullspace by SVD; if nullity is given use the last `nullity` right singular vectors."""
    t0 = time.time()
    U, s, Vt = sla.svd(M, full_matrices=True, lapack_driver=driver, check_finite=False)
    n = M.shape[1]
    s_full = np.zeros(n); s_full[:len(s)] = s
    if nullity is None:
        tol = max(M.shape) * np.finfo(float).eps * (s[0] if len(s) else 1.0) * 10
        nullity = int(np.sum(s_full <= tol))
    Z = Vt[n - nullity:, :].T
    info = {'nullity': nullity, 'smax': float(s_full[0]), 'last_nonzero': float(s_full[n - nullity - 1]) if n - nullity - 1 >= 0 else None,
            'first_zero': float(s_full[n - nullity]) if nullity > 0 else None, 't_svd': round(time.time() - t0, 1)}
    return Z, info


def e1_nullity_modp(F, A, d, pD, p=None):
    """Exact (mod p) nullity of the E1 matrix and the admissible dimension, via jacobian.mult_kernel."""
    out = jacobian.mult_kernel(F, A, d, p=p, parity=pD)
    return out


class Family:
    """Linear family of denominators with linear E1 lifts.  Dbasis: nD x K (orthonormal), N0: nN x K."""
    def __init__(self, d, pD, Dm, Ncols, Dbasis, N0, tag=''):
        self.d, self.pD, self.Dm, self.Ncols = d, pD, Dm, Ncols
        self.Dbasis = np.asarray(Dbasis, float)
        self.N0 = np.asarray(N0, float)
        self.K = self.Dbasis.shape[1]
        self.tag = tag

    def D_vec(self, a):
        return self.Dbasis @ a

    def N0_vec(self, a):
        return self.N0 @ a

    def D_poly(self, a, tol=1e-14):
        return self.Dm.poly(self.D_vec(a), tol)


def _orth_with_lifts(Dcols, Ncols_mat, tol=1e-10):
    """Orthonormalize the columns of Dcols (nD x k) and transport the lifts (nN x k) along.
    Returns (U, N0U) with U orthonormal spanning col(Dcols) and N0U lifts of U's columns."""
    if Dcols.shape[1] == 0:
        return Dcols, Ncols_mat
    U, s, Vt = sla.svd(Dcols, full_matrices=False, check_finite=False)
    r = int(np.sum(s > tol * max(1.0, s[0])))
    U = U[:, :r]
    # Dcols = U S Vt  ->  U = Dcols Vt^T S^{-1}  -> lifts: N0 Vt^T S^{-1}
    coeff = Vt[:r, :].T / s[:r]
    return U, Ncols_mat @ coeff


def j_family(F, A, d, pD, Dm=None, Ncols=None, exclude_F=True):
    """Jacobian-type family D = dF(G) = sum g_{i,m} m F_i with exact lifts N0 = A G (w-space)."""
    pA = core.parity_of(A)
    if Dm is None:
        Dm = Index(monomials(d, pD))
    parsN = [(pA - pD - PARF[i]) % 2 for i in range(NV)]
    if Ncols is None:
        Ncols = [(i, m) for i in range(NV) for m in monomials(d + 1, parsN[i])]
    Nidx = {c: j for j, c in enumerate(Ncols)}
    Fi = [to_float(pderiv(F, i)) for i in range(NV)]
    Af = to_float(A)
    cols = []; lifts = []
    for i in range(NV):
        for m in monomials(d - 5, (pD - PARF[i]) % 2):
            v = np.zeros(Dm.n)
            for e, c in Fi[i].items():
                v[Dm.idx[addexp(e, m)]] += c
            n = np.zeros(len(Ncols))
            for e, c in Af.items():
                n[Nidx[(i, addexp(e, m))]] += c
            cols.append(v); lifts.append(n)
    if not cols:
        return Family(d, pD, Dm, Ncols, np.zeros((Dm.n, 0)), np.zeros((len(Ncols), 0)), tag='J')
    Dc = np.array(cols).T; Nc = np.array(lifts).T
    U, N0 = _orth_with_lifts(Dc, Nc)
    fam = Family(d, pD, Dm, Ncols, U, N0, tag='J')
    if exclude_F and d >= 6:
        fam = remove_F_multiples(F, fam)
    return fam


def F_multiple_basis(F, d, pD, Dm):
    Ff = to_float(F)
    cols = []
    for m in monomials(d - 6, pD) if d >= 6 else []:
        v = np.zeros(Dm.n)
        for e, c in Ff.items():
            v[Dm.idx[addexp(e, m)]] += c
        cols.append(v)
    if not cols:
        return np.zeros((Dm.n, 0))
    W = np.array(cols).T
    Q, _ = np.linalg.qr(W)
    return Q


def remove_F_multiples(F, fam):
    """Restrict the family to the orthogonal complement of F * S_{d-6} (spurious witnesses)."""
    W = F_multiple_basis(F, fam.d, fam.pD, fam.Dm)
    if W.shape[1] == 0:
        return fam
    P = fam.Dbasis - W @ (W.T @ fam.Dbasis)
    U, N0 = _orth_with_lifts(P, fam.N0)
    out = Family(fam.d, fam.pD, fam.Dm, fam.Ncols, U, N0, tag=fam.tag + '-noF')
    out.removed_F_dim = fam.K - U.shape[1]
    return out


def full_family(F, A, d, pD, nullity=None, exclude_F=True, verbose=False, J=None):
    """Full E1-admissible family (kernel + J_d) with lifts from the float nullspace of the E1 matrix.
    The J-part keeps exact lifts A G; the non-J complement comes from the SVD nullspace."""
    t0 = time.time()
    M, Dm, Ncols, rows = e1_matrix_float(F, A, d, pD)
    if verbose:
        print(f'[full_family] d={d} pD={pD} E1 matrix {M.shape} built {time.time()-t0:.1f}s', file=sys.stderr, flush=True)
    Z, info = nullspace_float(M, nullity=nullity)
    if verbose:
        print(f'[full_family] nullspace {info}', file=sys.stderr, flush=True)
    nD = Dm.n
    ZD = Z[:nD, :]; ZN = Z[nD:, :]
    # admissible D space = col(ZD)
    UD, sD, VtD = sla.svd(ZD, full_matrices=False, check_finite=False)
    r = int(np.sum(sD > 1e-9 * max(1.0, sD[0])))
    UD = UD[:, :r]
    # lifts: minimal-norm z with ZD z = u  ->  z = VtD^T S^-1 UD^T u
    coeff = VtD[:r, :].T / sD[:r]
    N0D = ZN @ coeff            # lifts of the columns of UD
    if J is None:
        J = j_family(F, A, d, pD, Dm, Ncols, exclude_F=False)
    # complement of J inside the admissible space
    Pc = UD - J.Dbasis @ (J.Dbasis.T @ UD)
    Uc, sc, Vtc = sla.svd(Pc, full_matrices=False, check_finite=False)
    rc = int(np.sum(sc > 1e-7 * max(1.0, sc[0])))
    Uc = Uc[:, :rc]
    # lifts of Uc: Uc = Pc Vtc^T S^-1 = (UD - J J^T UD) Vtc^T S^-1 -> lift = (N0D - J.N0 (J^T UD)) Vtc^T S^-1
    coef2 = Vtc[:rc, :].T / sc[:rc]
    N0c = (N0D - J.N0 @ (J.Dbasis.T @ UD)) @ coef2
    Dbasis = np.hstack([J.Dbasis, Uc]); N0 = np.hstack([J.N0, N0c])
    fam = Family(d, pD, Dm, Ncols, Dbasis, N0, tag='full')
    fam.info = dict(info, adm_dim=r, J_dim=J.K, kernel_dim=rc, sD_gap=(float(sD[r-1]), float(sD[r]) if r < len(sD) else None),
                    sc_gap=(float(sc[rc-1]) if rc else None, float(sc[rc]) if rc < len(sc) else None), t=round(time.time()-t0, 1))
    fam.Ucomplement = Uc; fam.N0complement = N0c
    if exclude_F and d >= 6:
        fam = remove_F_multiples(F, fam)
        fam.info = dict(getattr(fam, 'info', {}), **{'adm_dim': r, 'J_dim': J.K, 'kernel_dim': rc})
    return fam


def subfamily(fam, cols, tag):
    out = Family(fam.d, fam.pD, fam.Dm, fam.Ncols, fam.Dbasis[:, cols], fam.N0[:, cols], tag=tag)
    return out


def check_family_e1(F, A, fam, nsamples=3, seed=0):
    """Residual of N0(a)(F) - A D(a) for random a (should be ~1e-12 relative)."""
    rng = np.random.default_rng(seed)
    Fi = [to_float(pderiv(F, i)) for i in range(NV)]
    Af = to_float(A)
    out = []
    for _ in range(nsamples):
        a = rng.standard_normal(fam.K); a /= np.linalg.norm(a)
        Dv = fam.D_vec(a); Nv = fam.N0_vec(a)
        D = fam.Dm.poly(Dv); N = [{} for _ in range(NV)]
        for (i, m), c in zip(fam.Ncols, Nv):
            if c:
                N[i][m] = c
        lhs = {}
        for i in range(NV):
            lhs = add(lhs, mul(N[i], Fi[i]))
        rhs = mul(Af, D)
        diff = sub(lhs, rhs)
        num = np.sqrt(sum(v * v for v in diff.values())); den = np.sqrt(sum(v * v for v in rhs.values()))
        out.append(num / den)
    return out


# ----------------------------------------------------------------------------
# Reduced E2 system
# ----------------------------------------------------------------------------
class ReducedE2:
    def __init__(self, F, A, d, pD, Dm=None, Ncols=None, build_koszul=True, verbose=False, reduce_var=0):
        t0 = time.time()
        self.F, self.A, self.d, self.pD = F, A, d, pD
        pA = core.parity_of(A); self.pA = pA
        if Dm is None:
            Dm = Index(monomials(d, pD))
        parsN = [(pA - pD - PARF[i]) % 2 for i in range(NV)]
        if Ncols is None:
            Ncols = [(i, m) for i in range(NV) for m in monomials(d + 1, parsN[i])]
        self.Dm, self.Ncols = Dm, Ncols
        self.nD, self.nN = Dm.n, len(Ncols)
        self.rows = Index(monomials(2 * d, pA))
        R = self.rows.n; self.R = R
        # bilinear table: B_m(e_{i,m'}) = (m'_i - m_i) z^{m+m'-e_i}
        Rtab = np.full((self.nD, self.nN), R, dtype=np.int32)   # R = trash row
        Ctab = np.zeros((self.nD, self.nN), dtype=np.float64)
        ridx = self.rows.idx
        for a, m in enumerate(Dm.mons):
            for b, (i, mp) in enumerate(Ncols):
                c = mp[i] - m[i]
                if c:
                    e = list(m)
                    for j in range(NV):
                        e[j] += mp[j]
                    e[i] -= 1
                    Rtab[a, b] = ridx[tuple(e)]
                    Ctab[a, b] = c
        self.Rtab, self.Ctab = Rtab, Ctab
        # reduction mod F (monomial X_{reduce_var}^6 -> rest)
        self.reduce_var = reduce_var
        self._build_reduction()
        if build_koszul:
            self._build_koszul()
        else:
            self.nc = 0
        self.build_time = time.time() - t0
        if verbose:
            print(f'[ReducedE2] d={d} pD={pD} nD={self.nD} nN={self.nN} rows={R} rows_red={self.Rred} nc={self.nc} '
                  f'nnzS={getattr(self, "nnzS", 0)} built in {self.build_time:.1f}s', file=sys.stderr, flush=True)

    # -- reduction modulo F ---------------------------------------------------
    def _build_reduction(self):
        F = to_float(self.F); v = self.reduce_var
        lead = tuple(6 if j == v else 0 for j in range(NV))
        assert abs(F.get(lead, 0.0) - 1.0) < 1e-12, 'leading coefficient of X_v^6 must be 1'
        rest = {e: -c for e, c in F.items() if e != lead}   # X_v^6 == rest (mod F)
        mons = self.rows.mons
        red_mons = [m for m in mons if m[v] < 6]
        self.red = Index(red_mons); self.Rred = self.red.n
        memo = {}

        def nf(m):
            if m[v] < 6:
                return {m: 1.0}
            if m in memo:
                return memo[m]
            base = list(m); base[v] -= 6; base = tuple(base)
            out = {}
            for e, c in rest.items():
                sub_ = nf(addexp(base, e))
                for k, cc in sub_.items():
                    out[k] = out.get(k, 0.0) + c * cc
            memo[m] = out
            return out
        ri, ci, vals = [], [], []
        for j, m in enumerate(mons):
            for k, c in nf(m).items():
                ri.append(self.red.idx[k]); ci.append(j); vals.append(c)
        self.RF = sp.csr_matrix((vals, (ri, ci)), shape=(self.Rred, self.R + 1))  # extra trash column (all zero)
        self.RFT = self.RF.T.tocsr()
        self.rf_maxabs = float(np.max(np.abs(vals)))

    # -- Koszul correction columns -------------------------------------------
    def _build_koszul(self):
        d = self.d
        Fi = [to_float(pderiv(self.F, i)) for i in range(NV)]
        Nidx = {c: j for j, c in enumerate(self.Ncols)}
        cols = []   # each: list of (nindex, value)
        if d - 4 >= 0:
            for i in range(NV):
                for j in range(i + 1, NV):
                    pc = (self.pA - self.pD - PARF[i] - PARF[j]) % 2
                    for m in monomials(d - 4, pc):
                        ent = {}
                        for e, c in Fi[j].items():
                            k = Nidx[(i, addexp(e, m))]; ent[k] = ent.get(k, 0.0) + c
                        for e, c in Fi[i].items():
                            k = Nidx[(j, addexp(e, m))]; ent[k] = ent.get(k, 0.0) - c
                        cols.append(ent)
        self.nc = len(cols)
        self.koszul_cols = cols
        # Kmat (nN x nc) sparse
        ri, ci, vals = [], [], []
        for c, ent in enumerate(cols):
            for k, v in ent.items():
                ri.append(k); ci.append(c); vals.append(v)
        self.Kmat = sp.csr_matrix((vals, (ri, ci)), shape=(self.nN, self.nc)) if self.nc else None
        # S: (nD, Rred*nc) sparse, so that M_red(a) = (D(a) @ S).reshape(Rred, nc)
        # entry for D-monomial m, koszul col c, N-entry (k, v): row RF-reduced ... we first build in the
        # unreduced row space then apply RF column-wise; cheaper: build T (nD, (R+1)*nc) then reduce.
        if self.nc == 0:
            self.S = None; self.nnzS = 0
            return
        # Build unreduced S0 (nD, (R+1)*nc) in chunks of columns, reduce each chunk by RF, concatenate.
        Rtab, Ctab = self.Rtab, self.Ctab
        nD, R, nc, Rred = self.nD, self.R, self.nc, self.Rred
        blocks = []
        chunk = max(1, min(nc, int(2e7 // max(1, nD * 8))))
        for c0 in range(0, nc, chunk):
            c1 = min(nc, c0 + chunk)
            ri, ci, vals = [], [], []
            for c in range(c0, c1):
                ent = cols[c]
                ks = np.fromiter(ent.keys(), dtype=np.int64, count=len(ent))
                vs = np.fromiter(ent.values(), dtype=np.float64, count=len(ent))
                rows_sub = Rtab[:, ks]            # nD x nk
                coef_sub = Ctab[:, ks] * vs[None, :]
                mask = coef_sub != 0.0
                aidx, kidx = np.nonzero(mask)
                ri.append(aidx)
                ci.append(rows_sub[aidx, kidx].astype(np.int64) * nc + c)
                vals.append(coef_sub[aidx, kidx])
            ri = np.concatenate(ri); ci = np.concatenate(ci); vals = np.concatenate(vals)
            S0 = sp.csr_matrix((vals, (ri, ci)), shape=(nD, (R + 1) * nc))
            blocks.append(S0)
        S0 = blocks[0] if len(blocks) == 1 else sum(blocks[1:], blocks[0])
        # reduce rows mod F: for each D-row a, reshape its (R+1) x nc block and apply RF.  Do it via kron:
        # S_red = S0 @ kron(RF^T, I_nc)  (shape (R+1)*nc x Rred*nc)
        Kron = sp.kron(self.RFT, sp.identity(nc, format='csr'), format='csr')
        self.S = (S0 @ Kron).tocsr()
        self.S.sum_duplicates()
        self.nnzS = self.S.nnz

    # -- evaluation ------------------------------------------------------------
    def B_vec(self, Dv, Nv):
        """B_D(N) in the unreduced row space (length R+1 incl. trash)."""
        w = (Dv[:, None] * Nv[None, :]) * self.Ctab
        return np.bincount(self.Rtab.ravel(), weights=w.ravel(), minlength=self.R + 1)

    def B_red(self, Dv, Nv):
        return self.RF @ self.B_vec(Dv, Nv)

    def adjoint_D(self, rt, Nv):
        """gradient wrt D-coefficients of  rt . B_D(N)   (rt in unreduced row space, length R+1)."""
        return (self.Ctab * rt[self.Rtab]) @ Nv

    def adjoint_N(self, rt, Dv):
        return Dv @ (self.Ctab * rt[self.Rtab])

    def M_red(self, Dv):
        """Reduced Koszul matrix T_D (Rred x nc)."""
        return np.asarray(self.S.T @ Dv).reshape(self.Rred, self.nc)

    def value_and_grad(self, u, fam, koszul=True, cond=1e-9, want_x=False, solver='lstsq'):
        """f(a) and gradient wrt u for a = u/|u| (fam.Dbasis columns)."""
        u = np.asarray(u, float); nu = np.linalg.norm(u); a = u / nu
        Dv = fam.D_vec(a); Nv = fam.N0_vec(a)
        b_full = -self.B_vec(Dv, Nv)                # unreduced, length R+1
        b = self.RF @ b_full
        if koszul and self.nc:
            M = self.M_red(Dv)
            if solver == 'lstsq':
                x, *_ = sla.lstsq(M, b, cond=cond, lapack_driver='gelsd', check_finite=False)
            else:
                G = M.T @ M; lam = cond * max(1.0, float(np.trace(G)) / max(1, self.nc))
                G[np.diag_indices_from(G)] += lam
                x = sla.solve(G, M.T @ b, assume_a='pos', check_finite=False)
            r = M @ x - b
            Nc = np.asarray(self.Kmat @ x).ravel()
        else:
            x = None; r = -b; Nc = np.zeros(self.nN)
        f = float(r @ r)
        # gradient: g_k = 2 r . (dM/da_k x - db/da_k) ;  M x - b = RF[ B_D(Nc) + B_D(N0) ] = RF B_D(Nc + N0)
        rt = self.RFT @ r                           # length R+1
        Ntot = Nc + Nv
        gD = self.adjoint_D(rt, Ntot)               # d/dD of rt.B_D(Ntot)   (nD)
        gN = self.adjoint_N(rt, Dv)                 # d/dN of rt.B_D(N) at N  (nN) -> acts on N0 part
        g = 2.0 * (fam.Dbasis.T @ gD + fam.N0.T @ gN)
        g_u = (g - float(g @ a) * a) / nu
        if want_x:
            return f, g_u, dict(a=a, x=x, r=r, bnorm=float(np.linalg.norm(b)), Dv=Dv)
        return f, g_u

    def value(self, a, fam, koszul=True, cond=1e-9):
        f, g, info = self.value_and_grad(a, fam, koszul=koszul, cond=cond, want_x=True)
        rel = f / max(info['bnorm'] ** 2, 1e-300)
        return f, rel, info

    def search(self, fam, n_restarts=10, seed=0, maxiter=300, gtol=1e-12, koszul=True, cond=1e-9,
               verbose=False, stop_rel=1e-13, init=None, callback=None, polish_rel=1e-6):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = np.asarray(init, float)
            else:
                u0 = rng.standard_normal(fam.K)
            u0 /= np.linalg.norm(u0)
            fun = lambda u: self.value_and_grad(u, fam, koszul=koszul, cond=cond)
            res = minimize(fun, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            f, rel, info = self.value(res.x, fam, koszul=koszul, cond=cond)
            if rel < polish_rel:
                # polish with more iterations and tighter tolerance
                res2 = minimize(fun, res.x / np.linalg.norm(res.x), jac=True, method='L-BFGS-B',
                                options={'maxiter': 4 * maxiter, 'gtol': 1e-16, 'ftol': 1e-20, 'maxcor': 30})
                f2, rel2, info2 = self.value(res2.x, fam, koszul=koszul, cond=cond)
                if rel2 < rel:
                    res, f, rel, info = res2, f2, rel2, info2
            cand = {'value': f, 'rel': rel, 'a': info['a'], 'iters': int(res.nit), 'restart': it, 'nfev': int(res.nfev),
                    'bnorm': info['bnorm'], 'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if verbose:
                print(f'[search {fam.tag} d={fam.d}] restart {it}: rel={rel:.3e} f={f:.3e} iters={res.nit} ({time.time()-t0:.1f}s)',
                      file=sys.stderr, flush=True)
            if callback is not None:
                callback(cand)
            if rel < stop_rel:
                break
        out.sort(key=lambda c: c['rel'])
        return out


# ----------------------------------------------------------------------------
# Exact checks
# ----------------------------------------------------------------------------
def rationalize_vec(v, max_den=1000, tol=1e-8):
    """Rational reconstruction of a float vector (normalized so that the entry of largest |.| is 1)."""
    v = np.asarray(v, float)
    j = int(np.argmax(np.abs(v)))
    v = v / v[j]
    out = []
    # try common denominators first
    for den in list(range(1, 200)) + [210, 240, 252, 280, 294, 300, 336, 343, 360, 420, 480, 490, 504, 540, 560, 588, 600, 630, 672, 686, 700, 720, 735, 756, 784, 840, 882, 900, 960, 980, 1000, 1008, 1029, 1050, 1080, 1120, 1176, 1200, 1260, 1344, 1372, 1400, 1440, 1470, 1512, 1568, 1680, 1764, 1800, 1960, 2000, 2016, 2058, 2100, 2160, 2205, 2240, 2352, 2401, 2520, 2688, 2744, 2800, 2880, 2940, 3024, 3087, 3136, 3360, 3430, 3528, 3600, 3920, 4116, 4200, 4320, 4410, 4480, 4704, 4802, 5040, 5488, 5880, 6174, 6860, 7056, 7203, 7350, 8232, 8820, 9604, 10080, 10290, 12348, 14406, 16807]:
        if den > max_den * 20:
            break
        w = v * den
        if np.max(np.abs(w - np.round(w))) < tol * den:
            out.append(tuple(Fraction(int(round(x)), den) for x in w))
            break
    fr = tuple(Fraction(float(x)).limit_denominator(max_den) for x in v)
    if max(abs(float(f) - x) for f, x in zip(fr, v)) < tol:
        if fr not in out:
            out.append(fr)
    return out


def exact_D_from_frac(Dm, frac):
    return {m: Fraction(c) for m, c in zip(Dm.mons, frac) if c != 0}


def exact_check_D(F, A, q, D, verbose=False, want_exact=True):
    """E1 (mod p) then the full fixed-D system (compressed mod p + exact lift).  Returns dict."""
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    res = {'e1': bool(ok), 'e1_info': {k: v for k, v in info.items() if k in ('rank', 'cols', 'rows')}}
    if not ok:
        res.update(consistent=False, seconds=round(time.time() - t0, 1))
        return res
    out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose, want_exact=want_exact)
    res.update(consistent=bool(out['consistent']), rank=int(out['rank']), nullity=int(out['nullity']),
               rows=out.get('rows'), cols=out.get('cols'), seconds=round(time.time() - t0, 1))
    if out.get('consistent') and out.get('exact'):
        res['N'] = [core.to_str(n) for n in out['N']]; res['Q'] = core.to_str(out['Q'])
        res['exact'] = True
    return res


def e2_modp_check(F, A, D, p=None, verbose=False):
    """Exact normalized E2 consistency test mod p for a rational D (E1 lift via nullspace, Koszul + F*C columns).
    Returns dict(consistent, rank, rank_aug, rows, cols)."""
    import flint
    from flint import nmod_mat
    p = p or core.PRIMES[0]
    t0 = time.time()
    q = 1
    d = core.degree(D); pD = core.parity_of(D); pA = core.parity_of(A)
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    # E1: solve sum N_i F_i = A D
    parsN = [(pA - pD - PARF[i]) % 2 for i in range(NV)]
    Ncols = [(i, m) for i in range(NV) for m in monomials(d + 1, parsN[i])]
    rows1 = Index(monomials(d + 6, (pA + pD) % 2))
    flat = [0] * (rows1.n * (len(Ncols) + 1)); W = len(Ncols) + 1
    for j, (i, m) in enumerate(Ncols):
        for e, c in Fi[i].items():
            flat[rows1.idx[addexp(e, m)] * W + j] = c
    AD = mul(Ap, Dp, p)
    for e, c in AD.items():
        flat[rows1.idx[e] * W + len(Ncols)] = c
    M1 = nmod_mat(rows1.n, W, flat, p)
    Mr, rk = M1.rref()
    pivots = []
    for i in range(rk):
        j = 0
        while j < W and int(Mr[i, j]) == 0:
            j += 1
        pivots.append(j)
    if pivots and pivots[-1] == W - 1:
        return {'e1': False, 'seconds': round(time.time() - t0, 1)}
    N = [{} for _ in range(NV)]
    for i, j in enumerate(pivots):
        c = int(Mr[i, W - 1])
        if c:
            ii, m = Ncols[j]; N[ii][m] = c
    # B_D(N)
    Di = [pderiv(Dp, i, p) for i in range(NV)]
    B = {}
    for i in range(NV):
        B = add(B, mul(Dp, pderiv(N[i], i, p), p), p)
        B = sub(B, mul(N[i], Di[i], p), p)
    rows2 = Index(monomials(2 * d, pA))
    cols = []
    if d - 4 >= 0:
        for i in range(NV):
            for j in range(i + 1, NV):
                pc = (pA - pD - PARF[i] - PARF[j]) % 2
                for m in monomials(d - 4, pc):
                    cm = shift({(0,) * NV: 1}, m, 1, p)
                    Ki = mul(cm, Fi[j], p); Kj = scale(mul(cm, Fi[i], p), p - 1, p)
                    # B_D(K) = D (d_i Ki + d_j Kj) - Ki D_i - Kj D_j
                    col = add(mul(Dp, add(pderiv(Ki, i, p), pderiv(Kj, j, p), p), p),
                              scale(add(mul(Ki, Di[i], p), mul(Kj, Di[j], p), p), p - 1, p), p)
                    cols.append(col)
    nc = len(cols)
    for m in monomials(2 * d - 6, pA) if 2 * d - 6 >= 0 else []:
        cols.append(scale(shift(Fp, m, 1, p), p - 1, p))
    W2 = len(cols) + 1
    flat = [0] * (rows2.n * W2)
    for j, col in enumerate(cols):
        for e, c in col.items():
            flat[rows2.idx[e] * W2 + j] = c
    for e, c in B.items():
        flat[rows2.idx[e] * W2 + len(cols)] = (p - c) % p
    M2 = nmod_mat(rows2.n, W2, flat, p)
    rank_aug = M2.rank()
    M2b = nmod_mat(rows2.n, W2 - 1, [flat[r * W2 + j] for r in range(rows2.n) for j in range(W2 - 1)], p)
    rank = M2b.rank()
    return {'e1': True, 'consistent': rank == rank_aug, 'rank': rank, 'rank_aug': rank_aug, 'rows': rows2.n,
            'cols': W2 - 1, 'nc': nc, 'p': p, 'seconds': round(time.time() - t0, 1)}
