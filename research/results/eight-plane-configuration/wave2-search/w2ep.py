"""w2ep -- wave-2 search library for the eight-plane-configuration lane (q = 1).

Koszul-reduced E2 variable projection.  For a linear family D(a) = sum_k a_k D_k with EXACT E1 lifts
N^(k) (sum_i N^(k)_i F_i = A D_k as polynomials) the certificate exists iff

      B_D(N) + T_D(c) = F Q,       B_D(N) = D div N - N(D),
      T_D(c) = sum_{i<j} [ D (F_j d_i c_ij - F_i d_j c_ij) - c_ij (D_i F_j - D_j F_i) ],   deg c_ij = d - 4,

(docs/05-computation.md, "full correction map"; the F-multiple freedom G z/6 in N drops out of E2 exactly).
Rows live in S_{2d}^{par A}; Q is eliminated by exact reduction modulo F (pivot X0^6, triangular), so the
per-evaluation least-squares problem is (|S_2d| - |S_{2d-6}|) x nc with nc = sum_{i<j} |S_{d-4}^{par c_ij}|.

    f(a) = min_c || Red( T_{D(a)} c + B_{D(a)}(N(a)) ) ||^2,     a on the unit sphere, F-multiples projected out.

Gradient by the variable-projection formula (exact at the least-squares optimum), implemented with adjoints of
product tables (no per-basis matrices are stored, so K ~ 1000 is fine).

Families: jacobian_family(d, parD) (D = w(F), N = A w) and kernel_family(d, parD) (float E1 kernel + lifts from
kernel_float.py).  HessianRelax: the necessary condition A (w^T Hess F w) in J, i.e. w^T H w in Ker_{2d-6}^+.

Engine: /home/kejace/devel/sextics-work/hodge (not modified).
"""
from __future__ import annotations
import sys, os, time, json, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power

LANE = 'eight-plane-configuration'
WD = f'/tmp/hodge/results/{LANE}'
SD_DIR = f'{WD}/wave2-search'
LOG = f'{WD}/wave2-search.log.jsonl'
BK = 64
POW = np.array([BK ** i for i in range(6)], dtype=np.int64)
VARPAR = [0, 0, 0, 1, 1, 1]


def load():
    F, A, q, rec = core.load_instance(LANE)
    return F, A, q


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%S')


def log_batch(entry):
    entry = dict(entry); entry.setdefault('ts', now()); entry.setdefault('lane', LANE); entry.setdefault('engine', 2)
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry) + '\n')


# ----------------------------------------------------------------------------------------------------------
# monomial spaces and packed keys
# ----------------------------------------------------------------------------------------------------------
def keys_of(mons):
    return (np.array(mons, dtype=np.int64) @ POW) if len(mons) else np.zeros(0, np.int64)


def unkey(k):
    e = []
    for _ in range(6):
        e.append(int(k % BK)); k //= BK
    return tuple(e)


def poly_arrays(P):
    ks = np.array([int(np.dot(e, POW)) for e in P], dtype=np.int64)
    vs = np.array([float(v) for v in P.values()], dtype=float)
    return ks, vs


class MonoSpace:
    def __init__(self, deg, par):
        self.deg, self.par = deg, par
        self.mons = monomials(deg, par) if deg >= 0 else []
        self.n = len(self.mons)
        self.exps = np.array(self.mons, dtype=np.int64).reshape(self.n, 6)
        self.keys = keys_of(self.mons)
        order = np.argsort(self.keys)
        self.sorted_keys = self.keys[order]; self.order = order

    def index(self, ks):
        ks = np.asarray(ks, dtype=np.int64)
        pos = np.searchsorted(self.sorted_keys, ks)
        pos = np.minimum(pos, max(self.n - 1, 0))
        if self.n == 0 or not np.all(self.sorted_keys[pos] == ks):
            raise RuntimeError(f'monomial outside space S_{self.deg}^{self.par}')
        return self.order[pos]

    def vec(self, P):
        v = np.zeros(self.n)
        if P:
            ks, vs = poly_arrays(P)
            v[self.index(ks)] = vs
        return v

    def poly(self, v, tol=0.0):
        return {self.mons[j]: float(v[j]) for j in range(self.n) if abs(v[j]) > tol}

    def __repr__(self):
        return f'S_{self.deg}^{self.par}({self.n})'


_SPACES = {}


def space(deg, par):
    k = (deg, par % 2)
    if k not in _SPACES:
        _SPACES[k] = MonoSpace(deg, par % 2)
    return _SPACES[k]


def mult_matrix(Pk, Pv, src: MonoSpace, dst: MonoSpace):
    """Sparse matrix (dst.n x src.n) of multiplication by the polynomial (Pk, Pv)."""
    if len(Pk) == 0 or src.n == 0:
        return sp.csr_matrix((dst.n, src.n))
    kk = (Pk[:, None] + src.keys[None, :]).ravel()
    rows = dst.index(kk)
    cols = np.broadcast_to(np.arange(src.n, dtype=np.int64)[None, :], (len(Pk), src.n)).ravel()
    vals = np.broadcast_to(Pv[:, None], (len(Pk), src.n)).ravel()
    return sp.csr_matrix((vals, (rows, cols)), shape=(dst.n, src.n))


def deriv_matrix(i, src: MonoSpace, dst: MonoSpace):
    mi = src.exps[:, i]
    sel = np.nonzero(mi)[0]
    if sel.size == 0:
        return sp.csr_matrix((dst.n, src.n))
    rows = dst.index(src.keys[sel] - BK ** i)
    return sp.csr_matrix((mi[sel].astype(float), (rows, sel)), shape=(dst.n, src.n))


class ProductTable:
    """All pairwise products src1 x src2 -> dst as index arrays (I, J, R)."""

    def __init__(self, s1: MonoSpace, s2: MonoSpace, dst: MonoSpace):
        self.s1, self.s2, self.dst = s1, s2, dst
        if s1.n == 0 or s2.n == 0:
            self.I = self.J = self.R = np.zeros(0, np.int64); return
        I, J = np.meshgrid(np.arange(s1.n), np.arange(s2.n), indexing='ij')
        self.I = I.ravel().astype(np.int32); self.J = J.ravel().astype(np.int32)
        self.R = dst.index(s1.keys[self.I] + s2.keys[self.J]).astype(np.int32)

    def apply(self, u, v):
        if self.I.size == 0:
            return np.zeros(self.dst.n)
        return np.bincount(self.R, weights=u[self.I] * v[self.J], minlength=self.dst.n)

    def adj_u(self, rho, v):
        if self.I.size == 0:
            return np.zeros(self.s1.n)
        return np.bincount(self.I, weights=rho[self.R] * v[self.J], minlength=self.s1.n)

    def adj_v(self, rho, u):
        if self.I.size == 0:
            return np.zeros(self.s2.n)
        return np.bincount(self.J, weights=rho[self.R] * u[self.I], minlength=self.s2.n)


# ----------------------------------------------------------------------------------------------------------
# exact reduction modulo F (pivot monomial X0^6)
# ----------------------------------------------------------------------------------------------------------
class ModFReducer:
    def __init__(self, F, SR: MonoSpace, SQ: MonoSpace):
        self.SR, self.SQ = SR, SQ
        nQ, R = SQ.n, SR.n
        self.nQ = nQ
        if nQ == 0:
            self.piv = np.zeros(0, np.int64); self.nonpiv = np.arange(R); self.P1inv = None; self.P2 = None
            self.nred = R
            return
        Fk, Fv = poly_arrays(F)
        c0 = float(F[(6, 0, 0, 0, 0, 0)])
        # order Q-monomials by X0-degree descending -> P1 lower triangular
        order = np.argsort(-SQ.exps[:, 0], kind='stable')
        self.qorder = order
        qkeys = SQ.keys[order]
        piv = SR.index(qkeys + 6)  # X0^6 * m  (key of X0 is 1)
        self.piv = piv
        mask = np.ones(R, bool); mask[piv] = False
        self.nonpiv = np.nonzero(mask)[0]
        self.nred = R - nQ
        MQ = mult_matrix(Fk, Fv, MonoSpaceView(qkeys), SR).tocsr()
        P1 = MQ[piv, :].toarray()
        # sanity: lower triangular with diagonal c0
        assert np.allclose(np.diag(P1), c0), 'pivot block diagonal mismatch'
        assert np.allclose(np.triu(P1, 1), 0.0), 'pivot block not lower triangular'
        self.P1inv = sla.solve_triangular(P1, np.eye(nQ), lower=True, check_finite=False)
        self.P2 = MQ[self.nonpiv, :].tocsr()
        self.P2T = self.P2.T.tocsr()

    def reduce(self, V):
        """V: dense (R x k) or (R,) -> reduced (nred x k)."""
        if self.nQ == 0:
            return V
        one = (V.ndim == 1)
        if one:
            V = V[:, None]
        V1 = V[self.piv, :]; V2 = V[self.nonpiv, :]
        out = V2 - self.P2 @ (self.P1inv @ V1)
        return out[:, 0] if one else out

    def reduce_sparse(self, M):
        """M: sparse (R x k) -> dense reduced (nred x k)."""
        if self.nQ == 0:
            return M.toarray()
        M = M.tocsr()
        V1 = M[self.piv, :].toarray(); V2 = M[self.nonpiv, :].toarray()
        return V2 - self.P2 @ (self.P1inv @ V1)

    def adjoint(self, r):
        """Red^T r: vector in the full row space."""
        out = np.zeros(self.SR.n)
        if self.nQ == 0:
            return r.copy()
        out[self.nonpiv] = r
        out[self.piv] = -(self.P1inv.T @ (self.P2T @ r))
        return out

    def quotient(self, V1):
        """Q coefficients (in qorder) for the pivot part V1 = V[piv]."""
        return self.P1inv @ V1


class MonoSpaceView:
    """Minimal src-space stand-in for mult_matrix with a given key order."""

    def __init__(self, keys):
        self.keys = np.asarray(keys, np.int64); self.n = len(keys)


# ----------------------------------------------------------------------------------------------------------
# families
# ----------------------------------------------------------------------------------------------------------
def numerator_parities(parD, parA):
    return [(parA + parD + VARPAR[i]) % 2 for i in range(NV)]


def jacobian_family(F, A, d, parD):
    """w-basis {z^m e_i}: D_k = z^m F_i, N^(k) = A z^m e_i.  Returns dict with Dmat, Nmats, wdesc."""
    parA = core.parity_of(A)
    SD = space(d, parD)
    pars = numerator_parities(parD, parA)
    SN = [space(d + 1, pars[i]) for i in range(NV)]
    Fi = [pderiv(F, i) for i in range(NV)]
    wdesc = []
    for i in range(NV):
        for m in monomials(d - 5, (parD + VARPAR[i]) % 2):
            wdesc.append((i, m))
    K = len(wdesc)
    Dmat = np.zeros((K, SD.n))
    Nmats = [np.zeros((K, SN[i].n)) for i in range(NV)]
    for k, (i, m) in enumerate(wdesc):
        Dmat[k] = SD.vec(shift(Fi[i], m, 1))
        Nmats[i][k] = SN[i].vec(shift(A, m, 1))
    return {'d': d, 'parD': parD, 'Dmat': Dmat, 'Nmats': Nmats, 'wdesc': wdesc, 'kind': 'jacobian'}


def kernel_family(d, parD, path=None):
    path = path or f'{SD_DIR}/kernel_float_d{d}_par{parD}.npz'
    z = np.load(path)
    Nmats = [z[f'N{i}'] for i in range(NV)]
    return {'d': d, 'parD': parD, 'Dmat': z['D'], 'Nmats': Nmats, 'kind': 'kernel', 'path': path}


def subfamily(fam, coeffs):
    """Linear subfamily: rows of coeffs (k x K) combine the family's basis."""
    C = np.asarray(coeffs, float)
    return {**fam, 'Dmat': C @ fam['Dmat'], 'Nmats': [C @ M for M in fam['Nmats']], 'kind': fam['kind'] + '-sub'}


# ----------------------------------------------------------------------------------------------------------
# reduced E2 variable projection
# ----------------------------------------------------------------------------------------------------------
class ReducedE2:
    def __init__(self, F, A, fam, cond=1e-9, verbose=False, exclude_F_multiples=True, dtype=np.float64):
        t0 = time.time()
        self.F, self.A = F, A
        self.q = 1
        d = fam['d']; parD = fam['parD']; self.d, self.parD = d, parD
        parA = core.parity_of(A); self.parA = parA
        self.cond = cond
        self.Dmat = np.asarray(fam['Dmat'], float); self.Nmats = [np.asarray(M, float) for M in fam['Nmats']]
        K = self.Dmat.shape[0]; self.K = K
        self.SD = space(d, parD)
        self.parsN = numerator_parities(parD, parA)
        self.SN = [space(d + 1, self.parsN[i]) for i in range(NV)]
        self.SdivN = space(d, (parA + parD) % 2)
        self.SdD = [space(d - 1, (parD + VARPAR[i]) % 2) for i in range(NV)]
        self.SR = space(2 * d, parA)
        self.SQ = space(2 * d - 6, parA)
        # derivative operators
        self.DerN = [deriv_matrix(i, self.SN[i], self.SdivN) for i in range(NV)]
        self.DerD = [deriv_matrix(i, self.SD, self.SdD[i]) for i in range(NV)]
        # product tables for B_D(N)
        self.PT_DdivN = ProductTable(self.SD, self.SdivN, self.SR)
        self.PT_NdD = [ProductTable(self.SN[i], self.SdD[i], self.SR) for i in range(NV)]
        # Koszul data
        Fi = [pderiv(F, i) for i in range(NV)]
        Fk = [poly_arrays(Fi[i]) for i in range(NV)]
        self.pairs = [(i, j) for i in range(NV) for j in range(i + 1, NV)]
        self.cspaces, self.Phi, self.Psi, self.PsiP, self.Sphi, self.Spsi, self.SpsiP = [], [], [], [], [], [], []
        ncol = 0; self.cslices = []
        for (i, j) in self.pairs:
            parc = (self.parsN[i] + VARPAR[j]) % 2
            Sc = space(d - 4, parc)
            Sdc_i = space(d - 5, parc)  # d_i c has same parity as c for X-derivatives? no: parity flips for Y
            Sdc_i = space(d - 5, (parc + VARPAR[i]) % 2)
            Sdc_j = space(d - 5, (parc + VARPAR[j]) % 2)
            Sphi = space(d, (parc + VARPAR[i] + VARPAR[j]) % 2)
            Spsi = space(d + 1, (parc + VARPAR[j]) % 2)    # F_j c
            SpsiP = space(d + 1, (parc + VARPAR[i]) % 2)   # F_i c
            Phi = (mult_matrix(Fk[j][0], Fk[j][1], Sdc_i, Sphi) @ deriv_matrix(i, Sc, Sdc_i)
                   - mult_matrix(Fk[i][0], Fk[i][1], Sdc_j, Sphi) @ deriv_matrix(j, Sc, Sdc_j)).tocsc()
            Psi = mult_matrix(Fk[j][0], Fk[j][1], Sc, Spsi).tocsc()
            PsiP = mult_matrix(Fk[i][0], Fk[i][1], Sc, SpsiP).tocsc()
            self.cspaces.append(Sc); self.Phi.append(Phi.toarray()); self.Psi.append(Psi.toarray()); self.PsiP.append(PsiP.toarray())
            self.Sphi.append(Sphi); self.Spsi.append(Spsi); self.SpsiP.append(SpsiP)
            self.cslices.append(slice(ncol, ncol + Sc.n)); ncol += Sc.n
        self.nc = ncol
        # product tables for D * (.) and dD_i * (.) onto rows (keyed by source space)
        self._pt_cache = {}
        # mod-F reducer
        self.red = ModFReducer(F, self.SR, self.SQ)
        self.nred = self.red.nred
        # F-multiples (spurious) subspace in a-coordinates
        self.P = np.eye(K); self.excluded_dim = 0
        if exclude_F_multiples and d >= 6:
            Gm = monomials(d - 6, (parD - core.parity_of(F)) % 2)
            Fcols = np.array([self.SD.vec(shift(F, m, 1)) for m in Gm]) if Gm else np.zeros((0, self.SD.n))
            Mb = np.vstack([self.Dmat, -Fcols]) if len(Gm) else self.Dmat
            ns = sla.null_space(Mb.T)   # vectors (a, g) with Dmat^T a = Fcols^T g
            S = ns[:K, :] if ns.size else np.zeros((K, 0))
            if S.shape[1]:
                U, s, _ = np.linalg.svd(S, full_matrices=True)
                r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
                self.P = U[:, r:]; self.excluded_dim = r
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        self._nfev = 0
        self.verbose = verbose
        if verbose:
            print(f'[reducedE2] d={d} parD={parD} K={K} dimP={self.dimP} rows={self.SR.n} reduced={self.nred} '
                  f'nc={self.nc} nQ={self.SQ.n} built in {self.build_time:.1f}s', file=sys.stderr)

    # --- polynomial-level pieces -----------------------------------------------------------------------
    def _pt(self, src: MonoSpace):
        k = (src.deg, src.par)
        if k not in self._pt_cache:
            self._pt_cache[k] = ProductTable(self.SD, src, self.SR) if src.deg == self.d else None
        return self._pt_cache[k]

    def D_of(self, a):
        return self.Dmat.T @ a

    def N_of(self, a):
        return [M.T @ a for M in self.Nmats]

    def b_vec(self, Dv, Nv):
        """-B_D(N) = -(D div N - sum N_i dD_i)  in SR coordinates."""
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.DerN[i] @ Nv[i]
        B = self.PT_DdivN.apply(Dv, divN)
        for i in range(NV):
            dDi = self.DerD[i] @ Dv
            B -= self.PT_NdD[i].apply(Nv[i], dDi)
        return -B

    def T_matrix(self, Dv):
        """Dense reduced matrix Red(T_D) of shape (nred x nc)."""
        Dk = self.SD.keys; nz = np.nonzero(Dv)[0]
        Dk, Dvv = Dk[nz], Dv[nz]
        dD = []
        for i in range(NV):
            v = self.DerD[i] @ Dv
            nzi = np.nonzero(v)[0]
            dD.append((self.SdD[i].keys[nzi], v[nzi]))
        T = np.zeros((self.SR.n, self.nc))
        for p, (i, j) in enumerate(self.pairs):
            sl = self.cslices[p]
            MD = mult_matrix(Dk, Dvv, self.Sphi[p], self.SR)
            T[:, sl] += MD @ self.Phi[p]
            Mi = mult_matrix(dD[i][0], dD[i][1], self.Spsi[p], self.SR)
            T[:, sl] -= Mi @ self.Psi[p]
            Mj = mult_matrix(dD[j][0], dD[j][1], self.SpsiP[p], self.SR)
            T[:, sl] += Mj @ self.PsiP[p]
        return self.red.reduce(T)

    def T_adjoint_D(self, rho, c):
        """Gradient wrt the D coefficient vector of rho . T_D(c) (rho in full row space)."""
        g = np.zeros(self.SD.n)
        for p, (i, j) in enumerate(self.pairs):
            cp = c[self.cslices[p]]
            phi = self.Phi[p] @ cp; psi = self.Psi[p] @ cp; psiP = self.PsiP[p] @ cp
            # rho . (D * phi) -> adjoint on D
            nz = np.nonzero(phi)[0]
            if nz.size:
                g += mult_matrix(self.Sphi[p].keys[nz], phi[nz], self.SD, self.SR).T @ rho
            nz = np.nonzero(psi)[0]
            if nz.size:
                g -= self.DerD[i].T @ (mult_matrix(self.Spsi[p].keys[nz], psi[nz], self.SdD[i], self.SR).T @ rho)
            nz = np.nonzero(psiP)[0]
            if nz.size:
                g += self.DerD[j].T @ (mult_matrix(self.SpsiP[p].keys[nz], psiP[nz], self.SdD[j], self.SR).T @ rho)
        return g

    def b_adjoint(self, rho, Dv, Nv):
        """Gradients wrt D and N_i of rho . b(D, N), b = -(D div N - sum N_i dD_i)."""
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.DerN[i] @ Nv[i]
        gD = -self.PT_DdivN.adj_u(rho, divN)
        gdivN = -self.PT_DdivN.adj_v(rho, Dv)
        gN = [self.DerN[i].T @ gdivN for i in range(NV)]
        for i in range(NV):
            dDi = self.DerD[i] @ Dv
            gN[i] += self.PT_NdD[i].adj_u(rho, dDi)
            gD += self.DerD[i].T @ self.PT_NdD[i].adj_v(rho, Nv[i])
        return gD, gN

    # --- objective --------------------------------------------------------------------------------------
    def _point(self, u):
        a_full = self.P @ np.asarray(u, float)
        na = np.linalg.norm(a_full)
        return a_full / na, na

    def value_full(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        Dv = self.D_of(a); Nv = self.N_of(a)
        b = self.red.reduce(self.b_vec(Dv, Nv))
        T = self.T_matrix(Dv)
        c, *_ = sla.lstsq(T, b, cond=self.cond, lapack_driver='gelsy', check_finite=False)
        r = T @ c - b
        self._nfev += 1
        return {'f': float(r @ r), 'bn': float(np.linalg.norm(b)), 'rel': float(np.linalg.norm(r) / max(np.linalg.norm(b), 1e-300)),
                'c': c, 'r': r, 'Dv': Dv, 'Nv': Nv, 'T': T, 'b': b, 'a': a}

    def value_and_grad(self, u):
        a, na = self._point(u)
        Dv = self.D_of(a); Nv = self.N_of(a)
        bfull = self.b_vec(Dv, Nv)
        b = self.red.reduce(bfull)
        T = self.T_matrix(Dv)
        c, *_ = sla.lstsq(T, b, cond=self.cond, lapack_driver='gelsy', check_finite=False)
        r = T @ c - b
        f = float(r @ r)
        rho = self.red.adjoint(2.0 * r)
        gD = self.T_adjoint_D(rho, c)
        gbD, gbN = self.b_adjoint(rho, Dv, Nv)
        gD = gD - gbD
        g = self.Dmat @ gD
        for i in range(NV):
            g -= self.Nmats[i] @ gbN[i]
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def fd_check(self, seed=0, h=1e-6):
        rng = np.random.default_rng(seed)
        u = rng.standard_normal(self.dimP); u /= np.linalg.norm(u)
        f0, g = self.value_and_grad(u)
        v = rng.standard_normal(self.dimP); v -= (v @ u) * u; v /= np.linalg.norm(v)
        fp, _ = self.value_and_grad(u + h * v); fm, _ = self.value_and_grad(u - h * v)
        return float(g @ v), (fp - fm) / (2 * h), f0

    def is_spurious(self, a, tol=1e-8):
        """Is D(a) (numerically) divisible by F?  Uses the excluded subspace complement."""
        if self.excluded_dim == 0:
            return False, 1.0
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        comp = self.P @ (self.P.T @ a)
        rel = float(np.linalg.norm(comp))      # fraction of a inside the allowed complement
        return rel < tol, 1.0 - rel

    # --- search -----------------------------------------------------------------------------------------
    def minimize_from(self, seed, maxiter=300, gtol=1e-10, init=None):
        rng = np.random.default_rng(seed)
        if init is not None:
            u0 = self.P.T @ np.asarray(init, float)
        else:
            u0 = rng.standard_normal(self.dimP)
        u0 /= np.linalg.norm(u0)
        t0 = time.time()
        res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                       options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-16})
        a, _ = self._point(res.x)
        out = self.value_full(a)
        spur, srel = self.is_spurious(a)
        return {'seed': seed, 'f': out['f'], 'rel': out['rel'], 'bn': out['bn'], 'iters': int(res.nit), 'nfev': int(res.nfev),
                'spurious': bool(spur), 'seconds': round(time.time() - t0, 2), 'a': a, 'c': out['c'], 'Dv': out['Dv'],
                'msg': str(res.message)[:60]}

    def rank_stats(self, a):
        out = self.value_full(a)
        s = np.linalg.svd(out['T'], compute_uv=False)
        rk = int(np.sum(s > self.cond * s[0]))
        return {'rows': self.nred, 'cols': self.nc, 'rank': rk, 'nullity': self.nc - rk, 'smax': float(s[0]), 'smin_nonzero': float(s[rk - 1]) if rk else 0.0}

    # --- reconstruction and exact confirmation ---------------------------------------------------------
    def D_poly_float(self, a):
        return self.SD.poly(self.D_of(np.asarray(a, float) / np.linalg.norm(a)))

    def reconstruct_vec(self, v, max_den=1000, tol=1e-7):
        """Rational candidates for the vector v (normalized by its largest entry)."""
        v = np.asarray(v, float)
        j = int(np.argmax(np.abs(v)))
        v = v / v[j]
        cands = []
        fr = [Fraction(float(x)).limit_denominator(max_den) for x in v]
        if max(abs(float(f) - x) for f, x in zip(fr, v)) < tol:
            cands.append(tuple(fr))
        for den in range(1, max_den + 1):
            w = v * den
            if np.max(np.abs(w - np.round(w))) < tol * den:
                t = tuple(Fraction(int(round(x)), den) for x in w)
                if t not in cands:
                    cands.append(t)
                break
        return cands

    def D_exact_from_a(self, frac_a, exact_basis=None):
        """Exact D from rational family coordinates (exact_basis: list of exact polys D_k)."""
        D = {}
        for k, c in enumerate(frac_a):
            if c:
                D = add(D, scale(exact_basis[k], Fraction(c)))
        return D

    def D_exact_from_vec(self, frac_v):
        return {m: Fraction(c) for m, c in zip(self.SD.mons, frac_v) if c}

    def confirm(self, D, verbose=False):
        if not D:
            return 'inconsistent', {'note': 'D = 0'}
        ok, info = fixed_d.e1_solvable(self.F, self.A, 1, D)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails', 'e1': {k: int(v) for k, v in info.items()}}
        try:
            out = fixed_d.solve_fixed_D_fast(self.F, self.A, 1, D, verbose=verbose)
            if not out['consistent']:
                return 'inconsistent', {'rank': out['rank'], 'nullity': out['nullity']}
            if out.get('exact'):
                return 'solved', {'N': out['N'], 'Q': out['Q'], 'rank': out['rank'], 'nullity': out['nullity']}
        except Exception as ex:
            if verbose:
                print(f'[reducedE2] fast solver failed ({ex}); dense fallback', file=sys.stderr)
        status, data = fixed_d.solve_fixed_D(self.F, self.A, 1, D, verbose=verbose)
        return status, data


# ----------------------------------------------------------------------------------------------------------
# Hessian relaxation (Jacobian-type degrees): w^T Hess(F) w in Ker_{2d-6}^+
# ----------------------------------------------------------------------------------------------------------
class HessianRelax:
    def __init__(self, F, A, d, parD, ker_basis_vecs=None, verbose=False):
        """ker_basis_vecs: array (k x |S_{2d-6}^+|) spanning the E1 kernel of degree 2d-6 (even); None -> J only."""
        t0 = time.time()
        self.F, self.A, self.d, self.parD = F, A, d, parD
        self.fam = jacobian_family(F, A, d, parD)
        self.wdesc = self.fam['wdesc']; K = len(self.wdesc); self.K = K
        self.Sw = [space(d - 5, (parD + VARPAR[i]) % 2) for i in range(NV)]
        self.wslices = []; self.widx = [[] for _ in range(NV)]
        for k, (i, m) in enumerate(self.wdesc):
            self.widx[i].append(k)
        self.Wsel = [np.array(self.widx[i], dtype=np.int64) for i in range(NV)]
        # component coefficient vectors: w_i = a[Wsel[i]] in the order of monomials(d-5, par) -- jacobian_family uses that order
        self.ST = space(2 * d - 6, 0)
        self.Sprod = {}
        Fi = [pderiv(F, i) for i in range(NV)]
        self.Hk = {}
        self.PT = {}
        for i in range(NV):
            for j in range(NV):
                Fij = pderiv(Fi[i], j)
                Sp = space(2 * d - 10, (self.Sw[i].par + self.Sw[j].par) % 2)
                self.PT[(i, j)] = ProductTable(self.Sw[i], self.Sw[j], Sp)
                ks, vs = poly_arrays(Fij)
                self.Hk[(i, j)] = mult_matrix(ks, vs, Sp, self.ST).tocsr()
        # kernel subspace of S_{2d-6}^+: J generators (+ extra kernel vectors)
        gens = []
        for i in range(NV):
            for m in monomials(2 * d - 11, VARPAR[i]):
                gens.append(self.ST.vec(shift(Fi[i], m, 1)))
        G = np.array(gens) if gens else np.zeros((0, self.ST.n))
        if ker_basis_vecs is not None and len(ker_basis_vecs):
            G = np.vstack([G, np.asarray(ker_basis_vecs, float)])
        if G.shape[0]:
            U, s, Vt = np.linalg.svd(G, full_matrices=False)
            r = int(np.sum(s > 1e-10 * s[0]))
            self.U = Vt[:r, :].T       # ST.n x r orthonormal basis of the kernel subspace
            self.kdim = r
        else:
            self.U = np.zeros((self.ST.n, 0)); self.kdim = 0
        # spurious (Euler cone) complement
        Gm = monomials(d - 6, parD)
        E = np.zeros((len(Gm), K))
        for g, m in enumerate(Gm):
            for i in range(NV):
                mm = tuple(m[t] + (1 if t == i else 0) for t in range(6))
                k = self.wdesc.index((i, mm))
                E[g, k] = 1.0
        if len(Gm):
            U, s, Vt = np.linalg.svd(E, full_matrices=True)
            r = int(np.sum(s > 1e-10 * s[0]))
            self.P = Vt[r:, :].T
        else:
            self.P = np.eye(K)
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        if verbose:
            print(f'[hessrelax] d={d} parD={parD} K={K} dimP={self.dimP} target S_{2*d-6}^+ dim {self.ST.n} kernel dim {self.kdim} '
                  f'-> {self.ST.n - self.kdim} quadrics; built {self.build_time:.1f}s', file=sys.stderr)

    def quad(self, a):
        """w^T H w as a vector in ST."""
        w = [a[self.Wsel[i]] for i in range(NV)]
        v = np.zeros(self.ST.n)
        for i in range(NV):
            for j in range(NV):
                v += self.Hk[(i, j)] @ self.PT[(i, j)].apply(w[i], w[j])
        return v

    def value_and_grad(self, u):
        a_full = self.P @ u; na = np.linalg.norm(a_full); a = a_full / na
        v = self.quad(a)
        res = v - self.U @ (self.U.T @ v)
        f = float(res @ res)
        rho = 2.0 * res
        w = [a[self.Wsel[i]] for i in range(NV)]
        g = np.zeros(self.K)
        for i in range(NV):
            for j in range(NV):
                rr = self.Hk[(i, j)].T @ rho
                g[self.Wsel[i]] += self.PT[(i, j)].adj_u(rr, w[j])
                g[self.Wsel[j]] += self.PT[(i, j)].adj_v(rr, w[i])
        g_full = (g - float(g @ a) * a) / na
        return f, self.P.T @ g_full

    def value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        v = self.quad(a)
        res = v - self.U @ (self.U.T @ v)
        return float(res @ res), float(np.linalg.norm(res) / max(np.linalg.norm(v), 1e-300))

    def minimize_from(self, seed, maxiter=500, gtol=1e-12):
        rng = np.random.default_rng(seed)
        u0 = rng.standard_normal(self.dimP); u0 /= np.linalg.norm(u0)
        t0 = time.time()
        res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B', options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18})
        a_full = self.P @ res.x; a = a_full / np.linalg.norm(a_full)
        f, rel = self.value(a)
        return {'seed': seed, 'f': f, 'rel': rel, 'iters': int(res.nit), 'seconds': round(time.time() - t0, 2), 'a': a}

    def exact_quadrics_modp(self, p):
        """Rows of the quadratic map modulo the kernel, as exact integer data mod p:
        returns (list of quadric coefficient dicts {(k,l): coeff}, number of quadrics)."""
        raise NotImplementedError


# ----------------------------------------------------------------------------------------------------------
# multiprocessing driver
# ----------------------------------------------------------------------------------------------------------
_G = {}


def _worker(args):
    seed, maxiter = args
    fam = _G['obj']
    r = fam.minimize_from(seed, maxiter=maxiter)
    out = {k: v for k, v in r.items() if k not in ('a', 'c', 'Dv')}
    keep = r['rel'] < _G.get('keep_rel', 1e-6)
    if keep:
        out['a'] = [float(x) for x in r['a']]
    with open(_G['restart_log'], 'a') as f:
        f.write(json.dumps(out) + '\n')
    return out, (r['a'] if keep else None)


def run_restarts(obj, name, seeds, maxiter=300, workers=4, keep_rel=1e-6, threads_per_worker=1):
    """Run restarts with a fork pool; per-restart lines in wave2-search/<name>.restarts.jsonl."""
    import multiprocessing as mp
    os.environ.setdefault('OMP_NUM_THREADS', str(threads_per_worker))
    _G['obj'] = obj; _G['restart_log'] = f'{SD_DIR}/{name}.restarts.jsonl'; _G['keep_rel'] = keep_rel
    done = set()
    if os.path.exists(_G['restart_log']):
        for line in open(_G['restart_log']):
            try:
                done.add(json.loads(line)['seed'])
            except Exception:
                pass
    todo = [(s, maxiter) for s in seeds if s not in done]
    results = []
    t0 = time.time()
    if workers <= 1:
        for a in todo:
            results.append(_worker(a))
    else:
        ctx = mp.get_context('fork')
        with ctx.Pool(workers) as pool:
            for r in pool.imap_unordered(_worker, todo, chunksize=1):
                results.append(r)
    return results, time.time() - t0


def summarize(results):
    rels = np.array([r[0]['rel'] for r in results]) if results else np.zeros(0)
    nonsp = np.array([not r[0]['spurious'] for r in results]) if results else np.zeros(0, bool)
    if rels.size == 0:
        return {'n': 0}
    rr = rels[nonsp] if nonsp.any() else rels
    hist = {}
    for lo, hi, tag in [(0, 1e-10, '<1e-10'), (1e-10, 1e-8, '1e-10..1e-8'), (1e-8, 1e-6, '1e-8..1e-6'), (1e-6, 1e-4, '1e-6..1e-4'),
                        (1e-4, 1e-2, '1e-4..1e-2'), (1e-2, 1e-1, '1e-2..1e-1'), (1e-1, 10, '>1e-1')]:
        hist[tag] = int(np.sum((rr >= lo) & (rr < hi)))
    return {'n': int(rels.size), 'n_spurious': int((~nonsp).sum()), 'best_rel': float(rr.min()), 'median_rel': float(np.median(rr)),
            'worst_rel': float(rr.max()), 'hist': hist,
            'mean_iters': float(np.mean([r[0]['iters'] for r in results])),
            'mean_seconds': float(np.mean([r[0]['seconds'] for r in results]))}
