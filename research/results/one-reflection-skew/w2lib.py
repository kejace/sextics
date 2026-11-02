"""w2lib -- wave-2 machinery for lane one-reflection-skew (does not modify the engine).

Exact five-variable reduction (fringewidth; re-derived in the wave-2 report):
    F = f(x) + c*Y2^6 (x = X0,X1,X2,Y0,Y1), A free of Y2, q = 1.
    A certificate exists  <=>  a rational field u(x) (5 components, homogeneous degree 1) with
        A - u(f) + (6/5) f div u = 0.
    Writing u = N/D (deg D = d, deg N_i = d+1) this is the strict polynomial identity
        A D^2 - D N(f) + (6/5) f (D div N - N(D)) = 0                     (*)
    and (D^2, (N_i D)_{i<5}, Y2 (D div N - N(D))/5, Q = 0) is a strict 6-variable certificate.
    Necessary E1-type condition: A D in J_f (f irreducible, f does not divide D).

Sectors: rho = sign change of (Y0,Y1) [weight e3+e4], sigma = sign change of (X1,Y1) [weight e1+e4];
D may be taken in any (rho,sigma) character; targets A = X0^3 Y0^3 (rho-odd, sigma-even).

Numerics: for a linear family D = D(c), c = B u (B orthonormal, nD x K, monomial coordinates), the
fixed-D system is M(c) x = b(c) with M linear and b = -T D^2 quadratic in c.  The reduced functional
phi(c) = min_x ||M x - b||^2 + lam ||x||^2 (Tikhonov, lam ~ 1e-12 tr(M^T M)/n) is evaluated through
the normal equations (dsyrk + Cholesky) and differentiated by the envelope theorem; the exact
least-squares residual (SVD) is recomputed at the end of every restart.  Candidates with a tiny
exact residual are rationalised (limit_denominator ladder + LLL) and confirmed by an exact modular
solve with CRT lifting and core.verify_witness on the 6-variable lift.
"""
from __future__ import annotations
import sys, os, time, json, math, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.linalg import blas as _blas
from scipy.optimize import minimize
from fractions import Fraction
from hodge import core
from hodge.core import mul, add, sub, scale, pderiv, shift
from flint import nmod_mat, fmpz_mat

LANE = 'one-reflection-skew'
WD = f'/tmp/hodge/results/{LANE}'
import socket as _socket
# Remote jobs log to a file that never exists on the local machine: the runner's rsync (`up`) pushes the local
# lane dir without --update, so a local copy of wave2-search.log.jsonl would clobber the remote one.
_ON_REMOTE = not _socket.gethostname().startswith('brunos')
W2LOG = f'{WD}/wave2-remote.log.jsonl' if _ON_REMOTE else f'{WD}/wave2-search.log.jsonl'
PR = core.PRIMES
SIXFIFTH = 1.2


def _is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    dd, r = n - 1, 0
    while dd % 2 == 0:
        dd //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17):
        x = pow(a, dd, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes31(count, start=2147483647):
    """`count` primes below 2^31, descending from `start` (deterministic)."""
    out = []
    n = start
    while len(out) < count:
        if _is_prime(n):
            out.append(n)
        n -= 2 if n % 2 else 1
    return out
DEFL_RUNS = 10
DEFL_PERTURB = 0.0


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%S')


def log_event(rec, path=W2LOG):
    rec = dict(rec)
    rec.setdefault('ts', now())
    rec.setdefault('lane', LANE)
    rec.setdefault('engine', 'w2lib-5var-v1')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'a') as fh:
        fh.write(json.dumps(rec, default=str) + '\n')
    return rec


# ----------------------------------------------------------------------------
# monomials / parities (6-tuples with e5 = 0 so that core helpers apply)
# ----------------------------------------------------------------------------

def rpar(e):
    return (e[3] + e[4]) & 1


def spar(e):
    return (e[1] + e[4]) & 1


def parity_of(P, fn):
    s = {fn(e) for e in P}
    return s.pop() if len(s) == 1 else None


def mons5(deg, rp=None, sp_=None):
    if deg < 0:
        return []
    out = []
    for comb in itertools.combinations(range(deg + 4), 4):
        e = []
        prev = -1
        for x in comb:
            e.append(x - prev - 1)
            prev = x
        e.append(deg + 3 - prev)
        e = tuple(e) + (0,)
        if rp is not None and rpar(e) != rp:
            continue
        if sp_ is not None and spar(e) != sp_:
            continue
        out.append(e)
    return out


def code(E):
    E = np.asarray(E, dtype=np.int64)
    return E[..., 0] + 64 * (E[..., 1] + 64 * (E[..., 2] + 64 * (E[..., 3] + 64 * E[..., 4])))


def poly_from_vec(mons, vec, tol=0.0):
    out = {}
    for m, v in zip(mons, vec):
        if (abs(v) > tol) if tol else v:
            out[m] = v
    return out


def vec_from_poly(mons, P, dtype=float):
    idx = {m: k for k, m in enumerate(mons)}
    v = np.zeros(len(mons), dtype=dtype)
    for e, c in P.items():
        v[idx[e]] = dtype(c) if dtype is float else c
    return v


# ----------------------------------------------------------------------------
# problem data
# ----------------------------------------------------------------------------

class Problem:
    def __init__(self, F, A, name='?'):
        c = Fraction(F.get((0, 0, 0, 0, 0, 6), 0))
        assert c != 0, 'F must contain Y2^6'
        f = {e: Fraction(v) for e, v in F.items() if e[5] == 0}
        assert len(f) == len(F) - 1, 'F must be f(x) + c*Y2^6'
        assert all(e[5] == 0 for e in A), 'A must be free of Y2'
        self.name = name
        self.F = {e: Fraction(v) for e, v in F.items()}
        self.A = {e: Fraction(v) for e, v in A.items()}
        self.f, self.c = f, c
        self.fi = [pderiv(f, i) for i in range(5)]
        self.rA, self.sA = parity_of(self.A, rpar), parity_of(self.A, spar)
        self.rf, self.sf = parity_of(f, rpar), parity_of(f, spar)
        assert self.rA is not None and self.rf == 0, 'need f rho-even and A rho-pure'
        self.use_sigma = (self.sf == 0) and (self.sA is not None)
        self.fterms = [(np.array(e[:5], dtype=np.int64), float(v)) for e, v in f.items()]

    @classmethod
    def from_instance(cls, iid):
        F, A, q, rec = core.load_instance(iid)
        assert q == 1
        return cls(F, A, iid)

    def sectors(self, pD, sD, T):
        """Row parities and per-component N parities for D in (pD, sD) and target T."""
        rT = parity_of(T, rpar)
        sT = parity_of(T, spar) if self.use_sigma else None
        use_s = self.use_sigma and sT is not None and sD is not None
        Npar = []
        for i in range(5):
            rp = (rT - pD - (1 if i >= 3 else 0)) % 2
            spp = ((sT - sD - (1 if i in (1, 4) else 0)) % 2) if use_s else None
            Npar.append((rp, spp))
        return rT, (sT if use_s else None), Npar

    def Dmons(self, d, pD, sD):
        return mons5(d, pD, sD if (self.use_sigma and sD is not None) else None)


# ----------------------------------------------------------------------------
# the fixed-D operator, linear in D (monomial coordinates)
# ----------------------------------------------------------------------------

class Operator5:
    """Columns (i, m): N_i = x^m; rows: monomials of degree 2d+6; entries linear in D.

    For D = x^eD the column (i, m) is  sum_t c_t [ -e_t[i] + (6/5)(m_i - eD_i) ] x^{eD + m + e_t - e_i}
    (t over the terms c_t x^{e_t} of f).  b(c) = -T D(c)^2.
    """

    def __init__(self, prob, d, pD, sD, T, verbose=False):
        t0 = time.time()
        self.prob, self.d, self.pD, self.sD = prob, d, pD, sD
        self.T = {e: Fraction(v) for e, v in T.items()}
        rT, sT, Npar = prob.sectors(pD, sD, T)
        self.rT, self.sT, self.Npar = rT, sT, Npar
        self.Dm = prob.Dmons(d, pD, sD)
        self.Ncols = [(i, m) for i in range(5) for m in mons5(d + 1, *Npar[i])]
        self.rows = mons5(2 * d + 6, rT, sT)
        nD, n, m = len(self.Dm), len(self.Ncols), len(self.rows)
        self.nD, self.n, self.m = nD, n, m
        rcode = code(np.array(self.rows))
        order = np.argsort(rcode)
        self._rc = rcode[order]
        self._ro = order
        ED = np.array([e[:5] for e in self.Dm], dtype=np.int64)
        self.ED = ED
        RI, CI, AI, VA = [], [], [], []
        for j, (i, mm) in enumerate(self.Ncols):
            mi = np.array(mm[:5], dtype=np.int64)
            for et, ct in prob.fterms:
                sh = mi + et
                sh[i] -= 1
                coef = ct * (-float(et[i]) + SIXFIFTH * (float(mm[i]) - ED[:, i].astype(float)))
                nz = np.nonzero(coef)[0]
                if nz.size == 0:
                    continue
                E = ED[nz] + sh[None, :]
                if (E < 0).any():
                    raise RuntimeError('negative exponent in operator (bug)')
                r = self._lookup(E)
                RI.append(r); CI.append(np.full(nz.size, j, dtype=np.int64)); AI.append(nz.astype(np.int64)); VA.append(coef[nz])
        self.g_rows = np.concatenate(RI); self.g_cols = np.concatenate(CI)
        self.g_a = np.concatenate(AI); self.g_val = np.concatenate(VA)
        self.T_op = sp.csr_matrix((self.g_val, (self.g_rows * n + self.g_cols, self.g_a)), shape=(m * n, nD))
        # target block
        PR_, PC_, PV_ = [], [], []
        ab = np.arange(nD * nD, dtype=np.int64)
        EA = (ED[:, None, :] + ED[None, :, :]).reshape(nD * nD, 5)
        for es, ts in self.T.items():
            E = EA + np.array(es[:5], dtype=np.int64)[None, :]
            PR_.append(self._lookup(E)); PC_.append(ab); PV_.append(np.full(nD * nD, -float(ts)))
        self.P_op = sp.csr_matrix((np.concatenate(PV_), (np.concatenate(PR_), np.concatenate(PC_))), shape=(m, nD * nD))
        self.nnz = self.g_val.size
        self.build_time = time.time() - t0
        if verbose:
            print(f'[op5 {prob.name}] d={d} sector=({pD},{sD}) rows={m} cols={n} nD={nD} nnz={self.nnz} '
                  f'build {self.build_time:.1f}s', file=sys.stderr)

    def _lookup(self, E):
        cc = code(E)
        pos = np.searchsorted(self._rc, cc)
        pos = np.minimum(pos, len(self._rc) - 1)
        if not (self._rc[pos] == cc).all():
            raise RuntimeError('exponent outside the row space (parity bookkeeping bug)')
        return self._ro[pos]

    def assemble(self, c):
        c = np.asarray(c, float)
        M = (self.T_op @ c).reshape(self.m, self.n)
        b = self.P_op @ np.outer(c, c).ravel()
        return M, b

    def grad_c(self, c, r, x):
        gT = np.bincount(self.g_a, weights=self.g_val * r[self.g_rows] * x[self.g_cols], minlength=self.nD)
        G = (self.P_op.T @ r).reshape(self.nD, self.nD)
        gb = (G + G.T) @ c
        return 2.0 * (gT - gb)

    def D_poly(self, c, tol=0.0):
        return poly_from_vec(self.Dm, c, tol)


# ----------------------------------------------------------------------------
# linear families of D (orthonormal bases in monomial coordinates)
# ----------------------------------------------------------------------------

def orthonormalize(V, tol=1e-9, rank=None):
    V = np.asarray(V, float)
    if V.size == 0 or V.shape[1] == 0:
        return np.zeros((V.shape[0], 0))
    U, s, _ = np.linalg.svd(V, full_matrices=False)
    if rank is not None:
        r = min(int(rank), U.shape[1])
    else:
        r = int(np.sum(s > tol * s[0])) if s.size else 0
    return U[:, :r]


def complement_within(B, W, tol=1e-9):
    """Orthonormal basis of span(B) ∩ span(W)^perp (W assumed inside span(B))."""
    if W.shape[1] == 0:
        return B
    Pm = B - W @ (W.T @ B)
    return orthonormalize(Pm, tol)


def intersect(U, V, rcond=1e-10):
    """Orthonormal basis of span(U) ∩ span(V)."""
    if U.shape[1] == 0 or V.shape[1] == 0:
        return np.zeros((U.shape[0], 0))
    ns = sla.null_space(np.hstack([U, -V]), rcond=rcond)
    if ns.shape[1] == 0:
        return np.zeros((U.shape[0], 0))
    return orthonormalize(U @ ns[:U.shape[1]])


def _dense(rows, polys):
    ridx = {e: k for k, e in enumerate(rows)}
    Mx = np.zeros((len(rows), len(polys)))
    for j, P in enumerate(polys):
        for e, v in P.items():
            Mx[ridx[e], j] = float(v)
    return Mx


def _modp_rank(rows, polys, p=PR[0]):
    ridx = {e: k for k, e in enumerate(rows)}
    R, C = len(rows), len(polys)
    if R == 0 or C == 0:
        return 0
    flat = [0] * (R * C)
    for j, P in enumerate(polys):
        Pp = core.reduce_mod(P, p)
        for e, v in Pp.items():
            flat[ridx[e] * C + j] = v
    return nmod_mat(R, C, flat, p).rank()


def e1_kernel(prob, d, pD, sD, p=PR[0], verbose=False, T=None):
    """Orthonormal float basis (nD x k) of {c : T D(c) in J_f} (T = A by default) in sector (pD, sD); dims cross-checked mod p."""
    t0 = time.time()
    A, f = (prob.A if T is None else {e: Fraction(v) for e, v in T.items()}), prob.f
    rT_, sT_ = parity_of(A, rpar), parity_of(A, spar)
    rAD = (rT_ + pD) % 2
    use_s = prob.use_sigma and sD is not None and sT_ is not None
    sAD = ((sT_ + sD) % 2) if use_s else None
    Dm = prob.Dmons(d, pD, sD)
    rows = mons5(d + 6, rAD, sAD)
    gens = []
    for i in range(5):
        rp = (rAD - (1 if i >= 3 else 0)) % 2
        spp = ((sAD - (1 if i in (1, 4) else 0)) % 2) if use_s else None
        for mm in mons5(d + 1, rp, spp):
            gens.append(shift(prob.fi[i], mm, 1))
    Acols = [shift(A, mm, 1) for mm in Dm]
    nD = len(Dm)
    rJ_p = _modp_rank(rows, gens, p)
    rAJ_p = _modp_rank(rows, Acols + gens, p)
    kerD_p = nD - (rAJ_p - rJ_p)
    info = {'d': d, 'pD': pD, 'sD': sD if use_s else None, 'nD': nD, 'rows_d6': len(rows), 'nJgens': len(gens),
            'dimJ_d6': rJ_p, 'rank_AD_mod_J': rAJ_p - rJ_p, 'kerD': kerD_p, 'p': p}
    # Jacobian part J_d in the sector
    Jd = []
    for i in range(5):
        rp = (pD - (1 if i >= 3 else 0)) % 2
        spp = ((sD - (1 if i in (1, 4) else 0)) % 2) if use_s else None
        for g in mons5(d - 5, rp, spp):
            Jd.append(shift(prob.fi[i], g, 1))
    info['dimJ_d'] = _modp_rank(Dm, Jd, p) if Jd else 0
    info['ker_mod_J'] = kerD_p - info['dimJ_d']
    # f-multiples
    fm = [shift(f, g, 1) for g in mons5(d - 6, pD, sD if use_s else None)]
    info['dim_fmult'] = _modp_rank(Dm, fm, p) if fm else 0
    if kerD_p == 0:
        info['seconds'] = round(time.time() - t0, 1)
        return np.zeros((nD, 0)), info
    MJ = _dense(rows, gens)
    MA = _dense(rows, Acols)
    if MJ.shape[1]:
        U, s, _ = np.linalg.svd(MJ, full_matrices=True)
        rJ = int(np.sum(s > 1e-9 * s[0]))
        if rJ != rJ_p:
            raise RuntimeError(f'float rank {rJ} != mod-p rank {rJ_p} for J_{d+6}')
        Qperp = U[:, rJ:]
        Kmat = Qperp.T @ MA
    else:
        Kmat = MA
    # the kernel dimension is known exactly (mod p); take exactly kerD_p right singular vectors and record the gap
    Uk, sk, Vk = np.linalg.svd(Kmat, full_matrices=True)
    nD_ = Kmat.shape[1]
    rank_k = nD_ - kerD_p
    ns = Vk[rank_k:].T.copy()
    s_small = float(sk[rank_k]) if rank_k < len(sk) else 0.0          # largest singular value treated as zero
    s_keep = float(sk[rank_k - 1]) if rank_k >= 1 else float('inf')   # smallest singular value kept as nonzero
    info['sv_gap'] = [s_keep, s_small]
    ns_auto = sla.null_space(Kmat, rcond=1e-9)
    info['float_kernel_dim_auto'] = int(ns_auto.shape[1])
    if s_small > 1e-6 * s_keep and s_small > 1e-6:
        raise RuntimeError(f'ambiguous kernel at d={d} sector ({pD},{sD}): mod-p dim {kerD_p}, singular values around the cut {s_keep:.2e} / {s_small:.2e}')
    info['seconds'] = round(time.time() - t0, 1)
    if verbose:
        print(f'[e1ker {prob.name}] {info}', file=sys.stderr)
    return ns, info


def jacobian_family(prob, d, pD, sD):
    Dm = prob.Dmons(d, pD, sD)
    use_s = prob.use_sigma and sD is not None
    vecs = []
    for i in range(5):
        rp = (pD - (1 if i >= 3 else 0)) % 2
        spp = ((sD - (1 if i in (1, 4) else 0)) % 2) if use_s else None
        for g in mons5(d - 5, rp, spp):
            vecs.append(vec_from_poly(Dm, shift(prob.fi[i], g, 1)))
    if not vecs:
        return np.zeros((len(Dm), 0))
    polys = []
    for i in range(5):
        rp = (pD - (1 if i >= 3 else 0)) % 2
        spp = ((sD - (1 if i in (1, 4) else 0)) % 2) if use_s else None
        for g in mons5(d - 5, rp, spp):
            polys.append(shift(prob.fi[i], g, 1))
    rk = _modp_rank(Dm, polys)
    return orthonormalize(np.array(vecs).T, rank=rk)


def fmult_family(prob, d, pD, sD):
    Dm = prob.Dmons(d, pD, sD)
    use_s = prob.use_sigma and sD is not None
    gs = mons5(d - 6, pD, sD if use_s else None)
    vecs = [vec_from_poly(Dm, shift(prob.f, g, 1)) for g in gs]
    if not vecs:
        return np.zeros((len(Dm), 0))
    rk = _modp_rank(Dm, [shift(prob.f, g, 1) for g in gs])
    return orthonormalize(np.array(vecs).T, rank=rk)


def monomial_subfamily(prob, d, pD, sD, pred):
    """Coordinate vectors of the D-monomials satisfying pred(exponent)."""
    Dm = prob.Dmons(d, pD, sD)
    cols = [k for k, m in enumerate(Dm) if pred(m)]
    B = np.zeros((len(Dm), len(cols)))
    for j, k in enumerate(cols):
        B[k, j] = 1.0
    return B


# ----------------------------------------------------------------------------
# variable projection
# ----------------------------------------------------------------------------

class VarPro5:
    def __init__(self, op, B, lam_rel=1e-12):
        self.op = op
        self.B = np.ascontiguousarray(np.asarray(B, float))
        self.K = self.B.shape[1]
        self.lam_rel = lam_rel
        self.nfev = 0

    def inner(self, M, b):
        n = M.shape[1]
        G = _blas.dsyrk(1.0, M, trans=1)      # upper triangle of M^T M
        tr = float(np.trace(G))
        lam = self.lam_rel * tr / n if tr > 0 else 1e-300
        G[np.diag_indices(n)] += lam
        Mtb = M.T @ b
        try:
            cf = sla.cho_factor(G, lower=False, check_finite=False)
            x = sla.cho_solve(cf, Mtb, check_finite=False)
        except np.linalg.LinAlgError:
            x, *_ = sla.lstsq(M, b, cond=1e-10, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        return x, r, lam

    def value_and_grad(self, u):
        nu = float(np.linalg.norm(u))
        c = self.B @ (u / nu)
        M, b = self.op.assemble(c)
        x, r, lam = self.inner(M, b)
        phi = float(r @ r + lam * (x @ x))
        gc = self.op.grad_c(c, r, x)
        g = self.B.T @ (gc - float(gc @ c) * c) / nu
        self.nfev += 1
        return phi, g

    def exact_residual(self, c):
        """Exact least-squares relative residual ||r||/||b|| (SVD) at the monomial vector c (normalised)."""
        c = np.asarray(c, float)
        c = c / np.linalg.norm(c)
        M, b = self.op.assemble(c)
        x, *_ = sla.lstsq(M, b, cond=1e-10, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        return float(np.linalg.norm(r) / max(np.linalg.norm(b), 1e-300)), x, r

    def sparse_start(self, rng, k):
        nD = self.op.nD
        c0 = np.zeros(nD)
        idx = rng.choice(nD, size=min(k, nD), replace=False)
        c0[idx] = rng.integers(1, 4, size=idx.size) * rng.choice([-1.0, 1.0], size=idx.size)
        u0 = self.B.T @ c0
        if np.linalg.norm(u0) < 1e-8:
            u0 = rng.standard_normal(self.K)
        return u0

    # -- Levenberg-Marquardt on the projected residual r(c) = M(c) x(c) - b(c), full Golub-Pereyra Jacobian --
    def _prep(self, u):
        key = np.asarray(u, float).tobytes()
        if getattr(self, '_ckey', None) == key:
            return self._cache
        nu = float(np.linalg.norm(u))
        c = self.B @ (u / nu)
        op = self.op
        M, b = op.assemble(c)
        n = M.shape[1]
        G = _blas.dsyrk(1.0, M, trans=1)
        tr = float(np.trace(G))
        lam = self.lam_rel * tr / n if tr > 0 else 1e-300
        G[np.diag_indices(n)] += lam
        cf = sla.cho_factor(G, lower=False, check_finite=False)
        x = sla.cho_solve(cf, M.T @ b, check_finite=False)
        r = M @ x - b
        self.nfev += 1
        self._ckey = key
        self._cache = (nu, c, M, b, cf, x, r)
        return self._cache

    def residual(self, u):
        return self._prep(u)[6]

    def jacobian(self, u, full=True):
        nu, c, M, b, cf, x, r = self._prep(u)
        op = self.op
        # W[:, a] = M_a x ;  Bb[:, a] = d b / d c_a = sum_b (P[:, (a,b)] + P[:, (b,a)]) c_b
        Wm = np.bincount(op.g_rows * op.nD + op.g_a, weights=op.g_val * x[op.g_cols],
                         minlength=op.m * op.nD).reshape(op.m, op.nD)
        cs = sp.csr_matrix(c[:, None])
        PC1 = (op.P_op @ sp.kron(sp.identity(op.nD, format='csr'), cs)).toarray()
        PC2 = (op.P_op @ sp.kron(cs, sp.identity(op.nD, format='csr'))).toarray()
        V = Wm - (PC1 + PC2)
        V = V - M @ sla.cho_solve(cf, M.T @ V, check_finite=False)          # P_perp (M_a x - d_a b)
        if full:
            Z = np.bincount(op.g_cols * op.nD + op.g_a, weights=op.g_val * r[op.g_rows],
                            minlength=op.n * op.nD).reshape(op.n, op.nD)      # Z[:, a] = M_a^T r
            V = V - M @ sla.cho_solve(cf, Z, check_finite=False)              # - (M^+)^T M_a^T r
        Jc = V @ self.B
        u = np.asarray(u, float)
        Ju = (Jc - np.outer(Jc @ u, u) / (nu * nu)) / nu
        return Ju

    def residual_and_jac(self, u, want_jac=True, full=True):
        r = self.residual(u)
        return r, (self.jacobian(u, full) if want_jac else None)

    # -- deflation (Farrell-Birkisson-Funke): r_defl = eta(u) r(u), eta = prod_i (s_i^{-p/2} + alpha), s_i = 1 - (û·û_i)^2 --
    def _deflation(self, u, defl, p=2.0, alpha=1.0):
        if not defl:
            return 1.0, np.zeros_like(u)
        nu = float(np.linalg.norm(u)); uh = u / nu
        eta = 1.0
        grad = np.zeros_like(u)
        for ui in defl:
            cth = float(uh @ ui)
            s = max(1.0 - cth * cth, 1e-14)
            term = s ** (-p / 2) + alpha
            dterm = -(p / 2) * s ** (-p / 2 - 1) * (-2.0 * cth) * (ui - cth * uh) / nu
            eta *= term
            grad = grad + dterm / term
        return eta, eta * grad

    def restart_lm(self, seed, maxiter=80, init='gauss', sparse_k=6, full=True, defl=None, u0=None):
        rng = np.random.default_rng(seed)
        if u0 is None:
            u0 = self.sparse_start(rng, sparse_k) if init == 'sparse' else rng.standard_normal(self.K)
        u0 = np.asarray(u0, float) / np.linalg.norm(u0)
        t0 = time.time()
        self.nfev = 0
        self._ckey = None
        from scipy.optimize import least_squares
        defl = list(defl) if defl else []

        def fun(u):
            r = self.residual(u)
            if not defl:
                return r
            eta, _ = self._deflation(u, defl)
            return eta * r

        def jac(u):
            J = self.jacobian(u, full)
            if not defl:
                return J
            r = self.residual(u)
            eta, geta = self._deflation(u, defl)
            return eta * J + np.outer(r, geta)
        res = least_squares(fun, u0, jac=jac, method='lm', max_nfev=maxiter, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        u = res.x / np.linalg.norm(res.x)
        c = self.B @ u
        rel, x, r = self.exact_residual(c)
        return {'seed': int(seed), 'init': init, 'method': ('lm' if full else 'lmk') + ('-defl' if defl else ''),
                'phi': float(2 * res.cost), 'rel': rel, 'iters': int(res.nfev), 'nfev': int(self.nfev),
                'secs': round(time.time() - t0, 2), 'c': c, 'u': u, 'msg': str(res.message)[:60], 'n_defl': len(defl)}

    def deflated_chain(self, seed, n_runs=10, maxiter=80, init='gauss', sparse_k=6, perturb=0.0):
        """Sequential LM runs; every local minimum found is deflated for the subsequent runs.
        perturb > 0: start the next run from a perturbation of the last minimum instead of a fresh random point."""
        rng = np.random.default_rng(seed)
        defl, outs = [], []
        last = None
        for k in range(n_runs):
            u0 = None
            if perturb > 0 and last is not None:
                u0 = last + perturb * rng.standard_normal(self.K)
            out = self.restart_lm(seed * 1000 + k, maxiter=maxiter, init=init if (k % 2 == 0) else 'sparse',
                                  sparse_k=sparse_k, defl=defl, u0=u0)
            outs.append(out)
            if out['rel'] > 1e-8:          # a genuine (non-zero) local minimum: deflate it
                defl.append(out['u'])
            last = out['u']
            if out['rel'] < 1e-8:
                break
        return outs

    def restart(self, seed, maxiter=300, init='gauss', sparse_k=6, gtol=1e-13, method='lbfgs'):
        if method in ('lm', 'lmk'):
            return self.restart_lm(seed, maxiter=maxiter, init=init, sparse_k=sparse_k, full=(method == 'lm'))
        if method == 'defl':
            return self.deflated_chain(seed, n_runs=DEFL_RUNS, maxiter=maxiter, init=init, sparse_k=sparse_k, perturb=DEFL_PERTURB)
        rng = np.random.default_rng(seed)
        if init == 'sparse':
            u0 = self.sparse_start(rng, sparse_k)
        else:
            u0 = rng.standard_normal(self.K)
        u0 = u0 / np.linalg.norm(u0)
        t0 = time.time()
        self.nfev = 0
        res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                       options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
        u = res.x / np.linalg.norm(res.x)
        c = self.B @ u
        rel, x, r = self.exact_residual(c)
        return {'seed': int(seed), 'init': init, 'method': 'lbfgs', 'phi': float(res.fun), 'rel': rel,
                'iters': int(res.nit), 'nfev': int(self.nfev), 'secs': round(time.time() - t0, 2), 'c': c,
                'msg': str(res.message)[:60]}


# ----------------------------------------------------------------------------
# rationalisation and exact confirmation
# ----------------------------------------------------------------------------

DENS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 24, 28, 30, 35, 36, 42, 48, 56, 60, 63, 70, 72,
        84, 90, 105, 120, 126, 140, 168, 180, 210, 252, 280, 315, 336, 360, 420, 504, 560, 630, 720, 840, 1008,
        1260, 1680, 2520, 5040, 10**4, 3 * 10**4, 10**5, 3 * 10**5, 10**6, 10**7)


def rationalize(c, tol=1e-6, dens=DENS, lll_dim=60):
    """Rational candidate vectors for c (normalised so that the largest entry is 1)."""
    c = np.asarray(c, float)
    j = int(np.argmax(np.abs(c)))
    c = c / c[j]
    out = []
    for den in dens:
        fr = [Fraction(float(v)).limit_denominator(den) for v in c]
        err = max(abs(float(x) - v) for x, v in zip(fr, c))
        if err < tol:
            cand = tuple(fr)
            if cand not in out:
                out.append(cand)
    # LLL simultaneous Diophantine approximation on the largest coordinates
    try:
        order = np.argsort(-np.abs(c))[:lll_dim]
        sub_ = c[order]
        S = 10**10
        n = len(sub_)
        rows = [[1] + [int(round(S * v)) for v in sub_]]
        for i in range(n):
            rows.append([0] * (1 + i) + [S] + [0] * (n - 1 - i))
        Lm = fmpz_mat(rows).lll()
        for i in range(min(3, n + 1)):
            qd = int(Lm[i, 0])
            if qd == 0:
                continue
            qd = abs(qd)
            fr = tuple(Fraction(int(round(qd * v)), qd) for v in c)
            err = max(abs(float(x) - v) for x, v in zip(fr, c))
            if err < tol and fr not in out:
                out.append(fr)
    except Exception:
        pass
    return out


def fixed_D_system(prob, D, pD, sD, T, p):
    """Mod-p columns (N) and rhs for the fixed-D 5-variable system (*). Returns (cols, colpolys, rows, rhs)."""
    rT, sT, Npar = prob.sectors(pD, sD, T)
    fp = core.reduce_mod(prob.f, p)
    Dp = core.reduce_mod(D, p)
    d = core.degree(D)
    c65 = (6 * pow(5, -1, p)) % p
    fD = scale(mul(fp, Dp, p), c65, p)
    G = []
    for i in range(5):
        Di = pderiv(Dp, i, p)
        G.append(scale(add(mul(Dp, pderiv(fp, i, p), p), scale(mul(fp, Di, p), c65, p), p), p - 1, p))
    cols, colpolys = [], []
    for i in range(5):
        for m in mons5(d + 1, *Npar[i]):
            col = shift(G[i], m, 1, p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(fD, tuple(m2), m[i], p), p)
            cols.append((i, m)); colpolys.append(col)
    rhs = scale(mul(core.reduce_mod(T, p), mul(Dp, Dp, p), p), p - 1, p)
    rows = {}
    for col in colpolys:
        for e in col:
            rows.setdefault(e, len(rows))
    for e in rhs:
        rows.setdefault(e, len(rows))
    return cols, colpolys, rows, rhs


def e1_exact_modp(prob, D, p=PR[0]):
    """A*D in J_f mod p?"""
    AD = core.reduce_mod(mul(prob.A, D), p)
    if not AD:
        return True
    deg = core.degree(AD)
    rp, spp = parity_of(AD, rpar), (parity_of(AD, spar) if prob.use_sigma else None)
    gens = []
    for i in range(5):
        for mm in mons5(deg - 5, None, None):
            P = shift(prob.fi[i], mm, 1)
            if P and parity_of(P, rpar) == rp and (spp is None or parity_of(P, spar) == spp):
                gens.append(P)
    rows = sorted(set(e for g in gens for e in g) | set(AD))
    r0 = _modp_rank(rows, gens, p)
    r1 = _modp_rank(rows, gens + [AD], p)
    return r0 == r1


def exact_confirm(prob, D, pD, sD, T, primes=None, verbose=False, max_primes=6):
    """E1 mod p, fixed-D solve mod p, CRT lift, exact 6-variable verification.

    Returns (status, info): status in {'e1-fails', 'inconsistent', 'solved', 'lift-failed'}."""
    primes = primes or PR
    D = {e: Fraction(v) for e, v in D.items() if v}
    if not D:
        return 'inconsistent', {'note': 'D = 0'}
    t0 = time.time()
    if not e1_exact_modp(prob, D, primes[0]):
        return 'e1-fails', {'secs': round(time.time() - t0, 2)}
    acc = None; piv0 = None; last = {}
    for k, p in enumerate(primes[:max_primes]):
        cols, colpolys, rows, rhs = fixed_D_system(prob, D, pD, sD, T, p)
        R, C = len(rows), len(cols)
        flat = [0] * (R * (C + 1))
        for j, col in enumerate(colpolys):
            for e, v in col.items():
                flat[rows[e] * (C + 1) + j] = v
        for e, v in rhs.items():
            flat[rows[e] * (C + 1) + C] = v
        M = nmod_mat(R, C + 1, flat, p)
        Mr, rank = M.rref()
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0:
                j += 1
            pivots.append(j)
        last = {'p': p, 'rank': rank, 'rows': R, 'cols': C, 'nullity': C - rank, 'secs': round(time.time() - t0, 2)}
        if pivots and pivots[-1] == C:
            return 'inconsistent', last
        x = [0] * C
        for i, j in enumerate(pivots):
            x[j] = int(Mr[i, C])
        if piv0 is None or piv0 != pivots:
            acc = (x, p); piv0 = pivots
        else:
            xs, mmod = acc
            xs = [core.crt_pair(a, mmod, b, p)[0] for a, b in zip(xs, x)]
            acc = (xs, mmod * p)
        xs, mmod = acc
        rat = [core.ratrecon(a, mmod) for a in xs]
        if all(v is not None for v in rat):
            N = [{} for _ in range(6)]
            for (i, mon), v in zip(cols, rat):
                if v:
                    N[i][mon] = v
            divN = {}; ND = {}
            for i in range(5):
                divN = add(divN, pderiv(N[i], i)); ND = add(ND, mul(N[i], pderiv(D, i)))
            BD = sub(mul(D, divN), ND)
            # exact 5-variable identity
            lhs = add(sub(mul(T, mul(D, D)), mul(D, sum_((mul(N[i], prob.fi[i]) for i in range(5))))),
                      scale(mul(prob.f, BD), Fraction(6, 5)))
            if lhs:
                if verbose:
                    print('[confirm] rational lift failed the 5-var identity; continuing', file=sys.stderr)
                continue
            Y2 = {(0, 0, 0, 0, 0, 1): Fraction(1)}
            N6 = [mul(N[i], D) for i in range(5)] + [scale(mul(Y2, BD), Fraction(1, 5))]
            D6 = mul(D, D)
            ok, vinfo = core.verify_witness(prob.F, T, 1, D6, N6, {})
            last.update({'N5': N, 'D6': D6, 'N6': N6, 'verified6': ok, 'primes_used': k + 1, 'secs': round(time.time() - t0, 2)})
            return ('solved' if ok else 'lift-failed'), last
    last['primes_used'] = max_primes
    return 'lift-failed', last


def sum_(it):
    out = {}
    for P in it:
        out = add(out, P)
    return out


def served_analysis(prob, D, pD, sD, p=PR[0]):
    """Dimension of the space of degree-6 targets (A's parities) served by D, and modulo J_f,6."""
    t0 = time.time()
    T0 = prob.A
    cols, colpolys, rows, _ = fixed_D_system(prob, D, pD, sD, T0, p)
    Dp = core.reduce_mod(D, p)
    D2 = mul(Dp, Dp, p)
    tm = mons5(6, prob.rA, prob.sA if prob.use_sigma else None)
    tcols = [shift(D2, m, 1, p) for m in tm]
    for col in tcols:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols) + len(tcols)
    flat = [0] * (R * C)
    for j, col in enumerate(colpolys + tcols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    nN = len(cols); nT = len(tm)
    if nullity == 0:
        return {'served_dim': 0, 'served_mod_J': 0, 'rank_N': nN, 'nullity_N': 0, 'secs': round(time.time() - t0, 1)}
    proj = nmod_mat(nT, nullity, [int(X[nN + j, k]) for j in range(nT) for k in range(nullity)], p)
    Tr, rk = proj.transpose().rref()
    basis = []
    for i in range(rk):
        basis.append({tm[j]: int(Tr[i, j]) for j in range(nT) if int(Tr[i, j])})
    Jb = []
    for i in range(5):
        for j in range(5):
            P = core.reduce_mod(mul(core.parse(core.VARS[j]), prob.fi[i]), p)
            if P and parity_of(P, rpar) == prob.rA and (not prob.use_sigma or parity_of(P, spar) == prob.sA):
                Jb.append(P)
    rJ = _modp_rank(tm, Jb, p)
    rJS = _modp_rank(tm, Jb + basis, p)
    rA = _modp_rank(tm, Jb + basis + [core.reduce_mod(prob.A, p)], p)
    return {'served_dim': rk, 'served_mod_J': rJS - rJ, 'serves_A_mod_J': rA == rJS, 'dimJ6': rJ,
            'rank_N': nN - (nullity - rk), 'nullity_N': nullity - rk, 'rows': R, 'cols': C, 'secs': round(time.time() - t0, 1)}


def save_witness_candidate(prob, D6, N6, Q, path=f'{WD}/witness-candidate.json', extra=None):
    ok, info = core.verify_witness(prob.F, prob.A, 1, D6, N6, Q)
    w = {'instance': prob.name, 'q': 1, 'd': core.degree(D6), 'D': core.to_str(D6), 'N': [core.to_str(n) for n in N6],
         'Q': core.to_str(Q) if Q else '0', 'verified': ok, 'info': str(info), 'ts': now()}
    if extra:
        w.update(extra)
    with open(path, 'w') as fh:
        json.dump(w, fh, indent=1)
    return ok


# ----------------------------------------------------------------------------
# snapping a numerical zero to a rational point (zero sets are often positive-dimensional)
# ----------------------------------------------------------------------------

def zero_set_dim(vp, c, tol=1e-6):
    """Estimated local dimension of the zero set at c: near-null directions of the LM Jacobian beyond scaling."""
    u = vp.B.T @ (np.asarray(c, float) / np.linalg.norm(c))
    J = vp.jacobian(u, True)
    s = np.linalg.svd(J, compute_uv=False)
    k = int(np.sum(s < tol * s[0]))
    return max(k - 1, 0), [float(x) for x in s[-6:]]


def snap_confirm(vp, prob, c, pD, sD, T, taus=(0.0, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1), maxiter=150,
                 rel_ok=1e-7, verbose=False, logger=None):
    """Try to turn a numerical zero into an exact witness: rationalise; if that fails, sparsify (drop monomials
    with |coefficient| < tau), refit by LM inside the sub-family, rationalise again, confirm exactly.
    Returns (status, info); status 'solved' carries D5/N5/D6/N6."""
    op, B = vp.op, vp.B
    c = np.asarray(c, float)
    c = c / c[int(np.argmax(np.abs(c)))]
    tried = set()
    attempts = []
    for tau in taus:
        if tau == 0.0:
            c2, rel2, kept = c, vp.exact_residual(c)[0], op.nD
        else:
            keep = [i for i in range(op.nD) if abs(c[i]) > tau]
            if not keep or len(keep) == op.nD:
                continue
            Bc = np.zeros((op.nD, len(keep)))
            for j, i in enumerate(keep):
                Bc[i, j] = 1.0
            Bs = intersect(B, Bc)
            if Bs.shape[1] == 0:
                continue
            vps = VarPro5(op, Bs)
            try:
                o2 = vps.restart_lm(1, maxiter=maxiter, u0=Bs.T @ c)
            except Exception as ex:
                attempts.append({'tau': tau, 'error': str(ex)[:80]})
                continue
            c2 = o2['c']; c2 = c2 / c2[int(np.argmax(np.abs(c2)))]
            rel2, kept = o2['rel'], len(keep)
            if rel2 > rel_ok:
                attempts.append({'tau': tau, 'kept': kept, 'rel': rel2, 'note': 'sparsification breaks the zero'})
                break
        rats = rationalize(c2)
        attempts.append({'tau': tau, 'kept': kept, 'rel': rel2, 'n_rational': len(rats)})
        for fr in rats[:6]:
            key = tuple(fr)
            if key in tried:
                continue
            tried.add(key)
            Dc = poly_from_vec(op.Dm, fr)
            st, inf = exact_confirm(prob, Dc, pD, sD, T)
            rec = {'tau': tau, 'kept': kept, 'rel': rel2, 'status': st, 'D': core.to_str(Dc)[:200], 'terms': len(Dc),
                   'rank': inf.get('rank'), 'nullity': inf.get('nullity')}
            attempts.append(rec)
            if logger:
                logger(rec)
            if verbose:
                print('[snap]', json.dumps(rec)[:300], file=sys.stderr)
            if st == 'solved':
                inf['D5'] = Dc
                inf['attempts'] = attempts
                return 'solved', inf
    return 'unconfirmed', {'attempts': attempts}


# ----------------------------------------------------------------------------
# exact E1 kernel (mod p nullspace -> canonical rref -> CRT -> rational reconstruction)
# ----------------------------------------------------------------------------

def _kernel_rref_modp(prob, d, pD, sD, T, p):
    A = {e: Fraction(v) for e, v in T.items()}
    rT_, sT_ = parity_of(A, rpar), parity_of(A, spar)
    rAD = (rT_ + pD) % 2
    use_s = prob.use_sigma and sD is not None and sT_ is not None
    sAD = ((sT_ + sD) % 2) if use_s else None
    Dm = prob.Dmons(d, pD, sD)
    rows = mons5(d + 6, rAD, sAD)
    gens = []
    for i in range(5):
        rp = (rAD - (1 if i >= 3 else 0)) % 2
        spp = ((sAD - (1 if i in (1, 4) else 0)) % 2) if use_s else None
        for mm in mons5(d + 1, rp, spp):
            gens.append(shift(prob.fi[i], mm, 1))
    Acols = [shift(A, mm, 1) for mm in Dm]
    ridx = {e: k for k, e in enumerate(rows)}
    R, nD, C = len(rows), len(Dm), len(Dm) + len(gens)
    flat = [0] * (R * C)
    for j, P in enumerate(Acols + gens):
        for e, v in core.reduce_mod(P, p).items():
            flat[ridx[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    if nullity == 0:
        return [], 0, Dm, {'rows': R, 'nJgens': len(gens)}
    proj = nmod_mat(nullity, nD, [int(X[j, k]) for k in range(nullity) for j in range(nD)], p)
    Pr, rk = proj.rref()
    basis = [[int(Pr[i, j]) for j in range(nD)] for i in range(rk)]
    return basis, rk, Dm, {'rows': R, 'nJgens': len(gens)}


def e1_kernel_exact(prob, d, pD, sD, T=None, primes=None, verbose=False, max_primes=400):
    """Exact rational rref basis of {c : T D(c) in J_f} (T = A by default) in sector (pD, sD), by CRT over 31-bit primes
    with rational reconstruction, verified against one further prime.  Returns (B_float_orthonormal (nD x k), info);
    info['exact_basis'] holds the rational rows (Fractions)."""
    t0 = time.time()
    T = prob.A if T is None else T
    # cache (exact bases are expensive at high degree: heights of hundreds/thousands of bits)
    import hashlib
    th = hashlib.sha1(core.to_str(T).encode()).hexdigest()[:8]
    cache_dir = f'{WD}/w2cache'
    os.makedirs(cache_dir, exist_ok=True)
    cpath = f'{cache_dir}/kernel-{prob.name}-d{d}-r{pD}s{sD}-{th}.json'
    if os.path.exists(cpath):
        try:
            cj = json.load(open(cpath))
            rat = [[Fraction(x) for x in row] for row in cj['basis']]
            Dm = prob.Dmons(d, pD, sD)
            k0 = len(rat)
            info = dict(cj['info']); info['cached'] = True
            if not rat:
                return np.zeros((len(Dm), 0)), dict(info, exact_basis=[])
            Bf = np.array([[float(v) for v in row] for row in rat]).T
            B = orthonormalize(Bf, rank=k0)
            info['exact_basis'] = rat
            return B, info
        except Exception:
            pass
    primes = primes or primes31(max_primes)
    acc = None; mod = None; k0 = None; Dm = None; meta = {}
    rat = None; verified = False; used = 0
    next_try = 2
    for p in primes[:max_primes]:
        basis, rk, Dm, meta = _kernel_rref_modp(prob, d, pD, sD, T, p)
        used += 1
        if k0 is None:
            k0 = rk
        elif rk != k0:
            raise RuntimeError(f'kernel dimension differs between primes ({k0} vs {rk}) at d={d} sector ({pD},{sD})')
        if rk == 0:
            rat = []; verified = True
            break
        if rat is not None:
            ok = True
            for ra, rb in zip(rat, basis):
                for a, b in zip(ra, rb):
                    if (a.numerator * pow(a.denominator, -1, p) - b) % p:
                        ok = False; break
                if not ok:
                    break
            if ok:
                verified = True
                break
        if acc is None:
            acc, mod = basis, p
        else:
            acc = [[core.crt_pair(a, mod, b, p)[0] for a, b in zip(ra, rb)] for ra, rb in zip(acc, basis)]
            mod *= p
        if used >= next_try:
            cand = []
            okall = True
            for row in acc:
                r2 = []
                for a in row:
                    v = core.ratrecon(a, mod)
                    if v is None:
                        okall = False; break
                    r2.append(v)
                if not okall:
                    break
                cand.append(r2)
            rat = cand if okall else None
            if rat is None:
                next_try = used + max(1, used // 4)   # geometric cadence of reconstruction attempts
        else:
            rat = None
    if not verified:
        raise RuntimeError(f'exact kernel not reconstructed/verified after {used} primes (d={d}, sector ({pD},{sD}), modulus {mod.bit_length() if mod else 0} bits)')
    nD = len(Dm)
    info = {'d': d, 'pD': pD, 'sD': sD, 'nD': nD, 'kerD': k0, 'primes_used': used, 'secs': round(time.time() - t0, 1),
            'max_height_bits': (max((max(abs(v.numerator), v.denominator) for row in rat for v in row), default=1)).bit_length() if rat else 0}
    info.update(meta)
    if not rat:
        try:
            with open(cpath, 'w') as fh:
                json.dump({'info': info, 'basis': []}, fh)
        except Exception:
            pass
        return np.zeros((nD, 0)), dict(info, exact_basis=[])
    Bf = np.array([[float(v) for v in row] for row in rat]).T          # nD x k
    B = orthonormalize(Bf, rank=k0)
    info['exact_basis'] = rat
    try:
        with open(cpath, 'w') as fh:
            json.dump({'info': {k: v for k, v in info.items() if k != 'exact_basis'}, 'basis': [[str(v) for v in row] for row in rat]}, fh)
    except Exception:
        pass
    if verbose:
        print(f'[e1ker-exact {prob.name}] d={d} sector=({pD},{sD}) kerD={k0} primes={used} heights<=2^{info["max_height_bits"]} ({info["secs"]}s)', file=sys.stderr)
    return B, info
