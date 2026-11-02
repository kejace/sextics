"""w2lib: wave-2 E2-reduced joint-search toolkit for eight-plane-configuration (q = 1, A = X0^3 Y0^3).

Normalization.  Any certificate can be normalized to N(F) = A*D exactly (subtract (G/6)*Euler where
N(F) = A D + F G; the E2 identity and Q are unchanged).  Since F is smooth, syzygies of (F_0..F_5) are
Koszul, so the E1 solutions are N = N_0(a) + K_c with K_c = sum_{i<j} c_ij (F_j d_i - F_i d_j),
deg c_ij = d - 4.  E2 then reads

    beta(D, N) := D div N - N(D) = F Q        (degree 2d, odd Y-parity, deg Q = 2d - 6).

Families.  Jacobian type (every E1-admissible D of degree <= 9 here): D = w(F), N_0 = A w, with w a
vector of forms of degree d-5 (parities forced by D's parity).  Kernel type (d >= 10): D = sum a_k D_k
over an E1-admissible basis with particular solutions N_k (N_k(F) = A D_k) computed in floats.

Reduced functional (variable projection).  For a on the unit sphere,
    f(a) = min_{c,Q} || Omega(D(a)) (N_0(a) + Kmat c) - T_F Q ||^2
where Omega(D) is the sparse matrix of N -> beta(D, N) (linear in D) and T_F is multiplication by F.
Q is eliminated with an orthonormal basis U of col(T_F); c by an SVD.  The Kaufman variable-projection
Jacobian is used in Levenberg-Marquardt.  The spurious subspace Z = {a : D(a) in F*S_{d-6}} (where f = 0
trivially) is projected out; this is lossless because (D + F g, N + (A g/6) Euler) is again a certificate.

Necessary condition ("Hessian relaxation").  Reducing E2 modulo J(F) kills the Koszul and F Q terms:
    [beta(D, N_0)] = 0 in R_{2d} = S/J(F);   for Jacobian type this is  A * Hess_F(w, w) in J(F)_{2d},
i.e. [Hess_F(w,w)] in Ann_R(A)_{2d-6}; a system of quadrics in w (exact mod p in hess_quadrics()).
"""
from __future__ import annotations
import sys, os, json, time, math, pickle
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import least_squares
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import monomials, mul, add, sub, scale, pderiv, shift, reduce_mod, parity_of, degree, NV

IID = 'eight-plane-configuration'
WD = f'/tmp/hodge/results/{IID}'
F, A, q, rec = core.load_instance(IID)
assert q == 1
pA = parity_of(A)
Fi = [pderiv(F, i) for i in range(NV)]
PRIME = core.PRIMES[1]   # 2147483629: same prime as the wave-1 remote kernel pickles; < 2^31 for msolve
BASE = 64
POW = np.array([BASE ** i for i in range(NV)], dtype=np.int64)


def jlog(path, obj):
    with open(path, 'a') as f:
        f.write(json.dumps(obj, default=str) + '\n')


def enc_rows(E):
    return (np.asarray(E, dtype=np.int64) @ POW)


class Space:
    """Monomials of one degree and Y-parity, with fast vectorized index lookup."""
    def __init__(self, deg, par):
        self.deg, self.par = deg, par
        self.mons = monomials(deg, par) if deg >= 0 else []
        self.n = len(self.mons)
        self.idx = {m: i for i, m in enumerate(self.mons)}
        self.E = np.array(self.mons, dtype=np.int64).reshape(-1, NV)
        self.codes = enc_rows(self.E) if self.n else np.zeros(0, dtype=np.int64)
        self.order = np.argsort(self.codes)
        self.sorted = self.codes[self.order]

    def lookup(self, codes):
        """indices of the monomials with the given codes (all must exist)."""
        pos = np.searchsorted(self.sorted, codes)
        pos = np.minimum(pos, max(self.n - 1, 0))
        ok = self.sorted[pos] == codes
        if not np.all(ok):
            raise KeyError('monomial not in space')
        return self.order[pos]

    def vec(self, P):
        v = np.zeros(self.n)
        for e, c in P.items():
            v[self.idx[e]] = float(c)
        return v

    def poly(self, v, tol=0.0):
        return {m: float(v[i]) for i, m in enumerate(self.mons) if abs(v[i]) > tol}


def mul_matrix(P, src: Space, dst: Space):
    """Sparse matrix (dst.n x src.n) of multiplication by the polynomial P."""
    if not P or src.n == 0:
        return sp.csr_matrix((dst.n, src.n))
    Pe = np.array(list(P.keys()), dtype=np.int64)
    Pv = np.array([float(v) for v in P.values()])
    codes = src.codes[None, :] + enc_rows(Pe)[:, None]       # nP x nsrc
    rows = dst.lookup(codes.ravel())
    cols = np.tile(np.arange(src.n), len(Pv))
    data = np.repeat(Pv, src.n)
    return sp.csr_matrix((data, (rows, cols)), shape=(dst.n, src.n))


class OmegaPattern:
    """Structure of Omega(D): N-coefficients (blocks N_0..N_5) -> coefficients of beta(D, N) in dst.

    beta(D, m e_i) = D d_i m - m D_i = sum_e D_e (m_i - e_i) z^{e + m - e_i}.
    Stores a CSR skeleton (rows x ncols) with per-entry (eidx, fac) so that Omega(D) for any D vector is
    data = D[eidx]*fac, and a second skeleton Psi for the transposed role (rows x nD) with
    data = fac*y[col] giving Psi_y with Psi_y D = Omega(D) y.
    """
    def __init__(self, Ds: Space, Ns: list, dst: Space):
        rows_l, cols_l, eidx_l, fac_l = [], [], [], []
        base = 0
        self.colbase = []
        for i in range(NV):
            self.colbase.append(base)
            Ni = Ns[i]
            if Ni.n == 0 or Ds.n == 0:
                base += Ni.n
                continue
            codes = Ds.codes[:, None] + Ni.codes[None, :] - POW[i]     # nD x nN
            fac = (Ni.E[None, :, i] - Ds.E[:, None, i])                # m_i - e_i
            mask = fac != 0
            r = dst.lookup(codes[mask])
            dI, nI = np.nonzero(mask)
            rows_l.append(r); cols_l.append(nI + base); eidx_l.append(dI); fac_l.append(fac[mask])
            base += Ni.n
        self.ncols = base
        self.nrows = dst.n
        self.nD = Ds.n
        rows = np.concatenate(rows_l); cols = np.concatenate(cols_l)
        eidx = np.concatenate(eidx_l); fac = np.concatenate(fac_l).astype(np.float64)
        self.nnz = len(rows)
        # skeleton 1: Omega (rows x ncols), sorted by (row, col)
        perm = np.lexsort((cols, rows))
        self.ind1 = cols[perm].astype(np.int32)
        self.indptr1 = np.zeros(self.nrows + 1, dtype=np.int64)
        np.add.at(self.indptr1, rows[perm] + 1, 1)
        self.indptr1 = np.cumsum(self.indptr1)
        self.eidx1 = eidx[perm].astype(np.int32)
        self.fac1 = fac[perm]
        # skeleton 2: Psi (rows x nD), sorted by (row, eidx)
        perm2 = np.lexsort((eidx, rows))
        self.ind2 = eidx[perm2].astype(np.int32)
        self.indptr2 = np.zeros(self.nrows + 1, dtype=np.int64)
        np.add.at(self.indptr2, rows[perm2] + 1, 1)
        self.indptr2 = np.cumsum(self.indptr2)
        self.col2 = cols[perm2].astype(np.int32)
        self.fac2 = fac[perm2]

    def omega(self, Dvec):
        data = Dvec[self.eidx1] * self.fac1
        return sp.csr_matrix((data, self.ind1, self.indptr1), shape=(self.nrows, self.ncols))

    def psi(self, yvec):
        data = self.fac2 * yvec[self.col2]
        return sp.csr_matrix((data, self.ind2, self.indptr2), shape=(self.nrows, self.nD))


class E2Family:
    """E2-reduced variable-projection family.

    Dfull: (nD x K) array/sparse, columns = basis denominators (coefficients over Ds).
    Nfull: (nNtot x K) array/sparse, columns = particular E1 solutions N_k (over the N blocks).
    """
    def __init__(self, d, pD, Dfull, Nfull, name='', verbose=True, cond=1e-10, want_U=True, want_C=False,
                 threads_note=''):
        t0 = time.time()
        self.d, self.pD, self.name, self.cond = d, pD, name, cond
        self.pars = [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(NV)]
        self.Ds = Space(d, pD)
        self.Ns = [Space(d + 1, self.pars[i]) for i in range(NV)]
        self.Rs = Space(2 * d, 1)
        self.Qs = Space(2 * d - 6, pA)
        self.rows = self.Rs.n
        self.pat = OmegaPattern(self.Ds, self.Ns, self.Rs)
        self.nN = self.pat.ncols
        # Koszul matrix
        blocks = []
        self.Cs = {}
        cbase = 0
        kr, kc, kv = [], [], []
        for i in range(NV):
            for j in range(i + 1, NV):
                pc = (self.pars[i] - (0 if j < 3 else 1)) % 2
                assert pc == (self.pars[j] - (0 if i < 3 else 1)) % 2
                Cs = Space(d - 4, pc)
                self.Cs[(i, j)] = (cbase, Cs)
                if Cs.n:
                    Bi = mul_matrix(Fi[j], Cs, self.Ns[i]).tocoo()
                    Bj = mul_matrix(Fi[i], Cs, self.Ns[j]).tocoo()
                    kr += [Bi.row + self.pat.colbase[i], Bj.row + self.pat.colbase[j]]
                    kc += [Bi.col + cbase, Bj.col + cbase]
                    kv += [Bi.data, -Bj.data]
                cbase += Cs.n
        self.nc = cbase
        self.Kmat = sp.csr_matrix((np.concatenate(kv), (np.concatenate(kr), np.concatenate(kc))),
                                  shape=(self.nN, self.nc)) if kv else sp.csr_matrix((self.nN, 0))
        self.TF = mul_matrix(F, self.Qs, self.Rs)
        self.nQ = self.Qs.n
        self.Dfull = Dfull
        self.Nfull = Nfull
        self.K = Dfull.shape[1]
        # spurious subspace Z = {a : Dfull a in F * S_{d-6}}
        self.Gs = Space(d - 6, pD)
        if self.Gs.n:
            TFd = mul_matrix(F, self.Gs, self.Ds)
            M = np.hstack([_dense(Dfull), -TFd.toarray()])
            ns = sla.null_space(M, rcond=1e-11)
            Za = ns[:self.K, :]
        else:
            Za = sla.null_space(_dense(Dfull), rcond=1e-11)
        if Za.size:
            Uz, sz, _ = np.linalg.svd(Za, full_matrices=True)
            rz = int(np.sum(sz > 1e-9 * max(1.0, sz[0])))
            self.P = Uz[:, rz:]
            self.dimZ = rz
        else:
            self.P = np.eye(self.K); self.dimZ = 0
        self.dimP = self.P.shape[1]
        # Q-elimination projector
        self.U = None
        if want_U and self.nQ:
            Qf, _ = sla.qr(self.TF.toarray(), mode='economic')
            self.U = np.ascontiguousarray(Qf)
        self.C = None
        if want_C:
            self.C = self.r_projector()
        self.build_time = time.time() - t0
        if verbose:
            print(f'[E2Family {name}] d={d} pD={pD} K={self.K} (spurious dim {self.dimZ}, free {self.dimP}) '
                  f'rows={self.rows} nN={self.nN} nc={self.nc} nQ={self.nQ} nnz(Omega)={self.pat.nnz} '
                  f'built in {self.build_time:.1f}s', file=sys.stderr, flush=True)

    # ------------------------------------------------------------------ exact-ish helpers
    def r_projector(self, cache=True):
        """Orthonormal basis C (rows x dim R_{2d}^-) of the orthogonal complement of J(F)_{2d}^- (floats)."""
        path = f'{WD}/w2_Cproj_d{self.d}.npy'
        if cache and os.path.exists(path):
            return np.load(path)
        t0 = time.time()
        gens = []
        for i in range(NV):
            pm = (1 - (0 if i < 3 else 1)) % 2
            Ms = Space(2 * self.d - 5, pm)
            gens.append(mul_matrix(Fi[i], Ms, self.Rs))
        G = sp.hstack(gens).tocsr()
        ns = sla.null_space(G.T.toarray(), rcond=1e-10)
        if cache:
            np.save(path, ns)
        print(f'[r_projector] d={self.d} dim={ns.shape[1]} ({time.time()-t0:.1f}s)', file=sys.stderr, flush=True)
        return ns

    # ------------------------------------------------------------------ core evaluation
    def a_of_u(self, u):
        af = self.P @ u
        na = np.linalg.norm(af)
        return af / na, na

    def linear_solve(self, a):
        """Given a (unit), returns dict with D, Omega, y0=N_0, b, c, r (full residual after optimal c,Q), and
        the SVD pieces needed for the Jacobian."""
        D = _col(self.Dfull @ a)
        Om = self.pat.omega(D)
        y0 = _col(self.Nfull @ a)
        b = Om @ y0
        out = {'D': D, 'Om': Om, 'y0': y0, 'b': b}
        bt = b - self.U @ (self.U.T @ b) if self.U is not None else b.copy()
        if self.nc:
            B = (Om @ self.Kmat)
            B = B.toarray() if sp.issparse(B) else B
            Bt = B - self.U @ (self.U.T @ B) if self.U is not None else B
            Ub, s, Vt = sla.svd(Bt, full_matrices=False, check_finite=False, lapack_driver='gesdd')
            k = int(np.sum(s > self.cond * s[0])) if s.size else 0
            Ubk = Ub[:, :k]
            coef = Ubk.T @ bt
            c = -(Vt[:k].T @ (coef / s[:k]))
            r = bt - Ubk @ coef
            out.update({'c': c, 'Ubk': Ubk, 'rankB': k, 'sB': s})
        else:
            c = np.zeros(0); r = bt
            out.update({'c': c, 'Ubk': np.zeros((self.rows, 0)), 'rankB': 0})
        out['r'] = r
        out['bnorm'] = float(np.linalg.norm(b))
        return out

    def value(self, a):
        o = self.linear_solve(a)
        f = float(o['r'] @ o['r'])
        return f, o

    def residual_u(self, u):
        a, _ = self.a_of_u(u)
        o = self.linear_solve(a)
        self._last = (np.array(u, copy=True), o)
        return o['r']

    def jac_a(self, a, o):
        """Kaufman Jacobian of the reduced residual w.r.t. a (rows x K)."""
        y = o['y0'] + (self.Kmat @ o['c'] if self.nc else 0.0)
        Psi = self.pat.psi(y)
        G1 = Psi @ _dense(self.Dfull)
        G2 = o['Om'] @ self.Nfull
        G2 = G2.toarray() if sp.issparse(G2) else G2
        J = G1 + G2
        if self.U is not None:
            J = J - self.U @ (self.U.T @ J)
        if o['rankB']:
            J = J - o['Ubk'] @ (o['Ubk'].T @ J)
        return J

    def jac_u(self, u):
        a, na = self.a_of_u(u)
        last = getattr(self, '_last', None)
        if last is not None and last[0].shape == np.shape(u) and np.array_equal(last[0], u):
            o = last[1]
        else:
            o = self.linear_solve(a)
        J = self.jac_a(a, o)
        Ju = (J - np.outer(J @ a, a)) @ self.P / na
        return Ju

    # ------------------------------------------------------------------ R-part (Hessian relaxation)
    def r_residual_u(self, u):
        a, _ = self.a_of_u(u)
        D = _col(self.Dfull @ a); Om = self.pat.omega(D); y0 = _col(self.Nfull @ a)
        return self.C.T @ (Om @ y0)

    def r_jac_u(self, u):
        a, na = self.a_of_u(u)
        D = _col(self.Dfull @ a); Om = self.pat.omega(D); y0 = _col(self.Nfull @ a)
        Psi = self.pat.psi(y0)
        G = Psi @ _dense(self.Dfull)
        G2 = Om @ self.Nfull
        G2 = G2.toarray() if sp.issparse(G2) else G2
        J = self.C.T @ (G + G2)
        return (J - np.outer(J @ a, a)) @ self.P / na

    # ------------------------------------------------------------------ searches
    def lm(self, u0, max_nfev=200, stage1=False, xtol=1e-15, ftol=1e-15, gtol=1e-15):
        fun = self.r_residual_u if stage1 else self.residual_u
        jac = self.r_jac_u if stage1 else self.jac_u
        m = (self.C.shape[1] if stage1 else self.rows)
        method = 'lm' if m >= self.dimP else 'trf'
        res = least_squares(fun, u0, jac=jac, method=method, max_nfev=max_nfev, xtol=xtol, ftol=ftol, gtol=gtol)
        return res

    def restart(self, seed, max_nfev=200, stage1_first=False, stage1_nfev=100, init=None):
        rng = np.random.default_rng(seed)
        u0 = rng.standard_normal(self.dimP) if init is None else np.asarray(init, float)
        u0 /= np.linalg.norm(u0)
        t0 = time.time()
        info = {'seed': int(seed)}
        if stage1_first and self.C is not None:
            r1 = self.lm(u0, max_nfev=stage1_nfev, stage1=True)
            u0 = r1.x / np.linalg.norm(r1.x)
            info['stage1_cost'] = float(2 * r1.cost); info['stage1_nfev'] = int(r1.nfev)
        res = self.lm(u0, max_nfev=max_nfev)
        a, _ = self.a_of_u(res.x)
        f, o = self.value(a)
        bn = o['bnorm']
        rel = math.sqrt(f) / max(bn, 1e-300)
        rpart = float(np.linalg.norm(self.C.T @ o['b'])) / max(bn, 1e-300) if self.C is not None else None
        info.update({'f': f, 'res': math.sqrt(f), 'bnorm': bn, 'rel': rel, 'r_part_rel': rpart,
                     'nfev': int(res.nfev), 'status': int(res.status), 'seconds': round(time.time() - t0, 2),
                     'rankB': int(o['rankB'])})
        return info, a, o

    # ------------------------------------------------------------------ exact confirmation
    def D_poly_float(self, a):
        return self.Ds.poly(_col(self.Dfull @ a))

    def confirm_rational(self, a, max_dens=(1, 10, 100, 1000, 10**4, 10**5, 10**6), use_w=True, log_path=None,
                         note=''):
        """Rational reconstruction of a (if use_w and the family is Jacobian-type the w-vector is reconstructed,
        else the D coefficients), then exact E1 + fixed-D solve.  Returns list of attempt records."""
        attempts = []
        vecs = []
        a = np.asarray(a, float)
        if use_w and getattr(self, 'is_jacobian', False):
            v = a / a[np.argmax(np.abs(a))]
            vecs.append(('w', v))
        D = _col(self.Dfull @ a)
        vecs.append(('D', D / D[np.argmax(np.abs(D))]))
        tried = set()
        for kind, v in vecs:
            for md in max_dens:
                fr = [Fraction(float(x)).limit_denominator(md) for x in v]
                err = max(abs(float(f) - x) for f, x in zip(fr, v))
                key = (kind, tuple(fr))
                if key in tried:
                    continue
                tried.add(key)
                if kind == 'w':
                    Dex = self.w_to_D(fr)
                else:
                    Dex = {m: fr[i] for i, m in enumerate(self.Ds.mons) if fr[i] != 0}
                recd = {'kind': kind, 'max_den': md, 'recon_err': err, 'note': note}
                if not Dex:
                    recd['status'] = 'zero'; attempts.append(recd); continue
                ok, info = fixed_d.e1_solvable(F, A, q, Dex)
                recd['e1'] = bool(ok)
                if ok:
                    divF, _ = jacobian.divisible_by(F, Dex)
                    recd['F_divides_D'] = bool(divF)
                    if not divF:
                        t0 = time.time()
                        out = fixed_d.solve_fixed_D_fast(F, A, q, Dex, want_exact=True)
                        recd.update({'consistent': bool(out['consistent']), 'rank': out['rank'], 'nullity': out['nullity'],
                                     'seconds': round(time.time() - t0, 1)})
                        if out['consistent'] and out.get('exact'):
                            okv, vinfo = core.verify_witness(F, A, q, Dex, out['N'], out['Q'])
                            recd['verified'] = bool(okv)
                            if okv:
                                w = {'instance': IID, 'd': core.degree(Dex), 'D': core.to_str(Dex),
                                     'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q']),
                                     'family': self.name, 'engine': fixed_d.ENGINE_VERSION}
                                with open(f'{WD}/witness-candidate.json', 'w') as fh:
                                    json.dump(w, fh, indent=1)
                                recd['WITNESS'] = True
                                print('!!!! WITNESS FOUND AND VERIFIED', file=sys.stderr, flush=True)
                attempts.append(recd)
                if log_path:
                    jlog(log_path, {'event': 'reconstruct', **recd})
                if recd.get('WITNESS'):
                    return attempts
        return attempts

    def w_to_D(self, fr):
        raise NotImplementedError


def _dense(M):
    return M.toarray() if sp.issparse(M) else np.asarray(M)


def _col(v):
    v = np.asarray(v)
    return v.ravel() if v.ndim > 1 else v


# ---------------------------------------------------------------------------------------------------
# Jacobian-type family: D = w(F), N_0 = A w
# ---------------------------------------------------------------------------------------------------
class JacobianFamily(E2Family):
    def __init__(self, d, pD, verbose=True, **kw):
        self.Ws = []
        wb = 0
        self.wbase = []
        for i in range(NV):
            pw = (pD - (0 if i < 3 else 1)) % 2
            Ws = Space(d - 5, pw)
            self.Ws.append(Ws); self.wbase.append(wb); wb += Ws.n
        self.nW = wb
        Ds = Space(d, pD)
        pars = [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(NV)]
        Ns = [Space(d + 1, pars[i]) for i in range(NV)]
        dcols, acols = [], []
        for i in range(NV):
            dcols.append(mul_matrix(Fi[i], self.Ws[i], Ds))
            acols.append(mul_matrix(A, self.Ws[i], Ns[i]))
        Dfull = sp.hstack(dcols).tocsr()
        Nfull = sp.block_diag(acols).tocsr()
        super().__init__(d, pD, Dfull, Nfull, name=f'jacobian-d{d}-p{pD}', verbose=verbose, **kw)
        self.is_jacobian = True

    def w_to_D(self, fr):
        D = {}
        for i in range(NV):
            for k, m in enumerate(self.Ws[i].mons):
                c = fr[self.wbase[i] + k]
                if c:
                    D = add(D, scale(shift(Fi[i], m), Fraction(c)))
        return D

    def w_polys(self, a):
        return [self.Ws[i].poly(a[self.wbase[i]:self.wbase[i] + self.Ws[i].n]) for i in range(NV)]


# ---------------------------------------------------------------------------------------------------
# Kernel family from a precomputed float basis (npz with Dfull, Nfull)
# ---------------------------------------------------------------------------------------------------
class KernelFamily(E2Family):
    def __init__(self, d, pD, npz_path, verbose=True, **kw):
        z = np.load(npz_path)
        super().__init__(d, pD, z['Dfull'], z['Nfull'], name=f'kernel-d{d}-p{pD}', verbose=verbose, **kw)
        self.is_jacobian = False


# ---------------------------------------------------------------------------------------------------
# Exact (mod p) Hessian-relaxation quadrics for the Jacobian family
# ---------------------------------------------------------------------------------------------------
def hess_quadrics(d, pD, p=PRIME, verbose=True):
    """Quadrics in the w-coefficients (monomial-field basis, ordered as JacobianFamily.wbase) cutting out
    {w : A*Hess_F(w,w) in J(F)_{2d}} mod p.  Returns dict with 'Qmat' (nquad x npairs) over F_p,
    'pairs' (list of (k,l)), 'nW', ranks, timings."""
    from flint import nmod_mat
    t0 = time.time()
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
    Fij = [[pderiv(pderiv(Fp, i, p), j, p) for j in range(NV)] for i in range(NV)]
    Wm = []  # (i, m)
    for i in range(NV):
        pw = None if pD is None else (pD - (0 if i < 3 else 1)) % 2
        for m in monomials(d - 5, pw):
            Wm.append((i, m))
    nW = len(Wm)
    Rs = Space(2 * d, None if pD is None else 1)
    # left null space of J_{2d} (odd part, or both parities for mixed D): c . (m F_i) = 0 for all generators
    gens = jacobian.jacobian_generators(Fp, 2 * d, p, parity=(None if pD is None else 1))
    ng = len(gens)
    flat = [0] * (ng * Rs.n)
    for r, (_, _, P) in enumerate(gens):
        base = r * Rs.n
        for e, v in P.items():
            flat[base + Rs.idx[e]] = v
    G = nmod_mat(ng, Rs.n, flat, p)        # generators as rows
    t1 = time.time()
    Cn, dimR = G.nullspace()               # nullspace of G: vectors x with G x = 0  -> functionals on R_{2d}
    t2 = time.time()
    # pair vectors v_kl = A m_k m_l F_{i_k i_l}, processed in chunks of pairs (d=9 has 71k pairs)
    pairs = [(k, l) for k in range(nW) for l in range(k, nW)]
    Cm = nmod_mat(dimR, Rs.n, [int(Cn[i, j]) for j in range(dimR) for i in range(Rs.n)], p)  # dimR x rows
    AF = {}
    Qcols = []
    chunk = max(1, min(len(pairs), 4000000 // max(Rs.n, 1)))
    for c0 in range(0, len(pairs), chunk):
        sub_pairs = pairs[c0:c0 + chunk]
        V = [0] * (Rs.n * len(sub_pairs))
        for j, (k, l) in enumerate(sub_pairs):
            (i1, m1), (i2, m2) = Wm[k], Wm[l]
            key = (i1, i2) if i1 <= i2 else (i2, i1)
            if key not in AF:
                AF[key] = mul(Ap, Fij[key[0]][key[1]], p)
            P = shift(AF[key], tuple(a + b for a, b in zip(m1, m2)), 1 if k == l else 2, p)   # w_k w_l appears twice for k<l
            for e, v in P.items():
                V[Rs.idx[e] * len(sub_pairs) + j] = v
        Vm = nmod_mat(Rs.n, len(sub_pairs), V, p)
        Qc = Cm * Vm
        Qcols.append([[int(Qc[i, j]) for j in range(len(sub_pairs))] for i in range(dimR)])
    flatQ = []
    for i in range(dimR):
        for part in Qcols:
            flatQ.extend(part[i])
    Qm = nmod_mat(dimR, len(pairs), flatQ, p)   # dimR x npairs : quadric coefficients (k<=l, coeff of w_k w_l)
    rank2 = Qm.rank()
    t3 = time.time()
    out = {'d': d, 'pD': pD, 'p': p, 'nW': nW, 'npairs': len(pairs), 'dimR_2d_odd': dimR, 'rank_quadrics': rank2,
           'Wm': Wm, 'pairs': pairs, 'Qm': Qm, 't_null': round(t2 - t1, 1), 't_quad': round(t3 - t2, 1)}
    if verbose:
        print(f'[hess_quadrics] d={d} pD={pD} nW={nW} pairs={len(pairs)} dim R_2d^-={dimR} rank(quadrics)={rank2} '
              f'({t3-t0:.1f}s)', file=sys.stderr, flush=True)
    return out


def quadric_rows(res):
    """Reduced row basis of the quadric space as list of dict {(k,l): coeff mod p}."""
    Qm = res['Qm']; pairs = res['pairs']
    Qr, rk = Qm.rref()
    rows = []
    for i in range(rk):
        row = {}
        for j in range(len(pairs)):
            v = int(Qr[i, j])
            if v:
                row[pairs[j]] = v
        rows.append(row)
    return rows


def write_msolve(res, path):
    """Write the (reduced) quadric system in msolve input format over F_p."""
    rows = quadric_rows(res)
    nW = res['nW']; p = res['p']
    names = [f'w{k}' for k in range(nW)]
    with open(path, 'w') as f:
        f.write(','.join(names) + '\n')
        f.write(f'{p}\n')
        polys = []
        for row in rows:
            terms = []
            for (k, l), v in row.items():
                mon = f'w{k}^2' if k == l else f'w{k}*w{l}'
                terms.append(f'{v}*{mon}')
            polys.append('+'.join(terms))
        f.write(',\n'.join(polys) + '\n')
    return len(rows)


def macaulay_rank(rows, nW, deg, p):
    """Rank mod p of the degree-`deg` Macaulay matrix of the quadrics (homogeneous ideal, deg >= 2)."""
    from flint import nmod_mat
    import itertools
    mons = list(itertools.combinations_with_replacement(range(nW), deg))
    idx = {m: i for i, m in enumerate(mons)}
    mult = list(itertools.combinations_with_replacement(range(nW), deg - 2))
    R, C = len(rows) * len(mult), len(mons)
    flat = [0] * (R * C)
    r = 0
    for row in rows:
        for mu in mult:
            for (k, l), v in row.items():
                m = tuple(sorted(mu + (k, l)))
                flat[r * C + idx[m]] = (flat[r * C + idx[m]] + v) % p
            r += 1
    M = nmod_mat(R, C, flat, p)
    return M.rank(), R, C


# ---------------------------------------------------------------------------------------------------
# exact E2 residual for self-tests (Fractions)
# ---------------------------------------------------------------------------------------------------
def beta_exact(D, N):
    divN = {}; ND = {}
    for i in range(NV):
        if N[i]:
            divN = add(divN, pderiv(N[i], i))
            ND = add(ND, mul(N[i], pderiv(D, i)))
    return sub(mul(D, divN), ND)
