"""w2lib — wave-2 search library for the off-paired lane.

Sector-restricted variable projection (VarPro) over EXACT rational bases of linear denominator families
(colon space (J(F):A)_{d,chi}, Jacobian part J_{d,chi}, sub-families), matrix-free gradient, multiprocessing
driver with per-seed JSONL checkpoints, rational reconstruction + exact sector-reduced re-solve.

Character bookkeeping comes from the wave-1 oplib (generators g = (0,2,0,1,3,3), h = (0,0,1,0,0,5) mod 6).
Engine: /home/kejace/devel/sextics-work/hodge (not modified).
"""
import sys, os, time, json, math, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, scipy.linalg as sla, scipy.sparse as sp
from scipy.optimize import minimize
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from flint import nmod_mat
import oplib

LANE = HERE
char, cadd, cneg, VARCHAR = oplib.char, oplib.cadd, oplib.cneg, oplib.VARCHAR
P0 = core.PRIMES[0]
BK = 64
POW = np.array([BK ** i for i in range(6)], dtype=np.int64)


def monos_char(deg, c):
    return oplib.monos_char(deg, c)


def key(e):
    return int(sum(int(a) * int(b) for a, b in zip(e, POW)))


def keys_of(mons):
    return (np.array(mons, dtype=np.int64) @ POW) if len(mons) else np.zeros(0, np.int64)


def poly_arrays(P):
    ks = np.array([key(e) for e in P], dtype=np.int64)
    vs = np.array([float(v) for v in P.values()], dtype=float)
    return ks, vs


def load():
    F, A, q, rec = core.load_instance('off-paired')
    return F, A, q


# ----------------------------------------------------------------------------------------------------------
# exact rational row-space bases via modular rref + CRT + rational reconstruction (+ third-prime verification)
# ----------------------------------------------------------------------------------------------------------
def _rref_rows_modp(vecs, ncols, p):
    """vecs: list of lists of ints (mod p). Returns (rank, rref rows as lists, pivots)."""
    if not vecs:
        return 0, [], []
    M = nmod_mat(len(vecs), ncols, [v % p for row in vecs for v in row], p)
    Mr, rk = M.rref()
    rows, piv = [], []
    for i in range(rk):
        row = [int(Mr[i, j]) for j in range(ncols)]
        j = 0
        while j < ncols and row[j] == 0:
            j += 1
        piv.append(j); rows.append(row)
    return rk, rows, piv


def exact_rowspace_basis(vec_fn, ncols, primes=None, verbose=False, tag=''):
    """vec_fn(p) -> list of integer row vectors mod p spanning the same Q-space reduced mod p.
    Returns (list of Fraction rows in rref, pivots)."""
    primes = primes or core.PRIMES[:6]
    acc = None; mod = None; piv0 = None
    used = 0
    for p in primes:
        rk, rows, piv = _rref_rows_modp(vec_fn(p), ncols, p)
        used += 1
        if acc is None:
            acc, mod, piv0 = rows, p, piv
        else:
            if piv != piv0:
                raise RuntimeError(f'{tag}: pivot mismatch between primes')
            acc = [[core.crt_pair(a, mod, b, p)[0] for a, b in zip(ra, rb)] for ra, rb in zip(acc, rows)]
            mod *= p
        rat = [[core.ratrecon(v, mod) for v in row] for row in acc]
        if all(v is not None for row in rat for v in row):
            # verify with the next prime
            pv = primes[used] if used < len(primes) else None
            if pv is None:
                return rat, piv0
            rkv, rowsv, pivv = _rref_rows_modp(vec_fn(pv), ncols, pv)
            if pivv != piv0:
                raise RuntimeError(f'{tag}: pivot mismatch at verification prime')
            ok = True
            for ra, rb in zip(rat, rowsv):
                for v, b in zip(ra, rb):
                    num, den = v.numerator, v.denominator
                    if (num * pow(den, -1, pv) - b) % pv:
                        ok = False; break
                if not ok:
                    break
            if ok:
                if verbose:
                    print(f'[exact basis {tag}] rank {len(rat)} with {used} primes (verified)', file=sys.stderr)
                return rat, piv0
            # else continue accumulating
    raise RuntimeError(f'{tag}: rational reconstruction failed with {len(primes)} primes')


def colon_vecs_modp(F, A, q, d, c, Dpolys_or_mons, p):
    """Rows (in the coordinates of the given D-basis) spanning {D in span : A*D in J(F)} mod p.
    Dpolys_or_mons: list of exact polys (the D-family basis)."""
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    cA = oplib.char_of(Ap)
    degN = d + 6 * q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    cols = []; kinds = []
    for i in range(NV):
        ci = cadd(cadd(c, VARCHAR[i]), cA)
        for m in monos_char(degN, ci):
            cols.append(shift(Fi[i], m, q, p)); kinds.append(('N', i, m))
    nW = len(Dpolys_or_mons)
    for k, Bk in enumerate(Dpolys_or_mons):
        Bp = core.reduce_mod(Bk, p)
        cols.append(scale(mul(Ap, Bp, p), p - 1, p)); kinds.append(('D', k, None))
    rows = {}
    for col in cols:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    if C == 0 or nW == 0:
        return []
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    if nullity == 0:
        return []
    didx = [j for j, k in enumerate(kinds) if k[0] == 'D']
    vecs = [[int(X[j, kk]) for j in didx] for kk in range(nullity)]
    return vecs


def exact_colon_basis(F, A, q, d, c, family=None, verbose=False):
    """Exact rational rref basis of (J:A)_{d,c} (family=None: all monomials of degree d, character c) or of
    {D in span(family) : A D in J}.  Returns (polys, coord_vectors, Dbasis) where polys are exact dict polys."""
    Dbasis = family if family is not None else [{m: Fraction(1)} for m in monos_char(d, c)]
    nW = len(Dbasis)
    if nW == 0:
        return [], [], Dbasis
    rat, piv = exact_rowspace_basis(lambda p: colon_vecs_modp(F, A, q, d, c, Dbasis, p), nW, verbose=verbose,
                                    tag=f'colon d={d} c={c}')
    polys = []
    for row in rat:
        P = {}
        for k, v in enumerate(row):
            if v:
                P = add(P, scale(Dbasis[k], v))
        polys.append(P)
    return polys, rat, Dbasis


def jacobian_generators_exact(F, d, c):
    """m*F_i with deg m = d-5, char(m) = c + VARCHAR[i] (so that the product has character c)."""
    Fi = [pderiv(F, i) for i in range(NV)]
    out = []
    if d < 5:
        return out
    for i in range(NV):
        for m in monos_char(d - 5, cadd(c, VARCHAR[i])):
            out.append(shift(Fi[i], m, 1))
    return out


def exact_span_basis(polys, mons, verbose=False, tag='span'):
    """Exact rref basis (as polys) of the span of exact polys in the monomial list mons."""
    idx = {m: i for i, m in enumerate(mons)}
    if not polys:
        return []

    def vec_fn(p):
        vecs = []
        for P in polys:
            Pp = core.reduce_mod(P, p)
            row = [0] * len(mons)
            for e, v in Pp.items():
                row[idx[e]] = v
            vecs.append(row)
        return vecs
    rat, piv = exact_rowspace_basis(vec_fn, len(mons), verbose=verbose, tag=tag)
    return [{mons[j]: v for j, v in enumerate(row) if v} for row in rat]


def gauge_polys(F, d, c):
    """F * m for monomials m of degree d-6 and character c (F has character (0,0))."""
    if d < 6:
        return []
    return [shift(F, m, 1) for m in monos_char(d - 6, c)]


# ----------------------------------------------------------------------------------------------------------
# sector-restricted matrix-free VarPro
# ----------------------------------------------------------------------------------------------------------
class SectorVarPro:
    """Family D(a) = sum_k a_k B_k (exact polys B_k of degree d and character c), fixed-D system restricted to the
    character sector; f(a) = min_{N,Q} ||M(a) x - b(a)||^2 with b = -A D^2, evaluated on the unit sphere of the
    complement of the gauge subspace span(F*R_{d-6}) (so D is never an F-multiple)."""

    def __init__(self, F, A, q, basis, c, gauge=None, verbose=False, cond=1e-9, remove_gauge=True):
        t0 = time.time()
        self.F, self.A, self.q, self.c = F, A, q, c
        self.basis = basis
        K = len(basis); self.K = K
        d = core.degree(basis[0]); self.d = d
        for Bk in basis:
            assert core.degree(Bk) == d and oplib.char_of(Bk) == c, 'basis must be homogeneous of one character'
        cA = oplib.char_of(A); self.cA = cA
        self.cond = cond
        Dm = monos_char(d, c); self.Dm = Dm; nD = len(Dm); self.nD = nD
        Didx = {m: i for i, m in enumerate(Dm)}
        Bex = np.zeros((K, nD))
        for k, Bk in enumerate(basis):
            for e, v in Bk.items():
                Bex[k, Didx[e]] = float(v)
        self.Bex = Bex
        # orthonormal basis of the row space of Bex (in D-monomial coordinates)
        U, s, Vt = np.linalg.svd(Bex, full_matrices=False)
        r = int(np.sum(s > 1e-10 * s[0]))
        self.Borth = Vt[:r, :]              # r x nD, orthonormal rows
        self.rankB = r
        # gauge complement
        gauge = gauge if gauge is not None else (gauge_polys(F, d, c) if remove_gauge else [])
        ng = 0
        if gauge:
            G = np.zeros((len(gauge), nD))
            for k, Gk in enumerate(gauge):
                for e, v in Gk.items():
                    G[k, Didx[e]] = float(v)
            Gc = G @ self.Borth.T          # coordinates in the orthonormal basis (gauge subset of the family)
            resid = np.linalg.norm(G - Gc @ self.Borth) / max(np.linalg.norm(G), 1e-300)
            if resid > 1e-8:
                raise ValueError(f'gauge vectors not inside the family (rel {resid:.1e})')
            Ug, sg, Vgt = np.linalg.svd(Gc, full_matrices=True)
            ng = int(np.sum(sg > 1e-10 * max(sg[0], 1e-300)))
            self.P = Vgt[ng:, :].T          # r x (r-ng): orthonormal complement of the gauge in orth coords
        else:
            self.P = np.eye(r)
        self.ng = ng
        self.dimP = self.P.shape[1]
        degN = d + 6 * q - 5; degQ = 2 * d - 6
        self.degN, self.degQ = degN, degQ
        crow = cadd(cA, cadd(c, c))
        rowmons = monos_char(2 * d + 6, crow)
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
        # N columns grouped by i
        Ncols = []; self.Ncols_by_i = []
        for i in range(NV):
            ci = cadd(cadd(c, VARCHAR[i]), cA)
            ms = monos_char(degN, ci)
            self.Ncols_by_i.append((len(Ncols), ms))
            Ncols += [(i, m) for m in ms]
        nN = len(Ncols); self.nN = nN; self.Ncols = Ncols
        # T_N: sparse (R*nN) x nD
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
                    rows = row_index(kk)
                    ri_l.append(rows * nN + np.broadcast_to(cols[None, :], (len(gk), len(ms))).ravel())
                    ci_l.append(np.full(kk.size, mu_idx, dtype=np.int64))
                    va_l.append(np.broadcast_to(gv[:, None], (len(gk), len(ms))).ravel().copy())
                # F*mu * m_i * (m / z_i)
                mi = np.array([m[i] for m in ms], dtype=np.int64)
                sel = np.nonzero(mi)[0]
                if sel.size:
                    mk2 = mk[sel] - BK ** i
                    kk = (fk[:, None] + mk2[None, :]).ravel()
                    rows = row_index(kk)
                    ri_l.append(rows * nN + np.broadcast_to(cols[sel][None, :], (len(fk), sel.size)).ravel())
                    ci_l.append(np.full(kk.size, mu_idx, dtype=np.int64))
                    va_l.append((fv[:, None] * mi[sel][None, :].astype(float)).ravel())
        ri = np.concatenate(ri_l); ci = np.concatenate(ci_l); va = np.concatenate(va_l)
        self.TN = sp.csr_matrix((va, (ri, ci)), shape=(R * nN, nD))
        self.TNT = self.TN.T.tocsr()
        # Q block (constant) and its orthonormal column basis
        Qm = monos_char(degQ, crow) if degQ >= 0 else []
        self.Qm = Qm; nQ = len(Qm); self.nQ = nQ
        if nQ:
            F2n = scale(mul(F, F), -1)
            f2k, f2v = poly_arrays(F2n)
            MQ = np.zeros((R, nQ))
            qk = keys_of(Qm)
            kk = (f2k[:, None] + qk[None, :]).ravel()
            rows = row_index(kk)
            cols = np.broadcast_to(np.arange(nQ)[None, :], (len(f2k), nQ)).ravel()
            np.add.at(MQ, (rows, cols), np.broadcast_to(f2v[:, None], (len(f2k), nQ)).ravel())
            self.MQ = MQ
            Qq, Rq = np.linalg.qr(MQ, mode='reduced')
            self.QQ = Qq
        else:
            self.MQ = np.zeros((R, 0)); self.QQ = np.zeros((R, 0))
        # b pairs: -A mu nu
        kA = key(next(iter(A)))
        assert len(A) == 1, 'A must be a monomial for the pair table (true for X0^3 Y0^3)'
        self.cAval = float(A[next(iter(A))])
        dk = keys_of(Dm)
        mu_i, nu_i = np.meshgrid(np.arange(nD), np.arange(nD), indexing='ij')
        mu_i = mu_i.ravel(); nu_i = nu_i.ravel()
        self.pair_mu, self.pair_nu = mu_i, nu_i
        self.pair_row = row_index(kA + dk[mu_i] + dk[nu_i])
        self.ncols = nN + nQ
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[sector varpro] d={d} c={c} K={K} rank={r} gauge={ng} dimP={self.dimP} rows={R} '
                  f'N={nN} Q={nQ} nnz(TN)={self.TN.nnz} built in {self.build_time:.1f}s', file=sys.stderr)

    # --- evaluation ---------------------------------------------------------------------------------------
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

    # --- search ---------------------------------------------------------------------------------------------
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

    # --- exact confirmation --------------------------------------------------------------------------------
    def D_exact_from_coef(self, Dc, max_den):
        j = int(np.argmax(np.abs(Dc)))
        v = Dc / Dc[j]
        fr = [Fraction(float(x)).limit_denominator(max_den) for x in v]
        err = max(abs(float(a) - b) for a, b in zip(fr, v))
        D = {self.Dm[i]: fr[i] for i in range(self.nD) if fr[i] != 0}
        return D, err

    def D_exact_from_basis_coords(self, Dc, max_den):
        """Reconstruct in the coordinates of the exact rref basis (often sparser)."""
        a, *_ = np.linalg.lstsq(self.Bex.T, Dc, rcond=None)
        j = int(np.argmax(np.abs(a)))
        v = a / a[j]
        fr = [Fraction(float(x)).limit_denominator(max_den) for x in v]
        err = max(abs(float(x) - y) for x, y in zip(fr, v))
        D = {}
        for k, Bk in enumerate(self.basis):
            if fr[k]:
                D = add(D, scale(Bk, fr[k]))
        return D, err


def confirm_exact(F, A, q, D, verbose=False, max_primes=12):
    """E1 (sector-reduced, mod p) then exact sector-reduced multi-prime solve; returns (status, data)."""
    if not D:
        return 'inconsistent', {'note': 'D=0'}
    divF, _ = jacobian.divisible_by(F, D)
    if divF:
        return 'spurious', {'note': 'F divides D'}
    ok, info = oplib.e1_reduced(F, A, q, D)
    if not ok:
        return 'inconsistent', {'note': 'e1 fails', 'e1': info}
    st, data = oplib.lift_solution(F, A, q, D, verbose=verbose, max_primes=max_primes)
    if st == 'solved':
        ok2, info2 = core.verify_witness(F, A, q, D, data['N'], data['Q'])
        data['verified'] = bool(ok2)
    return st, data


def try_reconstruct_and_confirm(fam, Dc, dens=(1, 2, 6, 12, 60, 120, 1000, 10000, 100000, 10 ** 6), tol=1e-5,
                                verbose=False):
    """Returns list of attempt records; a record with status 'solved' is a witness (N, Q included)."""
    out = []
    seen = set()
    for mode in ('coef', 'basis'):
        for md in dens:
            try:
                D, err = (fam.D_exact_from_coef(Dc, md) if mode == 'coef' else fam.D_exact_from_basis_coords(Dc, md))
            except Exception as ex:
                out.append({'mode': mode, 'max_den': md, 'error': str(ex)}); continue
            if err > tol or not D:
                out.append({'mode': mode, 'max_den': md, 'recon_err': err, 'status': 'no-fit'}); continue
            s = core.to_str(D)
            if s in seen:
                continue
            seen.add(s)
            t0 = time.time()
            st, data = confirm_exact(fam.F, fam.A, fam.q, D, verbose=verbose)
            rec = {'mode': mode, 'max_den': md, 'recon_err': err, 'status': st, 'D': s, 'time': round(time.time() - t0, 1)}
            for k in ('rank', 'nullity', 'note', 'verified'):
                if k in data:
                    rec[k] = data[k]
            if st == 'solved':
                rec['N'] = [core.to_str(n) for n in data['N']]; rec['Q'] = core.to_str(data['Q'])
                out.append(rec)
                return out
            out.append(rec)
    return out


# ----------------------------------------------------------------------------------------------------------
# families
# ----------------------------------------------------------------------------------------------------------
def var_subset_family(d, c, vars_):
    """All monomials of degree d, character c, in the given variable subset (cone over a lower-dim variety)."""
    vs = set(vars_)
    return [{m: Fraction(1)} for m in monos_char(d, c) if all(m[i] == 0 for i in range(NV) if i not in vs)]


def sigma2_symmetrize(polys, sign=+1):
    """Span of P +/- sigma2(P) (sigma2: X2 <-> Y2) — only meaningful for h-invariant characters (a,0)."""
    out = []
    for P in polys:
        Ps = {}
        for e, v in P.items():
            e2 = (e[0], e[1], e[5], e[3], e[4], e[2])
            Ps[e2] = Ps.get(e2, 0) + v
        Q = add(P, scale(Ps, sign))
        if Q:
            out.append(Q)
    return out


def build_family(F, A, q, d, c, kind, verbose=False):
    """kind: 'colon' | 'jac' | 'vars:X0,X1,Y0,Y1' | 'colon-s2even' | 'colon-s2odd' | 'jac-s2even' | 'colon-vars:...'
    Returns (basis polys, gauge polys, info)."""
    cA = oplib.char_of(A)
    gauge = gauge_polys(F, d, c)
    info = {'d': d, 'c': list(c), 'kind': kind}
    if kind == 'colon':
        basis, _, _ = exact_colon_basis(F, A, q, d, c, verbose=verbose)
    elif kind == 'jac':
        gens = jacobian_generators_exact(F, d, c)
        basis = exact_span_basis(gens, monos_char(d, c), verbose=verbose, tag=f'jac d={d} c={c}')
    elif kind.startswith('colon-vars:') or kind.startswith('vars:'):
        vs = [core.VARS.index(v) for v in kind.split(':', 1)[1].split(',')]
        fam = var_subset_family(d, c, vs)
        if kind.startswith('colon-vars:'):
            basis, _, _ = exact_colon_basis(F, A, q, d, c, family=fam, verbose=verbose)
        else:
            basis = fam
        # gauge inside this family: F*m with all variables in the subset -> only if F in the subset (never); so the
        # gauge intersection is {0} unless the family is the full sector; use gauge vectors that lie in the span
        gauge = _gauge_inside(basis, gauge)
    elif kind in ('colon-s2even', 'colon-s2odd', 'jac-s2even', 'jac-s2odd'):
        base_kind = 'colon' if kind.startswith('colon') else 'jac'
        sign = +1 if kind.endswith('even') else -1
        b0, g0, _ = build_family(F, A, q, d, c, base_kind, verbose=verbose)
        sym = sigma2_symmetrize(b0, sign)
        basis = exact_span_basis(sym, monos_char(d, c), verbose=verbose, tag=kind) if sym else []
        gauge = exact_span_basis(sigma2_symmetrize(gauge, sign), monos_char(d, c), tag='gauge-s2') if gauge else []
        gauge = _gauge_inside(basis, gauge)
    else:
        raise ValueError(kind)
    info['K'] = len(basis); info['gauge'] = len(gauge)
    return basis, gauge, info


def _gauge_inside(basis, gauge):
    """Keep only the part of span(gauge) that lies inside span(basis) (numerically)."""
    if not basis or not gauge:
        return []
    mons = sorted({e for P in basis + gauge for e in P})
    idx = {m: i for i, m in enumerate(mons)}
    Bm = np.zeros((len(basis), len(mons))); Gm = np.zeros((len(gauge), len(mons)))
    for k, P in enumerate(basis):
        for e, v in P.items(): Bm[k, idx[e]] = float(v)
    for k, P in enumerate(gauge):
        for e, v in P.items(): Gm[k, idx[e]] = float(v)
    U, s, Vt = np.linalg.svd(Bm, full_matrices=False)
    r = int(np.sum(s > 1e-10 * s[0])); Vr = Vt[:r]
    # intersection of span(G) with span(B): vectors g = sum c_j G_j with g (I - Vr^T Vr) = 0
    Proj = Gm - (Gm @ Vr.T) @ Vr
    Uc, sc, Vct = np.linalg.svd(Proj, full_matrices=True)
    rc = int(np.sum(sc > 1e-8 * max(sc[0], 1e-300))) if sc.size else 0
    C = Vct[rc:]  # combinations of gauge rows lying inside span(B)
    out = []
    for crow in C:
        v = crow @ Gm
        if np.linalg.norm(v) < 1e-12:
            continue
        # return as float polys (only used numerically)
        out.append({mons[j]: float(v[j]) for j in range(len(mons)) if abs(v[j]) > 1e-14})
    return out


# ----------------------------------------------------------------------------------------------------------
# logging helpers
# ----------------------------------------------------------------------------------------------------------
def append_json(path, rec):
    with open(path, 'a') as f:
        f.write(json.dumps(rec) + '\n')


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%S')
