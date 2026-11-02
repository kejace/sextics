"""w2k3 — wave-2 joint-search library for the lane polarized-k3-remainder.

Instance: q = 1, A = X0 X1 X2 Y0 Y1 Y2 (Y-odd), F = (f(X+Y) + f(X-Y))/2 (Y-even), diagonal group {1, tau}:
the only character is Y-parity, so the lossless restriction is D Y-even, N_X Y-odd, N_Y Y-even, C Y-odd.

Normalized certificate system (q = 1):   E1:  N(F) = A D            (degree d+6)
                                         E2:  D div N - N(D) = F C  (degree 2d, deg C = 2d-6)

Parametrization (engine-independent, numpy):
  * Jacobian part        D = sum_i w_i F_i  (deg w_i = d-5),  N_0 = A w            (E1 solved identically)
  * non-Jacobian part    D = sum_j a_j k_j   (E1-kernel representatives, d >= 10) with lifts L_j, L_j(F) = A k_j
  * every other E1 lift  N = N_0 + K(c),  K(c) = sum_{i<j} c_ij (F_j e_i - F_i e_j),  deg c_ij = d-4
  * residual             r(a; c, C) = B(a) + T(a) c - F C,  B(a) = D(a) div N_0(a) - N_0(a).grad D(a)
    (B quadratic in a, T(a) linear in a); f(a) = min_{c,C} |S r|^2 with S a sparse random row sketch that
    preserves exact consistency (k_s >= #linear unknowns + margin), minimized on the unit sphere of the
    gauge complement (F * R_{d-6} and the Koszul syzygies of w are removed) by Levenberg-Marquardt with the
    Kaufman (variable-projection) Jacobian.
  * relaxations          A1 'c = 0' (N = A w only, C free)      -- sufficient sub-search, very cheap
                         A2 'mod L' (B(a) in L_{2d}, L = (F) + J^2 + (F_k H_ij - F_j H_ik)) -- the necessary
                         quadric system of the discussion's DEGREE_EIGHT.md, here solved numerically
  * exact confirmation   reduced modular solve (flint nmod_mat) + CRT/rational lift + core.verify_witness.
Engine files are not modified.  Engine: /home/kejace/devel/sextics-work/hodge.
"""
from __future__ import annotations
import os, sys, time, json, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from flint import nmod_mat

LANE = '/tmp/hodge/results/polarized-k3-remainder'
INSTANCE = 'polarized-k3-remainder'
BK = 64
POW = np.array([BK ** i for i in range(NV)], dtype=np.int64)
PRIMES = list(core.PRIMES) if hasattr(core, 'PRIMES') else [2147483647, 2147483629, 2147483587, 2147483579, 2147483563, 2147483549]
PARZ = [0, 0, 0, 1, 1, 1]            # Y-parity of the variable z_i
PAR_F = [0, 0, 0, 1, 1, 1]           # parity of F_i = dF/dz_i (F Y-even)
PAR_N = [1, 1, 1, 0, 0, 0]           # parity of N_i for Y-even D, Y-odd A
PAR_W = [0, 0, 0, 1, 1, 1]           # parity of w_i (w_i F_i Y-even)
PAIRS = [(i, j) for i in range(NV) for j in range(i + 1, NV)]
PAR_C = {(i, j): (PAR_F[i] + PAR_F[j] + 1) % 2 for (i, j) in PAIRS}   # parity of c_ij


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%S')


def append_json(path, rec):
    with open(path, 'a') as f:
        f.write(json.dumps(rec) + '\n')


def log_batch(rec, path=os.path.join(LANE, 'wave2-search.log.jsonl')):
    rec = dict(rec); rec.setdefault('ts', now()); rec.setdefault('lane', INSTANCE); rec.setdefault('engine', 2)
    append_json(path, rec)
    return rec


def load():
    F, A, q, rec = core.load_instance(INSTANCE)
    return F, A, q


def key(e):
    return int(sum(int(a) * int(b) for a, b in zip(e, POW)))


def keys_of(mons):
    return (np.array(mons, dtype=np.int64) @ POW) if len(mons) else np.zeros(0, np.int64)


def unkey(k):
    e = []
    for _ in range(NV):
        e.append(int(k % BK)); k //= BK
    return tuple(e)


class Space:
    """Monomials of one degree and Y-parity with fast key -> index lookup."""

    def __init__(self, deg, par):
        self.deg, self.par = deg, par
        self.mons = monomials(deg, par) if deg >= 0 else []
        self.n = len(self.mons)
        self.keys = keys_of(self.mons)
        self.order = np.argsort(self.keys)
        self.sorted = self.keys[self.order]
        self.idx = {m: i for i, m in enumerate(self.mons)}

    def index(self, ks):
        ks = np.asarray(ks, dtype=np.int64)
        if ks.size == 0:
            return np.zeros(0, dtype=np.int64)
        pos = np.searchsorted(self.sorted, ks)
        pos = np.minimum(pos, max(self.n - 1, 0))
        if self.n == 0 or not np.all(self.sorted[pos] == ks):
            raise RuntimeError(f'monomial outside Space(deg={self.deg}, par={self.par})')
        return self.order[pos]

    def vec(self, P):
        v = np.zeros(self.n)
        if P:
            ks = np.array([key(e) for e in P], dtype=np.int64)
            v[self.index(ks)] = [float(c) for c in P.values()]
        return v

    def poly(self, v, tol=0.0):
        return {self.mons[i]: float(v[i]) for i in range(self.n) if abs(v[i]) > tol}

    def poly_exact(self, fracs):
        return {self.mons[i]: Fraction(fracs[i]) for i in range(self.n) if fracs[i] != 0}


def ptable(U, W, T):
    """All pairs (u, w) with target index in T.  Returns (iu, iw, it) int arrays of length U.n*W.n."""
    iu = np.repeat(np.arange(U.n, dtype=np.int64), W.n)
    iw = np.tile(np.arange(W.n, dtype=np.int64), U.n)
    it = T.index(U.keys[iu] + W.keys[iw])
    return iu, iw, it


def pmul(tab, x, y, nT):
    iu, iw, it = tab
    return np.bincount(it, weights=x[iu] * y[iw], minlength=nT)


def padj_u(tab, rho, y, nU):
    """Adjoint of x -> x*y:  returns vector in U with (result)_u = sum_w y_w rho_{u+w}."""
    iu, iw, it = tab
    return np.bincount(iu, weights=rho[it] * y[iw], minlength=nU)


def mult_matrix(tab, x, nW, nT):
    """Sparse matrix of multiplication by x (in U): W -> T."""
    iu, iw, it = tab
    return sp.csr_matrix((x[iu], (it, iw)), shape=(nT, nW))


def deriv_matrix(i, U, T):
    """Sparse matrix of d/dz_i : U -> T."""
    ri, ci, va = [], [], []
    for j, m in enumerate(U.mons):
        if m[i]:
            m2 = list(m); m2[i] -= 1
            ri.append(T.idx[tuple(m2)]); ci.append(j); va.append(float(m[i]))
    return sp.csr_matrix((va, (ri, ci)), shape=(T.n, U.n))


def poly_columns(polys, T):
    """Sparse matrix whose columns are the given polys expressed in T."""
    ri, ci, va = [], [], []
    for j, P in enumerate(polys):
        for e, c in P.items():
            ri.append(T.idx[e]); ci.append(j); va.append(float(c))
    return sp.csr_matrix((va, (ri, ci)), shape=(T.n, len(polys)))


def count_sketch(nrows, k, seed=0, per_row=2):
    """Sparse random sketch S (k x nrows): each original row is spread into `per_row` random buckets
    with Gaussian weights.  For generic S of this shape, S x = 0 for x in a fixed subspace of dimension
    < k happens with probability zero, so exact consistency of a linear system is preserved."""
    rng = np.random.default_rng(seed)
    rows = np.repeat(np.arange(nrows), per_row)
    cols = rng.integers(0, k, size=nrows * per_row)
    vals = rng.standard_normal(nrows * per_row)
    S = sp.csr_matrix((vals, (cols, rows)), shape=(k, nrows))
    return S


# ---------------------------------------------------------------------------------------------------------
# E1 kernel (non-Jacobian part) with lifts, mod p and in floating point
# ---------------------------------------------------------------------------------------------------------
def modp_rows_to_rational(rows_by_prime, primes, verify_prime_rows=None):
    """CRT + rational reconstruction of a list of integer rows given mod each prime (same rref pivots)."""
    acc = [list(r) for r in rows_by_prime[0]]
    mod = primes[0]
    for rows, p in zip(rows_by_prime[1:], primes[1:]):
        acc = [[core.crt_pair(a, mod, b, p)[0] for a, b in zip(ra, rb)] for ra, rb in zip(acc, rows)]
        mod *= p
    rat = [[core.ratrecon(v, mod) for v in row] for row in acc]
    ok = all(v is not None for row in rat for v in row)
    return rat, ok


# ---------------------------------------------------------------------------------------------------------
# the reduced family
# ---------------------------------------------------------------------------------------------------------
class Family:
    """Reduced (normalized, Koszul-parametrized) certificate system for Y-even denominators of degree d.

    kind: 'jac'  -> D = sum w_i F_i                           (all E1-admissible D for d <= 9)
          'full' -> jac + non-Jacobian E1 kernel (d >= 10), lifts solved in floating point
          'jac-sub' with wsub: restrict w to a subspace given by a list of (i, monomial) allowed pairs
    """

    def __init__(self, F, A, d, kind='jac', sketch_margin=300, seed=0, verbose=False, cond=1e-10,
                 kernel_data=None, wpred=None, with_relaxations=False, relax_only=False):
        t0 = time.time()
        self.F, self.A, self.d, self.kind = F, A, d, kind
        self.cond = cond
        self.q = 1
        self.verbose = verbose
        Fi = [pderiv(F, i) for i in range(NV)]
        self.Fi = Fi
        # spaces
        self.SD = Space(d, 0)
        self.SdD = [Space(d - 1, PARZ[i]) for i in range(NV)]
        self.SN = [Space(d + 1, PAR_N[i]) for i in range(NV)]
        self.SdivN = Space(d, 1)
        self.SR = Space(2 * d, 1)
        self.Sw = [Space(d - 5, PAR_W[i]) for i in range(NV)]
        self.Sc = {ij: Space(d - 4, PAR_C[ij]) for ij in PAIRS}
        self.SC = Space(2 * d - 6, 1)
        SD, SR = self.SD, self.SR
        # --- w-part: columns (i, m) ---
        wcols = []
        for i in range(NV):
            for m in self.Sw[i].mons:
                if wpred is None or wpred(i, m):
                    wcols.append((i, m))
        self.wcols = wcols
        Kw = len(wcols)
        Dcols = [shift(Fi[i], m, 1) for (i, m) in wcols]
        Dmat = poly_columns(Dcols, SD).toarray() if Kw else np.zeros((SD.n, 0))
        # N_0 = A w : for column (i,m) the N_i block gets A*m
        N0mats = []
        Am = next(iter(A)); Ac = float(A[Am])
        for i in range(NV):
            ri, ci, va = [], [], []
            for j, (ii, m) in enumerate(wcols):
                if ii == i:
                    e = tuple(a + b for a, b in zip(Am, m))
                    ri.append(self.SN[i].idx[e]); ci.append(j); va.append(Ac)
            N0mats.append(sp.csr_matrix((va, (ri, ci)), shape=(self.SN[i].n, Kw)).toarray())
        # --- kernel part ---
        self.Kker = 0
        self.kernel_info = None
        if kind == 'full':
            kd = kernel_data if kernel_data is not None else self.float_kernel_with_lifts(verbose=verbose)
            self.kernel_info = {k: v for k, v in kd.items() if k in ('kernel_dim', 'seconds', 'lift_resid')}
            if kd['kernel_dim']:
                Dk = np.array(kd['Dfloat']).T            # SD.n x Kk
                Dmat = np.hstack([Dmat, Dk])
                for i in range(NV):
                    N0mats[i] = np.hstack([N0mats[i], np.array(kd['Nfloat'][i]).T])
                self.Kker = Dk.shape[1]
        self.K = Dmat.shape[1]
        self.Dmat = Dmat
        self.N0mats = N0mats
        # derivative operators
        self.dD = [deriv_matrix(i, SD, self.SdD[i]) for i in range(NV)]          # SdD_i x SD
        self.dN = [deriv_matrix(i, self.SN[i], self.SdivN) for i in range(NV)]   # SdivN x SN_i
        # product tables
        self.PT_Ddiv = ptable(SD, self.SdivN, SR)
        self.PT_NdD = [ptable(self.SN[i], self.SdD[i], SR) for i in range(NV)]
        # Koszul correction columns c = (ij, m)
        ccols = []
        for ij in PAIRS:
            for m in self.Sc[ij].mons:
                ccols.append((ij, m))
        self.ccols = ccols
        self.nc = len(ccols)
        Kmats = [[[], [], []] for _ in range(NV)]
        for j, ((i1, i2), m) in enumerate(ccols):
            # K_{i1} += m F_{i2},  K_{i2} -= m F_{i1}
            for e, c in shift(Fi[i2], m, 1).items():
                Kmats[i1][0].append(self.SN[i1].idx[e]); Kmats[i1][1].append(j); Kmats[i1][2].append(float(c))
            for e, c in shift(Fi[i1], m, 1).items():
                Kmats[i2][0].append(self.SN[i2].idx[e]); Kmats[i2][1].append(j); Kmats[i2][2].append(-float(c))
        self.Kmat = [sp.csr_matrix((Kmats[i][2], (Kmats[i][0], Kmats[i][1])), shape=(self.SN[i].n, self.nc)).tocsc()
                     for i in range(NV)]
        DivK = None
        for i in range(NV):
            t = self.dN[i] @ self.Kmat[i]
            DivK = t if DivK is None else DivK + t
        self.DivK = DivK.tocsc()                                                  # SdivN x nc
        # F*C columns
        self.MC = poly_columns([shift(F, m, 1) for m in self.SC.mons], SR) if self.SC.n else sp.csr_matrix((SR.n, 0))
        self.nC = self.SC.n
        # gauge: w = z g / 6  (D = F g), g in S_{d-6} even; plus the numerical kernel of a -> D(a) (Koszul syzygies)
        gauge_vecs = []
        if d >= 6:
            Sg = Space(d - 6, 0)
            widx = {wc: j for j, wc in enumerate(wcols)}
            for g in Sg.mons:
                v = np.zeros(self.K)
                okg = True
                for i in range(NV):
                    m = list(g); m[i] += 1; m = tuple(m)
                    if (i, m) in widx:
                        v[widx[(i, m)]] = 1.0 / 6
                    else:
                        okg = False
                if okg:
                    gauge_vecs.append(v)
        self.ngauge = len(gauge_vecs)
        # kernel of Dmat (syzygies among the family) -- exclude
        if self.K:
            U, s, Vt = np.linalg.svd(Dmat, full_matrices=False)
            rD = int(np.sum(s > 1e-9 * s[0]))
            syz = Vt[rD:, :]                      # (K - rD) x K
        else:
            syz = np.zeros((0, 0)); rD = 0
        self.rankD = rD
        bad = np.vstack([np.array(gauge_vecs), syz]) if gauge_vecs else syz
        if bad.size:
            Ub, sb, Vbt = np.linalg.svd(bad, full_matrices=True)
            nb = int(np.sum(sb > 1e-9 * max(sb[0], 1e-300)))
            self.P = Vbt[nb:, :].T                # K x dimP orthonormal complement
            self.nbad = nb
        else:
            self.P = np.eye(self.K); self.nbad = 0
        self.dimP = self.P.shape[1]
        # sketch
        self.nlin = self.nc + self.nC
        self.ks = min(SR.n, self.nlin + sketch_margin)
        self.S = count_sketch(SR.n, self.ks, seed=seed)
        self.ST = self.S.T.tocsr()
        SMC = (self.S @ self.MC).toarray() if self.nC else np.zeros((self.ks, 0))
        if self.nC:
            Qc, Rc = np.linalg.qr(SMC, mode='reduced')
            self.QC = Qc
        else:
            self.QC = np.zeros((self.ks, 0))
        self.SMC = SMC
        self._nfev = 0
        self.relax = None
        if with_relaxations:
            self.build_relaxation_L(verbose=verbose)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[family d={d} kind={kind}] K={self.K} (w {Kw}, ker {self.Kker}) rankD={rD} gauge={self.ngauge} '
                  f'bad={self.nbad} dimP={self.dimP} rows={SR.n} nc={self.nc} nC={self.nC} ks={self.ks} '
                  f'built {self.build_time:.1f}s', file=sys.stderr)

    # ------------------------------------------------------------------------------------------------
    def float_kernel_with_lifts(self, p=None, verbose=False):
        """Non-Jacobian part of the E1 kernel {D : A D in J(F)} in floating point, with lifts.
        The exact dimension (mod p) comes from jacobian.mult_kernel; the vectors come from the float null
        space of [A*m | m F_i] (SVD), projected to a complement of J(F)_d."""
        d = self.d
        p = p or PRIMES[0]
        t0 = time.time()
        mk = jacobian.mult_kernel(self.F, self.A, d, p, parity=0)
        kdim = mk['dim']
        info = {'kernel_dim': kdim, 'modp': {k: mk[k] for k in ('dim', 'kernel_dim', 'jacobian_dim', 'ring_dim')}}
        if kdim == 0:
            info.update({'Dfloat': [], 'Nfloat': [[] for _ in range(NV)], 'seconds': round(time.time() - t0, 1)})
            return info
        SD, SN, ST = self.SD, self.SN, Space(d + 6, 1)
        Acols = poly_columns([shift(self.A, m, 1) for m in SD.mons], ST)
        Jcols = poly_columns([shift(self.Fi[i], m, 1) for i in range(NV) for m in SN[i].mons], ST)
        M = sp.hstack([Acols, Jcols]).toarray()
        nD = SD.n
        U, s, Vt = np.linalg.svd(M, full_matrices=True)
        r = int(np.sum(s > 1e-9 * s[0]))
        NS = Vt[r:, :]                                   # nullity x (nD + nN), orthonormal rows
        del U, Vt
        PD = NS[:, :nD]
        # complement of J(F)_d: orthonormal basis Vw of range(Dmat_w) in SD
        Dw = np.array(self.Dmat[:, :len(self.wcols)]) if len(self.wcols) else np.zeros((nD, 0))
        if Dw.shape[1]:
            Uw, sw, Vwt = np.linalg.svd(Dw, full_matrices=False)
            rw = int(np.sum(sw > 1e-9 * sw[0]))
            Vw = Uw[:, :rw]                              # nD x rw
            PDperp = PD - (PD @ Vw) @ Vw.T
        else:
            PDperp = PD
        Up, sp_, Vpt = np.linalg.svd(PDperp, full_matrices=False)
        rp = int(np.sum(sp_ > 1e-8 * sp_[0]))
        if rp != kdim and verbose:
            print(f'[float kernel d={d}] WARNING float non-J rank {rp} vs mod-p {kdim}', file=sys.stderr)
        take = min(rp, kdim) if rp != kdim else kdim
        W = Up[:, :take].T @ NS                          # take x (nD + nN)
        Dfloat = W[:, :nD]
        Nparts = W[:, nD:]
        Nfloat = []
        off = 0
        for i in range(NV):
            n_i = SN[i].n
            Nfloat.append((-Nparts[:, off:off + n_i]).tolist())
            off += n_i
        res = np.linalg.norm(M @ W.T, axis=0) / np.maximum(np.linalg.norm(W, axis=1), 1e-300)
        info.update({'kernel_dim': take, 'Dfloat': Dfloat.tolist(), 'Nfloat': Nfloat, 'float_rank': rp,
                     'lift_resid': float(res.max()), 'seconds': round(time.time() - t0, 1),
                     'E1_rows': ST.n, 'E1_cols': M.shape[1], 'E1_rank': r})
        if verbose:
            print(f'[float kernel d={d}] non-J dim {take} (mod p {kdim}, float {rp}) lift resid {res.max():.1e} '
                  f'({info["seconds"]}s)', file=sys.stderr)
        return info

    # ------------------------------------------------------------------------------------------------
    # evaluation
    # ------------------------------------------------------------------------------------------------
    def coords(self, u):
        a_full = self.P @ np.asarray(u, float)
        na = np.linalg.norm(a_full)
        return a_full / na, na

    def D_of(self, a):
        return self.Dmat @ a

    def N0_of(self, a):
        return [self.N0mats[i] @ a for i in range(NV)]

    def B_of(self, a, Dv=None, N0=None):
        """Quadratic part B(a) = D div N_0 - sum_i N_0,i d_i D  in SR."""
        Dv = self.D_of(a) if Dv is None else Dv
        N0 = self.N0_of(a) if N0 is None else N0
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.dN[i] @ N0[i]
        B = pmul(self.PT_Ddiv, Dv, divN, self.SR.n)
        for i in range(NV):
            B -= pmul(self.PT_NdD[i], N0[i], self.dD[i] @ Dv, self.SR.n)
        return B

    def ST_of(self, Dv):
        """Sketched correction matrix S T(D) (ks x nc)."""
        MD = mult_matrix(self.PT_Ddiv, Dv, self.SdivN.n, self.SR.n)     # SR x SdivN
        SMD = (self.S @ MD).toarray()                                     # ks x SdivN
        T = np.asarray(self.DivK.T @ SMD.T).T                              # ks x nc  (dense @ sparse)
        for i in range(NV):
            dDi = self.dD[i] @ Dv
            iu, iw, it = self.PT_NdD[i]
            Mi = sp.csr_matrix((dDi[iw], (it, iu)), shape=(self.SR.n, self.SN[i].n))   # mult by d_i D: SN_i -> SR
            SMi = (self.S @ Mi).toarray()
            T -= np.asarray(self.Kmat[i].T @ SMi.T).T
        return T

    def solve_linear(self, a, want_x=True):
        """At a (unit): sketched residual after optimal (c, C).  Returns r (ks), c, Q1 (orthonormal basis of
        the projected correction range), B, Dv, N0."""
        Dv = self.D_of(a); N0 = self.N0_of(a)
        B = self.B_of(a, Dv, N0)
        SB = self.S @ B
        T = self.ST_of(Dv)
        # project out F*C
        if self.nC:
            Tp = T - self.QC @ (self.QC.T @ T)
            bp = SB - self.QC @ (self.QC.T @ SB)
        else:
            Tp, bp = T, SB
        Qf, Rf, piv = sla.qr(Tp, mode='economic', pivoting=True, check_finite=False)
        dg = np.abs(np.diag(Rf))
        rank = int(np.sum(dg > self.cond * dg[0])) if dg.size and dg[0] > 0 else 0
        Q1 = Qf[:, :rank]
        y = Q1.T @ bp
        c = np.zeros(self.nc)
        if rank:
            z = sla.solve_triangular(Rf[:rank, :rank], y, check_finite=False)
            c[piv[:rank]] = -z          # r = bp + Tp c  minimal  => c = -Tp^+ bp
        r = bp - Q1 @ y
        self._nfev += 1
        return r, c, Q1, B, Dv, N0, SB

    def value(self, u):
        a, na = self.coords(u)
        r, c, Q1, B, Dv, N0, SB = self.solve_linear(a)
        return float(r @ r), float(SB @ SB), r, c, a

    def relres(self, u):
        f, bn, r, c, a = self.value(u)
        return math.sqrt(f / max(bn, 1e-300))

    def full_residual(self, a):
        """Unsketched relative residual at a (unit) with the sketched optimum (c, C) re-solved exactly in the
        full row space: min over C of |B + T c - F C| using the full rows."""
        Dv = self.D_of(a); N0 = self.N0_of(a)
        B = self.B_of(a, Dv, N0)
        # full T(D) c for the sketched c
        r, c, Q1, Bs, Dv2, N02, SB = self.solve_linear(a)
        MD = mult_matrix(self.PT_Ddiv, Dv, self.SdivN.n, self.SR.n)
        Tc = MD @ (self.DivK @ c)
        for i in range(NV):
            iu, iw, it = self.PT_NdD[i]
            dDi = self.dD[i] @ Dv
            Mi = sp.csr_matrix((dDi[iw], (it, iu)), shape=(self.SR.n, self.SN[i].n))
            Tc -= Mi @ (self.Kmat[i] @ c)
        v = B + Tc
        if self.nC:
            x, *_ = sla.lstsq(self.MC.toarray(), v, cond=1e-12, lapack_driver='gelsd', check_finite=False)
            v = v - self.MC @ x
        return float(np.linalg.norm(v) / max(np.linalg.norm(B), 1e-300))

    # Jacobian ---------------------------------------------------------------------------------------
    def resid_jac(self, u):
        """Sketched residual r(u) and the Kaufman Jacobian dr/du (ks x dimP)."""
        u = np.asarray(u, float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        r, c, Q1, B, Dv, N0, SB = self.solve_linear(a)
        SR, K = self.SR, self.K
        # dB/da_k : B = beta(a,a), beta(a,a') = D(a) divN0(a') - sum_i N0_i(a') d_iD(a)
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.dN[i] @ N0[i]
        # term 1: D(e_k) * divN(a)   -> Mult_{divN} @ Dmat
        Mdiv = mult_matrix((self.PT_Ddiv[1], self.PT_Ddiv[0], self.PT_Ddiv[2]), divN, self.SD.n, SR.n)  # SR x SD
        JB = Mdiv @ self.Dmat                                                   # SR x K
        # term 2: D(a) * divN0(e_k)  -> Mult_D @ (sum_i dN_i N0mats_i)
        MD = mult_matrix(self.PT_Ddiv, Dv, self.SdivN.n, SR.n)
        divN0mat = None
        for i in range(NV):
            t = self.dN[i] @ self.N0mats[i]
            divN0mat = t if divN0mat is None else divN0mat + t
        JB = JB + MD @ divN0mat
        # term 3: - N0_i(e_k) d_iD(a)  and  - N0_i(a) d_iD(e_k)
        for i in range(NV):
            iu, iw, it = self.PT_NdD[i]
            dDi = self.dD[i] @ Dv
            Mi = sp.csr_matrix((dDi[iw], (it, iu)), shape=(SR.n, self.SN[i].n))      # mult by d_iD: SN_i -> SR
            JB = JB - Mi @ self.N0mats[i]
            Mn = sp.csr_matrix((N0[i][iu], (it, iw)), shape=(SR.n, self.SdD[i].n))    # mult by N0_i: SdD_i -> SR
            JB = JB - Mn @ (self.dD[i] @ self.Dmat)
        # T_k c for all k: D(e_k) divK(c) - sum_i K_i(c) d_iD(e_k)
        divKc = self.DivK @ c
        MdivK = mult_matrix((self.PT_Ddiv[1], self.PT_Ddiv[0], self.PT_Ddiv[2]), divKc, self.SD.n, SR.n)
        JT = MdivK @ self.Dmat
        for i in range(NV):
            iu, iw, it = self.PT_NdD[i]
            Kic = self.Kmat[i] @ c
            Mk = sp.csr_matrix((Kic[iu], (it, iw)), shape=(SR.n, self.SdD[i].n))
            JT = JT - Mk @ (self.dD[i] @ self.Dmat)
        J = np.asarray(self.S @ (JB + JT))                                        # ks x K
        if self.nC:
            J = J - self.QC @ (self.QC.T @ J)
        J = J - Q1 @ (Q1.T @ J)
        J = (J - np.outer(J @ a, a)) / na
        return r, J @ self.P, c, a, SB

    def lm(self, seed=None, init=None, maxiter=100, tol=1e-13, lam0=1e-3, verbose=False, time_limit=None):
        rng = np.random.default_rng(seed)
        u = rng.standard_normal(self.dimP) if init is None else np.asarray(init, float)
        u = u / np.linalg.norm(u)
        t0 = time.time()
        r, J, c, a, SB = self.resid_jac(u)
        f = float(r @ r); bn = float(SB @ SB)
        lam = lam0; it = 0; status = 'maxiter'
        hist = []
        while it < maxiter:
            it += 1
            g = J.T @ r; H = J.T @ J
            dH = np.diag(H).copy(); dH[dH < 1e-300] = 1e-300
            accepted = False
            for _ in range(40):
                try:
                    delta = np.linalg.solve(H + lam * np.diag(dH), -g)
                except np.linalg.LinAlgError:
                    lam *= 10; continue
                u_new = u + delta; u_new /= np.linalg.norm(u_new)
                f_new, bn_new, r_new, c_new, a_new = self.value(u_new)
                if f_new < f:
                    accepted = True; break
                lam *= 4.0
                if lam > 1e16:
                    break
            if not accepted:
                status = 'stalled'; break
            impr = (f - f_new) / max(f, 1e-300)
            u, f, bn = u_new, f_new, bn_new
            lam = max(lam / 3.0, 1e-15)
            hist.append(math.sqrt(f / max(bn, 1e-300)))
            if verbose:
                print(f'  it {it} relres {hist[-1]:.3e} lam {lam:.1e}', file=sys.stderr)
            if hist[-1] < tol:
                status = 'converged'; break
            if time_limit and time.time() - t0 > time_limit:
                status = 'timeout'; break
            r, J, c, a, SB = self.resid_jac(u)
            f = float(r @ r); bn = float(SB @ SB)
            if impr < 1e-10 and lam > 1e4:
                status = 'flat'; break
        a, na = self.coords(u)
        rel = math.sqrt(f / max(bn, 1e-300))
        return {'seed': seed, 'relres': rel, 'iters': it, 'status': status, 'time': round(time.time() - t0, 2),
                'u': u, 'a': a, 'hist': hist[-5:]}

    # ------------------------------------------------------------------------------------------------
    # relaxations
    # ------------------------------------------------------------------------------------------------
    def build_relaxation_L(self, verbose=False, col_sketch=None, seed=1):
        """Orthonormal basis W of the complement of L_{2d} (Y-odd part) in SR, where
        L = (F) + J(F)^2 + (F_k H_ij - F_j H_ik).  Stores self.W (SR.n x codim) and self.WB = projected pair
        tensors are NOT stored; use relax_value()."""
        t0 = time.time()
        d, F, Fi, SR = self.d, self.F, self.Fi, self.SR
        cols = [shift(F, m, 1) for m in Space(2 * d - 6, 1).mons]
        # J^2 * S_{2d-10}
        for i in range(NV):
            for j in range(i, NV):
                P = mul(Fi[i], Fi[j])
                par = (1 - PAR_F[i] - PAR_F[j]) % 2
                for m in Space(2 * d - 10, par).mons:
                    cols.append(shift(P, m, 1))
        # R_ijk = F_k H_ij - F_j H_ik,  deg 9
        H = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
        for i in range(NV):
            for j in range(NV):
                for k in range(j + 1, NV):
                    P = sub(mul(Fi[k], H[i][j]), mul(Fi[j], H[i][k]))
                    if not P:
                        continue
                    par = core.parity_of(P)
                    if par is None:
                        raise RuntimeError('mixed parity R_ijk')
                    for m in Space(2 * d - 9, (1 - par) % 2).mons:
                        cols.append(shift(P, m, 1))
        Lm = poly_columns(cols, SR)
        ncol = Lm.shape[1]
        if verbose:
            print(f'[relax L d={d}] rows {SR.n} spanning columns {ncol}', file=sys.stderr)
        # orthonormal basis of the column space: sketch columns if too many
        rng = np.random.default_rng(seed)
        target = col_sketch or min(ncol, SR.n + 500)
        if ncol > target:
            G = rng.standard_normal((ncol, target))
            Ld = np.asarray(Lm @ G)
        else:
            Ld = Lm.toarray()
        U, s, Vt = np.linalg.svd(Ld, full_matrices=True)
        rk = int(np.sum(s > 1e-9 * s[0]))
        self.W = U[:, rk:]                        # SR.n x codim
        self.rankL = rk
        self.codimL = SR.n - rk
        # also the plain (F)-complement for A1
        MCd = self.MC.toarray()
        Qf, Rf = np.linalg.qr(MCd, mode='reduced')
        self.QF = Qf
        self.relax = True
        if verbose:
            print(f'[relax L d={d}] rank L_{2*d} = {rk}, codim {self.codimL}, (F)-rank {Qf.shape[1]} '
                  f'({time.time()-t0:.1f}s)', file=sys.stderr)
        return self.codimL

    def relax_value_grad(self, u, mode='L'):
        """f(u) = |W^T B(a)|^2 (mode 'L') or |P_F^perp B(a)|^2 (mode 'F' = c = 0 relaxation A1), with gradient."""
        u = np.asarray(u, float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full); a = a_full / na
        Dv = self.D_of(a); N0 = self.N0_of(a)
        B = self.B_of(a, Dv, N0)
        if mode == 'L':
            y = self.W.T @ B
            rho = self.W @ y                      # gradient direction in SR
        else:
            rho = B - self.QF @ (self.QF.T @ B)
            y = rho
        f = float(y @ y)
        # dB/da_k . rho  for all k  (adjoint products)
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.dN[i] @ N0[i]
        g = self.Dmat.T @ padj_u(self.PT_Ddiv, rho, divN, self.SD.n)                       # D(e_k) divN
        adjD = padj_u((self.PT_Ddiv[1], self.PT_Ddiv[0], self.PT_Ddiv[2]), rho, Dv, self.SdivN.n)  # D * divN0(e_k)
        for i in range(NV):
            g += self.N0mats[i].T @ (self.dN[i].T @ adjD)
        for i in range(NV):
            iu, iw, it = self.PT_NdD[i]
            dDi = self.dD[i] @ Dv
            g -= self.N0mats[i].T @ padj_u(self.PT_NdD[i], rho, dDi, self.SN[i].n)          # N0_i(e_k) d_iD
            adjN = padj_u((iw, iu, it), rho, N0[i], self.SdD[i].n)                           # N0_i(a) d_iD(e_k)
            g -= self.Dmat.T @ (self.dD[i].T @ adjN)
        g = 2.0 * g
        g_full = (g - float(g @ a) * a) / na
        return f, self.P.T @ g_full, float(B @ B)

    def relax_minimize(self, seed=None, init=None, mode='L', maxiter=500, gtol=1e-12):
        from scipy.optimize import minimize
        rng = np.random.default_rng(seed)
        u0 = rng.standard_normal(self.dimP) if init is None else np.asarray(init, float)
        u0 /= np.linalg.norm(u0)
        t0 = time.time()
        res = minimize(lambda u: self.relax_value_grad(u, mode)[:2], u0, jac=True, method='L-BFGS-B',
                       options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 30})
        u = res.x / np.linalg.norm(res.x)
        f, g, bn = self.relax_value_grad(u, mode)
        a, na = self.coords(u)
        return {'seed': seed, 'relres': math.sqrt(f / max(bn, 1e-300)), 'iters': int(res.nit), 'time': round(time.time() - t0, 2),
                'u': u, 'a': a, 'status': int(res.status)}

    # ------------------------------------------------------------------------------------------------
    # diagnostics on a candidate D
    # ------------------------------------------------------------------------------------------------
    def D_poly(self, a, tol=1e-12):
        Dv = self.D_of(a)
        s = np.max(np.abs(Dv))
        return self.SD.poly(Dv / s, tol)

    def divisibility(self, a, G, deg_G):
        """Relative residual of D(a) in G * S_{d - deg G} (numerically)."""
        Dv = self.D_of(a); Dv = Dv / np.linalg.norm(Dv)
        parG = core.parity_of(G)
        Sg = Space(self.d - deg_G, (0 - parG) % 2)
        if Sg.n == 0:
            return 1.0
        M = poly_columns([shift(G, m, 1) for m in Sg.mons], self.SD).toarray()
        x, *_ = sla.lstsq(M, Dv, cond=1e-12, check_finite=False)
        return float(np.linalg.norm(M @ x - Dv))

    def reconstruct(self, a, max_den=1000, tol=1e-7, mode='w'):
        """Rational candidates (Fractions) for the family coordinates a (mode 'w') or the D coefficients (mode 'D')."""
        v = np.asarray(a, float) if mode == 'w' else self.D_of(a)
        j = int(np.argmax(np.abs(v)))
        v = v / v[j]
        fr = [Fraction(float(x)).limit_denominator(max_den) for x in v]
        err = max(abs(float(x) - y) for x, y in zip(fr, v))
        return fr, err


# ---------------------------------------------------------------------------------------------------------
# exact reduced solve for D in J(F) (given exact w) -- modular + rational lift
# ---------------------------------------------------------------------------------------------------------
def exact_reduced_system_modp(F, A, d, w, p):
    """Build [T_D | -F C | rhs] mod p for D = sum w_i F_i with w = list of 6 exact polys (Fractions).
    Returns (nmod_mat M (rows x (nc+nC+1)), column descriptors, SR, SC, D poly mod p, N0 polys mod p)."""
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    wp = [core.reduce_mod(wi, p) for wi in w]
    D = {}
    for i in range(NV):
        D = add(D, mul(wp[i], Fi[i], p), p)
    Di = [pderiv(D, i, p) for i in range(NV)]
    N0 = [mul(Ap, wp[i], p) for i in range(NV)]
    SR = Space(2 * d, 1); SC = Space(2 * d - 6, 1)
    Sc = {ij: Space(d - 4, PAR_C[ij]) for ij in PAIRS}
    cols = []; desc = []
    for (i, j) in PAIRS:
        for m in Sc[(i, j)].mons:
            # K_i += m F_j, K_j -= m F_i ;  T_D(c) = D div K - K . grad D
            Ki = shift(Fi[j], m, 1, p); Kj = scale(shift(Fi[i], m, 1, p), p - 1, p)
            divK = add(pderiv(Ki, i, p), pderiv(Kj, j, p), p)
            col = sub(mul(D, divK, p), add(mul(Ki, Di[i], p), mul(Kj, Di[j], p), p), p)
            cols.append(col); desc.append(('c', (i, j), m))
    for m in SC.mons:
        cols.append(scale(shift(Fp, m, 1, p), p - 1, p)); desc.append(('C', None, m))
    # rhs = -B_D(N0)
    divN0 = {}
    for i in range(NV):
        divN0 = add(divN0, pderiv(N0[i], i, p), p)
    B = mul(D, divN0, p)
    for i in range(NV):
        B = sub(B, mul(N0[i], Di[i], p), p)
    rhs = scale(B, p - 1, p)
    R, C = SR.n, len(cols)
    flat = [0] * (R * (C + 1))
    W = C + 1
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[SR.idx[e] * W + j] = v
    for e, v in rhs.items():
        flat[SR.idx[e] * W + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    return M, desc, SR, SC, D, N0


def exact_reduced_solve(F, A, d, w, primes=None, verbose=False, max_primes=10, want_exact=True):
    """Decide the reduced system for D = sum w_i F_i (exact w).  Returns (status, data) with status in
    {'inconsistent', 'solved', 'lift-failed', 'modp-consistent'}.  On 'solved', data has exact D, N, Q verified."""
    primes = primes or PRIMES
    t0 = time.time()
    sols = []; pivots0 = None; rank0 = None
    for p in primes[:max_primes]:
        M, desc, SR, SC, Dp, N0p = exact_reduced_system_modp(F, A, d, w, p)
        C = M.ncols() - 1
        Mr, rank = M.rref()
        piv = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0:
                j += 1
            piv.append(j)
        consistent = not (piv and piv[-1] == C)
        if verbose:
            print(f'[exact-reduced d={d} p={p}] rows={SR.n} cols={C} rank={rank} consistent={consistent} '
                  f'({time.time()-t0:.1f}s)', file=sys.stderr)
        if not consistent:
            return 'inconsistent', {'rank': rank, 'nullity': C - rank, 'rows': SR.n, 'cols': C, 'p': p,
                                    'seconds': round(time.time() - t0, 1)}
        if not want_exact:
            return 'modp-consistent', {'rank': rank, 'nullity': C - rank, 'rows': SR.n, 'cols': C, 'p': p}
        if pivots0 is None:
            pivots0, rank0 = piv, rank
        elif piv != pivots0:
            return 'lift-failed', {'note': 'pivot mismatch between primes'}
        x = [0] * C
        for i, j in enumerate(piv):
            x[j] = int(Mr[i, C])
        sols.append((p, x))
        # try to lift
        acc = list(sols[0][1]); mod = sols[0][0]
        for pp, xx in sols[1:]:
            acc = [core.crt_pair(a_, mod, b_, pp)[0] for a_, b_ in zip(acc, xx)]
            mod *= pp
        rat = [core.ratrecon(v, mod) for v in acc]
        if all(v is not None for v in rat):
            # assemble exact N, Q and verify
            Fi = [pderiv(F, i) for i in range(NV)]
            D = {}
            for i in range(NV):
                D = add(D, mul(w[i], Fi[i]))
            N = [mul(A, w[i]) for i in range(NV)]
            Q = {}
            for (kind, ij, m), v in zip(desc, rat):
                if v == 0:
                    continue
                if kind == 'c':
                    i, j = ij
                    N[i] = add(N[i], scale(shift(Fi[j], m, 1), v))
                    N[j] = sub(N[j], scale(shift(Fi[i], m, 1), v))
                else:
                    Q[m] = Q.get(m, 0) + v
            Q = {e: c for e, c in Q.items() if c}
            ok, info = core.verify_witness(F, A, 1, D, N, Q)
            if ok:
                return 'solved', {'D': D, 'N': N, 'Q': Q, 'rank': rank, 'nullity': C - rank, 'primes': len(sols),
                                  'seconds': round(time.time() - t0, 1)}
            # else continue with more primes
    return 'lift-failed', {'note': f'no rational lift after {len(sols)} primes', 'rank': rank0}


def w_from_D_exact(F, d, D):
    """Solve D = sum w_i F_i exactly (Fractions) if possible: returns list of 6 polys or None."""
    from flint import fmpq_mat, fmpq
    Fi = [pderiv(F, i) for i in range(NV)]
    SD = Space(d, core.parity_of(D))
    wcols = [(i, m) for i in range(NV) for m in Space(d - 5, (core.parity_of(D) - PAR_F[i]) % 2).mons]
    cols = [shift(Fi[i], m, 1) for (i, m) in wcols]
    R, C = SD.n, len(cols)
    M = fmpq_mat(R, C + 1)
    for j, col in enumerate(cols):
        for e, v in col.items():
            M[SD.idx[e], j] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
    for e, v in D.items():
        M[SD.idx[e], C] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
    Mr, rank = M.rref()
    piv = []
    for i in range(rank):
        j = 0
        while j <= C and Mr[i, j] == 0:
            j += 1
        piv.append(j)
    if piv and piv[-1] == C:
        return None
    w = [{} for _ in range(NV)]
    for i, j in enumerate(piv):
        v = Mr[i, C]
        if v != 0:
            ii, m = wcols[j]
            w[ii][m] = Fraction(int(v.p), int(v.q))
    return w


def coprime_to_F(F, D):
    """F irreducible (smooth) => gcd(D,F) != 1 iff F | D."""
    if core.degree(D) < 6:
        return True
    divF, _ = jacobian.divisible_by(F, D)
    return not divF


# ---------------------------------------------------------------------------------------------------------
# self-test against the dict arithmetic of the engine
# ---------------------------------------------------------------------------------------------------------
def selftest(d=6, seed=0, verbose=True):
    F, A, q = load()
    fam = Family(F, A, d, kind='jac', verbose=verbose, sketch_margin=50)
    rng = np.random.default_rng(seed)
    a = rng.standard_normal(fam.K); a /= np.linalg.norm(a)
    # dict version
    Fi = fam.Fi
    w = [{} for _ in range(NV)]
    for j, (i, m) in enumerate(fam.wcols):
        w[i][m] = float(a[j])
    D = {}
    for i in range(NV):
        D = add(D, mul(w[i], Fi[i]))
    N0 = [mul(A, w[i]) for i in range(NV)]
    c = rng.standard_normal(fam.nc)
    N = [dict(n) for n in N0]
    for j, ((i1, i2), m) in enumerate(fam.ccols):
        N[i1] = add(N[i1], scale(shift(Fi[i2], m, 1), float(c[j])))
        N[i2] = sub(N[i2], scale(shift(Fi[i1], m, 1), float(c[j])))
    # E1 check: N(F) = A D
    NF = {}
    for i in range(NV):
        NF = add(NF, mul(N[i], Fi[i]))
    e1 = sub(NF, mul(A, D))
    e1n = max((abs(v) for v in e1.values()), default=0.0)
    # B_D(N) dict
    divN = {}
    for i in range(NV):
        divN = add(divN, pderiv(N[i], i))
    Bd = mul(D, divN)
    for i in range(NV):
        Bd = sub(Bd, mul(N[i], pderiv(D, i)))
    Bvec_dict = fam.SR.vec({e: v for e, v in Bd.items() if abs(v) > 0})
    # numpy version: B(a) + T(a) c  (unsketched, assemble via mult matrices)
    Dv = fam.D_of(a)
    Bnum = fam.B_of(a)
    MD = mult_matrix(fam.PT_Ddiv, Dv, fam.SdivN.n, fam.SR.n)
    Tc = MD @ (fam.DivK @ c)
    for i in range(NV):
        iu, iw, it = fam.PT_NdD[i]
        dDi = fam.dD[i] @ Dv
        Mi = sp.csr_matrix((dDi[iw], (it, iu)), shape=(fam.SR.n, fam.SN[i].n))
        Tc -= Mi @ (fam.Kmat[i] @ c)
    err = np.linalg.norm(Bvec_dict - (Bnum + Tc)) / np.linalg.norm(Bvec_dict)
    # sketched T consistency
    STc = fam.ST_of(Dv) @ c
    err2 = np.linalg.norm(STc - fam.S @ Tc) / np.linalg.norm(STc)
    # Jacobian finite-difference check
    u = rng.standard_normal(fam.dimP); u /= np.linalg.norm(u)
    r0, J, c0, a0, SB = fam.resid_jac(u)
    v = rng.standard_normal(fam.dimP); v -= (v @ u) * u; v /= np.linalg.norm(v)
    h = 1e-6
    fp = fam.value(u + h * v)[0]; fm = fam.value(u - h * v)[0]
    fd = (fp - fm) / (2 * h); an = 2 * float(r0 @ (J @ v))
    # relaxation gradient check
    out = {'E1_dict_err': e1n, 'B_rel_err': float(err), 'sketch_err': float(err2), 'grad_fd': fd, 'grad_an': an}
    if verbose:
        print('[selftest]', out, file=sys.stderr)
    return out, fam


if __name__ == '__main__':
    selftest(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
