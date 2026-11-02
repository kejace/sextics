"""w2rcr -- wave-2 search library for the lane reflection-coupled-residual (q = 1, A = X0^3 Y0^3).

Klein sector bookkeeping: sector(m) = (m5 % 2, (m3+m4) % 2)  [Y2-parity, (Y0+Y1)-parity];
Y-parity = sum of the two.  F has sector (0,0); F_i has sector VARSEC[i]; A has sector (0,1).

Part 1: exact (mod p) E1-admissible spaces  (J:A)_{d,chi} = {D in S_d^chi : A D in J(F)},
        Jacobian part J_{d,chi}, canonical representatives of the quotient, CRT rational lift.
Part 2: E1-reduced variable projection.  For a linear family D(c) = sum c_k M_k with E1 solutions
        N_k (sum_i N_k,i F_i = A M_k), E2 reads   B_D(N) + Kos_D(K) - F Q = 0  with
        B_D(N) = D div N - N(D),  Kos_D(K) = sum_{i<j} D(d_i K_ij F_j - d_j K_ij F_i) - K_ij (F_j D_i - F_i D_j),
        K antisymmetric (deg d-4), Q (deg 2d-6).  Unknowns (K,Q) are linear, c nonlinear (quadratic rhs).
        f(c) = min_{K,Q} || ... ||^2 on the unit sphere; gradient by Golub-Pereyra.  'ob' mode drops (K,Q) and
        projects onto the complement of J_{2d} (necessary condition ob(D) = [B_D(N)] = 0 in R_{2d}).
Engine (/home/kejace/devel/sextics-work/hodge) is used read-only.
"""
import sys, os, time, json, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np, scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
from hodge import core, fixed_d, jacobian, symmetry
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from flint import nmod_mat

LANE = 'reflection-coupled-residual'
LDIR = '/tmp/hodge/results/reflection-coupled-residual'
PRIMES = core.PRIMES
VARSEC = [(0, 0), (0, 0), (0, 0), (0, 1), (0, 1), (1, 0)]
CHI_A = (0, 1)
POW = np.array([64 ** i for i in range(6)], dtype=np.int64)


def sector(m):
    return (m[5] % 2, (m[3] + m[4]) % 2)


def sadd(a, b):
    return ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2)


def monos_sec(deg, s):
    return [m for m in monomials(deg) if sector(m) == s]


def sector_of(P):
    secs = {sector(e) for e in P}
    return secs.pop() if len(secs) == 1 else None


def key(e):
    return int(sum(int(a) * int(b) for a, b in zip(e, POW)))


def unkey(k):
    e = []
    for _ in range(6):
        e.append(int(k % 64)); k //= 64
    return tuple(e)


def load():
    F, A, q, rec = core.load_instance(LANE)
    return F, A, q


def log_batch(rec, path=os.path.join(LDIR, 'wave2-search.log.jsonl')):
    rec = dict(rec); rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    with open(path, 'a') as f:
        f.write(json.dumps(rec) + '\n')


# =====================================================================================================
# Part 1: exact E1-admissible spaces mod p
# =====================================================================================================
def _flat_cols(cols, ridx, p):
    """cols: list of polys mod p (dicts).  Returns flat row-major list for nmod_mat(R, C)."""
    R, C = len(ridx), len(cols)
    flat = [0] * (R * C)
    for j, P in enumerate(cols):
        for e, v in P.items():
            flat[ridx[e] * C + j] = v % p
    return flat


def _rows_to_polys(M, nrows, mons, p):
    out = []
    for r in range(nrows):
        P = {}
        for j, m in enumerate(mons):
            v = int(M[r, j])
            if v:
                P[m] = v
        out.append(P)
    return out


def _rref_rows(vecs, ncols, p):
    if not vecs:
        return 0, None
    M = nmod_mat(len(vecs), ncols, [int(v) % p for row in vecs for v in row], p)
    Mr, rk = M.rref()
    return rk, Mr


def e1_space_modp(F, A, d, chi, p, want_basis=True, verbose=False):
    """Returns dict: dims and (optionally) canonical mod-p representatives of (J:A)_{d,chi} / J_{d,chi}
    expressed in standard monomials of S_d^chi w.r.t. J_{d,chi}, plus a basis of J_{d,chi} (as (i,m) pairs)."""
    t0 = time.time()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    tsec = sadd(chi, CHI_A)
    rows = monos_sec(d + 6, tsec); ridx = {m: i for i, m in enumerate(rows)}; R = len(rows)
    Dm = monos_sec(d, chi); nD = len(Dm)
    gens = [(i, m) for i in range(NV) for m in monos_sec(d + 1, sadd(tsec, VARSEC[i]))]
    cols = [mul(Ap, {m: 1}, p) for m in Dm] + [shift(Fi[i], m, 1, p) for (i, m) in gens]
    C = len(cols)
    M = nmod_mat(R, C, _flat_cols(cols, ridx, p), p)
    del cols
    X, nul = M.nullspace()
    del M
    # D-parts of the nullspace vectors
    Xd = [[int(X[j, k]) for j in range(nD)] for k in range(nul)]
    del X
    rkD, Dr = _rref_rows(Xd, nD, p)   # basis of {D : A D in J}
    out = {'d': d, 'chi': list(chi), 'p': p, 'dim_S': nD, 'rows': R, 'ncols': C, 'nullity': nul,
           'dim_colon': rkD}
    # Jacobian part in degree d, sector chi
    if d >= 5:
        gens_d = [(i, m) for i in range(NV) for m in monos_sec(d - 5, sadd(chi, VARSEC[i]))]
        didx = {m: i for i, m in enumerate(Dm)}
        Jvecs = []
        for (i, m) in gens_d:
            P = shift(Fi[i], m, 1, p)
            row = [0] * nD
            for e, v in P.items():
                row[didx[e]] = v
            Jvecs.append(row)
        rkJ, Jr = _rref_rows(Jvecs, nD, p)
    else:
        gens_d, rkJ, Jr = [], 0, None
    out['dim_J'] = rkJ
    out['dim_nonjac'] = rkD - rkJ
    out['seconds'] = round(time.time() - t0, 1)
    if want_basis and rkD > rkJ:
        # normal forms of the colon basis modulo J_d (pivot elimination), then rref -> canonical reps
        piv = []
        if Jr is not None:
            for r in range(rkJ):
                j = 0
                while int(Jr[r, j]) == 0:
                    j += 1
                piv.append(j)
        reps = []
        for r in range(rkD):
            v = [int(Dr[r, j]) for j in range(nD)]
            for r2, j in enumerate(piv):
                if v[j]:
                    c = v[j]
                    for jj in range(nD):
                        v[jj] = (v[jj] - c * int(Jr[r2, jj])) % p
            if any(v):
                reps.append(v)
        rk2, Rr = _rref_rows(reps, nD, p)
        assert rk2 == rkD - rkJ, (rk2, rkD, rkJ)
        out['reps_modp'] = [[int(Rr[r, j]) for j in range(nD)] for r in range(rk2)]
        out['Dmons'] = Dm
        out['Jpivots'] = piv
    if verbose:
        print(f"[e1] d={d} chi={chi} p={p}: dim S={nD} rows={R} cols={C} colon={rkD} J={rkJ} nonjac={rkD-rkJ} ({out['seconds']}s)",
              file=sys.stderr, flush=True)
    return out


def crt_lift_reps(res1, res2, res3=None):
    """CRT + rational reconstruction of canonical representatives computed mod two primes; optional check mod a third."""
    p1, p2 = res1['p'], res2['p']
    assert res1['Jpivots'] == res2['Jpivots'] and len(res1['reps_modp']) == len(res2['reps_modp'])
    Dm = res1['Dmons']
    polys = []
    for v1, v2 in zip(res1['reps_modp'], res2['reps_modp']):
        P = {}
        for j, (a1, a2) in enumerate(zip(v1, v2)):
            if a1 == 0 and a2 == 0:
                continue
            a, _ = core.crt_pair(a1, p1, a2, p2)
            fr = core.ratrecon(a, p1 * p2)
            if fr is None:
                return None
            P[Dm[j]] = fr
        polys.append(P)
    if res3 is not None:
        p3 = res3['p']
        for P, v3 in zip(polys, res3['reps_modp']):
            for j, m in enumerate(Dm):
                fr = P.get(m, Fraction(0))
                val = (fr.numerator * pow(fr.denominator, -1, p3)) % p3
                if val != v3[j] % p3:
                    return None
    return polys


# =====================================================================================================
# Part 2: dense-vector polynomial spaces and the reduced variable-projection family
# =====================================================================================================
class Space:
    """Monomials of one degree and Klein sector, sorted by packed key; dense vectors live here."""
    _cache = {}

    def __new__(cls, deg, sec):
        k = (deg, tuple(sec))
        if k in cls._cache:
            return cls._cache[k]
        self = super().__new__(cls)
        self.deg, self.sec = deg, tuple(sec)
        mons = monos_sec(deg, sec) if deg >= 0 else []
        keys = np.array([key(m) for m in mons], dtype=np.int64)
        order = np.argsort(keys)
        self.keys = keys[order]
        self.mons = [mons[i] for i in order]
        self.n = len(mons)
        self.exps = np.array(self.mons, dtype=np.int64).reshape(self.n, 6) if self.n else np.zeros((0, 6), np.int64)
        cls._cache[k] = self
        return self

    def index(self, keys):
        idx = np.searchsorted(self.keys, keys)
        if idx.size and (idx.max() >= self.n or not np.array_equal(self.keys[idx], keys)):
            raise KeyError('monomial not in space')
        return idx

    def from_poly(self, P):
        v = np.zeros(self.n)
        if P:
            ks = np.array([key(e) for e in P], dtype=np.int64)
            v[self.index(ks)] = [float(c) for c in P.values()]
        return v

    def to_poly(self, v, tol=0.0):
        return {self.mons[i]: float(v[i]) for i in range(self.n) if abs(v[i]) > tol}


class Bilin:
    """Dense product S_a x S_b -> S_out with precomputed index grid (|a| x |b|)."""
    def __init__(self, sa, sb, sout):
        self.sa, self.sb, self.sout = sa, sb, sout
        grid = sa.keys[:, None] + sb.keys[None, :]
        self.idx = sout.index(grid.ravel()).reshape(sa.n, sb.n)

    def mul(self, u, v):
        return np.bincount(self.idx.ravel(), weights=np.outer(u, v).ravel(), minlength=self.sout.n)

    def corr_a(self, r, v):
        """w[a] = sum_b v[b] r[out(a,b)]  (adjoint in the first slot)."""
        return (r[self.idx] * v[None, :]).sum(axis=1)

    def corr_b(self, r, u):
        return (r[self.idx] * u[:, None]).sum(axis=0)


class SparseMul:
    """Multiplication by a fixed sparse polynomial: S_in -> S_out as a list of index maps."""
    def __init__(self, sin, P, sout):
        self.sin, self.sout = sin, sout
        self.terms = []
        for e, c in P.items():
            idx = sout.index(sin.keys + key(e))
            self.terms.append((idx, float(c)))

    def apply(self, v):
        out = np.zeros(self.sout.n)
        for idx, c in self.terms:
            out[idx] += c * v
        return out

    def adjoint(self, r):
        out = np.zeros(self.sin.n)
        for idx, c in self.terms:
            out += c * r[idx]
        return out


class Deriv:
    """d/dz_i : S_in -> S_out (degree-1, sector + VARSEC[i])."""
    def __init__(self, sin, i, sout):
        self.sin, self.sout, self.i = sin, sout, i
        ex = sin.exps[:, i]
        self.mask = ex > 0
        self.fac = ex[self.mask].astype(float)
        self.idx = sout.index(sin.keys[self.mask] - 64 ** i)

    def apply(self, v):
        out = np.zeros(self.sout.n)
        out[self.idx] += self.fac * v[self.mask]
        return out

    def adjoint(self, r):
        out = np.zeros(self.sin.n)
        out[self.mask] = self.fac * r[self.idx]
        return out


def poly_float(P):
    return {e: float(c) for e, c in P.items()}


class ReducedFamily:
    """E1-reduced variable projection for a linear family of denominators in one Klein sector.

    basis: list of polys (Fraction or float coefficients), degree d, sector chi.
    Nbasis: list of 6-tuples of polys (E1 solutions), or None -> computed by least squares.
    mode: 'full' (Kos + Q unknowns) or 'ob' (projection onto complement of J_{2d}).
    """
    def __init__(self, F, A, basis, Nbasis=None, mode='full', exclude_F_multiples=True, cond=1e-9,
                 verbose=False, nthreads=None):
        t0 = time.time()
        self.F, self.A = F, A
        self.mode = mode; self.cond = cond
        self.basis = [dict(B) for B in basis]
        K = len(basis); self.K = K
        d = core.degree(basis[0]); self.d = d
        chi = sector_of(basis[0]); assert chi is not None
        for B in basis:
            assert core.degree(B) == d and sector_of(B) == chi, 'family must be homogeneous in one sector'
        self.chi = chi
        chiB = CHI_A                      # sector of A D^2 and of B_D(N)
        self.SD = Space(d, chi)
        self.SB = Space(2 * d, chiB)      # rows
        self.SdivN = Space(d, sadd(chiB, chi))
        self.SDi = [Space(d - 1, sadd(chi, VARSEC[i])) for i in range(NV)]
        self.SN = [Space(d + 1, sadd(sadd(chiB, chi), VARSEC[i])) for i in range(NV)]
        Fi = [pderiv(F, i) for i in range(NV)]
        self.Fi = Fi
        # operators
        self.mul_D_divN = Bilin(self.SD, self.SdivN, self.SB)
        self.mul_N_Di = [Bilin(self.SN[i], self.SDi[i], self.SB) for i in range(NV)]
        self.dN = [Deriv(self.SN[i], i, self.SdivN) for i in range(NV)]
        self.dD = [Deriv(self.SD, i, self.SDi[i]) for i in range(NV)]
        # dense basis arrays
        self.Mk = np.array([self.SD.from_poly(B) for B in self.basis])          # K x nD
        # E1 solutions
        if Nbasis is None:
            Nbasis = self.solve_e1_all(verbose=verbose)
        self.Nk = [np.array([self.SN[i].from_poly(N[i]) for N in Nbasis]) for i in range(NV)]  # 6 x (K x nN_i)
        # F-multiple exclusion
        self.P = np.eye(K); self.excluded_dim = 0
        if exclude_F_multiples and d >= 6:
            self.P, self.excluded_dim = self._complement_of_F_multiples()
        self.dimP = self.P.shape[1]
        if mode == 'full':
            self._build_full(verbose)
        else:
            self._build_ob(verbose)
        self._nfev = 0
        self.build_time = time.time() - t0
        if verbose:
            extra = f'nK={self.nK} nQ={self.nQ}' if mode == 'full' else f'dim R_2d sector={self.rdim}'
            print(f'[rf] d={d} chi={chi} K={K} rows={self.SB.n} {extra} excludedF={self.excluded_dim} '
                  f'built in {self.build_time:.1f}s', file=sys.stderr, flush=True)

    # --- E1 ------------------------------------------------------------------------------------------
    def solve_e1_all(self, verbose=False):
        """Least-squares E1 solutions N_k for every basis element (exact for Jacobian members)."""
        d, chi = self.d, self.chi
        ST = Space(d + 6, sadd(chi, CHI_A))
        ops = [SparseMul(self.SN[i], self.Fi[i], ST) for i in range(NV)]
        ncols = sum(S.n for S in self.SN)
        M = np.zeros((ST.n, ncols)); off = [0]
        for i in range(NV):
            for idx, c in ops[i].terms:
                M[idx, off[-1] + np.arange(self.SN[i].n)] += c
            off.append(off[-1] + self.SN[i].n)
        mulA = SparseMul(self.SD, self.A, ST)
        rhs = np.array([mulA.apply(self.Mk[k]) for k in range(self.K)]).T   # ST.n x K
        t0 = time.time()
        X, res, rk, sv = sla.lstsq(M, rhs, cond=1e-12, lapack_driver='gelsd', check_finite=False)
        R = M @ X - rhs
        self.e1_relres = float(np.linalg.norm(R) / max(np.linalg.norm(rhs), 1e-300))
        self.e1_rank = int(rk)
        if verbose:
            print(f'[rf] E1 lstsq {M.shape} rank {rk}: relres {self.e1_relres:.2e} ({time.time()-t0:.1f}s)',
                  file=sys.stderr, flush=True)
        Nb = []
        for k in range(self.K):
            Nb.append([self.SN[i].to_poly(X[off[i]:off[i + 1], k]) for i in range(NV)])
        return Nb

    def _complement_of_F_multiples(self):
        dG = self.d - 6
        SG = Space(dG, self.chi)
        if SG.n == 0:
            return np.eye(self.K), 0
        mulF = SparseMul(SG, self.F, self.SD)
        FM = np.array([mulF.apply(np.eye(SG.n)[j]) for j in range(SG.n)])   # nG x nD
        # a with sum a_k M_k in span(FM): solve [Mk^T | -FM^T] z = 0
        Mb = np.hstack([self.Mk.T, -FM.T])
        ns = sla.null_space(Mb)
        S = ns[:self.K, :]
        if S.shape[1] == 0:
            return np.eye(self.K), 0
        U, s, _ = np.linalg.svd(S, full_matrices=True)
        r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
        return U[:, r:], r

    # --- structure -------------------------------------------------------------------------------------
    def _build_full(self, verbose):
        d, chi = self.d, self.chi
        chiB = CHI_A
        # Q columns: F * m, m in S_{2d-6}^{chiB}
        self.SQ = Space(2 * d - 6, chiB)
        self.nQ = self.SQ.n
        if self.nQ:
            mulF = SparseMul(self.SQ, self.F, self.SB)
            MQ = np.zeros((self.SB.n, self.nQ))
            for idx, c in mulF.terms:
                MQ[idx, np.arange(self.nQ)] += c
            self.MQ = -MQ
        else:
            self.MQ = np.zeros((self.SB.n, 0))
        # Kos columns: pairs (i<j), K_ij in S_{d-4}^{chiB + chi + VARSEC[i] + VARSEC[j]}
        self.SU = [Space(d + 5, sadd(chi, VARSEC[j])) for j in range(NV)]       # D F_j
        self.mulF_D = [SparseMul(self.SD, self.Fi[j], self.SU[j]) for j in range(NV)]
        self.SP = {}
        self.mulF_Di = {}
        self.kcols = []   # list of (i, j, m_key, idxU_j, fac_i, idxU_i, fac_j, idxP)
        for i in range(NV):
            for j in range(i + 1, NV):
                sK = sadd(sadd(chiB, chi), sadd(VARSEC[i], VARSEC[j]))
                SK = Space(d - 4, sK)
                sP = sadd(sadd(chi, VARSEC[i]), VARSEC[j])     # F_j D_i
                SPij = Space(d + 4, sP)
                self.SP[(i, j)] = SPij
                self.mulF_Di[(i, j)] = (SparseMul(self.SDi[i], self.Fi[j], SPij), SparseMul(self.SDi[j], self.Fi[i], SPij))
                for m in SK.mons:
                    km = key(m)
                    idxUj = self.SB.index(self.SU[j].keys + (km - 64 ** i)) if m[i] > 0 else None
                    idxUi = self.SB.index(self.SU[i].keys + (km - 64 ** j)) if m[j] > 0 else None
                    idxP = self.SB.index(SPij.keys + km)
                    self.kcols.append((i, j, idxUj, float(m[i]), idxUi, float(m[j]), idxP))
        self.nK = len(self.kcols)
        self.ncols = self.nK + self.nQ
        self.SK = {(i, j): Space(d - 4, sadd(sadd(chiB, chi), sadd(VARSEC[i], VARSEC[j]))) for i in range(NV) for j in range(i + 1, NV)}

    def _build_ob(self, verbose):
        """Orthonormal basis W of the complement of J_{2d} in the row sector (ob(D) = W^T B_D(N))."""
        d = self.d
        t0 = time.time()
        gens = [(i, m) for i in range(NV) for m in monos_sec(2 * d - 5, sadd(CHI_A, VARSEC[i]))]
        G = np.zeros((self.SB.n, len(gens)))
        for j, (i, m) in enumerate(gens):
            P = shift(self.Fi[i], m, 1)
            G[self.SB.index(np.array([key(e) for e in P], dtype=np.int64)), j] = [float(c) for c in P.values()]
        # left null space of G: vectors v with G^T v = 0
        W = sla.null_space(G.T, rcond=1e-11)
        self.W = W
        self.rdim = W.shape[1]
        if verbose:
            print(f'[rf] ob projector: rows {self.SB.n} gens {len(gens)} -> dim R_2d sector = {self.rdim} '
                  f'({time.time()-t0:.1f}s)', file=sys.stderr, flush=True)

    # --- evaluation ------------------------------------------------------------------------------------
    def _DN(self, a):
        D = a @ self.Mk
        N = [a @ self.Nk[i] for i in range(NV)]
        return D, N

    def B(self, D, N):
        """B_D(N) = D div N - sum N_i D_i as a dense vector in SB."""
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.dN[i].apply(N[i])
        out = self.mul_D_divN.mul(D, divN)
        for i in range(NV):
            out -= self.mul_N_Di[i].mul(N[i], self.dD[i].apply(D))
        return out

    def B_adjoint_D(self, r, N):
        """gradient of r . B_D(N) with respect to the dense vector D."""
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.dN[i].apply(N[i])
        g = self.mul_D_divN.corr_a(r, divN)
        for i in range(NV):
            g -= self.dD[i].adjoint(self.mul_N_Di[i].corr_b(r, N[i]))
        return g

    def B_adjoint_N(self, r, D):
        """gradient of r . B_D(N) with respect to N_i (list of dense vectors)."""
        cD = self.mul_D_divN.corr_b(r, D)      # over SdivN
        out = []
        for i in range(NV):
            gi = self.dN[i].adjoint(cD) - self.mul_N_Di[i].corr_a(r, self.dD[i].apply(D))
            out.append(gi)
        return out

    def kos_matrix(self, D):
        U = [self.mulF_D[j].apply(D) for j in range(NV)]
        Di = [self.dD[i].apply(D) for i in range(NV)]
        Pij = {}
        for (i, j), (m1, m2) in self.mulF_Di.items():
            Pij[(i, j)] = m1.apply(Di[i]) - m2.apply(Di[j])
        M = np.zeros((self.SB.n, self.nK))
        for c, (i, j, idxUj, fi, idxUi, fj, idxP) in enumerate(self.kcols):
            col = M[:, c]
            if idxUj is not None:
                col[idxUj] += fi * U[j]
            if idxUi is not None:
                col[idxUi] -= fj * U[i]
            col[idxP] -= Pij[(i, j)]
        return M

    def kos_to_N(self, xK):
        """N^K_i = sum_j K_ij F_j as dense vectors in SN (from the K coefficient vector)."""
        # assemble K_ij dense, then multiply
        Kd = {}
        c = 0
        for (i, j), SK in self.SK.items():
            Kd[(i, j)] = xK[c:c + SK.n]; c += SK.n
        assert c == self.nK
        NK = [np.zeros(self.SN[i].n) for i in range(NV)]
        for (i, j), v in Kd.items():
            SK = self.SK[(i, j)]
            if not hasattr(self, '_kmul'):
                self._kmul = {}
            if (i, j) not in self._kmul:
                self._kmul[(i, j)] = (SparseMul(SK, self.Fi[j], self.SN[i]), SparseMul(SK, self.Fi[i], self.SN[j]))
            m1, m2 = self._kmul[(i, j)]
            NK[i] += m1.apply(v)
            NK[j] -= m2.apply(v)
        return NK

    def value(self, a, want_x=False):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        D, N = self._DN(a)
        b = -self.B(D, N)
        if self.mode == 'ob':
            ob = self.W.T @ b
            f = float(ob @ ob)
            self._nfev += 1
            return (f, None, ob, b) if want_x else (f, None, ob)
        M = self.kos_matrix(D)
        if self.nQ:
            M = np.hstack([M, self.MQ])
        x, res, rk, sv = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        self._nfev += 1
        self._last_rank = int(rk)
        return (f := float(r @ r), x, r, b) if want_x else (float(r @ r), x, r)

    def value_and_grad(self, u):
        u = np.asarray(u, float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full); a = a_full / na
        D, N = self._DN(a)
        b = -self.B(D, N)
        if self.mode == 'ob':
            ob = self.W.T @ b
            f = float(ob @ ob)
            r = self.W @ ob            # rows-space residual direction (r . B derivative)
            # f = |W^T B|^2 (sign of b irrelevant); d f = 2 ob . W^T dB = 2 r . dB
            gD = self.B_adjoint_D(r, N)
            gN = self.B_adjoint_N(r, D)
            g = -2.0 * (self.Mk @ gD + sum(self.Nk[i] @ gN[i] for i in range(NV)))   # b = -B
        else:
            M = self.kos_matrix(D)
            if self.nQ:
                M = np.hstack([M, self.MQ])
            x, res, rk, sv = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
            r = M @ x - b
            f = float(r @ r)
            # r = Kos_D(K) - F Q + B_D(N);  d r / d a_k = B_{M_k}(N + N^K) + B_D(N_k)
            NK = self.kos_to_N(x[:self.nK])
            Ntot = [N[i] + NK[i] for i in range(NV)]
            gD = self.B_adjoint_D(r, Ntot)
            gN = self.B_adjoint_N(r, D)
            g = 2.0 * (self.Mk @ gD + sum(self.Nk[i] @ gN[i] for i in range(NV)))
            self._last_rank = int(rk)
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def bnorm(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        D, N = self._DN(a)
        b = self.B(D, N)
        return float(np.linalg.norm(b))      # rel = |residual| / |B_D(N)| in both modes

    def is_spurious(self, a, tol=1e-8):
        if self.d < 6:
            return False, 1.0
        a = np.asarray(a, float)
        D = a @ self.Mk
        SG = Space(self.d - 6, self.chi)
        if SG.n == 0:
            return False, 1.0
        mulF = SparseMul(SG, self.F, self.SD)
        FM = np.array([mulF.apply(np.eye(SG.n)[j]) for j in range(SG.n)]).T
        x, *_ = sla.lstsq(FM, D, check_finite=False)
        rel = float(np.linalg.norm(FM @ x - D) / max(np.linalg.norm(D), 1e-300))
        return rel < tol, rel

    def check_grad(self, seed=0, eps=1e-6):
        rng = np.random.default_rng(seed)
        u = rng.standard_normal(self.dimP); u /= np.linalg.norm(u)
        f, g = self.value_and_grad(u)
        v = rng.standard_normal(self.dimP); v /= np.linalg.norm(v)
        fp, _ = self.value_and_grad(u + eps * v); fm, _ = self.value_and_grad(u - eps * v)
        num = (fp - fm) / (2 * eps)
        return f, float(g @ v), num

    def search(self, n_restarts=10, seed=0, maxiter=300, gtol=1e-12, verbose=False, init=None, stop_rel=1e-13,
               callback=None):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            f, x, r = self.value(a)
            bn = self.bnorm(a)
            rel = math.sqrt(f) / max(bn, 1e-300)
            spur, srel = self.is_spurious(a)
            cand = {'value': f, 'rel': rel, 'a': a, 'spurious': spur, 'div_rel': srel, 'iters': int(res.nit),
                    'restart': it, 'nfev': int(res.nfev), 'rank': getattr(self, '_last_rank', None),
                    'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if callback:
                callback(cand)
            if verbose:
                print(f'[rf] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} spurious={spur} ({time.time()-t0:.1f}s)',
                      file=sys.stderr, flush=True)
            if rel < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    # --- reconstruction / confirmation -----------------------------------------------------------------
    def D_float(self, a):
        a = np.asarray(a, float)
        return a @ self.Mk

    def reconstruct_D(self, a, max_den=1000, tol=1e-7):
        """Rational candidates for D(a) in the monomial basis (normalized so the largest |coef| is 1)."""
        v = self.D_float(a)
        j = int(np.argmax(np.abs(v)))
        v = v / v[j]
        out = []
        for md in sorted({10, 100, max_den}):
            fr = [Fraction(float(x)).limit_denominator(md) for x in v]
            err = max(abs(float(f) - x) for f, x in zip(fr, v))
            P = {self.SD.mons[i]: fr[i] for i in range(self.SD.n) if fr[i] != 0}
            out.append((md, err, P))
        return out

    def confirm(self, Dexact, verbose=False, symmetric=True):
        """Exact fixed-D solve via the engine.  Returns (status, data)."""
        if not Dexact:
            return 'inconsistent', {'note': 'D = 0'}
        ok, info = fixed_d.e1_solvable(self.F, self.A, 1, Dexact)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails'}
        try:
            if symmetric:
                status, data = symmetry.solve_fixed_D_symmetric(self.F, self.A, 1, Dexact, verbose=verbose)
            else:
                out = fixed_d.solve_fixed_D_fast(self.F, self.A, 1, Dexact, verbose=verbose)
                if not out['consistent']:
                    return 'inconsistent', {'rank': out['rank'], 'nullity': out['nullity']}
                status, data = ('solved' if out.get('exact') else 'lift-failed'), out
        except Exception as ex:
            status, data = fixed_d.solve_fixed_D(self.F, self.A, 1, Dexact, verbose=verbose)
        return status, data


def jacobian_family(F, d, chi, p=None, independent=True):
    """Exact basis {m F_i} of J_{d,chi} (independent subset chosen mod p when syzygies exist)."""
    Fi = [pderiv(F, i) for i in range(NV)]
    gens = [(i, m) for i in range(NV) for m in monos_sec(d - 5, sadd(chi, VARSEC[i]))]
    polys = [shift(Fi[i], m, 1) for (i, m) in gens]
    if not independent or d < 10:
        return polys, gens
    p = p or PRIMES[0]
    Dm = monos_sec(d, chi); didx = {m: i for i, m in enumerate(Dm)}
    keep, rows = [], []
    rk = 0
    for g, P in zip(gens, polys):
        row = [0] * len(Dm)
        for e, v in core.reduce_mod(P, p).items():
            row[didx[e]] = v
        M = nmod_mat(len(rows) + 1, len(Dm), [v for r in rows + [row] for v in r], p)
        r2 = M.rank()
        if r2 > rk:
            rk = r2; rows.append(row); keep.append((g, P))
    return [P for g, P in keep], [g for g, P in keep]


# =====================================================================================================
# Part 3: E1 annihilator W (left null space of J_{d+6} in a sector) and fast E1 tests
# =====================================================================================================
def e1_annihilator(F, d, chi, p, verbose=False):
    """W: list of row vectors (mod p) over S_{d+6}^{chi+CHI_A} spanning the annihilator of J_{d+6} there.
    Returns (rows monomial list, W as nmod_mat (nul x R))."""
    t0 = time.time()
    Fp = core.reduce_mod(F, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    tsec = sadd(chi, CHI_A)
    rows = monos_sec(d + 6, tsec); ridx = {m: i for i, m in enumerate(rows)}; R = len(rows)
    gens = [(i, m) for i in range(NV) for m in monos_sec(d + 1, sadd(tsec, VARSEC[i]))]
    flat = [0] * (len(gens) * R)
    for g, (i, m) in enumerate(gens):
        base = g * R
        for e, v in shift(Fi[i], m, 1, p).items():
            flat[base + ridx[e]] = v
    GT = nmod_mat(len(gens), R, flat, p)
    del flat
    X, nul = GT.nullspace()          # R x nul, GT X = 0
    del GT
    W = nmod_mat(nul, R, [int(X[i, k]) for k in range(nul) for i in range(R)], p)
    del X
    if verbose:
        print(f'[ann] d={d} chi={chi} p={p}: rows={R} gens={len(gens)} dim R_(d+6) sector={nul} ({time.time()-t0:.0f}s)',
              file=sys.stderr, flush=True)
    return rows, ridx, W


def e1_space_modp_v2(F, A, d, chi, p, want_basis=True, verbose=False, ann=None):
    """Same output as e1_space_modp, via the annihilator W (much smaller nullspace)."""
    t0 = time.time()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    rows, ridx, W = ann if ann is not None else e1_annihilator(F, d, chi, p, verbose)
    R = len(rows); nul = W.nrows()
    Dm = monos_sec(d, chi); nD = len(Dm)
    flat = [0] * (R * nD)
    for j, m in enumerate(Dm):
        for e, v in mul(Ap, {m: 1}, p).items():
            flat[ridx[e] * nD + j] = v
    AM = nmod_mat(R, nD, flat, p)
    del flat
    T = W * AM                        # nul x nD
    X, k = T.nullspace()              # nD x k
    Xd = [[int(X[j, c]) for j in range(nD)] for c in range(k)]
    rkD, Dr = _rref_rows(Xd, nD, p)
    out = {'d': d, 'chi': list(chi), 'p': p, 'dim_S': nD, 'rows': R, 'dim_R_target': nul, 'dim_colon': rkD}
    if d >= 5:
        gens_d = [(i, m) for i in range(NV) for m in monos_sec(d - 5, sadd(chi, VARSEC[i]))]
        didx = {m: i for i, m in enumerate(Dm)}
        Jvecs = []
        for (i, m) in gens_d:
            row = [0] * nD
            for e, v in shift(Fi[i], m, 1, p).items():
                row[didx[e]] = v
            Jvecs.append(row)
        rkJ, Jr = _rref_rows(Jvecs, nD, p)
    else:
        rkJ, Jr = 0, None
    out['dim_J'] = rkJ; out['dim_nonjac'] = rkD - rkJ
    out['seconds'] = round(time.time() - t0, 1)
    if want_basis and rkD > rkJ:
        piv = []
        if Jr is not None:
            for r in range(rkJ):
                j = 0
                while int(Jr[r, j]) == 0:
                    j += 1
                piv.append(j)
        reps = []
        for r in range(rkD):
            v = [int(Dr[r, j]) for j in range(nD)]
            for r2, j in enumerate(piv):
                if v[j]:
                    c = v[j]
                    for jj in range(nD):
                        v[jj] = (v[jj] - c * int(Jr[r2, jj])) % p
            if any(v):
                reps.append(v)
        rk2, Rr = _rref_rows(reps, nD, p)
        assert rk2 == rkD - rkJ, (rk2, rkD, rkJ)
        out['reps_modp'] = [[int(Rr[r, j]) for j in range(nD)] for r in range(rk2)]
        out['Dmons'] = Dm
        out['Jpivots'] = piv
    if verbose:
        print(f"[e1v2] d={d} chi={chi} p={p}: dim S={nD} dim R_target={nul} colon={rkD} J={rkJ} nonjac={rkD-rkJ} ({out['seconds']}s)",
              file=sys.stderr, flush=True)
    return out


def e1_test_modp(F, A, D, p, ann_cache={}):
    """Fast exact E1 test (A D in J mod p) for a sector-pure D using a cached annihilator."""
    d = core.degree(D); chi = sector_of(D)
    if chi is None:
        return None
    k = (d, chi, p)
    if k not in ann_cache:
        ann_cache[k] = e1_annihilator(F, d, chi, p)
    rows, ridx, W = ann_cache[k]
    Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
    v = [0] * len(rows)
    for e, c in mul(Ap, Dp, p).items():
        v[ridx[e]] = c
    r = W * nmod_mat(len(rows), 1, v, p)
    return all(int(r[i, 0]) == 0 for i in range(r.nrows()))
