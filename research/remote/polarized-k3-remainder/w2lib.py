"""Wave-2 search library for lane polarized-k3-remainder (q = 1, A = X0X1X2Y0Y1Y2).

Reduced ("Koszul") certificate system.  For an E1-admissible denominator D (A*D in J(F))
with one exact lift N0 (dF(N0) = A*D), every other lift is N0 + K(c),
    K(c)_i = sum_j c_ij F_j,   c antisymmetric, deg c_ij = d - 4,
and E2 reads   B_D(N0) + T_D(c) = F * Q,   B_D(N) = D div N - N . grad D,
    T_D(c) = sum_{i<j} ( D [(d_i c_ij) F_j - (d_j c_ij) F_i] - c_ij (D_i F_j - D_j F_i) )
(book chapter 5.2).  For Jacobian-type D = sum_i a_i F_i (every E1-admissible D for d <= 9
here, since [A] is injective on (R/J)_d for d <= 9) one may take N0_i = A a_i.

Everything is in the (X,Y) frame with pure Y-parity for D (required by the parity bookkeeping).
Vectors live on monomial index spaces MonSpace(deg, parity).  Multiplication by a polynomial is a
sparse matrix built from a precomputed index pattern (vectorised), so float and mod-p systems share code.
"""
import sys, os, json, time, math, itertools
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')   # the remote launcher sets these to --cpus; locally stay single-threaded
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from fractions import Fraction as Fr
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, mul, add, sub, scale, shift, pderiv, power

LANE = 'polarized-k3-remainder'
RES = f'/tmp/hodge/results/{LANE}'
P0 = 2147483647
PRIMES = core.PRIMES

# ----------------------------------------------------------------------------- instance
_INST = None
def load():
    global _INST
    if _INST is None:
        F, A, q, rec = core.load_instance(LANE)
        _INST = (F, A, q)
    return _INST

def Fparts(F):
    Fi = [pderiv(F, i) for i in range(NV)]
    Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
    return Fi, Fij

# ----------------------------------------------------------------------------- UV frame helpers
def var(i):
    e = [0] * 6; e[i] = 1
    return {tuple(e): Fr(1)}

def uv_to_xy(P):
    """P given in slots (U0,U1,U2,V0,V1,V2) -> polynomial in (X,Y) with U = X+Y, V = X-Y."""
    X = [var(i) for i in range(3)]; Y = [var(3 + i) for i in range(3)]
    U = [add(X[i], Y[i]) for i in range(3)]; V = [sub(X[i], Y[i]) for i in range(3)]
    lin = U + V
    cache = {}
    def pw(i, k):
        if (i, k) not in cache:
            cache[(i, k)] = power(lin[i], k)
        return cache[(i, k)]
    out = {}
    for e, c in P.items():
        t = {(0,) * 6: Fr(c)}
        for i in range(6):
            if e[i]:
                t = mul(t, pw(i, e[i]))
        out = add(out, t)
    return out

def xy_to_uv(P):
    """inverse: X = (U+V)/2, Y = (U-V)/2 (slots 0..2 = U, 3..5 = V)."""
    U = [var(i) for i in range(3)]; V = [var(3 + i) for i in range(3)]
    X = [scale(add(U[i], V[i]), Fr(1, 2)) for i in range(3)]
    Y = [scale(sub(U[i], V[i]), Fr(1, 2)) for i in range(3)]
    lin = X + Y
    cache = {}
    def pw(i, k):
        if (i, k) not in cache:
            cache[(i, k)] = power(lin[i], k)
        return cache[(i, k)]
    out = {}
    for e, c in P.items():
        t = {(0,) * 6: Fr(c)}
        for i in range(6):
            if e[i]:
                t = mul(t, pw(i, e[i]))
        out = add(out, t)
    return out

def swapUV(P):
    return {(e[3], e[4], e[5], e[0], e[1], e[2]): c for e, c in P.items()}

def f_ternary():
    """f with f(U) = F(U/2... ) i.e. F = (f(U)+f(V))/2; f in slots 0..2."""
    F, A, q = load()
    Fuv = xy_to_uv(F)
    f = {}
    for e, c in Fuv.items():
        if e[3] == e[4] == e[5] == 0:
            f[e] = 2 * c
    return f

# ----------------------------------------------------------------------------- monomial spaces
_SPACES = {}
class MonSpace:
    def __init__(self, deg, par):
        self.deg, self.par = deg, par
        self.mons = monomials(deg, par) if deg >= 0 else []
        self.n = len(self.mons)
        self.idx = {m: i for i, m in enumerate(self.mons)}
        self.E = np.array(self.mons, dtype=np.int64).reshape(-1, 6)
        self.keys = self.E @ (32 ** np.arange(6, dtype=np.int64))
        order = np.argsort(self.keys)
        self.skeys = self.keys[order]; self.sorder = order

    def lookup(self, Emat):
        k = Emat @ (32 ** np.arange(6, dtype=np.int64))
        pos = np.searchsorted(self.skeys, k)
        pos = np.clip(pos, 0, len(self.skeys) - 1)
        ok = self.skeys[pos] == k
        return self.sorder[pos], ok

    def vec(self, P, dtype=float, p=None):
        v = np.zeros(self.n, dtype=dtype)
        for e, c in P.items():
            if p is not None:
                c = Fr(c); c = (c.numerator * pow(c.denominator, -1, p)) % p
            else:
                c = float(c)
            v[self.idx[e]] = c
        return v

    def poly(self, v, tol=0.0):
        return {m: v[i] for i, m in enumerate(self.mons) if abs(v[i]) > tol}

def space(deg, par):
    key = (deg, par)
    if key not in _SPACES:
        _SPACES[key] = MonSpace(deg, par)
    return _SPACES[key]

def parsum(*ps):
    """parity sum with None = mixed."""
    if any(p is None for p in ps):
        return None
    return sum(ps) % 2

_PATTERNS = {}
def mult_pattern(dp, pp, d2, p2):
    """Index pattern for multiplication R_{d2,p2} -> R_{dp+d2, pp+p2} by a polynomial in R_{dp,pp}.
    Returns (rows, cols, eidx, shape)."""
    key = (dp, pp, d2, p2)
    if key in _PATTERNS:
        return _PATTERNS[key]
    SP, S2, SO = space(dp, pp), space(d2, p2), space(dp + d2, parsum(pp, p2))
    ne, nm = SP.n, S2.n
    Esum = (SP.E[:, None, :] + S2.E[None, :, :]).reshape(-1, 6)
    rows, ok = SO.lookup(Esum)
    assert ok.all()
    eidx = np.repeat(np.arange(ne), nm)
    cols = np.tile(np.arange(nm), ne)
    pat = (rows.astype(np.int64), cols.astype(np.int64), eidx.astype(np.int64), (SO.n, nm))
    _PATTERNS[key] = pat
    return pat

def mult_matrix(vecP, pat, dtype=None):
    rows, cols, eidx, shape = pat
    data = vecP[eidx]
    nz = data != 0
    M = sp.csr_matrix((data[nz], (rows[nz], cols[nz])), shape=shape, dtype=dtype or vecP.dtype)
    return M

def spmod(A, B, p):
    """(A @ B) mod p for sparse/dense int64 matrices with entries in [0,p), p < 2^31, via 16-bit split."""
    if sp.issparse(B):
        Blo = B.copy(); Blo.data = Blo.data & 0xFFFF
        Bhi = B.copy(); Bhi.data = Bhi.data >> 16
    else:
        Blo = B & 0xFFFF; Bhi = B >> 16
    lo = A @ Blo; hi = A @ Bhi
    if sp.issparse(lo):
        lo = lo.toarray(); hi = hi.toarray()
    lo = np.asarray(lo) % p; hi = np.asarray(hi) % p
    return (lo + (hi * (1 << 16)) % p) % p

# ----------------------------------------------------------------------------- Jacobian-type family
def jac_basis(d, pD):
    """Basis (i0, m, poly m*F_i0) of J(F)_d with Y-parity pD (independent for d <= 9)."""
    F, A, q = load()
    Fi, _ = Fparts(F)
    out = []
    for i0 in range(NV):
        pm = (pD - (1 if i0 >= 3 else 0)) % 2
        for m in monomials(d - 5, pm):
            out.append((i0, m, shift(Fi[i0], m, 1)))
    return out

def f_multiple_vectors(basis, d, pD):
    """Vectors (in a-space) of the F-multiples  F*h = sum_i (z_i h/6) F_i."""
    pos = {(i0, m): k for k, (i0, m, _) in enumerate(basis)}
    vecs = []
    for h in monomials(d - 6, pD):
        v = np.zeros(len(basis))
        for i in range(NV):
            m = list(h); m[i] += 1; m = tuple(m)
            v[pos[(i, m)]] = 1.0 / 6
        vecs.append(v)
    return np.array(vecs).T if vecs else np.zeros((len(basis), 0))

def complement(S, K):
    """Orthonormal basis of the orthogonal complement of the column space of S in R^K."""
    if S.shape[1] == 0:
        return np.eye(K)
    U, s, _ = np.linalg.svd(S, full_matrices=True)
    r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
    return U[:, r:]

# ----------------------------------------------------------------------------- the reduced E2 system
def c_pairs(d, pD):
    """Unknown c_ij monomials: list of (i, j, n) for i<j, deg n = d-4, parity 1 - pD - [i>=3] - [j>=3]."""
    out = []
    for i in range(NV):
        for j in range(i + 1, NV):
            par = (1 - pD - (1 if i >= 3 else 0) - (1 if j >= 3 else 0)) % 2
            for n in monomials(d - 4, par):
                out.append((i, j, n))
    return out

class ReducedE2:
    """Reduced E2 system for a linear family of E1-admissible denominators of degree d and parity pD.

    family: list of (Dpoly, N0list) where N0list are the six exact-lift numerators (dF(N0) = A D).
    For Jacobian-type families use ReducedE2.jacobian(d, pD).
    Provides float assembly for variable projection (with the Q-block eliminated by the X0-remainder
    map, F being monic in X0) and an exact mod-p consistency test for one D."""

    def __init__(self, d, pD, family, verbose=True, build_float=True, use_gauge=True):
        t0 = time.time()
        F, A, q = load()
        self.F, self.A, self.d, self.pD = F, A, d, pD
        Fi, Fij = Fparts(F)
        self.Fi = Fi
        self.pars = [(1 - pD - (1 if i >= 3 else 0)) % 2 for i in range(NV)]   # parity of N_i
        self.S_D = space(d, pD)
        self.S_W = space(d, (1 - pD) % 2)                  # W vectors, div N0
        self.S_Di = [space(d - 1, (pD + (1 if i >= 3 else 0)) % 2) for i in range(NV)]
        self.S_N = [space(d + 1, self.pars[i]) for i in range(NV)]
        self.S_B = space(2 * d, 1)
        self.S_Q = space(2 * d - 6, 1)
        self.K = len(family)
        # family matrices
        self.Mbasis = sp.csc_matrix(np.array([self.S_D.vec(D) for D, _ in family]).T)
        self.DMbasis = [sp.csc_matrix(np.array([self.S_Di[i].vec(pderiv(D, i)) for D, _ in family]).T)
                        for i in range(NV)]
        self.N0mat = [sp.csc_matrix(np.array([self.S_N[i].vec(N0[i]) for _, N0 in family]).T) for i in range(NV)]
        divs = []
        for _, N0 in family:
            dv = {}
            for i in range(NV):
                dv = add(dv, pderiv(N0[i], i))
            divs.append(self.S_W.vec(dv))
        self.DivN0mat = sp.csc_matrix(np.array(divs).T)
        self.family = family
        # c unknowns
        self.cpairs = c_pairs(d, pD)
        nc = len(self.cpairs); self.nc = nc
        Wr, Wc, Wv = [], [], []
        Vr = [[] for _ in range(NV)]; Vc = [[] for _ in range(NV)]; Vv = [[] for _ in range(NV)]
        for col, (i, j, n) in enumerate(self.cpairs):
            npoly = {n: Fr(1)}
            W = sub(mul(pderiv(npoly, i), Fi[j]), mul(pderiv(npoly, j), Fi[i]))
            for e, c in W.items():
                Wr.append(self.S_W.idx[e]); Wc.append(col); Wv.append(float(c))
            for e, c in mul(npoly, Fi[j]).items():      # slot i gets  - n F_j ... see T_D: - c_ij (D_i F_j - D_j F_i)
                Vr[i].append(self.S_N[i].idx[e]); Vc[i].append(col); Vv[i].append(float(c))
            for e, c in mul(npoly, Fi[i]).items():      # slot j gets  + n F_i
                Vr[j].append(self.S_N[j].idx[e]); Vc[j].append(col); Vv[j].append(-float(c))
        self.W = sp.csr_matrix((Wv, (Wr, Wc)), shape=(self.S_W.n, nc))
        # V_i: column (ij,n) holds the degree-(d+1) poly  n F_j  in slot i  and  -n F_i in slot j, so that
        #   T_D(c) = D * (W c) - sum_i D_i * (V_i c)
        self.V = [sp.csr_matrix((Vv[i], (Vr[i], Vc[i])), shape=(self.S_N[i].n, nc)) for i in range(NV)]
        # Q block (sparse): F * m
        qr, qc, qv = [], [], []
        for col, m in enumerate(self.S_Q.mons):
            for e, c in shift(F, m, 1).items():
                qr.append(self.S_B.idx[e]); qc.append(col); qv.append(float(c))
        self.Qmat = sp.csr_matrix((qv, (qr, qc)), shape=(self.S_B.n, self.S_Q.n))
        # patterns
        self.pat_D = mult_pattern(d, pD, d, (1 - pD) % 2)                  # D * (S_W)
        self.pat_Di = [mult_pattern(d - 1, self.S_Di[i].par, d + 1, self.pars[i]) for i in range(NV)]
        self.pat_W = mult_pattern(d, (1 - pD) % 2, d, pD)                   # (S_W elt) * (S_D)
        self.pat_Ni = [mult_pattern(d + 1, self.pars[i], d - 1, self.S_Di[i].par) for i in range(NV)]
        self.build_time = time.time() - t0
        self.Cp = None; self.nc_eff = self.nc
        if build_float:
            self._build_reduction()
            if use_gauge:
                self._gauge()
        if verbose:
            print(f'[ReducedE2] d={d} pD={pD} K={self.K} rows={self.S_B.n} nc={nc} nQ={self.S_Q.n} '
                  f'red_rows={getattr(self, "nred", None)} nc_eff={getattr(self, "nc_eff", None)} '
                  f'built in {time.time()-t0:.1f}s', file=sys.stderr, flush=True)

    @classmethod
    def jacobian(cls, d, pD, verbose=True, build_float=True):
        F, A, q = load()
        basis = jac_basis(d, pD)
        fam = []
        for (i0, m, Dk) in basis:
            N0 = [{} for _ in range(NV)]
            N0[i0] = mul(A, {m: Fr(1)})
            fam.append((Dk, N0))
        obj = cls(d, pD, fam, verbose=verbose, build_float=build_float)
        obj.basis = basis
        obj.SF = f_multiple_vectors(basis, d, pD)
        obj.P = complement(obj.SF, obj.K)          # a = P u
        obj.dimP = obj.P.shape[1]
        return obj

    # ---- X0-remainder reduction: kills exactly F * S_Q inside S_B
    def _build_reduction(self):
        SB = self.S_B
        e0 = SB.E[:, 0]
        self.top = np.where(e0 <= 5)[0]; self.bot = np.where(e0 >= 6)[0]
        self.nred = len(self.top)
        ntop, nbot = len(self.top), len(self.bot)
        # remainder of each bottom monomial modulo F (as a monic polynomial in X0), expressed on top monomials
        # X0^6 = X0^6 - F + F  =>  X0^6 * m' = (X0^6 - F) m'  mod F.   G := X0^6 - F has X0-degree <= 5.
        F = self.F
        G = sub({(6, 0, 0, 0, 0, 0): Fr(1)}, F)
        Gl = [(np.array(e), float(c)) for e, c in G.items()]
        top_pos = {int(r): k for k, r in enumerate(self.top)}
        bot_pos = {int(r): k for k, r in enumerate(self.bot)}
        Red = np.zeros((ntop, nbot))
        # process bottom monomials in increasing X0-degree
        order = sorted(range(nbot), key=lambda k: e0[self.bot[k]])
        E = SB.E
        for k in order:
            r = self.bot[k]
            m = E[r].copy(); m[0] -= 6
            vec = np.zeros(ntop)
            for ge, gc in Gl:
                e2 = m + ge
                j, ok = SB.lookup(e2[None, :])
                j = int(j[0]); assert ok[0]
                if e2[0] <= 5:
                    vec[top_pos[j]] += gc
                else:
                    vec += gc * Red[:, bot_pos[j]]
            Red[:, k] = vec
        self.Red = Red

    def subsample_rows(self, m, seed=0):
        """Keep a random subset of m reduced rows (necessary-condition relaxation of E2; cheaper evaluations)."""
        rng = np.random.default_rng(seed)
        if m >= len(self.top):
            return
        S = np.sort(rng.choice(len(self.top), size=m, replace=False))
        self.top = self.top[S]; self.Red = np.ascontiguousarray(self.Red[S, :]); self.nred = m
        self.row_subsample = m

    def reduce_rows(self, M):
        """Apply the remainder map to a dense matrix/vector with rows indexed by S_B."""
        if M.ndim == 1:
            return M[self.top] + self.Red @ M[self.bot]
        return M[self.top, :] + self.Red @ M[self.bot, :]

    def _gauge(self):
        """Remove c-directions with K(c) = 0 identically (second Koszul syzygies, d >= 9)."""
        if self.d - 9 < 0:
            self.Cp = None; self.nc_eff = self.nc; return
        Vst = sp.vstack(self.V).toarray()
        ns = sla.null_space(Vst)
        if ns.shape[1] == 0:
            self.Cp = None; self.nc_eff = self.nc; return
        self.Cp = complement(ns, self.nc)
        self.nc_eff = self.Cp.shape[1]

    # ---- float assembly
    def assemble(self, a, want_jac=True):
        """Returns dict with reduced dense Mc (nred x nc_eff), b0 (nred), and pieces for the Jacobian."""
        a = np.asarray(a, float)
        Dv = self.Mbasis @ a
        Div = [self.DMbasis[i] @ a for i in range(NV)]
        N0v = [self.N0mat[i] @ a for i in range(NV)]
        divN0 = self.DivN0mat @ a
        MD = mult_matrix(Dv, self.pat_D)
        MDi = [mult_matrix(Div[i], self.pat_Di[i]) for i in range(NV)]
        Mc = MD @ self.W
        for i in range(NV):
            Mc = Mc - MDi[i] @ self.V[i]
        Mc = np.asarray(Mc.todense())
        b0 = MD @ divN0
        for i in range(NV):
            b0 = b0 - MDi[i] @ N0v[i]
        out = {'Mc_full': Mc, 'b0_full': b0, 'MD': MD, 'MDi': MDi, 'divN0': divN0, 'N0v': N0v, 'Dv': Dv, 'Div': Div}
        Mcr = self.reduce_rows(Mc)
        if self.Cp is not None:
            Mcr = Mcr @ self.Cp
        out['Mc'] = Mcr
        out['b0'] = self.reduce_rows(b0)
        return out

    def residual(self, a, asm=None, with_jac=False):
        """Variable projection: min_c |Mc c + b0|.  Returns dict(f, rel, c, r, [J])."""
        if asm is None:
            asm = self.assemble(a)
        Mc, b0 = asm['Mc'], asm['b0']
        if getattr(self, 'b_offset', None) is not None:
            b0 = b0 + self.b_offset   # synthetic-consistency tests only
        Qm, Rm, piv = sla.qr(Mc, mode='economic', pivoting=True, check_finite=False)
        dg = np.abs(np.diag(Rm))
        rank = int(np.sum(dg > 1e-10 * max(dg[0], 1e-300))) if dg.size else 0
        Q1 = Qm[:, :rank]
        qb = Q1.T @ b0
        r = b0 - Q1 @ qb             # residual of  Mc c + b0 = 0  with optimal c : r = Mc c + b0
        # particular c (minimum over pivot columns)
        c = np.zeros(Mc.shape[1])
        if rank:
            y = sla.solve_triangular(Rm[:rank, :rank], -qb, check_finite=False)
            c[piv[:rank]] = y
        f = float(r @ r)
        bn = float(b0 @ b0)
        out = {'f': f, 'rel': f / max(bn, 1e-300), 'c': c, 'r': r, 'rank': rank, 'bnorm2': bn}
        if with_jac:
            out['J'] = self.kaufman_jacobian(a, asm, c, Q1)
        return out

    def kaufman_jacobian(self, a, asm, c, Q1):
        """d r / d a  (reduced rows, projected off the column space of Mc)."""
        MD, MDi, divN0, N0v = asm['MD'], asm['MDi'], asm['divN0'], asm['N0v']
        cfull = self.Cp @ c if self.Cp is not None else c
        # derivative of b0(a) = D divN0 - sum_i D_i N0_i  (bilinear) w.r.t. a:
        JB = mult_matrix(divN0, self.pat_W) @ self.Mbasis + MD @ self.DivN0mat
        for i in range(NV):
            JB = JB - mult_matrix(N0v[i], self.pat_Ni[i]) @ self.DMbasis[i] - MDi[i] @ self.N0mat[i]
        # derivative of Mc(a) c = D (W c) - sum_i D_i (V_i c):
        Wc = self.W @ cfull
        JT = mult_matrix(Wc, self.pat_W) @ self.Mbasis
        for i in range(NV):
            Vc = self.V[i] @ cfull
            JT = JT - mult_matrix(Vc, self.pat_Ni[i]) @ self.DMbasis[i]
        J = np.asarray((JB + JT).todense())
        J = self.reduce_rows(J)
        J = J - Q1 @ (Q1.T @ J)
        return J

    # ---- exact mod-p consistency for one member a (rational) or an explicit (D, N0)
    def exact_system_modp(self, D, N0, p, verbose=False):
        """Dense nmod_mat [Mc | Qmat | -b0] for the given exact D, N0 (Fraction polys), mod p."""
        from flint import nmod_mat
        t0 = time.time()
        Dv = self.S_D.vec(D, dtype=np.int64, p=p)
        Div = [self.S_Di[i].vec(pderiv(D, i), dtype=np.int64, p=p) for i in range(NV)]
        N0v = [self.S_N[i].vec(N0[i], dtype=np.int64, p=p) for i in range(NV)]
        dv = {}
        for i in range(NV):
            dv = add(dv, pderiv(N0[i], i))
        divN0 = self.S_W.vec(dv, dtype=np.int64, p=p)
        def modmat(M):
            # float sparse matrices with small integer/rational entries -> int64 mod p (entries are exact rationals
            # with denominator dividing 7^k * 6^k; rebuild from Fractions to be safe)
            return M
        # exact int64 versions of W, V, Q (rebuild with Fractions mod p)
        Wm = self._exact_sparse('W', p); Vm = [self._exact_sparse(('V', i), p) for i in range(NV)]
        Qm = self._exact_sparse('Q', p)
        MD = mult_matrix(Dv, self.pat_D, dtype=np.int64)
        MDi = [mult_matrix(Div[i], self.pat_Di[i], dtype=np.int64) for i in range(NV)]
        Mc = spmod(MD, Wm, p)
        for i in range(NV):
            Mc = (Mc - spmod(MDi[i], Vm[i], p)) % p
        b0 = spmod(MD, divN0.reshape(-1, 1), p)
        for i in range(NV):
            b0 = (b0 - spmod(MDi[i], N0v[i].reshape(-1, 1), p)) % p
        Qd = Qm.toarray() % p
        rhs = (-b0) % p
        Aug = np.concatenate([Mc, Qd, rhs], axis=1).astype(np.int64)
        R, C = Aug.shape
        t1 = time.time()
        M = nmod_mat(R, C, Aug.ravel().tolist(), p)
        if verbose:
            print(f'[exact] built {R}x{C} mod {p} in {time.time()-t0:.1f}s', file=sys.stderr, flush=True)
        return M, R, C - 1

    def _exact_sparse(self, which, p):
        key = ('_ex', which, p)
        if hasattr(self, '_excache') and key in self._excache:
            return self._excache[key]
        if not hasattr(self, '_excache'):
            self._excache = {}
        F = self.F; Fi = self.Fi
        r, c_, v = [], [], []
        def red(x):
            x = Fr(x); return (x.numerator * pow(x.denominator, -1, p)) % p
        if which == 'W':
            for col, (i, j, n) in enumerate(self.cpairs):
                npoly = {n: Fr(1)}
                W = sub(mul(pderiv(npoly, i), Fi[j]), mul(pderiv(npoly, j), Fi[i]))
                for e, cc in W.items():
                    r.append(self.S_W.idx[e]); c_.append(col); v.append(red(cc))
            M = sp.csr_matrix((v, (r, c_)), shape=(self.S_W.n, self.nc), dtype=np.int64)
        elif which == 'Q':
            for col, m in enumerate(self.S_Q.mons):
                for e, cc in shift(F, m, 1).items():
                    r.append(self.S_B.idx[e]); c_.append(col); v.append(red(cc))
            M = sp.csr_matrix((v, (r, c_)), shape=(self.S_B.n, self.S_Q.n), dtype=np.int64)
        else:
            _, i = which
            for col, (ii, jj, n) in enumerate(self.cpairs):
                npoly = {n: Fr(1)}
                if ii == i:
                    for e, cc in mul(npoly, Fi[jj]).items():
                        r.append(self.S_N[i].idx[e]); c_.append(col); v.append(red(cc))
                if jj == i:
                    for e, cc in mul(npoly, Fi[ii]).items():
                        r.append(self.S_N[i].idx[e]); c_.append(col); v.append(red(-cc))
            M = sp.csr_matrix((v, (r, c_)), shape=(self.S_N[i].n, self.nc), dtype=np.int64)
        self._excache[key] = M
        return M

    def exact_consistent(self, D, N0, p=P0, verbose=False, compress=None):
        """Exact mod-p consistency of the reduced E2 system for (D, N0).  Returns dict."""
        from flint import nmod_mat
        t0 = time.time()
        M, R, C = self.exact_system_modp(D, N0, p, verbose=verbose)
        if compress and R > C + compress:
            # random row compression (k = C + compress rows)
            k = C + compress
            rng = np.random.default_rng(1)
            Rm = rng.integers(0, p, size=(k, R), dtype=np.int64)
            Aug = np.array(M.entries(), dtype=np.int64).reshape(R, C + 1) if hasattr(M, 'entries') else None
            raise NotImplementedError
        Mr, rank = M.rref()
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0:
                j += 1
            pivots.append(j)
        consistent = not (pivots and pivots[-1] == C)
        x = None
        if consistent:
            x = [0] * C
            for i, j in enumerate(pivots):
                x[j] = int(Mr[i, C])
        return {'consistent': consistent, 'rank': rank, 'nullity': C - rank, 'rows': R, 'cols': C,
                'x': x, 'pivots': pivots, 'p': p, 'seconds': round(time.time() - t0, 1)}

    def unpack_exact(self, x, a_frac=None, D=None, N0=None):
        """Mod-p solution vector -> (c polys, Q poly) mod p."""
        c = {}
        for col, (i, j, n) in enumerate(self.cpairs):
            if x[col]:
                c[(i, j, n)] = x[col]
        Q = {}
        for k, m in enumerate(self.S_Q.mons):
            v = x[self.nc + k]
            if v:
                Q[m] = v
        return c, Q

    def numerators_from_c(self, N0, c, p=None):
        """N = N0 + K(c), K(c)_i = sum_j c_ij F_j (c antisymmetric)."""
        N = [dict(N0[i]) for i in range(NV)]
        for (i, j, n), v in c.items():
            N[i] = add(N[i], shift(self.Fi[j], n, v, p), p)
            N[j] = sub(N[j], shift(self.Fi[i], n, v, p), p)
        return N

# ----------------------------------------------------------------------------- rational lifting of an exact solution
def lift_solution(sysobj, D, N0, primes=PRIMES, max_primes=6, verbose=False):
    """CRT + rational reconstruction over several primes.  Returns (status, data)."""
    acc = None; piv0 = None
    for k, p in enumerate(primes[:max_primes]):
        out = sysobj.exact_consistent(D, N0, p=p, verbose=verbose)
        if not out['consistent']:
            return 'inconsistent', out
        if piv0 is None:
            piv0 = out['pivots']
        elif piv0 != out['pivots']:
            acc = None; piv0 = out['pivots']
        x = out['x']
        if acc is None:
            acc = (x, p)
        else:
            xs, m = acc
            xs = [core.crt_pair(a, m, b, p)[0] for a, b in zip(xs, x)]
            acc = (xs, m * p)
        xs, m = acc
        rat = [core.ratrecon(a, m) for a in xs]
        if all(v is not None for v in rat):
            c, Q = sysobj.unpack_exact(rat)
            N = sysobj.numerators_from_c(N0, c)
            ok, info = core.verify_witness(sysobj.F, sysobj.A, 1, D, N, Q)
            if ok:
                return 'solved', {'N': N, 'Q': Q, 'primes_used': k + 1, 'rank': out['rank'], 'nullity': out['nullity']}
    return 'lift-failed', {'primes_used': max_primes}

# ----------------------------------------------------------------------------- Hessian relaxation (Jacobian-type degrees)
class HessianRelaxation:
    """r(a) = Z^T vec(A Hess_F(a,a)) where Z spans the complement of J(F)_{2d}^{odd}: necessary for E2.
    Pair tensor W[:, k, l] = -Z^T vec(A m_k m_l F_{i_k i_l})  (symmetrised)."""

    def __init__(self, d, pD, basis=None, verbose=True, par_free=False, cofactors=None):
        t0 = time.time()
        F, A, q = load()
        self.F, self.A, self.d, self.pD = F, A, d, pD
        Fi, Fij = Fparts(F)
        if basis is None and cofactors is None:
            basis = jac_basis(d, pD) if not par_free else jac_basis(d, 0) + jac_basis(d, 1)
        if cofactors is None:
            cofactors = []
            for (i0, m, _) in basis:
                v = [{} for _ in range(NV)]; v[i0] = {m: Fr(1)}
                cofactors.append(v)
        self.basis = basis
        self.cofactors = cofactors
        K = len(cofactors); self.K = K
        self.par_free = par_free
        SB = space(2 * d, None if par_free else 1)
        self.S_B = SB
        # Jacobian generators of degree 2d (same parity) -> left null space Z
        gens = []
        for i in range(NV):
            pm = None if par_free else (1 - (1 if i >= 3 else 0)) % 2
            for m in monomials(2 * d - 5, pm):
                gens.append(shift(Fi[i], m, 1))
        Jm = np.zeros((SB.n, len(gens)))
        for col, g in enumerate(gens):
            for e, c in g.items():
                Jm[SB.idx[e], col] = float(c)
        U, s, Vt = np.linalg.svd(Jm, full_matrices=True)
        rank = int(np.sum(s > 1e-9 * s[0]))
        self.Z = U[:, rank:]                          # SB.n x nquot
        self.nquot = self.Z.shape[1]
        self.jac_rank = rank
        # pair tensor
        Wt = np.zeros((self.nquot, K, K))
        for k in range(K):
            ak = cofactors[k]
            for l in range(k, K):
                al = cofactors[l]
                Pk = {}
                for i in range(NV):
                    if not ak[i]: continue
                    for j in range(NV):
                        if not al[j]: continue
                        Pk = add(Pk, mul(A, mul(mul(ak[i], al[j]), Fij[i][j])))
                if not Pk:
                    continue
                v = np.zeros(SB.n)
                for e, c in Pk.items():
                    v[SB.idx[e]] = -float(c)
                w = self.Z.T @ v
                if k == l:
                    Wt[:, k, k] = w
                else:
                    Wt[:, k, l] = w; Wt[:, l, k] = w      # both (k,l) and (l,k) occur in sum_{IJ} a_I a_J F_IJ
        self.Wt = Wt
        if basis is None:
            SF = np.zeros((K, 0))
        elif par_free:
            SF = self._fmult_parfree(basis, d)
        else:
            SF = f_multiple_vectors(basis, d, pD)
        self.SF = SF
        self.P = complement(SF, K)
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        if verbose:
            print(f'[HessRelax] d={d} pD={pD} par_free={par_free} K={K} dimP={self.dimP} quotient={self.nquot} '
                  f'(rank J_{2*d} = {rank}) built in {self.build_time:.1f}s', file=sys.stderr, flush=True)

    def _fmult_parfree(self, basis, d):
        pos = {(i0, m): k for k, (i0, m, _) in enumerate(basis)}
        vecs = []
        for h in monomials(d - 6, None):
            v = np.zeros(len(basis))
            for i in range(NV):
                m = list(h); m[i] += 1; m = tuple(m)
                v[pos[(i, m)]] = 1.0 / 6
            vecs.append(v)
        return np.array(vecs).T

    def r_and_J(self, a):
        Wa = np.einsum('ikl,l->ik', self.Wt, a)       # nquot x K
        r = Wa @ a
        J = 2 * Wa
        return r, J

    def value(self, u):
        a = self.P @ u
        r, J = self.r_and_J(a)
        return float(r @ r), r, J @ self.P

def lm_sphere(resfun, u0, maxiter=200, ftol=1e-28, lam0=1e-3, verbose=False, deflate=None, step_tol=1e-14):
    """Levenberg-Marquardt on the unit sphere for a homogeneous residual.
    resfun(u) -> (f, r, J) with J = dr/du (unprojected).  deflate(u) -> (rho, grad rho) optional:
    minimises |r|^2 / rho^2 instead (rho = distance to known solution components)."""
    u = np.asarray(u0, float); u = u / np.linalg.norm(u)
    lam = lam0
    f, r, J = resfun(u)
    if deflate is not None:
        rho, g = deflate(u)
        r_eff = r / rho; J_eff = J / rho - np.outer(r, g) / rho ** 2
        f_eff = float(r_eff @ r_eff)
    else:
        r_eff, J_eff, f_eff = r, J, f
    hist = []
    nfev = 1
    for it in range(maxiter):
        # project Jacobian off the radial direction
        Jp = J_eff - np.outer(J_eff @ u, u)
        g = Jp.T @ r_eff
        H = Jp.T @ Jp
        dn = np.diag(H).copy(); dn[dn < 1e-30] = 1e-30
        improved = False
        for _ in range(12):
            try:
                delta = np.linalg.solve(H + lam * np.diag(dn), -g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            delta -= (delta @ u) * u
            un = u + delta; un /= np.linalg.norm(un)
            fn, rn, Jn = resfun(un); nfev += 1
            if deflate is not None:
                rho, gg = deflate(un)
                rn_eff = rn / rho; Jn_eff = Jn / rho - np.outer(rn, gg) / rho ** 2
                fn_eff = float(rn_eff @ rn_eff)
            else:
                rn_eff, Jn_eff, fn_eff = rn, Jn, fn
            if fn_eff < f_eff:
                step = float(np.linalg.norm(delta))
                u, f, r, J = un, fn, rn, Jn
                r_eff, J_eff, f_eff = rn_eff, Jn_eff, fn_eff
                lam = max(lam / 3, 1e-15)
                improved = True
                break
            lam *= 4
            if lam > 1e12:
                break
        hist.append(f)
        if verbose:
            print(f'  it {it} f={f:.3e} f_eff={f_eff:.3e} lam={lam:.1e}', file=sys.stderr, flush=True)
        if not improved or f < ftol or (improved and step < step_tol):
            break
    return {'u': u, 'f': f, 'f_eff': f_eff, 'iters': it + 1, 'nfev': nfev, 'hist': hist}

# ----------------------------------------------------------------------------- logging
def log_line(path, rec):
    rec = dict(rec); rec.setdefault('time', time.strftime('%Y-%m-%dT%H:%M:%S')); rec.setdefault('engine', 2)
    rec.setdefault('lane', LANE); rec.setdefault('wave', 2)
    with open(path, 'a') as f:
        f.write(json.dumps(rec, default=_json_default) + '\n')

def _json_default(o):
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, Fr): return str(o)
    return str(o)

def a_to_D(basis, a, tol=0.0):
    D = {}
    for k, (i0, m, Dk) in enumerate(basis):
        if abs(a[k]) > tol:
            D = add(D, scale(Dk, a[k]))
    return D

def rational_candidates(a, max_den=1000, tol=1e-7):
    """Rational vectors approximating a (normalised to max |a_k| = 1), engine-style."""
    a = np.asarray(a, float)
    j = int(np.argmax(np.abs(a)))
    a = a / a[j]
    cands = []
    fr = [Fr(float(v)).limit_denominator(max_den) for v in a]
    if max(abs(float(f) - v) for f, v in zip(fr, a)) < tol:
        cands.append(tuple(fr))
    for den in range(1, max_den + 1):
        v = a * den
        if np.max(np.abs(v - np.round(v))) < tol * den:
            cands.append(tuple(Fr(int(round(x)), den) for x in v))
            break
    uniq = []
    for c in cands:
        if c not in uniq:
            uniq.append(c)
    return uniq

# ----------------------------------------------------------------------------- general exact E1 lift mod p
def e1_particular_modp(D, p=P0, verbose=False):
    """Solve dF(N) = A*D mod p for N (degree d+1, parities forced by D).  Returns (N0 list mod p, info) or (None, info)."""
    from flint import nmod_mat
    t0 = time.time()
    F, A, q = load()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
    d = core.degree(D)
    pars = fixed_d.numerator_parities(D, A)
    degN = d + 1
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    cols = []
    colpolys = []
    for i in range(NV):
        for m in monomials(degN, pars[i]):
            cols.append((i, m)); colpolys.append(shift(Fi[i], m, 1, p))
    rhs = mul(Ap, Dp, p)
    SB = space(d + 6, 1)
    R, C = SB.n, len(cols)
    flat = np.zeros((R, C + 1), dtype=np.int64)
    for j, col in enumerate(colpolys):
        for e, v in col.items():
            flat[SB.idx[e], j] = v
    for e, v in rhs.items():
        flat[SB.idx[e], C] = v
    M = nmod_mat(R, C + 1, flat.ravel().tolist(), p)
    Mr, rank = M.rref()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0:
            j += 1
        pivots.append(j)
    info = {'rank': rank, 'rows': R, 'cols': C, 'seconds': round(time.time() - t0, 1), 'p': p}
    if pivots and pivots[-1] == C:
        info['e1'] = False
        return None, info
    info['e1'] = True
    N0 = [{} for _ in range(NV)]
    for i, j in enumerate(pivots):
        v = int(Mr[i, C])
        if v:
            ii, m = cols[j]
            N0[ii][m] = v
    if verbose:
        print(f'[e1] d={d} rows={R} cols={C} rank={rank} e1=True ({info["seconds"]}s)', file=sys.stderr, flush=True)
    return N0, info

def exact_fixed_D(D, p=P0, verbose=False, lift=False, e1_only=False):
    """Full exact (mod p) decision for one denominator via the reduced system.  Returns record dict."""
    F, A, q = load()
    t0 = time.time()
    rec = {'d': core.degree(D), 'parity': core.parity_of(D), 'terms': len(D), 'p': p}
    if rec['parity'] is None:
        rec['note'] = 'mixed Y-parity'; return rec
    divF, _ = jacobian.divisible_by(F, D, p)
    rec['F_divides_D'] = bool(divF)
    if divF:
        rec['note'] = 'F | D (illegal)'; return rec
    N0, info = e1_particular_modp(D, p, verbose=verbose)
    rec['e1'] = info['e1']; rec['e1_rank'] = info['rank']; rec['e1_cols'] = info['cols']; rec['e1_seconds'] = info['seconds']
    if N0 is None:
        rec['consistent'] = False; rec['note'] = 'E1 fails: A*D not in J(F) mod p'
        rec['seconds'] = round(time.time() - t0, 1); return rec
    if e1_only:
        rec['note'] = 'E1 only'; rec['seconds'] = round(time.time() - t0, 1); return rec
    S = ReducedE2(rec['d'], rec['parity'], [(D, N0)], verbose=verbose, build_float=False)
    out = S.exact_consistent(D, N0, p=p, verbose=verbose)
    rec.update({'consistent': out['consistent'], 'rank': out['rank'], 'nullity': out['nullity'],
                'rows': out['rows'], 'cols': out['cols'], 'e2_seconds': out['seconds'], 'solver': 'reduced-koszul-e2'})
    if out['consistent'] and lift:
        # lift: need exact rational N0 -> recompute N0 over Q is expensive; instead CRT over primes on the full
        # reduced system with mod-p N0 per prime is not coherent.  Use the engine's dense solver for the lift.
        status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
        rec['lift_status'] = status
        if status == 'solved':
            rec['witness'] = {'d': rec['d'], 'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}
    rec['seconds'] = round(time.time() - t0, 1)
    return rec

# ----------------------------------------------------------------------------- UV-bidegree pieces of the Jacobian-type family
def uv_bidegree_family(d, pD, a):
    """Swap-symmetric (pD=0) / antisymmetric (pD=1) Jacobian-type denominators of (U,V)-bidegree (a, d-a) + (d-a, a):
    D = f_i(U) * m  ±  swap,  m a UV-monomial of bidegree (a-5, d-a), i = 0..2  (plus the mirror pieces when a-5 < 0
    is handled by the caller choosing a >= 5).  Returns list of (label, D_xy, N0_xy) with N0 = A * cofactor vector."""
    F, A, q = load()
    f = f_ternary()
    fi = [pderiv(f, i) for i in range(3)]
    fam = []
    sign = 1 if pD == 0 else -1
    for i in range(3):
        for mU in monomials(a - 5, None):
            if mU[3] or mU[4] or mU[5]:
                continue
            for mV in monomials(d - a, None):
                if mV[0] or mV[1] or mV[2]:
                    continue
                m = (mU[0], mU[1], mU[2], mV[3], mV[4], mV[5])
                m_xy = uv_to_xy({m: Fr(1)})
                ms_xy = uv_to_xy({(m[3], m[4], m[5], m[0], m[1], m[2]): Fr(1)})
                term = mul(uv_to_xy(fi[i]), m_xy)
                D = add(term, scale(mul(uv_to_xy(swapUV(fi[i])), ms_xy), sign))
                # cofactors: f_i(U) = F_i + F_{i+3} ;  f_i(V) = F_i - F_{i+3}
                avec = [{} for _ in range(NV)]
                avec[i] = add(avec[i], m_xy); avec[i + 3] = add(avec[i + 3], m_xy)
                avec[i] = add(avec[i], scale(ms_xy, sign)); avec[i + 3] = sub(avec[i + 3], scale(ms_xy, sign))
                N0 = [mul(A, avec[j]) if avec[j] else {} for j in range(NV)]
                fam.append((f'f_{i}(U)*{core.to_str({m: Fr(1)})}{"+" if sign > 0 else "-"}swap', D, N0, avec))
    return fam

def check_family_e1(fam):
    """exact check dF(N0) = A D for each member (Fractions)."""
    F, A, q = load()
    Fi = [pderiv(F, i) for i in range(NV)]
    bad = 0
    for lab, D, N0 in fam:
        s = {}
        for i in range(NV):
            s = add(s, mul(N0[i], Fi[i]))
        if sub(s, mul(A, D)):
            bad += 1
    return bad


def uv_bidegree_profile(D_xy, d):
    """Fraction of the squared coefficient norm of D (XY frame) in each (U,V)-bidegree a + mirror (a <= d/2)."""
    Duv = xy_to_uv({e: Fr(c).limit_denominator(10**9) if not isinstance(c, Fr) else c for e, c in D_xy.items()})
    tot = 0.0; by = {}
    for e, c in Duv.items():
        a = e[0] + e[1] + e[2]; key = min(a, d - a)
        v = float(c) ** 2; by[key] = by.get(key, 0.0) + v; tot += v
    return {str(k): round(v / tot, 4) for k, v in sorted(by.items())} if tot else {}
