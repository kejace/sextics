"""w2lib -- wave-2 EXCLUSION/RIGOR utilities for lane `diagonal-permutation` (q = 1).

Everything is finite linear algebra over GF(p) (python-flint nmod_mat + numpy) unless stated.

Notation.  F, A from the instance; J = J(F) = (F_0..F_5); R = S/J (Jacobian ring, Koszul Hilbert
function since F is smooth).  A vector field w = (w_0..w_5) acts by w(P) = sum_i w_i dP/dz_i.
Hessian form  H(u,v) = sum_{i,j} u_i v_j F_{ij}.

Lemma (E2 modulo J).  If (D,N,Q) is a certificate for A at q = 1 then N(F) = A D + F H and
N = u + K + (H/6) E  with u(F) = A D, K Koszul (K(F) = 0), E the Euler field.  Substituting,
    A D^2 - D N(F) + F (D div N - N(D)) = F ( D div u - u(D) + D div K - K(D) ),
so a certificate exists iff   Psi(D;u) + (D div K - K(D)) = F Q   is solvable in (K, Q), where
    Psi(D;u) := D div u - u(D).
For Koszul K = sum_{i<j} c_ij (F_i d_j - F_j d_i):  D div K - K(D) in J (antisymmetry of c, symmetry
of the Hessian), hence the NECESSARY condition   Psi(D;u) in J(F)   ("relaxation").  Psi mod J does not
depend on the choice of u (two choices differ by a Koszul field) and is invariant under D -> D + g F.
For D = w(F):  u = A w  and  Psi = -A H(w,w) mod J.
For D = lam D8 + w(F) with A D8 = u8(F):  u = lam u8 + A w and
    Psi = lam^2 (D8 div u8 - u8(D8)) + lam (D8 div(A w) - H(u8,w) - A w(D8)) - A H(w,w)   (mod J).

Reduced fixed-D system (E2').  Unknowns c_ij (deg d-4, 15 of them) and Q (deg 2d-6); equation
    sum_{i<j} [ c_ij (F_i D_j - F_j D_i) + D (F_j d_i c_ij - F_i d_j c_ij) ] - F Q = -Psi(D;u).
It is equivalent to the full fixed-D system (N = u + K) and about 3x smaller.

Rigor conventions.  Kernel dimension mod p is >= the rational one (a modular rank is a lower bound
for the rational rank): zero kernel mod p is a proof over Q.  Inconsistency mod p of an inhomogeneous
system proves inconsistency over Q only together with rank_p(M) = rank_Q(M), which we certify by lifting
a basis of the mod-p nullspace of M to exact rational vectors (see certify_rank).  Groebner results:
if the mod-p ideal is m-primary then so is the rational ideal (Hilbert-function argument), provided the
mod-p equations are the reduction of the rational equation lattice (true when dim J_n mod p has the
Koszul value, which is checked).
"""
from __future__ import annotations
import sys, os, time, json, itertools, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
from fractions import Fraction
from flint import nmod_mat, fmpq_mat, fmpq, fmpz_mat
from hodge import core, jacobian, fixed_d
from hodge.core import monomials, mul, add, sub, scale, shift, pderiv, reduce_mod, parity_of, degree, NV

IID = os.environ.get('W2IID', 'diagonal-permutation')
LANE = '/tmp/hodge/results/diagonal-permutation'
PKT = f'{LANE}/wave2-excl'
LOG = f'{LANE}/wave2-excl.log.jsonl'
P1, P2, P3 = 1000000007, 998244353, 2147483647   # P3 = wave-1 prime; P1/P2 fresh
F, A, q, rec = core.load_instance(IID)
assert q == 1
JOB = os.environ.get('W2JOB', '')


def log(exp, **kw):
    entry = {'exp': exp, 'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'job': JOB}
    entry.update(kw)
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry, default=str) + '\n')
    return entry


def koszul_hilbert(n, par=None):
    """dim R_n (Jacobian ring of a smooth sextic in 6 vars) by parity: product of (1+t+..+t^4)^3 (1+ty+..)^3."""
    one = [1, 1, 1, 1, 1]
    # polynomial in t with y-grading: dict (deg, ypar) -> count
    cnt = {(0, 0): 1}
    for k in range(6):
        new = {}
        for (dg, yp), c in cnt.items():
            for e in range(5):
                key = (dg + e, (yp + (e if k >= 3 else 0)) % 2)
                new[key] = new.get(key, 0) + c
        cnt = new
    if par is None:
        return cnt.get((n, 0), 0) + cnt.get((n, 1), 0)
    return cnt.get((n, par), 0)


# ----------------------------------------------------------------------------- mod-p numpy helpers
def matmul_mod(Aa, Bb, p):
    """(Aa @ Bb) mod p for int64 arrays with entries in [0,p), p < 2^31, using a 16-bit split."""
    Aa = np.ascontiguousarray(Aa, dtype=np.int64); Bb = np.ascontiguousarray(Bb, dtype=np.int64)
    if Aa.shape[1] == 0:
        return np.zeros((Aa.shape[0], Bb.shape[1]), dtype=np.int64)
    out = np.zeros((Aa.shape[0], Bb.shape[1]), dtype=np.int64)
    B_lo = Bb & 0xFFFF
    B_hi = Bb >> 16
    # chunk rows of Aa so that intermediate sums stay < 2^63: each product < 2^31 * 2^16 = 2^47, so
    # up to 2^16 terms are safe; chunk the inner dimension accordingly.
    K = Aa.shape[1]
    step = 1 << 15
    for k0 in range(0, K, step):
        a = Aa[:, k0:k0 + step]
        lo = a @ B_lo[k0:k0 + step]
        hi = a @ B_hi[k0:k0 + step]
        out = (out + lo % p + ((hi % p) << 16) % p) % p
    return out


def to_nmod(arr, p):
    arr = np.ascontiguousarray(arr, dtype=np.int64)
    return nmod_mat(arr.shape[0], arr.shape[1], arr.ravel().tolist(), p)


def from_nmod(M):
    R, C = M.nrows(), M.ncols()
    ents = M.entries()
    return np.fromiter((int(x) for x in ents), dtype=np.int64, count=R * C).reshape(R, C)


def rref_np(arr, p, cols=None):
    """rref over GF(p) via flint.  Returns (reduced rows restricted to `cols` (default all), pivots).

    Pivots are found by a monotone scan (O(C) element accesses); the submatrix is extracted by element
    access, avoiding the conversion of the whole (possibly huge) rref matrix."""
    R, C = arr.shape
    if R == 0:
        return np.zeros((0, C if cols is None else len(cols)), dtype=np.int64), []
    M = to_nmod(arr, p)
    return _rref_extract(M, R, C, cols)


def _rref_extract(M, R, C, cols=None):
    Mr, rank = M.rref()
    piv = []
    j = 0
    for i in range(rank):
        while j < C and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
        j += 1
    if cols is None:
        red = from_nmod(Mr)[:rank]
    else:
        red = np.zeros((rank, len(cols)), dtype=np.int64)
        for a, c in enumerate(cols):
            for i in range(rank):
                red[i, a] = int(Mr[i, c])
    return red, piv


def rref_sparse(polys, idx, C, p, cols=None):
    """rref of the matrix whose rows are the sparse polys (coordinates idx), built by element assignment."""
    R = len(polys)
    if R == 0:
        return np.zeros((0, C if cols is None else len(cols)), dtype=np.int64), []
    M = nmod_mat(R, C, p)
    for r, P in enumerate(polys):
        for e, c in P.items():
            j = idx.get(e)
            if j is not None:
                M[r, j] = int(c) % p
    return _rref_extract(M, R, C, cols)


def rank_np(arr, p):
    if arr.shape[0] == 0 or arr.shape[1] == 0:
        return 0
    return to_nmod(arr, p).rank()


def nullspace_np(arr, p):
    """Basis of the right nullspace {x : arr x = 0} as rows of an array (nullity x C)."""
    R, C = arr.shape
    if R == 0:
        return np.eye(C, dtype=np.int64)
    X, nullity = to_nmod(arr, p).nullspace()
    Xa = from_nmod(X)[:, :nullity]
    return np.ascontiguousarray(Xa.T)


# ----------------------------------------------------------------------------- polynomial <-> vector
def poly_to_vec(P, idx, C, p):
    v = np.zeros(C, dtype=np.int64)
    for e, c in P.items():
        v[idx[e]] = int(c) % p
    return v


def polys_to_mat(polys, idx, C, p):
    M = np.zeros((len(polys), C), dtype=np.int64)
    for r, P in enumerate(polys):
        for e, c in P.items():
            M[r, idx[e]] = int(c) % p
    return M


def vec_to_poly(v, mons):
    return {mons[j]: int(v[j]) for j in np.nonzero(v)[0]}


# ----------------------------------------------------------------------------- symmetry (01) and tau
def swap01(e):
    return (e[1], e[0], e[2], e[4], e[3], e[5])


def act_swap01(P):
    return {swap01(e): c for e, c in P.items()}


def act_cyc(e):
    """three-cycle (012): X0->X1->X2->X0 simultaneously with Y."""
    return (e[2], e[0], e[1], e[5], e[3], e[4])


def act_perm(P, perm):
    """Apply the simultaneous permutation perm (tuple of images of 0,1,2) to P: X_i -> X_perm[i]."""
    out = {}
    for e, c in P.items():
        ne = [0] * 6
        for i in range(3):
            ne[perm[i]] += e[i]
            ne[3 + perm[i]] += e[3 + i]
        out[tuple(ne)] = out.get(tuple(ne), 0) + c
    return {e: c for e, c in out.items() if c}


PERMS = [(0, 1, 2), (1, 0, 2), (2, 1, 0), (0, 2, 1), (1, 2, 0), (2, 0, 1)]
PERM_SIGN = [1, -1, -1, -1, 1, 1]


def reynolds_S3(P, sign=False, p=None):
    """S3-average (sign-twisted if sign) of P."""
    out = {}
    for perm, sg in zip(PERMS, PERM_SIGN):
        out = add(out, scale(act_perm(P, perm), sg if sign else 1, p), p)
    return out


class Sector:
    """Monomial basis of S_n^{par} restricted to a (01)-eigen-sector (sigma = +1 / -1 / None).

    Coordinates of a sigma-eigen polynomial P are its coefficients on orbit representatives."""

    def __init__(self, n, par=None, sigma=None):
        self.n, self.par, self.sigma = n, par, sigma
        mons = monomials(n, par)
        if sigma is None:
            self.reps = mons
        else:
            reps = []
            for m in mons:
                sm = swap01(m)
                if sm == m:
                    if sigma == 1:
                        reps.append(m)
                elif m < sm:
                    reps.append(m)
            self.reps = reps
        self.idx = {m: i for i, m in enumerate(self.reps)}
        self.C = len(self.reps)

    def coords(self, P, p):
        """Coordinates of P (assumed in the sector) -- coefficients on representatives."""
        v = np.zeros(self.C, dtype=np.int64)
        for e, c in P.items():
            j = self.idx.get(e)
            if j is not None:
                v[j] = int(c) % p
        return v

    def project(self, P, p):
        """Project an arbitrary P onto the sector (P + sigma*swap(P))/2 (or identity if sigma None), as a poly."""
        if self.sigma is None:
            return reduce_mod(P, p) if p else dict(P)
        Q = act_swap01(P)
        if self.sigma == -1:
            Q = scale(Q, -1, p)
        S = add(P, Q, p)
        inv2 = pow(2, -1, p)
        return scale(S, inv2, p)

    def poly(self, v, p):
        """Polynomial from sector coordinates."""
        P = {}
        for j in np.nonzero(v)[0]:
            m = self.reps[j]
            c = int(v[j]) % p
            P[m] = c
            if self.sigma is not None:
                sm = swap01(m)
                if sm != m:
                    P[sm] = (c * self.sigma) % p
        return P


# ----------------------------------------------------------------------------- Jacobian reducer
class Red:
    """Row-reduced J(F)_n in a sector (parity par, (01)-sign sigma) mod p: normal forms, standard monomials."""

    def __init__(self, n, p, par=None, sigma=None, verbose=False):
        t0 = time.time()
        self.n, self.p, self.par, self.sigma = n, p, par, sigma
        self.sec = Sector(n, par, sigma)
        C = self.sec.C
        gens = jacobian.jacobian_generators(F, n, p, par)
        rows = []
        for (_, _, P) in gens:
            Pp = self.sec.project(P, p)
            if Pp:
                rows.append(Pp)
        # first pass: pivots only (cheap scan), second: extract the std columns
        if rows:
            M = nmod_mat(len(rows), C, p)
            for r, P in enumerate(rows):
                for e, c in P.items():
                    j = self.sec.idx.get(e)
                    if j is not None:
                        M[r, j] = int(c) % p
            Mr, rank = M.rref()
            piv = []
            j = 0
            for i in range(rank):
                while j < C and int(Mr[i, j]) == 0:
                    j += 1
                piv.append(j)
                j += 1
            self.piv = piv
            self.rank = rank
            pivset = set(piv)
            self.std = [j for j in range(C) if j not in pivset]
            red = np.zeros((rank, len(self.std)), dtype=np.int64)
            for a, c in enumerate(self.std):
                for i in range(rank):
                    red[i, a] = int(Mr[i, c])
            self.Jr_std = red
            del Mr, M
        else:
            self.piv, self.rank = [], 0
            self.std = list(range(C))
            self.Jr_std = np.zeros((0, C), dtype=np.int64)
        self.std_mons = [self.sec.reps[j] for j in self.std]
        self.ring_dim = len(self.std)
        self.t = time.time() - t0
        if verbose:
            print(f'[Red n={n} par={par} sigma={sigma} p={p}] C={C} gens={len(rows)} rank(J)={self.rank} '
                  f'dim R={self.ring_dim} ({self.t:.1f}s)', file=sys.stderr, flush=True)

    def nf_rows(self, V):
        """Normal forms (standard-monomial coordinates) of the rows of V (array m x C in sector coords)."""
        V = np.ascontiguousarray(V, dtype=np.int64)
        if self.rank == 0:
            return np.ascontiguousarray(V[:, self.std])
        Vp = V[:, self.piv]
        return (V[:, self.std] - matmul_mod(Vp, self.Jr_std, self.p)) % self.p

    def nf_polys(self, polys):
        V = polys_to_mat([self.sec.project(P, self.p) for P in polys], self.sec.idx, self.sec.C, self.p)
        return self.nf_rows(V)

    def nf(self, P):
        return self.nf_polys([P])[0]

    def in_J(self, P):
        return not np.any(self.nf(P))


# ----------------------------------------------------------------------------- E1 kernel in a sector
def e1_kernel(d, p, par=0, sigma=None, verbose=False, red=None, redA=None):
    """Kernel of [A]: R_d^{par,sigma} -> R_{d+6}^{1-par,-sigma} (A is Y-odd and (01)-anti-invariant).

    Returns dict(dim, basis (polys in standard monomials of R_d), ring_dim, ...)."""
    t0 = time.time()
    if red is None:
        red = Red(d, p, par, sigma, verbose)
    sig2 = None if sigma is None else -sigma
    if redA is None:
        redA = Red(d + 6, p, 1 - par, sig2, verbose)
    Ap = reduce_mod(A, p)
    prods = [shift(Ap, m, 1, p) for m in red.std_mons]
    M = redA.nf_polys(prods)                 # rows: images of standard monomials
    rank = rank_np(M, p)
    ker = nullspace_np(np.ascontiguousarray(M.T), p)   # left kernel: rows = coefficient vectors over red.std_mons
    basis = [vec_to_poly(v, red.std_mons) for v in ker]
    out = {'d': d, 'p': p, 'par': par, 'sigma': sigma, 'ring_dim': red.ring_dim, 'target_dim': redA.ring_dim,
           'rank': int(rank), 'dim': int(ker.shape[0]), 'basis': basis, 'seconds': round(time.time() - t0, 1),
           'jac_dim_d': red.rank, 'jac_dim_d6': redA.rank}
    if verbose:
        print(f'[e1_kernel d={d} par={par} sigma={sigma} p={p}] R_d={red.ring_dim} -> R_{d+6}={redA.ring_dim} '
              f'rank={rank} ker={ker.shape[0]} ({out["seconds"]}s)', file=sys.stderr, flush=True)
    return out


# ----------------------------------------------------------------------------- vector fields
def field_basis(d, eps, euler_complement=True):
    """Basis fields w = m d_i with deg m = d-5 such that w(F) has Y-parity eps (eps in {0,1}).

    w_{X_i} must have parity eps (F_{X_i} even), w_{Y_i} parity 1-eps (F_{Y_i} odd).
    If euler_complement, drop the fields X0^a.. m d_{X0} with X0 | m (a complement of {g E})."""
    out = []
    for i in range(NV):
        par_i = eps if i < 3 else (1 - eps) % 2
        for m in monomials(d - 5, par_i):
            if euler_complement and i == 0 and m[0] > 0:
                continue
            out.append((i, m))
    return out


def field_apply(w, P, p):
    """w(P) for a field given as list of 6 polys."""
    out = {}
    for i in range(NV):
        if w[i]:
            out = add(out, mul(w[i], pderiv(P, i, p), p), p)
    return out


def field_div(w, p):
    out = {}
    for i in range(NV):
        if w[i]:
            out = add(out, pderiv(w[i], i, p), p)
    return out


def hessian(p):
    Fp = reduce_mod(F, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    Fij = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
    return Fp, Fi, Fij


def solve_E1(D, p, verbose=False):
    """Find u (6 polys of degree d+1) with u(F) = A D mod p, or None.  Unknown parities follow D's parity
    (if pure); mixed D: all monomials."""
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); Dp = reduce_mod(D, p)
    d = degree(D)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    pD = parity_of(Dp)
    cols, colpolys = [], []
    for i in range(NV):
        if pD is None:
            pars = [0, 1]
        else:
            pars = [(1 + pD - (0 if i < 3 else 1)) % 2]   # parity(u_i F_i) = parity(A D) = 1 + pD
        for par in pars:
            for m in monomials(d + 1, par):
                cols.append((i, m)); colpolys.append(shift(Fi[i], m, 1, p))
    rhs = mul(Ap, Dp, p)
    rows = {}
    for P in colpolys + [rhs]:
        for e in P:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    M = np.zeros((R, C + 1), dtype=np.int64)
    for j, P in enumerate(colpolys):
        for e, v in P.items():
            M[rows[e], j] = v
    for e, v in rhs.items():
        M[rows[e], C] = v
    red, piv = rref_np(M, p)
    if piv and piv[-1] == C:
        return None
    u = [{} for _ in range(NV)]
    for r, j in enumerate(piv):
        v = int(red[r, C])
        if v:
            i, m = cols[j]
            u[i][m] = v
    return u


def psi(D, u, p):
    """Psi(D;u) = D div u - u(D) mod p (as a polynomial, not reduced mod J)."""
    Dp = reduce_mod(D, p)
    return sub(mul(Dp, field_div(u, p), p), field_apply(u, Dp, p), p)


# ----------------------------------------------------------------------------- relaxation systems
class Relaxation:
    """Quadratic system  Psi(x) = 0 in R_{2d}  for  D = x_0 * D8 (optional) + w(F),  w = sum_k x_k w^k.

    unknowns: optional D8 coefficient first, then the basis fields of the requested parities
    (eps list, e.g. [0], [1] or [0,1] for mixed).  Equations: coordinates of NF(Psi) in R_{2d}^{even} and
    R_{2d}^{odd} (both reducers), each a quadratic form in x.  Stored as a dense matrix  Cq  of shape
    (#pairs, #equations) over the pair list [(k,l), k<=l]."""

    def __init__(self, d, eps_list, p, with_D8=None, euler_complement=True, verbose=False, D8=None, u8=None):
        t0 = time.time()
        self.d, self.p, self.eps_list = d, p, list(eps_list)
        Fp, Fi, Fij = hessian(p)
        Ap = reduce_mod(A, p)
        self.fields = []
        for eps in self.eps_list:
            self.fields += [(eps, i, m) for (i, m) in field_basis(d, eps, euler_complement)]
        self.with_D8 = with_D8 is not None
        names = []
        if self.with_D8:
            self.D8 = reduce_mod(with_D8, p)
            self.u8 = u8 if u8 is not None else solve_E1(with_D8, p)
            assert self.u8 is not None, 'D8 does not satisfy E1'
            names.append('lam')
        for (eps, i, m) in self.fields:
            names.append(f'w{i}_' + ''.join(f'{e}' for e in m) + f'_e{eps}')
        self.names = names
        n = len(names)
        self.n = n
        # unknown k -> (kind, data)
        unk = []
        if self.with_D8:
            unk.append(('D8', None))
        for (eps, i, m) in self.fields:
            unk.append(('w', (i, m)))
        self.unk = unk
        # quadratic contributions P_{kl} (k<=l):  Psi = sum_{k<=l} x_k x_l P_{kl}
        pairs = []
        polys = []
        # precompute pieces for the D8 part
        if self.with_D8:
            D8, u8 = self.D8, self.u8
            Psi8 = psi(D8, u8, p)
            D8i = [pderiv(D8, i, p) for i in range(NV)]
        for k in range(n):
            for l in range(k, n):
                kk, kd = unk[k]; ll, ld = unk[l]
                if kk == 'D8' and ll == 'D8':
                    P = Psi8
                elif kk == 'D8' and ll == 'w':
                    i, m = ld
                    # C(w) = D8 div(A w) - H(u8, w) - A w(D8)   for basis field w = m d_i
                    Aw = mul(Ap, {m: 1}, p)                      # (A w)_i
                    divAw = pderiv(Aw, i, p)
                    term1 = mul(D8, divAw, p)
                    Hu8w = {}
                    for j in range(NV):
                        if u8[j]:
                            Hu8w = add(Hu8w, mul(mul(u8[j], {m: 1}, p), Fij[j][i], p), p)
                    term3 = mul(Aw, D8i[i], p)
                    P = sub(sub(term1, Hu8w, p), term3, p)
                else:
                    i, m = kd; j, m2 = ld
                    mm = mul({m: 1}, {m2: 1}, p)
                    P = scale(mul(mul(Ap, mm, p), Fij[i][j], p), p - 1, p)   # -A m m' F_ij
                    if k != l:
                        P = scale(P, 2, p)
                pairs.append((k, l)); polys.append(P)
        self.pairs = pairs
        # reduce modulo J in degree 2d, both parities
        n2 = 2 * d
        self.reds = {}
        eqs = []
        self.eq_parity = []
        for par in (0, 1):
            red = Red(n2, p, par, None, verbose)
            self.reds[par] = red
            # split polys by parity (a poly may have both parts in the mixed case)
            parts = []
            for P in polys:
                parts.append({e: c for e, c in P.items() if core.yparity(e) == par})
            if red.ring_dim == 0:
                continue
            V = polys_to_mat(parts, red.sec.idx, red.sec.C, p)
            NFm = red.nf_rows(V)                 # (#pairs x ring_dim)
            if np.any(NFm):
                eqs.append(NFm)
                self.eq_parity += [par] * NFm.shape[1]
        self.Cq = np.concatenate(eqs, axis=1) if eqs else np.zeros((len(pairs), 0), dtype=np.int64)
        # drop zero columns, and reduce to independent equations
        nz = np.any(self.Cq != 0, axis=0)
        self.Cq = np.ascontiguousarray(self.Cq[:, nz])
        self.eq_parity = [pp for pp, keep in zip(self.eq_parity, nz) if keep]
        self.neq_raw = self.Cq.shape[1]
        self.t = time.time() - t0
        if verbose:
            print(f'[Relaxation d={d} eps={self.eps_list} D8={self.with_D8}] unknowns={n} pairs={len(pairs)} '
                  f'equations={self.neq_raw} ({self.t:.1f}s)', file=sys.stderr, flush=True)

    def independent_equations(self):
        """Row-reduce the equations (as vectors over pairs): returns array (neq_indep x #pairs)."""
        red, piv = rref_np(np.ascontiguousarray(self.Cq.T), self.p)
        return red

    def write_msolve(self, path, independent=True):
        E = self.independent_equations() if independent else self.Cq.T
        p = self.p
        with open(path, 'w') as f:
            f.write(','.join(self.names) + '\n')
            f.write(f'{p}\n')
            lines = []
            for r in range(E.shape[0]):
                terms = []
                for t in np.nonzero(E[r])[0]:
                    k, l = self.pairs[t]
                    c = int(E[r, t])
                    mon = f'{self.names[k]}^2' if k == l else f'{self.names[k]}*{self.names[l]}'
                    terms.append(f'{c}*{mon}')
                lines.append('+'.join(terms))
            f.write(',\n'.join(lines) + '\n')
        return E.shape[0]

    def evaluate(self, x):
        """Evaluate all raw equations at the integer vector x (mod p): returns residual vector."""
        x = np.array(x, dtype=np.int64) % self.p
        vals = np.array([(x[k] * x[l]) % self.p for (k, l) in self.pairs], dtype=np.int64)
        return matmul_mod(vals.reshape(1, -1), self.Cq, self.p)[0]

    def field_from_x(self, x):
        """The field w (6 polys) and lam from an unknown vector."""
        w = [{} for _ in range(NV)]
        lam = 0
        for k, (kind, data) in enumerate(self.unk):
            v = int(x[k]) % self.p
            if not v:
                continue
            if kind == 'D8':
                lam = v
            else:
                i, m = data
                w[i][m] = (w[i].get(m, 0) + v) % self.p
        return lam, w

    def D_from_x(self, x):
        lam, w = self.field_from_x(x)
        Fp = reduce_mod(F, self.p)
        D = field_apply(w, Fp, self.p)
        if lam:
            D = add(D, scale(self.D8, lam, self.p), self.p)
        return D


# ----------------------------------------------------------------------------- reduced fixed-D system E2'
class ReducedE2:
    """Linear system in (c_ij, Q):  sum_{i<j} [c_ij (F_i D_j - F_j D_i) + D (F_j d_i c_ij - F_i d_j c_ij)] - F Q = -Psi(D;u).

    Optional (01)-sector restriction sigma (D must be a (01)-eigenvector with sign sigma_D; then the
    unknowns are restricted to the compatible sector and the equations projected).  Unknown parities
    are fixed by D's parity when D is pure."""

    def __init__(self, D, p, u=None, sigma=None, verbose=False):
        t0 = time.time()
        self.p = p
        Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); Dp = reduce_mod(D, p)
        self.D = Dp
        d = degree(D); self.d = d
        Fi = [pderiv(Fp, i, p) for i in range(NV)]
        Di = [pderiv(Dp, i, p) for i in range(NV)]
        if u is None:
            u = solve_E1(D, p)
            assert u is not None, 'E1 fails for this D'
        self.u = u
        rhs = scale(psi(Dp, u, p), p - 1, p)       # -Psi
        pD = parity_of(Dp)
        self.pD = pD
        # For pure D the equation parity is 1 + pD (Psi has the parity of A D^2 / F); unknown columns of the
        # wrong parity are dropped below.
        cols, colpolys = [], []
        for i in range(NV):
            for j in range(i + 1, NV):
                base = sub(mul(Fi[i], Di[j], p), mul(Fi[j], Di[i], p), p)   # F_i D_j - F_j D_i
                DFj = mul(Dp, Fi[j], p); DFi = mul(Dp, Fi[i], p)
                for m in monomials(d - 4):
                    P = shift(base, m, 1, p)
                    if m[i]:
                        m2 = list(m); m2[i] -= 1
                        P = add(P, shift(DFj, tuple(m2), m[i], p), p)
                    if m[j]:
                        m2 = list(m); m2[j] -= 1
                        P = sub(P, shift(DFi, tuple(m2), m[j], p), p)
                    if P:
                        cols.append(('c', (i, j), m)); colpolys.append(P)
        for m in monomials(2 * d - 6):
            cols.append(('Q', None, m)); colpolys.append(shift(Fp, m, p - 1, p))
        # parity filter: keep columns whose polynomial parity matches the rhs parity (if pure)
        eqpar = parity_of(rhs) if rhs else None
        if pD is not None:
            eqpar = (1 + pD) % 2
        if eqpar is not None:
            keep = [k for k, P in enumerate(colpolys) if parity_of(P) == eqpar]
            cols = [cols[k] for k in keep]; colpolys = [colpolys[k] for k in keep]
        self.eqpar = eqpar
        # (01)-sector restriction
        self.sigma = sigma
        if sigma is not None:
            # D is sigma_D-eigen; Psi has sign sigma_D * sign(A)= -sigma_D ... we project equations onto the
            # sector of the rhs and symmetrize the unknown columns accordingly (orbit sums with sign).
            sD = self._sign01(Dp)
            assert sD is not None, 'D is not a (01)-eigenvector'
            srhs = self._sign01(rhs) if rhs else -sD
            self.srhs = srhs
            # group columns into (01)-orbits with the sign making the column polynomial srhs-eigen
            seen = {}
            newcols, newpolys = [], []
            for (kind, ij, m), P in zip(cols, colpolys):
                key = self._col_key(kind, ij, m)
                skey = self._col_key(*self._swap_col(kind, ij, m))
                if key in seen or skey in seen:
                    continue
                seen[key] = True
                Ps = act_swap01(P)
                Pc = add(P, scale(Ps, srhs, p), p)
                if Pc:
                    newcols.append((kind, ij, m, 'orb')); newpolys.append(Pc)
            cols, colpolys = newcols, newpolys
            sec = Sector(2 * d, eqpar, srhs)
        else:
            sec = Sector(2 * d, eqpar, None)
        self.sec = sec
        self.cols = cols
        C = len(cols)
        M = np.zeros((sec.C, C + 1), dtype=np.int64)
        for k, P in enumerate(colpolys):
            for e, v in P.items():
                j = sec.idx.get(e)
                if j is not None:
                    M[j, k] = v
        for e, v in rhs.items():
            j = sec.idx.get(e)
            if j is not None:
                M[j, C] = v
        # drop zero rows
        M = M[np.any(M != 0, axis=1)]
        self.M = M
        self.nrows, self.ncols = M.shape[0], C
        self.t_build = time.time() - t0
        if verbose:
            print(f'[ReducedE2 d={d} sigma={sigma}] rows={self.nrows} cols={C} (c: {sum(1 for c in cols if c[0]=="c")}, '
                  f'Q: {sum(1 for c in cols if c[0]=="Q")}) built in {self.t_build:.1f}s', file=sys.stderr, flush=True)

    def _sign01(self, P):
        Q = act_swap01(P)
        if Q == P:
            return 1
        if len(Q) == len(P) and all(Q.get(e, 0) == (-c) % self.p for e, c in P.items()):
            return -1
        return None

    @staticmethod
    def _col_key(kind, ij, m):
        return (kind, ij, m)

    @staticmethod
    def _swap_col(kind, ij, m):
        sm = swap01(m)
        if kind == 'Q':
            return ('Q', None, sm)
        i, j = ij
        perm = {0: 1, 1: 0, 2: 2, 3: 4, 4: 3, 5: 5}
        a, b = perm[i], perm[j]
        if a > b:
            a, b = b, a
        return ('c', (a, b), sm)

    def solve(self, compressed=True, k_extra=40, seed=1, verbose=False):
        """Consistency and ranks mod p.  Returns dict(consistent, rank, rank_aug, nullity, k)."""
        t0 = time.time()
        M = self.M; p = self.p; C = self.ncols
        if compressed and M.shape[0] > C + 1 + k_extra:
            k = C + 1 + k_extra
            rng = np.random.default_rng(seed)
            Rm = rng.integers(0, p, size=(k, M.shape[0]), dtype=np.int64)
            Mc = matmul_mod(Rm, M, p)
        else:
            Mc = M; k = M.shape[0]
        rank_aug = rank_np(Mc, p)
        rank0 = rank_np(np.ascontiguousarray(Mc[:, :C]), p)
        out = {'consistent': rank_aug == rank0, 'rank': int(rank0), 'rank_aug': int(rank_aug), 'nullity': int(C - rank0),
               'rows': int(M.shape[0]), 'cols': int(C), 'k': int(k), 'p': p, 'seconds': round(time.time() - t0, 1)}
        if verbose:
            print(f'[ReducedE2.solve p={p}] {out}', file=sys.stderr, flush=True)
        return out


# ----------------------------------------------------------------------------- served-target space (homogeneous)
def served_space(D, p, k_extra=40, seed=1, verbose=False):
    """Space of Y-odd degree-6 targets served by D mod p (compressed homogeneous system), as in wave-1 hlib.
    Returns (dim, basis polys, contains_A).  NOTE: mod p this is an upper bound for the rational served space
    only when the compression is exact; it is used as a screen and for the sum principle (see README)."""
    A0 = {monomials(6, 1)[0]: 1}
    S = fixed_d.FixedDSystem(F, A0, q, D, p, with_target_space=True, verbose=verbose)
    M = fixed_d._sparse_matrix(S, with_rhs=False)
    C = S.ncols
    k = min(S.nrows, C + k_extra)
    Cm = fixed_d._compress(M, p, k, seed)
    Mn = nmod_mat(k, C, Cm.ravel().tolist(), p)
    Xn, nullity = Mn.nullspace()
    a_idx = [j for j, c in enumerate(S.cols) if c[0] == 'A']
    if nullity == 0:
        return 0, [], False
    proj = nmod_mat(len(a_idx), nullity, [int(Xn[j, kk]) for j in a_idx for kk in range(nullity)], p)
    T = proj.transpose(); Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, j in enumerate(a_idx):
            v = int(Tr[i, jj])
            if v: poly[S.cols[j][2]] = v
        basis.append(poly)
    Ap = reduce_mod(A, p)
    mons = {m: i for i, m in enumerate(monomials(6, 1))}
    flat = np.zeros((rk + 1, len(mons)), dtype=np.int64)
    for i, poly in enumerate(basis):
        for e, v in poly.items(): flat[i, mons[e]] = v
    for e, v in Ap.items(): flat[rk, mons[e]] = v
    rk2 = rank_np(flat, p)
    return rk, basis, (rk2 == rk)


# ----------------------------------------------------------------------------- exact helpers
def lift_vector(residues, primes):
    """CRT + rational reconstruction of a list of per-prime residue vectors.  Returns list of Fractions or None."""
    acc = [int(a) for a in residues[0]]; m = primes[0]
    for res, p in zip(residues[1:], primes[1:]):
        acc = [core.crt_pair(int(a), m, int(b), p)[0] for a, b in zip(acc, res)]
        m *= p
    out = []
    for a in acc:
        r = core.ratrecon(a, m)
        if r is None:
            return None
        out.append(r)
    return out


def exact_apply_field(u, P):
    out = {}
    for i in range(NV):
        if u[i]:
            out = add(out, mul(u[i], pderiv(P, i)))
    return out
