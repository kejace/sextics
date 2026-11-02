"""Shared exact helpers for the polarized-k3-remainder exclusion packet (wave 2, lens excl).
Frames: XY = instance coordinates (X0,X1,X2,Y0,Y1,Y2); UV = (U,V) with U = X+Y, V = X-Y, slots 0..2 = U, 3..5 = V.
Fhat(U,V) = F((U+V)/2,(U-V)/2) = (f(U)+f(V))/2,  Ahat(U,V) = A((U+V)/2,(U-V)/2) = 2^-6 prod_i (U_i^2 - V_i^2).
All polynomials are core-style dicts {6-exponent tuple: Fraction} unless a modulus is given."""
import sys, json, time, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction as Fr
from hodge import core
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, power, shift
from flint import fmpq_mat, fmpz_mat, nmod_mat, fmpq, fmpz

IID = 'polarized-k3-remainder'

def load_xy():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q

def var(i):
    e = [0] * NV; e[i] = 1
    return {tuple(e): Fr(1)}

def substitute_linear(P, images, p=None):
    """P(z) with z_i := images[i] (polys)."""
    out = {}
    cache = {}
    for e, c in P.items():
        term = {(0,) * NV: Fr(c) if p is None else c % p}
        for i in range(NV):
            if e[i]:
                key = (i, e[i])
                if key not in cache:
                    cache[key] = power(images[i], e[i], p)
                term = mul(term, cache[key], p)
        out = add(out, term, p)
    return out

def xy_to_uv(P, p=None):
    """X = (U+V)/2, Y = (U-V)/2."""
    half = Fr(1, 2) if p is None else pow(2, -1, p)
    U = [var(i) for i in range(3)]; V = [var(3 + i) for i in range(3)]
    if p is not None:
        U = [core.reduce_mod(u, p) for u in U]; V = [core.reduce_mod(v, p) for v in V]
    imgs = [scale(add(U[i], V[i], p), half, p) for i in range(3)] + [scale(sub(U[i], V[i], p), half, p) for i in range(3)]
    return substitute_linear(P, imgs, p)

def uv_to_xy(P, p=None):
    """U = X + Y, V = X - Y."""
    X = [var(i) for i in range(3)]; Y = [var(3 + i) for i in range(3)]
    if p is not None:
        X = [core.reduce_mod(u, p) for u in X]; Y = [core.reduce_mod(v, p) for v in Y]
    imgs = [add(X[i], Y[i], p) for i in range(3)] + [sub(X[i], Y[i], p) for i in range(3)]
    return substitute_linear(P, imgs, p)

def load_uv():
    F, A, q = load_xy()
    return xy_to_uv(F), xy_to_uv(A), q

def ternary_f():
    """f(z0,z1,z2) with F = (f(X+Y)+f(X-Y))/2, as a dict over 3-tuples."""
    Fh, Ah, q = load_uv()
    f = {}
    for e, c in Fh.items():
        if e[3] == e[4] == e[5] == 0:
            f[(e[0], e[1], e[2])] = 2 * c
    return f

def vdeg(e):
    return e[3] + e[4] + e[5]

def split_vdeg_mod(P, n=6):
    """Split a UV-frame poly into components by V-degree mod n."""
    comps = {}
    for e, c in P.items():
        comps.setdefault(vdeg(e) % n, {})[e] = c
    return comps

def swap(P):
    return {(e[3], e[4], e[5], e[0], e[1], e[2]): c for e, c in P.items()}

# ---------------- exact linear algebra helpers ----------------

def polys_to_fmpq_mat(polys, rows=None):
    """Columns = polys.  Returns (fmpq_mat (nrows x ncols), row index dict)."""
    if rows is None:
        rows = {}
        for P in polys:
            for e in P:
                rows.setdefault(e, len(rows))
    R, C = len(rows), len(polys)
    M = fmpq_mat(R, C)
    for j, P in enumerate(polys):
        for e, c in P.items():
            c = Fr(c)
            M[rows[e], j] = fmpq(c.numerator, c.denominator)
    return M, rows

def polys_to_int_rows(polys, rows=None):
    """Columns = polys; returns (list of row-lists of python ints after scaling each column to integers, rows)."""
    if rows is None:
        rows = {}
        for P in polys:
            for e in P:
                rows.setdefault(e, len(rows))
    cols = []
    for P in polys:
        den = 1
        for c in P.values():
            c = Fr(c)
            from math import gcd
            den = den * c.denominator // gcd(den, c.denominator)
        cols.append({rows[e]: int(Fr(c) * den) for e, c in P.items()})
    return cols, rows

def int_cols_to_fmpz_mat(cols, nrows):
    C = len(cols)
    M = fmpz_mat(nrows, C)
    for j, col in enumerate(cols):
        for r, v in col.items():
            M[r, j] = v
    return M

def exact_rank(polys):
    """Exact rank over Q of the span of the given polys."""
    cols, rows = polys_to_int_rows(polys)
    M = int_cols_to_fmpz_mat(cols, len(rows))
    return M.rank()

# ---------------- exact ternary Jacobian ring and tensor ring over Q ----------------

def tern_monomials(k):
    return [(a, b, k - a - b) for a in range(k + 1) for b in range(k + 1 - a)]

def tmul(P, Q):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2])
            out[e] = out.get(e, 0) + c1 * c2
    return {e: v for e, v in out.items() if v}

def tderiv(P, i):
    out = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            out[tuple(e2)] = c * e[i]
    return out

class TernaryRingQ:
    """R_f = Q[z]/(f_0,f_1,f_2) degree by degree, exact (fmpq_mat rref)."""
    def __init__(self, f, maxdeg=18):
        self.f = f
        self.partials = [tderiv(f, i) for i in range(3)]
        self.data = {}
        for k in range(maxdeg + 1):
            self._build(k)
    def _build(self, k):
        mons = tern_monomials(k); idx = {m: i for i, m in enumerate(mons)}
        gens = []
        if k >= 5:
            for m in tern_monomials(k - 5):
                for Pi in self.partials:
                    gens.append(tmul({m: Fr(1)}, Pi))
        if gens:
            M = fmpq_mat(len(gens), len(mons))
            for r, g in enumerate(gens):
                for e, v in g.items():
                    v = Fr(v); M[r, idx[e]] = fmpq(v.numerator, v.denominator)
            Mr, rank = M.rref()
            pivots, rows = [], []
            for i in range(rank):
                j = 0
                while Mr[i, j] == 0: j += 1
                pivots.append(j)
                rows.append({jj: Fr(int(Mr[i, jj].p), int(Mr[i, jj].q)) for jj in range(len(mons)) if Mr[i, jj] != 0})
        else:
            rank, pivots, rows = 0, [], []
        std = [m for j, m in enumerate(mons) if j not in set(pivots)]
        self.data[k] = dict(mons=mons, idx=idx, pivots=pivots, rows=rows, std=std, dim=len(std),
                            stdidx={m: i for i, m in enumerate(std)})
    def dim(self, k):
        return self.data[k]['dim'] if k in self.data else 0
    def nf(self, P, k):
        d = self.data[k]
        vec = {}
        for e, v in P.items():
            vec[d['idx'][e]] = vec.get(d['idx'][e], 0) + Fr(v)
        for row, pc in zip(d['rows'], d['pivots']):
            c = vec.get(pc, 0)
            if c:
                for j, rv in row.items():
                    vec[j] = vec.get(j, 0) - c * rv
                vec.pop(pc, None)
        return {d['stdidx'][d['mons'][j]]: v for j, v in vec.items() if v}

class TensorRingQ:
    """R_F = R_f (x) R_f in the (U,V) frame, exact."""
    def __init__(self, R, maxdeg=24):
        self.R = R
        self.bases = {}
        for d in range(maxdeg + 1):
            basis = []
            for a in range(d + 1):
                b = d - a
                if R.dim(a) == 0 or R.dim(b) == 0: continue
                for i in range(R.dim(a)):
                    for j in range(R.dim(b)):
                        basis.append((a, b, i, j))
            self.bases[d] = basis
        self.index = {d: {t: i for i, t in enumerate(self.bases[d])} for d in self.bases}
        self._nfcache = {}
    def dim(self, d):
        return len(self.bases.get(d, []))
    def nfm(self, mon, k):
        key = (mon, k)
        if key not in self._nfcache:
            self._nfcache[key] = self.R.nf({mon: Fr(1)}, k)
        return self._nfcache[key]
    def vector_of(self, Q6, d):
        """Class of a UV-frame poly of degree d as a sparse dict {basis index: Fraction}."""
        vec = {}
        for e, v in Q6.items():
            a = e[0] + e[1] + e[2]; b = d - a
            if self.R.dim(a) == 0 or self.R.dim(b) == 0: continue
            nu = self.nfm((e[0], e[1], e[2]), a); nv = self.nfm((e[3], e[4], e[5]), b)
            for i, x in nu.items():
                for j, y in nv.items():
                    k = self.index[d][(a, b, i, j)]
                    vec[k] = vec.get(k, 0) + Fr(v) * x * y
        return {k: v for k, v in vec.items() if v}
    def mult_columns(self, Q6, dq, d):
        """Columns (sparse dicts over basis of degree d+dq) of multiplication by Q6 on basis of degree d."""
        R = self.R
        terms = [(Fr(v), (e[0], e[1], e[2]), e[0] + e[1] + e[2], (e[3], e[4], e[5]), e[3] + e[4] + e[5]) for e, v in Q6.items()]
        cols = []
        for (a, b, i, j) in self.bases[d]:
            mu = R.data[a]['std'][i]; mv = R.data[b]['std'][j]
            col = {}
            for v, eu, du, ev, dv in terms:
                a2, b2 = a + du, b + dv
                if R.dim(a2) == 0 or R.dim(b2) == 0: continue
                nu = self.nfm((mu[0] + eu[0], mu[1] + eu[1], mu[2] + eu[2]), a2)
                nv = self.nfm((mv[0] + ev[0], mv[1] + ev[1], mv[2] + ev[2]), b2)
                for ii, x in nu.items():
                    for jj, y in nv.items():
                        r = self.index[d + dq][(a2, b2, ii, jj)]
                        col[r] = col.get(r, 0) + v * x * y
            cols.append({r: c for r, c in col.items() if c})
        return cols

def sparse_cols_to_fmpz(cols, nrows):
    """Scale each sparse rational column to integers and return fmpz_mat (nrows x ncols)."""
    from math import gcd
    M = fmpz_mat(nrows, len(cols))
    for j, col in enumerate(cols):
        den = 1
        for v in col.values():
            v = Fr(v); den = den * v.denominator // gcd(den, v.denominator)
        for r, v in col.items():
            M[r, j] = int(Fr(v) * den)
    return M
