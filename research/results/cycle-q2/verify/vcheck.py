"""Independent (engine-free) verification library for the cycle-q2 lane.

Polynomials: dict {exponent tuple (6 ints): int or Fraction}.  Variables X0,X1,X2,Y0,Y1,Y2.
Identity (formal/README.md):  E = A D^2 - q D sum_i N_i F_i + F (D sum_i d_i N_i - sum_i N_i d_i D) = F^(q+1) Q.
Integer form used here (F2 = 2F has integer coefficients; Q rescaled, which does not change solvability):
   2^(q+1) A D^2 - 2^q q D sum_i N_i d_i(F2) + 2^q F2 (D div N - N(D)) - F2^(q+1) Q = 0.
Everything is written from scratch (sympy for parsing, python-flint for linear algebra only).
"""
import sys, json, time, itertools, re
from fractions import Fraction
import numpy as np
from flint import nmod_mat, fmpz_mat
import sympy

NV = 6
VARS = sympy.symbols('X0 X1 X2 Y0 Y1 Y2')
LOG = '/tmp/hodge/results/cycle-q2/wave2-verify.log.jsonl'
INSTANCE = '/home/kejace/devel/sextics/data/cycle-q2/instance.rec'


def log(**kw):
    kw.setdefault('t', time.strftime('%Y-%m-%dT%H:%M:%S'))
    kw.setdefault('tool', 'vcheck(own code, flint linalg)')
    with open(LOG, 'a') as f:
        f.write(json.dumps(kw, default=str) + '\n')


# ---------------------------------------------------------------- polynomials
def parse(s):
    e = sympy.sympify(s.replace('^', '**'), locals={str(v): v for v in VARS})
    P = sympy.Poly(sympy.expand(e), *VARS)
    out = {}
    for m, c in zip(P.monoms(), P.coeffs()):
        c = sympy.Rational(c)
        out[tuple(int(x) for x in m)] = Fraction(int(c.p), int(c.q))
    return out


def to_str(P):
    terms = []
    for m, c in sorted(P.items(), reverse=True):
        mon = '*'.join(f'{VARS[i]}**{e}' if e > 1 else str(VARS[i]) for i, e in enumerate(m) if e)
        terms.append(f'({c})*{mon}' if mon else f'({c})')
    return ' + '.join(terms) if terms else '0'


def load_instance(path=INSTANCE):
    txt = open(path).read()
    g = lambda k: re.search(r'^' + k + r':\s*(.*)$', txt, re.M).group(1).strip()
    return parse(g('F')), parse(g('A')), int(g('q'))


def add(P, Q, c=1):
    R = dict(P)
    for m, v in Q.items():
        R[m] = R.get(m, 0) + c * v
        if R[m] == 0:
            del R[m]
    return R


def scale(P, c):
    return {m: c * v for m, v in P.items()} if c else {}


def mul(P, Q):
    R = {}
    for m1, v1 in P.items():
        for m2, v2 in Q.items():
            m = tuple(a + b for a, b in zip(m1, m2))
            R[m] = R.get(m, 0) + v1 * v2
    return {m: v for m, v in R.items() if v != 0}


def power(P, k):
    R = {(0,) * NV: 1}
    for _ in range(k):
        R = mul(R, P)
    return R


def pderiv(P, i):
    R = {}
    for m, v in P.items():
        if m[i]:
            mm = list(m); mm[i] -= 1
            R[tuple(mm)] = v * m[i]
    return R


def degree(P):
    return max(sum(m) for m in P) if P else -1


def monomials(d, n=NV):
    """all exponent tuples of total degree d in n variables"""
    if n == 1:
        return [(d,)]
    out = []
    for e in range(d, -1, -1):
        for rest in monomials(d - e, n - 1):
            out.append((e,) + rest)
    return out


def mono(**kw):
    e = [0] * NV
    for k, v in kw.items():
        e[[str(x) for x in VARS].index(k)] = v
    return {tuple(e): Fraction(1)}


def is_integral(P):
    return all(Fraction(v).denominator == 1 for v in P.values())


# ---------------------------------------------------------------- diagonal group
def diag_group(F):
    """all w in (Z/6)^6 (including scalars) with equal weights on supp F (F must contain all sixth powers,
    which forces every diagonal symmetry of the line of F into mu_6^6 up to scalars)."""
    exps = np.array(list(F.keys()), dtype=np.int64)
    W = np.array(list(itertools.product(range(6), repeat=NV)), dtype=np.int64)
    wts = (exps @ W.T) % 6
    ok = np.all(wts == wts[0:1, :], axis=0)
    return W[ok]


def weight(P, w):
    """weight of eigenvector P under w (mod 6); None if not an eigenvector"""
    ws = {sum(a * b for a, b in zip(m, w)) % 6 for m in P}
    return ws.pop() if len(ws) == 1 else None


def stabilizer(G, P):
    exps = np.array(list(P.keys()), dtype=np.int64)
    wts = (exps @ G.T) % 6
    ok = np.all(wts == wts[0:1, :], axis=0)
    return G[ok]


def sector_monomials(d, H, target_weights):
    """monomials m of degree d with m.w = target_weights[k] (mod 6) for every row w_k of H"""
    allm = np.array(monomials(d), dtype=np.int64)
    if len(H) == 0:
        return [tuple(int(x) for x in m) for m in allm]
    wts = (allm @ H.T) % 6
    ok = np.all(wts == np.array(target_weights, dtype=np.int64)[None, :] % 6, axis=1)
    return [tuple(int(x) for x in m) for m in allm[ok]]


# ---------------------------------------------------------------- the fixed-D system
class FixedD:
    """integer sector-restricted linear system for (N, Q) given D.  Columns: N_i monomials then Q monomials.
    rhs b.  System: M x = b  (x = coefficients), consistent <=> certificate exists (over the field)."""

    def __init__(self, F, A, q, D, H=None, verbose=False):
        t0 = time.time()
        assert all(Fraction(v).denominator in (1, 2) for v in F.values())
        F2 = scale(F, 2)
        assert is_integral(F2) and is_integral(A) and is_integral(D), 'need integral 2F, A, D'
        F2 = {m: int(v) for m, v in F2.items()}
        A = {m: int(v) for m, v in A.items()}
        D = {m: int(v) for m, v in D.items()}
        d = degree(D); self.d = d; self.q = q
        G = diag_group(F) if H is None else H
        H = stabilizer(stabilizer(G, D), A)
        self.H = H; self.order = len(H)
        wA = [weight(A, w) for w in H]; wD = [weight(D, w) for w in H]; wF = [weight(F2, w) for w in H]
        degN = d + 6 * q - 5; degQ = 2 * d - 6
        F2i = [pderiv(F2, i) for i in range(NV)]
        Di = [pderiv(D, i) for i in range(NV)]
        c = 2 ** q
        # column generator for N_i, monomial m:  -c q D F2_i m + c F2 (D d_i m - m D_i)
        Gi = [add(scale(mul(D, F2i[i]), -c * q), mul(F2, Di[i]), -c) for i in range(NV)]   # -(c q D F2_i + c F2 D_i)
        cF2D = scale(mul(F2, D), c)
        Fq1 = scale(power(F2, q + 1), -1)
        rhs = scale(mul(A, mul(D, D)), -(2 ** (q + 1)))
        cols = []; meta = []
        for i in range(NV):
            tw = [(wA[k] + wD[k] - wF[k] + H[k][i]) % 6 for k in range(len(H))]
            for m in sector_monomials(degN, H, tw):
                col = add(mul(Gi[i], {m: 1}), mul(cF2D, pderiv({m: 1}, i)))
                cols.append(col); meta.append(('N', i, m))
        self.nN = len(cols)
        if degQ >= 0:
            tw = [(wA[k] + 2 * wD[k] - (q + 1) * wF[k]) % 6 for k in range(len(H))]
            for m in sector_monomials(degQ, H, tw):
                cols.append(mul(Fq1, {m: 1})); meta.append(('Q', -1, m))
        self.nQ = len(cols) - self.nN
        rows = {}
        for col in cols + [rhs]:
            for m in col:
                rows.setdefault(m, len(rows))
        self.rows = rows; self.cols = cols; self.rhs = rhs; self.meta = meta
        self.nrows = len(rows); self.ncols = len(cols)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[FixedD] d={d} |H|={self.order} rows={self.nrows} cols={self.ncols} (N {self.nN}, Q {self.nQ}) '
                  f'build {self.build_time:.1f}s', file=sys.stderr, flush=True)

    def dense(self, with_rhs=True):
        C = self.ncols + (1 if with_rhs else 0)
        M = np.zeros((self.nrows, C), dtype=object)
        M = [[0] * C for _ in range(self.nrows)]
        for j, col in enumerate(self.cols):
            for m, v in col.items():
                M[self.rows[m]][j] = int(v)
        if with_rhs:
            for m, v in self.rhs.items():
                M[self.rows[m]][self.ncols] = int(v)
        return M

    def ranks_mod(self, p):
        M = self.dense(True)
        flat = [x % p for r in M for x in r]
        Ma = nmod_mat(self.nrows, self.ncols + 1, flat, p)
        ra = Ma.rank()
        flat0 = [x % p for r in M for x in r[:-1]]
        r0 = nmod_mat(self.nrows, self.ncols, flat0, p).rank()
        return r0, ra

    def ranks_exact(self):
        M = self.dense(True)
        Ma = fmpz_mat(M)
        ra = Ma.rank()
        M0 = fmpz_mat([r[:-1] for r in M])
        r0 = M0.rank()
        return r0, ra

    def compressed(self, k_extra=40, seed=1, bound=2 ** 16):
        """random integer row compression (k x (ncols+1)); exact over Z.  Inconsistency of the compressed
        system over Q implies inconsistency of the original over Q (Farkas dual pulls back)."""
        rng = np.random.default_rng(seed)
        k = min(self.nrows, self.ncols + 1 + k_extra)
        if k == self.nrows:
            return self.dense(True)
        R = rng.integers(0, bound, size=(k, self.nrows), dtype=np.int64)
        M = np.zeros((self.nrows, self.ncols + 1), dtype=np.int64)
        for j, col in enumerate(self.cols):
            for m, v in col.items():
                M[self.rows[m], j] = int(v)
        for m, v in self.rhs.items():
            M[self.rows[m], self.ncols] = int(v)
        # exact integer product via object dtype would be slow; split R into 8-bit chunks to stay < 2^63
        out = np.zeros((k, self.ncols + 1), dtype=object)
        acc = None
        for s in range(0, 16, 8):
            Rc = (R >> s) & 0xFF
            part = Rc @ M            # |entries| <= 255 * nnz * max|M| : safe in int64 for our sizes
            part = part.astype(object) * (1 << s)
            acc = part if acc is None else acc + part
        return acc.tolist()


# ---------------------------------------------------------------- Jacobian-ideal helpers
def jacobian_cols(F2, deg, H=None, tw=None):
    """generators m*F2_i of J_deg (optionally restricted to a sector with target weights tw under H)"""
    F2i = [pderiv(F2, i) for i in range(NV)]
    out = []; meta = []
    for i in range(NV):
        if H is None:
            ms = monomials(deg - 5)
        else:
            wi = [(tw[k] - weight(F2i[i], H[k])) % 6 for k in range(len(H))]
            ms = sector_monomials(deg - 5, H, wi)
        for m in ms:
            out.append(mul(F2i[i], {m: 1})); meta.append((i, m))
    return out, meta


def rank_of_cols(cols, p=None):
    rows = {}
    for c in cols:
        for m in c:
            rows.setdefault(m, len(rows))
    if not cols or not rows:
        return 0, len(rows)
    if p is None:
        M = [[0] * len(cols) for _ in rows]
        for j, c in enumerate(cols):
            for m, v in c.items():
                M[rows[m]][j] = int(v)
        return fmpz_mat(M).rank(), len(rows)
    flat = [0] * (len(rows) * len(cols))
    nc = len(cols)
    for j, c in enumerate(cols):
        for m, v in c.items():
            flat[rows[m] * nc + j] = int(v) % p
    return nmod_mat(len(rows), nc, flat, p).rank(), len(rows)


def is_prime(n):
    return sympy.isprime(n)


if __name__ == '__main__':
    F, A, q = load_instance()
    print('F =', to_str(F)); print('A =', to_str(A)); print('q =', q)
    G = diag_group(F)
    print('|diag group incl. scalars| =', len(G))
