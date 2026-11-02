"""Independent referee library for reflection-coupled-residual (no dependency on the lane engine).
Polynomials are dicts {exponent-tuple: Fraction}; modular linear algebra via python-flint nmod_mat.
Variables X0,X1,X2,Y0,Y1,Y2 (indices 0..5). Klein character of a monomial e: (e5 mod 2, (e3+e4) mod 2)."""
import itertools, json, time, sys, re
from fractions import Fraction
from flint import nmod_mat, fmpz

NAMES = ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2'); NV = 6
CHI_A = (0, 1)
CHARS = [(0, 0), (0, 1), (1, 0), (1, 1)]
INSTANCE = '/home/kejace/devel/sextics/data/reflection-coupled-residual/instance.rec'

class Poly:
    __slots__ = ('t',)
    def __init__(self, t=None): self.t = dict(t) if t else {}
    @staticmethod
    def var(i):
        e = [0] * NV; e[i] = 1; return Poly({tuple(e): Fraction(1)})
    @staticmethod
    def const(c):
        c = Fraction(c); return Poly({(0,) * NV: c}) if c else Poly()
    def _p(self, o): return o if isinstance(o, Poly) else Poly.const(o)
    def __add__(self, o):
        o = self._p(o); r = dict(self.t)
        for e, c in o.t.items():
            v = r.get(e, 0) + c
            if v: r[e] = v
            else: r.pop(e, None)
        return Poly(r)
    __radd__ = __add__
    def __neg__(self): return Poly({e: -c for e, c in self.t.items()})
    def __sub__(self, o): return self + (-self._p(o))
    def __rsub__(self, o): return (-self) + o
    def __mul__(self, o):
        if not isinstance(o, Poly):
            o = Fraction(o)
            return Poly({e: c * o for e, c in self.t.items()}) if o else Poly()
        r = {}
        for e1, c1 in self.t.items():
            for e2, c2 in o.t.items():
                e = tuple(a + b for a, b in zip(e1, e2))
                v = r.get(e, 0) + c1 * c2
                if v: r[e] = v
                else: r.pop(e, None)
        return Poly(r)
    __rmul__ = __mul__
    def __truediv__(self, o): return self * Fraction(1, 1) * (Fraction(1) / Fraction(o))
    def __pow__(self, n):
        r = Poly.const(1)
        for _ in range(n): r = r * self
        return r
    def deriv(self, i):
        r = {}
        for e, c in self.t.items():
            if e[i]:
                f = list(e); f[i] -= 1; r[tuple(f)] = c * e[i]
        return Poly(r)
    def degree(self): return max((sum(e) for e in self.t), default=-1)
    def is_zero(self): return not self.t
    def chars(self): return sorted({chi(e) for e in self.t})
    def modp(self, p):
        return {e: (c.numerator * pow(c.denominator, -1, p)) % p for e, c in self.t.items()}

def parse(expr):
    env = {n: Poly.var(i) for i, n in enumerate(NAMES)}
    return eval(expr.replace('\n', ' '), {'__builtins__': {}}, env)

def load_instance(path=INSTANCE):
    txt = open(path).read()
    m = re.search(r'^F:\s*(.*?)(?=^\w+:)', txt, re.S | re.M)
    F = parse(m.group(1)); A = parse(re.search(r'^A:\s*(.*)$', txt, re.M).group(1))
    q = int(re.search(r'^q:\s*(\d+)', txt, re.M).group(1))
    return F, A, q

def chi(e): return (e[5] % 2, (e[3] + e[4]) % 2)
def cadd(a, b): return ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2)
CHI_Z = [chi(tuple(1 if j == i else 0 for j in range(NV))) for i in range(NV)]

_moncache = {}
def monomials(n, c=None):
    key = (n, c)
    if key not in _moncache:
        out = []
        if n >= 0:
            for cuts in itertools.combinations(range(n + NV - 1), NV - 1):
                e = []; prev = -1
                for x in cuts:
                    e.append(x - prev - 1); prev = x
                e.append(n + NV - 2 - prev)
                e = tuple(e)
                if c is None or chi(e) == c: out.append(e)
        _moncache[key] = out
    return _moncache[key]

def fermat_count(n, c=None):
    return sum(1 for e in monomials(n, c) if max(e) <= 4)

def shift(Pm, m, p):
    """monomial m times mod-p dict Pm."""
    return {tuple(a + b for a, b in zip(e, m)): v for e, v in Pm.items()}

def mulp(P, Q, p):
    r = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2)); r[e] = (r.get(e, 0) + c1 * c2) % p
    return {e: v for e, v in r.items() if v}

def addp(P, Q, p, s=1):
    r = dict(P)
    for e, v in Q.items():
        r[e] = (r.get(e, 0) + s * v) % p
    return {e: v for e, v in r.items() if v}

class Setup:
    def __init__(self, p):
        self.p = p
        assert fmpz(p).is_prime()
        self.F, self.A, self.q = load_instance()
        self.Fp = self.F.modp(p); self.Ap = self.A.modp(p)
        self.Fi = [self.F.deriv(i).modp(p) for i in range(NV)]
        self.Fij = {(i, j): self.F.deriv(i).deriv(j).modp(p) for i in range(NV) for j in range(NV)}

    def jac_rows(self, n, c):
        """generators m*F_i of J_n in character c (as mod-p dicts)."""
        rows = []
        if n >= 5:
            for i in range(NV):
                for m in monomials(n - 5, cadd(c, CHI_Z[i])):
                    rows.append(shift(self.Fi[i], m, self.p))
        return rows

def to_mat(rows, idx, C, p):
    flat = [0] * (len(rows) * C)
    for r, P in enumerate(rows):
        base = r * C
        for e, v in P.items(): flat[base + idx[e]] = v % p
    return nmod_mat(len(rows), C, flat, p)

class Reducer:
    """rref of J_n^c mod p: ranks, good-reduction flag, normal forms, standard-monomial coordinates."""
    def __init__(self, S, n, c):
        self.S, self.n, self.c, self.p = S, n, c, S.p
        self.mons = monomials(n, c); self.idx = {e: i for i, e in enumerate(self.mons)}
        C = len(self.mons); self.C = C
        rows = S.jac_rows(n, c)
        self.nrows = len(rows)
        if rows:
            M = to_mat(rows, self.idx, C, self.p)
            Jr, rank = M.rref()
            self.rank = rank
            piv = []
            for i in range(rank):
                j = 0
                while int(Jr[i, j]) == 0: j += 1
                piv.append(j)
            self.pivots = piv
            self.Jr = nmod_mat(rank, C, [int(Jr[i, j]) for i in range(rank) for j in range(C)], self.p) if rank else None
        else:
            self.rank = 0; self.pivots = []; self.Jr = None
        pivset = set(self.pivots)
        self.standard = [j for j in range(C) if j not in pivset]
        self.dimS = C; self.dimRJ = C - self.rank
        self.fermat = fermat_count(n, c)
        self.good = (self.dimRJ == self.fermat)

    def nf_mat(self, polys):
        """normal forms of polys (mod-p dicts, degree n, character c) as an nmod_mat (len x C)."""
        p, C = self.p, self.C
        V = to_mat(polys, self.idx, C, p)
        if self.rank == 0: return V
        k = self.rank
        pivmons = [self.mons[j] for j in self.pivots]
        flatp = [0] * (len(polys) * k)
        for r, P in enumerate(polys):
            base = r * k
            for t, m in enumerate(pivmons):
                v = P.get(m)
                if v: flatp[base + t] = v % p
        Vp = nmod_mat(len(polys), k, flatp, p)
        return V - Vp * self.Jr

    def coords(self, polys):
        """standard-monomial coordinates (list of lists) of the classes of polys in (S/J)_n^c."""
        NF = self.nf_mat(polys)
        return [[int(NF[r, j]) for j in self.standard] for r in range(len(polys))]

    def rank_mod_J(self, polys):
        """dimension of the span of the classes of polys in (S/J)_n^c."""
        if not polys: return 0
        NF = self.nf_mat(polys)
        return NF.rank()

def receipt(path, rec):
    rec = dict(rec); rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S'); rec['referee'] = True
    with open(path, 'a') as f: f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec)[:300], flush=True)
