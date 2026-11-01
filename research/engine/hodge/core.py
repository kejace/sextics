"""Exact polynomial and certificate machinery for the hodge.computer challenge.

Polynomials are dicts {exponent tuple (6 ints): coefficient}.  Coefficients are
python Fractions/ints (exact) or ints mod p when a modulus is given.

Variables: X0,X1,X2,Y0,Y1,Y2 (indices 0..5).  Y-parity of a monomial = (e3+e4+e5) % 2.

Certificate identity (cleared form):
  E = A*D^2 - q*D*sum_i N_i*F_i + F*(D*sum_i d_i N_i - sum_i N_i*D_i) = F^(q+1) * Q
with deg D = d, deg N_i = d+6q-5, deg Q = 2d-6.
"""
from __future__ import annotations
import ast, itertools, math, random, sys, time
from fractions import Fraction
from functools import lru_cache

NV = 6
VARS = ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2')

# ----------------------------------------------------------------------------
# parsing / printing
# ----------------------------------------------------------------------------

def parse(src: str, mod: int | None = None) -> dict:
    """Parse an exact polynomial expression in the instance.rec / witness.json syntax."""
    src = src.replace('\n', ' ')
    tree = ast.parse(src, mode='eval').body
    idx = {v: i for i, v in enumerate(VARS)}

    def t(n):
        if isinstance(n, ast.Name):
            e = [0] * NV; e[idx[n.id]] = 1
            return {tuple(e): Fraction(1)}
        if isinstance(n, ast.Constant):
            if isinstance(n.value, bool) or not isinstance(n.value, int):
                raise ValueError('only integer literals')
            return {(0,) * NV: Fraction(n.value)} if n.value else {}
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return scale(t(n.operand), Fraction(-1))
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.UAdd):
            return t(n.operand)
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, ast.Pow):
                assert isinstance(n.right, ast.Constant) and isinstance(n.right.value, int)
                return power(t(n.left), n.right.value)
            if isinstance(n.op, ast.Div):
                r = t(n.right)
                assert len(r) == 1 and (0,) * NV in r, 'only constant denominators'
                return scale(t(n.left), 1 / r[(0,) * NV])
            l, r = t(n.left), t(n.right)
            if isinstance(n.op, ast.Add): return add(l, r)
            if isinstance(n.op, ast.Sub): return sub(l, r)
            if isinstance(n.op, ast.Mult): return mul(l, r)
        raise ValueError(f'unsupported node {ast.dump(n)}')
    P = t(tree)
    return reduce_mod(P, mod) if mod else P


def to_str(P: dict) -> str:
    """Render in the witness.json syntax (integer coefficients over integer denominators)."""
    if not P:
        return '0'
    terms = []
    for e, c in sorted(P.items(), reverse=True):
        c = Fraction(c)
        mon = '*'.join(f'{VARS[i]}**{k}' if k > 1 else VARS[i] for i, k in enumerate(e) if k)
        num, den = c.numerator, c.denominator
        if mon:
            coef = '' if abs(num) == 1 else f'{abs(num)}*'
            s = f'{coef}{mon}'
        else:
            s = f'{abs(num)}'
        if den != 1:
            s = f'{s}/{den}'
        terms.append(('-' if num < 0 else '+', s))
    out = ''
    for i, (sg, s) in enumerate(terms):
        if i == 0:
            out = ('-' if sg == '-' else '') + s
        else:
            out += f' {sg} {s}'
    return out

# ----------------------------------------------------------------------------
# arithmetic
# ----------------------------------------------------------------------------

def reduce_mod(P, p):
    out = {}
    for e, c in P.items():
        c = Fraction(c)
        v = (c.numerator * pow(c.denominator, -1, p)) % p
        if v:
            out[e] = v
    return out


def add(P, Q, p=None):
    out = dict(P)
    for e, c in Q.items():
        v = out.get(e, 0) + c
        if p: v %= p
        if v: out[e] = v
        else: out.pop(e, None)
    return out


def sub(P, Q, p=None):
    out = dict(P)
    for e, c in Q.items():
        v = out.get(e, 0) - c
        if p: v %= p
        if v: out[e] = v
        else: out.pop(e, None)
    return out


def scale(P, c, p=None):
    if not c:
        return {}
    if p:
        return {e: (v * c) % p for e, v in P.items() if (v * c) % p}
    return {e: v * c for e, v in P.items()}


def mul(P, Q, p=None):
    out = {}
    if len(P) > len(Q):
        P, Q = Q, P
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = (e1[0]+e2[0], e1[1]+e2[1], e1[2]+e2[2], e1[3]+e2[3], e1[4]+e2[4], e1[5]+e2[5])
            out[e] = out.get(e, 0) + c1 * c2
    if p:
        return {e: v % p for e, v in out.items() if v % p}
    return {e: v for e, v in out.items() if v}


def power(P, k, p=None):
    R = {(0,) * NV: 1}
    B = P
    while k:
        if k & 1: R = mul(R, B, p)
        k >>= 1
        if k: B = mul(B, B, p)
    return R


def shift(P, m, c=1, p=None):
    """c * monomial(m) * P"""
    out = {}
    for e, v in P.items():
        e2 = (e[0]+m[0], e[1]+m[1], e[2]+m[2], e[3]+m[3], e[4]+m[4], e[5]+m[5])
        out[e2] = (v * c) % p if p else v * c
    return out


def pderiv(P, i, p=None):
    out = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            v = c * e[i]
            if p: v %= p
            if v: out[tuple(e2)] = v
    return out


def degree(P):
    return max((sum(e) for e in P), default=-1)


def is_homogeneous(P):
    degs = {sum(e) for e in P}
    return len(degs) <= 1


def yparity(e):
    return (e[3] + e[4] + e[5]) & 1


def parity_of(P):
    """0 = Y-even, 1 = Y-odd, None = mixed."""
    ps = {yparity(e) for e in P}
    return ps.pop() if len(ps) == 1 else None


def evaluate(P, pt, p=None):
    s = 0
    for e, c in P.items():
        t = c
        for i in range(NV):
            if e[i]:
                t *= pt[i] ** e[i]
        s += t
    return s % p if p else s


def monomials(deg, par=None):
    """All exponent tuples of total degree deg (optionally with given Y-parity)."""
    if deg < 0:
        return []
    out = []
    for c in itertools.combinations(range(deg + NV - 1), NV - 1):
        e = []
        prev = -1
        for x in c:
            e.append(x - prev - 1); prev = x
        e.append(deg + NV - 2 - prev)
        e = tuple(e)
        if par is None or yparity(e) == par:
            out.append(e)
    return out


def substitute_scalar(P, i, val, p=None):
    """Substitute variable i := val (a scalar)."""
    out = {}
    for e, c in P.items():
        v = c * (val ** e[i]) if e[i] else c
        if p: v %= p
        e2 = list(e); e2[i] = 0; e2 = tuple(e2)
        v2 = out.get(e2, 0) + v
        if p: v2 %= p
        if v2: out[e2] = v2
        else: out.pop(e2, None)
    return out

# ----------------------------------------------------------------------------
# certificate defect
# ----------------------------------------------------------------------------

def defect(F, A, D, q, N, p=None):
    """E = A D^2 - q D sum N_i F_i + F (D sum d_i N_i - sum N_i D_i)."""
    D2 = mul(D, D, p)
    E = mul(A, D2, p)
    dirF = {}
    divN = {}
    dirD = {}
    for i in range(NV):
        if not N[i]:
            continue
        dirF = add(dirF, mul(pderiv(F, i, p), N[i], p), p)
        divN = add(divN, pderiv(N[i], i, p), p)
        dirD = add(dirD, mul(pderiv(D, i, p), N[i], p), p)
    E = sub(E, scale(mul(D, dirF, p), q, p), p)
    E = add(E, mul(F, sub(mul(D, divN, p), dirD, p), p), p)
    return E


def verify_witness(F, A, q, D, N, Q, p=None):
    """Return (ok, info).  Exact check of homogeneity, degrees and identity (not coprimality)."""
    d = degree(D)
    info = {'d': d}
    if not D or not is_homogeneous(D):
        return False, 'D zero or not homogeneous'
    for i in range(NV):
        if N[i] and (not is_homogeneous(N[i]) or degree(N[i]) != d + 6*q - 5):
            return False, f'N{i} degree {degree(N[i])} != {d+6*q-5}'
    E = defect(F, A, D, q, N, p)
    R = mul(power(F, q + 1, p), Q, p)
    diff = sub(E, R, p)
    info['defect_terms'] = len(E)
    info['residual_terms'] = len(diff)
    return (not diff), info


# ----------------------------------------------------------------------------
# instance loading
# ----------------------------------------------------------------------------

ROOT = '/home/kejace/devel/sextics'


def load_instance(iid, mod=None):
    import yaml
    with open(f'{ROOT}/data/{iid}/instance.rec') as f:
        rec = yaml.safe_load(f)
    F = parse(rec['F'], mod)
    A = parse(rec['A'], mod)
    return F, A, int(rec['q']), rec


def load_witness(iid, mod=None):
    import json
    with open(f'{ROOT}/data/{iid}/witness.json') as f:
        w = json.load(f)
    D = parse(w['D'], mod)
    N = [parse(s, mod) for s in w['N']]
    Q = parse(w['Q'], mod)
    return w['d'], D, N, Q


# ----------------------------------------------------------------------------
# number theory helpers
# ----------------------------------------------------------------------------

PRIMES = [2147483647, 2147483629, 2147483587, 2147483579, 2147483563, 2147483549,
          2147483543, 2147483497, 2147483489, 2147483477, 2147483423, 2147483399]


def ratrecon(a, m):
    """Rational reconstruction of a mod m: returns Fraction or None."""
    a %= m
    if a == 0:
        return Fraction(0)
    bound = math.isqrt(m // 2)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    if math.gcd(r1, s1) != 1:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def crt_pair(a1, m1, a2, m2):
    inv = pow(m1, -1, m2)
    t = ((a2 - a1) * inv) % m2
    return (a1 + m1 * t) % (m1 * m2), m1 * m2
