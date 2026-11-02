"""Referee library (written from scratch, independent of ~/devel/sextics-work/hodge).
Exact parsing of data/eight-plane-configuration/instance.rec with sympy; polynomials are dicts {exponent tuple: coeff}
(Fraction for exact work, int for mod-p work). Variables X0,X1,X2,Y0,Y1,Y2 in this order. Y-parity = parity of total Y-degree.
Matrices: python-flint nmod_mat / fmpq_mat, filled sparsely (the local box has little free memory)."""
import sys, json, time, re
from fractions import Fraction
from math import comb
from itertools import combinations_with_replacement
from flint import nmod_mat, fmpq_mat, fmpq

IID = 'eight-plane-configuration'
REC = f'/home/kejace/devel/sextics/data/{IID}/instance.rec'
REF = f'/tmp/hodge/results/{IID}/wave2-excl-referee'
LOG = f'{REF}/referee.log.jsonl'
NV = 6
VARS = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']

def log(entry):
    entry.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    entry.setdefault('agent', 'referee')
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry) + '\n')
    print(json.dumps(entry), flush=True)

def _read_rec():
    txt = open(REC).read().splitlines()
    vals = {}; key = None
    for line in txt:
        m = re.match(r'^([A-Za-z]+):\s*(.*)$', line)
        if m and not line.startswith(' '):
            key = m.group(1); vals[key] = m.group(2)
        elif key and line.startswith(' '):
            vals[key] += ' ' + line.strip()
    return vals

def parse_poly(expr):
    import sympy
    syms = sympy.symbols(VARS)
    e = sympy.sympify(expr, locals=dict(zip(VARS, syms)))
    P = sympy.Poly(sympy.expand(e), *syms)
    out = {}
    for mon, c in P.as_dict().items():
        c = sympy.Rational(c)
        out[tuple(int(x) for x in mon)] = Fraction(int(c.p), int(c.q))
    return out

_cache = {}
def load():
    if 'F' not in _cache:
        v = _read_rec()
        _cache['F'] = parse_poly(v['F']); _cache['A'] = parse_poly(v['A']); _cache['q'] = int(v['q'])
    return _cache['F'], _cache['A'], _cache['q']

def ydeg(e): return e[3] + e[4] + e[5]
def deg(P): return sum(next(iter(P))) if P else -1
def parity(P):
    ps = {ydeg(e) % 2 for e in P}
    assert len(ps) <= 1, 'mixed parity'
    return ps.pop() if ps else None

def monomials(d, par=None):
    """all exponent tuples of total degree d (optionally with Y-parity par), in a fixed order."""
    out = []
    for c in combinations_with_replacement(range(NV), d):
        e = [0] * NV
        for i in c: e[i] += 1
        e = tuple(e)
        if par is None or ydeg(e) % 2 == par: out.append(e)
    return out

def S_dim(k, par=None):
    tot = 0
    for j in range(k + 1):
        if par is not None and j % 2 != par: continue
        tot += comb(j + 2, 2) * comb(k - j + 2, 2)
    return tot

def R_dim(k, par):
    """Koszul: R = S/(F_X0,F_X1,F_X2 (deg 5, even), F_Y0,F_Y1,F_Y2 (deg 5, odd)), regular sequence."""
    v = 0
    for a in range(4):
        for b in range(4):
            kk = k - 5 * (a + b)
            if kk < 0: continue
            v += (-1) ** (a + b) * comb(3, a) * comb(3, b) * S_dim(kk, (par - b) % 2)
    return v

def J_dim(k, par): return S_dim(k, par) - R_dim(k, par)

def red(P, p):
    """reduce an exact polynomial mod p (Fractions -> ints)."""
    out = {}
    for e, c in P.items():
        if isinstance(c, Fraction):
            v = (c.numerator % p) * pow(c.denominator % p, -1, p) % p
        else:
            v = c % p
        if v: out[e] = v
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

def mul(P, Q, p=None):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            v = out.get(e, 0) + c1 * c2
            if p: v %= p
            if v: out[e] = v
            elif e in out: del out[e]
    return out

def shift(P, m, c=1, p=None):
    out = {}
    for e, v in P.items():
        vv = v * c
        if p: vv %= p
        if vv: out[tuple(a + b for a, b in zip(e, m))] = vv
    return out

def add(P, Q, p=None, s=1):
    out = dict(P)
    for e, c in Q.items():
        v = out.get(e, 0) + s * c
        if p: v %= p
        if v: out[e] = v
        elif e in out: del out[e]
    return out

def nmod_from_rows(rows, cols_idx, ncols, p):
    """rows: list of dict polys; cols_idx: {monomial: column}. Sparse fill (memory-light)."""
    M = nmod_mat(len(rows), ncols, p)
    for r, P in enumerate(rows):
        for e, v in P.items():
            M[r, cols_idx[e]] = v % p
    return M

def fmpq_from_rows(rows, cols_idx, ncols):
    M = fmpq_mat(len(rows), ncols)
    for r, P in enumerate(rows):
        for e, v in P.items():
            v = Fraction(v)
            M[r, cols_idx[e]] = fmpq(v.numerator, v.denominator)
    return M

def J_generators(Fi, k, par, p=None):
    """generators of J_k^par: m * F_i with deg m = k-5 and parity(m) = par - parity(F_i)."""
    gens = []
    for i in range(NV):
        pm = (par - (0 if i < 3 else 1)) % 2
        for m in monomials(k - 5, pm):
            gens.append(shift(Fi[i], m, 1, p))
    return gens
