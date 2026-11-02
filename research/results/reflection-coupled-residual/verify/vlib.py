"""Independent (from-scratch) polynomial helpers for the adversarial verification of lane
reflection-coupled-residual.  No dependency on the hodge engine.  Polys = {exp tuple: coeff}."""
import itertools, re, sys, json, time
from fractions import Fraction
import sympy as sp

VARS = ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2')
NV = 6
REC = '/home/kejace/devel/sextics/data/reflection-coupled-residual/instance.rec'

def load_rec():
    txt = open(REC).read()
    # F spans two lines (continuation indented); A one line
    m = re.search(r'^F:\s*(.*?)(?=^\w+:)', txt, re.S | re.M)
    Fs = ' '.join(m.group(1).split())
    As = re.search(r'^A:\s*(.*)$', txt, re.M).group(1).strip()
    qv = int(re.search(r'^q:\s*(\d+)', txt, re.M).group(1))
    return Fs, As, qv

def parse(s):
    gens = sp.symbols(VARS)
    P = sp.Poly(sp.sympify(s, locals=dict(zip(VARS, gens))), *gens)
    return {tuple(int(x) for x in e): Fraction(int(c.p), int(c.q)) for e, c in P.as_dict().items()}

def modp(P, p):
    out = {}
    for e, c in P.items():
        c = Fraction(c); v = (c.numerator * pow(c.denominator, -1, p)) % p
        if v: out[e] = v
    return out

def deriv(P, i, p=None):
    out = {}
    for e, c in P.items():
        if e[i]:
            e2 = e[:i] + (e[i]-1,) + e[i+1:]
            v = c * e[i]
            if p: v %= p
            if v: out[e2] = v
    return out

def mul(P, Q, p=None):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a+b for a, b in zip(e1, e2))
            out[e] = out.get(e, 0) + c1*c2
    if p: return {e: v % p for e, v in out.items() if v % p}
    return {e: v for e, v in out.items() if v}

def add(P, Q, p=None):
    out = dict(P)
    for e, c in Q.items():
        v = out.get(e, 0) + c
        if p: v %= p
        if v: out[e] = v
        else: out.pop(e, None)
    return out

def scale(P, c, p=None):
    if p: return {e: (v*c) % p for e, v in P.items() if (v*c) % p}
    return {e: v*c for e, v in P.items() if v*c}

def shift(P, m, p=None):
    return {tuple(a+b for a, b in zip(e, m)): v for e, v in P.items()}

def monos(deg):
    """all exponent tuples of total degree deg in 6 vars"""
    out = []
    for c in itertools.combinations(range(deg + NV - 1), NV - 1):
        e = []; prev = -1
        for x in c:
            e.append(x - prev - 1); prev = x
        e.append(deg + NV - 2 - prev)
        out.append(tuple(e))
    return out

def chi(e):
    """Klein character of a monomial: (Y2-parity, (Y0+Y1)-parity)"""
    return (e[5] & 1, (e[3] + e[4]) & 1)

def chi_mul(a, b):
    return ((a[0]+b[0]) & 1, (a[1]+b[1]) & 1)

def char_of(P):
    cs = {chi(e) for e in P}
    assert len(cs) == 1, cs
    return cs.pop()

def monos_chi(deg, c):
    return [e for e in monos(deg) if chi(e) == c]

def fermat_count(deg, c=None):
    """dim of (S/(Z_i^5))_deg with given character = number of monomials with all exponents <= 4"""
    return sum(1 for e in monos(deg) if max(e) <= 4 and (c is None or chi(e) == c))

def nmod_from_rows(rows, cols_index, ncols, p):
    from flint import nmod_mat
    M = nmod_mat(len(rows), ncols, p)
    for r, poly in enumerate(rows):
        for e, v in poly.items():
            M[r, cols_index[e]] = v % p
    return M

def logjson(path, obj):
    obj = dict(obj); obj.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    with open(path, 'a') as f: f.write(json.dumps(obj) + '\n')

VLOG = '/tmp/hodge/results/reflection-coupled-residual/wave2-verify.log.jsonl'
