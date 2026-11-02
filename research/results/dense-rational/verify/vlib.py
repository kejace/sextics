"""Independent helpers for the dense-rational VERIFY lens.  Does NOT import the lane's dr.py / hodge.fixed_d.
Polynomials: dict {exp tuple (6,): int or Fraction}.  Variables X0,X1,X2,Y0,Y1,Y2 = indices 0..5."""
import re, json, time, itertools
from fractions import Fraction
from math import comb
from flint import nmod_mat, fmpz_mat, fmpq_mat, fmpz, fmpq

REC = '/home/kejace/devel/sextics/data/dense-rational/instance.rec'
VARS = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
NV = 6
LOGF = '/tmp/hodge/results/dense-rational/wave2-verify.log.jsonl'

def log(claim, verdict, evidence, notes='', **kw):
    rec = {'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'lane': 'dense-rational', 'lens': 'verify',
           'claim': claim, 'verdict': verdict, 'evidence': evidence, 'notes': notes}
    rec.update(kw)
    with open(LOGF, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    print('[LOG]', verdict, '|', claim, '|', evidence, '|', notes, flush=True)

def parse_rec():
    """From-scratch parser of instance.rec: returns (F, A, q) with Fraction coefficients."""
    txt = open(REC).read()
    def field(name):
        m = re.search(r'^' + name + r':(.*?)(?=^\w+:|\Z)', txt, re.S | re.M)
        return m.group(1)
    Fs = field('F').replace('\n', ' ').replace(' ', '')
    As = field('A').replace('\n', ' ').replace(' ', '')
    q = int(field('q').strip())
    return parse_poly(Fs), parse_poly(As), q

def parse_poly(s):
    s = s.replace('-', '+-')
    terms = [t for t in s.split('+') if t]
    P = {}
    for t in terms:
        sign = 1
        if t.startswith('-'):
            sign = -1; t = t[1:]
        num = Fraction(1); den = Fraction(1)
        e = [0] * NV
        # split on '*' but handle '**'
        parts = re.split(r'(?<!\*)\*(?!\*)', t)
        for part in parts:
            if '/' in part:
                a, b = part.split('/')
                part = a
                den *= Fraction(int(b))
            m = re.fullmatch(r'([XY])([012])(?:\*\*(\d+))?', part)
            if m:
                i = (0 if m.group(1) == 'X' else 3) + int(m.group(2))
                e[i] += int(m.group(3)) if m.group(3) else 1
            else:
                num *= Fraction(int(part))
        P[tuple(e)] = P.get(tuple(e), 0) + sign * num / den
    return {k: v for k, v in P.items() if v != 0}

def monomials(deg, par=None, nv=NV):
    out = []
    for c in itertools.combinations_with_replacement(range(nv), deg):
        e = [0] * nv
        for i in c: e[i] += 1
        if par is None or (e[3] + e[4] + e[5]) % 2 == par:
            out.append(tuple(e))
    out.sort(reverse=True)
    return out

def yparity(e): return (e[3] + e[4] + e[5]) % 2
def degree(P): return sum(next(iter(P))) if P else -1

def pmul(P, Q, p=None):
    R = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            R[e] = R.get(e, 0) + c1 * c2
    if p: R = {e: c % p for e, c in R.items() if c % p}
    else: R = {e: c for e, c in R.items() if c != 0}
    return R

def padd(P, Q, p=None, s=1):
    R = dict(P)
    for e, c in Q.items():
        R[e] = R.get(e, 0) + s * c
    if p: R = {e: c % p for e, c in R.items() if c % p}
    else: R = {e: c for e, c in R.items() if c != 0}
    return R

def pscale(P, c, p=None):
    if p: return {e: (v * c) % p for e, v in P.items() if (v * c) % p}
    return {e: v * c for e, v in P.items() if v * c != 0}

def pshift(P, m, c=1, p=None):
    R = {}
    for e, v in P.items():
        R[tuple(a + b for a, b in zip(e, m))] = (v * c) % p if p else v * c
    return {e: v for e, v in R.items() if v}

def pderiv(P, i, p=None):
    R = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            v = c * e[i]
            if p: v %= p
            if v: R[tuple(e2)] = v
    return R

def reduce_mod(P, p):
    R = {}
    for e, c in P.items():
        c = Fraction(c)
        v = (c.numerator * pow(c.denominator, -1, p)) % p
        if v: R[e] = v
    return R

def hilbert():
    """dim S_d^{+/-} and dim R_d^{+/-} (R = Jacobian ring of a smooth sextic in 6 vars, socle degree 24)."""
    def polymul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return out
    H = [1]
    for _ in range(6): H = polymul(H, [1] * 5)
    # tau-character: trace of tau on R_d: tau acts by (-1)^{Y-degree}; the Y-variables contribute (1 - t + t^2 - t^3 + t^4) each
    T = [1]
    for _ in range(3): T = polymul(T, [1] * 5)
    for _ in range(3): T = polymul(T, [1, -1, 1, -1, 1])
    H += [0] * 40; T += [0] * 40
    R = {d: ((H[d] + T[d]) // 2, (H[d] - T[d]) // 2) for d in range(0, 40)}
    S = {d: (len(monomials(d, 0)), len(monomials(d, 1))) for d in range(0, 30)}
    return S, R

def build_mat(cols, rowmons, p):
    """cols: list of polys mod p; rowmons: list of exponent tuples.  Returns nmod_mat."""
    idx = {m: i for i, m in enumerate(rowmons)}
    M = nmod_mat(len(rowmons), len(cols), p)
    for j, c in enumerate(cols):
        for e, v in c.items():
            M[idx[e], j] = v
    return M
