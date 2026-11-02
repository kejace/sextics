"""Shared helpers for the wave-2 exclusion packet of quadratic-cubic-pullback.
Usage: sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl'); from qx import *
"""
import sys, os, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from fractions import Fraction
from flint import nmod_mat, fmpq_mat, fmpq

IID = 'quadratic-cubic-pullback'
LANE = '/tmp/hodge/results/quadratic-cubic-pullback'
PK = LANE + '/wave2-excl'
LOG = LANE + '/wave2-excl.log.jsonl'
P1 = 2147483647      # engine prime (wave 1)
P2 = 1073741827      # second prime (this lens)
P3 = 1000003         # third prime

QSTR = [
    "X0**2 + 5*X1*X2/7 + 5*Y0*Y1/7",
    "5*X0*X2/11 + X1**2 + 5*Y1*Y2/11",
    "5*X0*X1/13 + X2**2 + 5*Y0*Y2/13",
    "5*X0*X1/17 + Y0**2 + 5*Y1*Y2/17",
    "5*X1*X2/19 + 5*Y0*Y2/19 + Y1**2",
    "5*X0*X2/23 + 5*Y0*Y1/23 + Y2**2",
]

def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q

def quadrics():
    return [core.parse(s) for s in QSTR]

def log(rec):
    rec = dict(rec); rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S')); rec.setdefault('engine', 2)
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec

def receipt(name, obj):
    fn = os.path.join(PK, name)
    with open(fn, 'w') as f:
        json.dump(obj, f, indent=1)
    return fn

def hilbert():
    """(S[d] = (even, odd) monomial counts, R[d] = (even, odd) Jacobian-ring dims) for a smooth sextic in 6 vars."""
    def pm(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): out[i + j] += x * y
        return out
    H = [1]
    for _ in range(6): H = pm(H, [1, 1, 1, 1, 1])
    Th = [1]
    for _ in range(3): Th = pm(Th, [1, 0, 1, 0, 1, 0, 1, 0, 1])
    H += [0] * 40; Th += [0] * 40
    R = {d: ((H[d] + Th[d]) // 2, (H[d] - Th[d]) // 2) for d in range(30)}
    S = {d: (len(monomials(d, 0)), len(monomials(d, 1))) for d in range(40)}
    return S, R

def det_poly(M, p=None):
    n = len(M)
    if n == 1: return M[0][0]
    out = {}
    for j in range(n):
        if not M[0][j]: continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = mul(M[0][j], det_poly(minor, p), p)
        if j % 2: term = scale(term, -1, p)
        out = add(out, term, p)
    return out

def delta(p=None):
    Q = quadrics()
    if p: Q = [core.reduce_mod(q_, p) for q_ in Q]
    M = [[pderiv(Q[i], j, p) for j in range(NV)] for i in range(NV)]
    return det_poly(M, p)

def mat_from_polys(polys, mons, p):
    idx = {m: j for j, m in enumerate(mons)}
    R, C = len(polys), len(mons)
    flat = [0] * (R * C)
    for r, P_ in enumerate(polys):
        base = r * C
        for e, v in P_.items():
            flat[base + idx[e]] = v % p
    return nmod_mat(R, C, flat, p)

def fmpq_from_polys(polys, mons):
    idx = {m: j for j, m in enumerate(mons)}
    R, C = len(polys), len(mons)
    M = fmpq_mat(R, C)
    for r, P_ in enumerate(polys):
        for e, v in P_.items():
            v = Fraction(v)
            M[r, idx[e]] = fmpq(v.numerator, v.denominator)
    return M
