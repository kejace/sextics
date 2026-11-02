"""Referee helpers (own implementation; only core.load_instance is borrowed, for parsing instance.rec)."""
import sys, itertools, json, os, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction
from flint import nmod_mat, fmpq_mat, fmpq
HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, 'referee.log.jsonl')
NV = 6
def log(**rec):
    rec.setdefault('lens', 'referee'); rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(LOG, 'a') as f: f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True); return rec
def load():
    F, A, q, rec = core.load_instance('dense-rational'); assert q == 1
    return F, A
def yp(e): return (e[3] + e[4] + e[5]) & 1
def pi(i): return 0 if i < 3 else 1
def monos(deg, par=None):
    if deg < 0: return []
    out = []
    for c in itertools.combinations_with_replacement(range(NV), deg):
        e = [0] * NV
        for i in c: e[i] += 1
        e = tuple(e)
        if par is None or yp(e) == par: out.append(e)
    return out
def modp(P, p):
    out = {}
    for e, c in P.items():
        c = Fraction(c); v = c.numerator * pow(c.denominator, -1, p) % p
        if v: out[e] = v
    return out
def pd(P, i, p=None):
    out = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1; v = c * e[i]
            if p: v %= p
            if v: out[tuple(e2)] = v
    return out
def shift(P, m, p=None, c=1):
    out = {}
    for e, v in P.items():
        e2 = tuple(a + b for a, b in zip(e, m)); vv = v * c
        out[e2] = vv % p if p else vv
    return out
def add(P, Q, p=None):
    out = dict(P)
    for e, c in Q.items():
        v = out.get(e, 0) + c
        if p: v %= p
        if v: out[e] = v
        else: out.pop(e, None)
    return out
def sub(P, Q, p=None): return add(P, {e: (-c % p if p else -c) for e, c in Q.items()}, p)
def mul(P, Q, p=None):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2)); v = out.get(e, 0) + c1 * c2
            if p: v %= p
            if v: out[e] = v
            else: out.pop(e, None)
    return out
def hilbert():
    def pm(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): out[i + j] += x * y
        return out
    Hs = [1]
    for _ in range(6): Hs = pm(Hs, [1] * 5)
    Th = [1]
    for _ in range(3): Th = pm(Th, [1, 0, 1, 0, 1, 0, 1, 0, 1])
    Hs += [0] * 40; Th += [0] * 40
    R = {d: ((Hs[d] + Th[d]) // 2, (Hs[d] - Th[d]) // 2) for d in range(30)}
    S = {d: (len(monos(d, 0)), len(monos(d, 1))) for d in range(26)}
    return S, R
def fill(M, cols, ridx):
    for j, col in enumerate(cols):
        for e, v in col.items(): M[ridx[e], j] = v
def pivots_of(Mr, rk):
    piv = []
    for r in range(rk):
        j = 0
        while int(Mr[r, j]) == 0: j += 1
        piv.append(j)
    return piv
def fill_rows(M, rows, cidx):
    for r, row in enumerate(rows):
        for e, v in row.items(): M[r, cidx[e]] = v
