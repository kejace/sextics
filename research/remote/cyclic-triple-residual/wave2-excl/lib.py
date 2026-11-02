"""Shared helpers for the wave-2 exclusion packet of lane cyclic-triple-residual.

Polynomials are dicts {exponent 6-tuple: coefficient}; coefficients are ints (after clearing
denominators) or ints mod p.  Variables X0,X1,X2,Y0,Y1,Y2 = indices 0..5.
"""
import sys, time, json, itertools, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction
from math import lcm
from itertools import product
from flint import fmpz_mat, nmod_mat, fmpq_mat, fmpq, fmpz

LANE = 'cyclic-triple-residual'
LOG = f'/tmp/hodge/results/{LANE}/wave2-excl.log.jsonl'
PACKET = f'/tmp/hodge/results/{LANE}/wave2-excl'
P1 = 2147483647          # wave-1 prime
P2 = 1000003             # fqlx prime
P3 = 1073741789          # new prime (largest prime < 2^30)
P4 = 998244353

def log(rec):
    rec = dict(rec); rec.setdefault('engine', 2); rec.setdefault('lens', 'excl'); rec.setdefault('time', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(LOG, 'a') as f: f.write(json.dumps(rec, default=str) + '\n')

def receipt(name, obj):
    os.makedirs(f'{PACKET}/receipts', exist_ok=True)
    with open(f'{PACKET}/receipts/{name}.json', 'w') as f: json.dump(obj, f, indent=1, default=str)

def integral(P):
    L = 1
    for c in P.values(): L = lcm(L, Fraction(c).denominator)
    return {m: int(Fraction(c) * L) for m, c in P.items() if c}

def load():
    F, A, q, rec = core.load_instance(LANE)
    return F, A, q, integral(F), integral(A)

def mons(deg):
    if deg < 0: return
    for c in itertools.combinations(range(deg+5), 5):
        prev = -1; e = []
        for x in c: e.append(x-prev-1); prev = x
        e.append(deg+4-prev); yield tuple(e)

def imul(p, q, mod=None):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(a+b for a, b in zip(m1, m2)); r[m] = r.get(m, 0) + c1*c2
    if mod: return {m: c % mod for m, c in r.items() if c % mod}
    return {m: c for m, c in r.items() if c}

def iadd(p, q, s=1):
    r = dict(p)
    for m, c in q.items(): r[m] = r.get(m, 0) + s*c
    return {m: c for m, c in r.items() if c}

def ideriv(p, i):
    r = {}
    for m, c in p.items():
        if m[i]:
            mm = list(m); mm[i] -= 1; r[tuple(mm)] = r.get(tuple(mm), 0) + c*m[i]
    return r

def ipow(p, k):
    r = {(0,)*6: 1}
    for _ in range(k): r = imul(r, p)
    return r

E6 = [tuple(int(j == i) for j in range(6)) for i in range(6)]

def group_elements(Fz):
    """All g in (Z/6)^6 with g5 = 0 and g.m = 0 mod 6 for every monomial m of F (projective diagonal group)."""
    ms = list(Fz); out = []
    for g in product(range(6), repeat=5):
        g = g + (0,)
        if all(sum(a*b for a, b in zip(m, g)) % 6 == 0 for m in ms): out.append(g)
    return out

def semi_invariant_subgroup(G, polys):
    out = []
    for g in G:
        ok = True
        for P in polys:
            ms = list(P)
            c0 = sum(a*b for a, b in zip(ms[0], g)) % 6
            if any(sum(a*b for a, b in zip(m, g)) % 6 != c0 for m in ms[1:]): ok = False; break
        if ok: out.append(g)
    return out

def generators(G):
    gens = []; span = {tuple([0]*6)}
    for g in G:
        if g in span: continue
        gens.append(g); new = set()
        for s in span:
            for k in range(6): new.add(tuple((a + k*b) % 6 for a, b in zip(s, g)))
        span = new
    assert len(span) == len(G)
    return gens

def make_chi(gens):
    def chi(m): return tuple(sum(a*b for a, b in zip(m, g)) % 6 for g in gens)
    return chi

def cadd(*cs): return tuple(sum(x) % 6 for x in zip(*cs))
def cneg(c): return tuple((-x) % 6 for x in c)

def chi_of(P, chi):
    cs = set(chi(m) for m in P); assert len(cs) == 1, 'not semi-invariant'; return cs.pop()

def rank_of_columns(cols, p=None):
    """Exact (p None) or modular rank of a matrix given as a list of sparse columns {rowkey: value}."""
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    if R == 0 or C == 0: return 0
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e]*C + j] = v if p is None else v % p
    M = fmpz_mat(R, C, flat) if p is None else nmod_mat(R, C, flat, p)
    return M.rank()

def ci_dim_J(d):
    """dim J(F)_d for a smooth sextic in 6 variables (complete intersection of six quintics):
    dim S_d - coefficient of t^d in (1-t^5)^6/(1-t)^6."""
    from math import comb
    def dimS(k): return comb(k+5, 5) if k >= 0 else 0
    r = 0
    for k in range(7):
        r += (-1)**k * comb(6, k) * dimS(d - 5*k)
    return dimS(d) - r
