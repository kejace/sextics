"""Independent verifier library (does NOT import hodge.*).  Polys = dict{exp6: coeff}.
Coefficients: Fraction (exact) or int mod p."""
import re, json, time, itertools, os
from fractions import Fraction
from math import comb
import sympy as sp

REC = '/home/kejace/devel/sextics/data/short-power-presentation/instance.rec'
VARS = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
LOG = '/tmp/hodge/results/short-power-presentation/wave2-verify.log.jsonl'

def logjson(**kw):
    kw['t'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    with open(LOG, 'a') as f:
        f.write(json.dumps(kw) + '\n')

def read_rec():
    txt = open(REC).read()
    fields = {}
    key = None
    for line in txt.splitlines():
        if line.startswith('#'): continue
        m = re.match(r'^([A-Za-z]+):\s*(.*)$', line)
        if m and not line.startswith(' '):
            key = m.group(1); fields[key] = m.group(2)
        elif key is not None:
            fields[key] += ' ' + line.strip()
    return fields

def sympy_F():
    f = read_rec()
    return sp.sympify(f['F'], locals=dict(zip([str(v) for v in VARS], VARS))), sp.sympify(f['A'], locals=dict(zip([str(v) for v in VARS], VARS))), int(f['q'])

def poly_from_sympy(expr):
    P = sp.Poly(sp.expand(expr), *VARS)
    out = {}
    for mon, c in P.terms():
        out[tuple(mon)] = Fraction(int(c.p), int(c.q))
    return out

def degree(P):
    return max(sum(e) for e in P) if P else -1

def ypar(e): return (e[3] + e[4] + e[5]) % 2

def parity_of(P):
    s = {ypar(e) for e in P}
    return s.pop() if len(s) == 1 else None

def monomials(d, par=None):
    out = []
    for e in itertools.product(range(d + 1), repeat=5):
        s = sum(e)
        if s > d: continue
        m = e + (d - s,)
        if par is None or ypar(m) == par: out.append(m)
    return out

def reduce_mod(P, p):
    out = {}
    for e, c in P.items():
        v = (c.numerator * pow(c.denominator, -1, p)) % p if isinstance(c, Fraction) else c % p
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

def shift(P, m, c=1, p=None):
    out = {}
    for e, v in P.items():
        w = v * c
        if p: w %= p
        if w: out[tuple(a + b for a, b in zip(e, m))] = w
    return out

def add(P, Q, p=None):
    out = dict(P)
    for e, v in Q.items():
        w = out.get(e, 0) + v
        if p: w %= p
        if w: out[e] = w
        else: out.pop(e, None)
    return out

def mul(P, Q, p=None):
    out = {}
    for e, v in P.items():
        for f, w in Q.items():
            g = tuple(a + b for a, b in zip(e, f))
            x = out.get(g, 0) + v * w
            if p: x %= p
            if x: out[g] = x
            else: out.pop(g, None)
    return out

def scale(P, c, p=None):
    return shift(P, (0,)*6, c, p)

# ---- Hilbert series of the Jacobian ring of a smooth sextic in 6 vars, tau-refined
def hilb():
    t = sp.symbols('t')
    H = sp.Poly(sp.expand((1 + t + t**2 + t**3 + t**4)**6), t)
    Th = sp.Poly(sp.expand((1 + t**2 + t**4 + t**6 + t**8)**3), t)
    def coeff(P, e):
        return int(P.coeff_monomial(t**e)) if e <= P.degree() else 0
    out = {}
    for e in range(0, 30):
        h = coeff(H, e); th = coeff(Th, e)
        S = comb(e + 5, 5)
        sig = comb(e // 2 + 2, 2) if e % 2 == 0 else 0
        Sp, Sm = (S + sig) // 2, (S - sig) // 2
        Rp, Rm = (h + th) // 2, (h - th) // 2
        out[e] = dict(S_plus=Sp, S_minus=Sm, R_plus=Rp, R_minus=Rm, J_plus=Sp - Rp, J_minus=Sm - Rm)
    return out

def jac_gens(Fp, deg, par, p):
    """Generators m*F_i of J_deg^{par} mod p."""
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    cols = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(deg - 5, (par - pi) % 2):
            cols.append(shift(Fi[i], m, 1, p))
    return cols

def nmod_from_polys(polys, p, rows=None):
    from flint import nmod_mat
    if rows is None:
        rows = {}
        for P in polys:
            for e in P: rows.setdefault(e, len(rows))
    R, C = len(rows), len(polys)
    flat = [0] * (R * C)
    for j, P in enumerate(polys):
        for e, v in P.items(): flat[rows[e] * C + j] = v
    return nmod_mat(R, C, flat, p), rows
