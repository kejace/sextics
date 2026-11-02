"""Independent (engine-free) polynomial helpers for the verify lens.
Polys = dict {exponent 6-tuple: coeff}; coeff = Fraction (exact) or int mod p.
Variables X0,X1,X2,Y0,Y1,Y2 = indices 0..5. Y-parity of monomial = (e3+e4+e5) % 2."""
import itertools, re
from fractions import Fraction
from flint import nmod_mat

NV = 6
VARS = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']

def parse(s):
    """parse a sympy-style polynomial string (sum of c*monomial terms, rational c) into a dict exactly."""
    import sympy
    syms = sympy.symbols(VARS)
    expr = sympy.sympify(s, locals=dict(zip(VARS, syms)))
    P = sympy.Poly(sympy.expand(expr), *syms)
    out = {}
    for mon, c in P.terms():
        c = sympy.Rational(c)
        out[tuple(int(x) for x in mon)] = Fraction(int(c.p), int(c.q))
    return out

def load_instance_F_A(path='/home/kejace/devel/sextics/data/quadratic-cubic-pullback/instance.rec'):
    txt = open(path).read()
    # F spans several lines (continuation lines start with whitespace)
    m = re.search(r'^F:\s*(.*?)(?=^\S)', txt, re.S | re.M)
    Fs = ' '.join(m.group(1).split())
    As = re.search(r'^A:\s*(.*)$', txt, re.M).group(1).strip()
    qs = int(re.search(r'^q:\s*(\d+)', txt, re.M).group(1))
    return parse(Fs), parse(As), qs

def ypar(e): return (e[3] + e[4] + e[5]) % 2
def deg(P): return max(sum(e) for e in P) if P else -1
def parity(P):
    s = {ypar(e) for e in P}
    return s.pop() if len(s) == 1 else None

def monomials(d, par=None):
    if d < 0: return []
    out = []
    for c in itertools.combinations_with_replacement(range(NV), d):
        e = [0]*NV
        for i in c: e[i] += 1
        e = tuple(e)
        if par is None or ypar(e) == par: out.append(e)
    return sorted(out)

def add(P, Q, p=None):
    out = dict(P)
    for e, v in Q.items():
        w = out.get(e, 0) + v
        if p: w %= p
        if w: out[e] = w
        else: out.pop(e, None)
    return out

def scale(P, c, p=None):
    return {e: (v*c % p if p else v*c) for e, v in P.items() if (v*c % p if p else v*c)}

def mul(P, Q, p=None):
    out = {}
    for e1, v1 in P.items():
        for e2, v2 in Q.items():
            e = tuple(a+b for a, b in zip(e1, e2))
            w = out.get(e, 0) + v1*v2
            if p: w %= p
            out[e] = w
    return {e: v for e, v in out.items() if v}

def shift(P, m, c=1, p=None):
    out = {}
    for e, v in P.items():
        w = v*c
        if p: w %= p
        if w: out[tuple(a+b for a, b in zip(e, m))] = w
    return out

def pderiv(P, i, p=None):
    out = {}
    for e, v in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            w = v*e[i]
            if p: w %= p
            if w: out[tuple(e2)] = w
    return out

def power(P, k, p=None):
    R = {(0,)*NV: 1}
    for _ in range(k): R = mul(R, P, p)
    return R

def modp(P, p):
    out = {}
    for e, v in P.items():
        if isinstance(v, Fraction):
            w = v.numerator * pow(v.denominator, -1, p) % p
        else:
            w = v % p
        if w: out[e] = w
    return out

def to_matrix_rows(polys, mons, p):
    col = {e: j for j, e in enumerate(mons)}
    flat = [0]*(len(polys)*len(mons))
    for r, P in enumerate(polys):
        base = r*len(mons)
        for e, v in P.items():
            flat[base + col[e]] = v % p
    return nmod_mat(len(polys), len(mons), flat, p)

def rank_mod_p(polys, mons, p):
    if not polys: return 0
    return to_matrix_rows(polys, mons, p).rank()

def rank_exact(polys, mons):
    """exact rank over Q via sympy (small matrices only)."""
    import sympy
    col = {e: j for j, e in enumerate(mons)}
    M = sympy.zeros(len(polys), len(mons))
    for r, P in enumerate(polys):
        for e, v in P.items():
            M[r, col[e]] = sympy.Rational(v.numerator, v.denominator)
    return M.rank()

def hilbert_R_parity(kmax=24):
    """bigraded Hilbert function of R = S/(six quintics, 3 Y-even + 3 Y-odd) assuming regular sequence:
    H(t,u) = (1-t^5)^3 (1-u t^5)^3 / ((1-t)^3 (1-u t)^3), u^2 = 1. Returns dict k -> (even, odd)."""
    # represent series as dict (k, par) -> coeff, truncate at kmax
    def smul(A, B):
        C = {}
        for (k1, p1), a in A.items():
            for (k2, p2), b in B.items():
                k = k1 + k2
                if k > kmax: continue
                key = (k, (p1 + p2) % 2)
                C[key] = C.get(key, 0) + a*b
        return C
    one = {(0, 0): 1}
    inv_even = {(k, 0): 1 for k in range(kmax+1)}          # 1/(1-t)
    inv_odd = {(k, k % 2): 1 for k in range(kmax+1)}       # 1/(1-ut)
    num_even = {(0, 0): 1, (5, 0): -1}                      # 1 - t^5
    num_odd = {(0, 0): 1, (5, 1): -1}                       # 1 - u t^5
    H = one
    for _ in range(3):
        H = smul(H, inv_even); H = smul(H, inv_odd); H = smul(H, num_even); H = smul(H, num_odd)
    return {k: (H.get((k, 0), 0), H.get((k, 1), 0)) for k in range(kmax+1)}

def hilbert_S_parity(k):
    return (len(monomials(k, 0)), len(monomials(k, 1)))
