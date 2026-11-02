# Referee's own minimal polynomial toolkit (dict exps->int mod p); engine used only to load (F, A).
import sys
from fractions import Fraction
from itertools import combinations_with_replacement
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
def load(p):
    F, A, q, _ = core.load_instance('symmetric-sextic-residual')
    def modp(P):
        out = {}
        for e, c in P.items():
            c = Fraction(c); v = (c.numerator * pow(c.denominator, -1, p)) % p
            if v: out[tuple(int(x) for x in e)] = v
        return out
    return modp(F), modp(A), q
def ypar(e): return (e[3] + e[4] + e[5]) % 2
_cache = {}
def monomials(n, par=None):
    key = (n, par)
    if key not in _cache:
        out = []
        for comb in combinations_with_replacement(range(6), n):
            e = [0]*6
            for i in comb: e[i] += 1
            e = tuple(e)
            if par is None or ypar(e) == par: out.append(e)
        _cache[key] = out
    return _cache[key]
def pderiv(P, i, p):
    out = {}
    for e, c in P.items():
        if e[i]:
            f = list(e); f[i] -= 1
            out[tuple(f)] = (c * e[i]) % p
    return {e: c for e, c in out.items() if c}
def mulmono(P, m, c0, p):
    return {tuple(a + b for a, b in zip(e, m)): (c * c0) % p for e, c in P.items()}
def mul(P, Q, p):
    out = {}
    for e, c in P.items():
        for f, d in Q.items():
            g = tuple(a + b for a, b in zip(e, f)); out[g] = (out.get(g, 0) + c * d) % p
    return {e: c for e, c in out.items() if c}
def add(P, Q, p, s=1):
    out = dict(P)
    for e, c in Q.items(): out[e] = (out.get(e, 0) + s * c) % p
    return {e: c for e, c in out.items() if c}
def deg(P): return sum(next(iter(P))) if P else None
def series_coeff(n):
    """[t^n] of ((1-t^5)/(1-t))^6 and ((1-t^10)/(1-t^2))^3 -> CI dims of R_n^even, R_n^odd"""
    def poly_pow(base, k, nmax):
        res = [1] + [0]*nmax
        for _ in range(k):
            new = [0]*(nmax+1)
            for i, a in enumerate(res):
                if a == 0: continue
                for j, b in enumerate(base):
                    if i + j <= nmax: new[i+j] += a*b
            res = new
        return res
    full = poly_pow([1,1,1,1,1], 6, n)[n]
    alt = poly_pow([1,0,1,0,1,0,1,0,1], 3, n)[n]
    return (full + alt)//2, (full - alt)//2   # (even, odd)
def dimR(n, par): return series_coeff(n)[par]
def dimJ_CI(n, par): return len(monomials(n, par)) - dimR(n, par)
