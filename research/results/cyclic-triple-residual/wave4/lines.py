"""Lines of W = {G=0} ⊂ P(1,1,1,1,1,3), G(x0,x1,z,y0,y1,y2) = F(x0,x1,X2,y0,y1,y2) with z = X2^3.
A weighted line: (s,t) -> (l(s,t), c(s,t)), l linear in C^5 (coords x0,x1,y0,y1,y2), c a binary cubic.
Incidence: the 7 coefficients of G(l(s,t), c(s,t)) in s^k t^(6-k) vanish."""
import sympy as sp
from fractions import Fraction
s, t = sp.symbols('s t')
x0, x1, z, y0, y1, y2 = sp.symbols('x0 x1 z y0 y1 y2')
c57 = sp.Rational(5, 7)
G = (x0**6 + c57*x0**3*x1*y0*y1 + c57*x0**3*z + c57*x0**2*y0*y1*y2**2 + c57*x0*x1*y1**3*y2
     + c57*x0*z*y1*y2 + x1**6 + c57*x1*z*y0*y1 + z**2 + y0**6 + y1**6 + y2**6)
Lcub = c57*(x0**3 + x0*y1*y2 + x1*y0*y1)          # G = z^2 + L z + M
Msex = sp.expand(G - z**2 - Lcub*z)
assert sp.expand(G - (z**2 + Lcub*z + Msex)) == 0

def incidence(lin, cub):
    """lin: 5 sympy exprs linear in s,t (x0,x1,y0,y1,y2); cub: binary cubic in s,t. Returns 7 coefficient polys."""
    e = sp.expand(G.subs({x0: lin[0], x1: lin[1], y0: lin[2], y1: lin[3], y2: lin[4], z: cub}, simultaneous=True))
    P = sp.Poly(e, s, t)
    out = []
    for k in range(7):
        out.append(sp.expand(P.coeff_monomial(s**k * t**(6-k))))
    return out

def msolve_file(path, polys, vars_, p):
    """Write an msolve input file over F_p (coefficients reduced mod p, 1/7 -> inverse)."""
    lines = [','.join(str(v) for v in vars_), str(p)]
    strs = []
    for f in polys:
        P = sp.Poly(sp.expand(f), *vars_)
        terms = []
        for mon, co in P.terms():
            co = Fraction(int(co.p), int(co.q)) if isinstance(co, sp.Rational) else Fraction(co)
            cm = (co.numerator * pow(co.denominator, -1, p)) % p
            if cm == 0: continue
            m = '*'.join(f'{v}^{e}' if e > 1 else f'{v}' for v, e in zip(vars_, mon) if e > 0)
            terms.append(f'{cm}*{m}' if m else f'{cm}')
        strs.append('+'.join(terms) if terms else '0')
    lines.append(',\n'.join(strs))
    open(path, 'w').write('\n'.join(lines) + '\n')

def singular_file(path, polys, vars_, p, extra=''):
    strs = []
    for f in polys:
        P = sp.Poly(sp.expand(f), *vars_)
        terms = []
        for mon, co in P.terms():
            co = Fraction(int(co.p), int(co.q))
            cm = (co.numerator * pow(co.denominator, -1, p)) % p
            if cm == 0: continue
            m = '*'.join(f'{v}^{e}' if e > 1 else f'{v}' for v, e in zip(vars_, mon) if e > 0)
            terms.append(f'{cm}*{m}' if m else f'{cm}')
        strs.append('+'.join(terms) if terms else '0')
    txt = f'ring r = {p}, ({",".join(map(str,vars_))}), dp;\nideal I = ' + ',\n'.join(strs) + ';\n' + extra
    open(path, 'w').write(txt)
