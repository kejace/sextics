import sys, ast
import sympy as sp
from parse_msolve import parse
path = sys.argv[1]
p, nvars, deg, vars_, lin, w, den, params = parse(path)
A = sp.Symbol('A')
W = sp.Poly(list(reversed(w)), A, modulus=p)
Den = sp.Poly(list(reversed(den)), A, modulus=p)
fac = sp.factor_list(W.as_expr(), modulus=p)[1]
for f, e in fac:
    q = sp.Poly(f, A, modulus=p)
    # work in F_p[A]/(q): represent elements as polys reduced mod q
    dinv = sp.invert(Den.as_expr(), q.as_expr(), modulus=p) if q.degree() > 1 else None
    coords = []
    for item in params:
        coeffs = item[0][1] if isinstance(item[0], list) else item[1]
        P = sp.Poly(list(reversed(coeffs)), A, modulus=p)
        if q.degree() == 1:
            r = (-q.all_coeffs()[-1]) % p
            val = (-P.eval(r) * pow(int(Den.eval(r)), -1, p)) % p
            coords.append(val)
        else:
            val = sp.rem(sp.Poly(-P.as_expr()*dinv, A, modulus=p), q)
            coords.append(val.as_expr())
    print('factor', q.as_expr(), '-> point:', dict(zip(vars_[:len(coords)], coords)))
