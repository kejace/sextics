import sys, ast, json
import sympy as sp
def parse(path):
    txt = open(path).read().strip().rstrip(':')
    data = ast.literal_eval(txt)
    dim, body = data[0], data[1]
    p, nvars, deg, vars_, lin, par = body[0], body[1], body[2], body[3], body[4], body[5]
    _, [[dw, w], [dd, den], params] = par
    return p, nvars, deg, vars_, lin, w, den, params
if __name__ == '__main__':
    p, nvars, deg, vars_, lin, w, den, params = parse(sys.argv[1])
    A = sp.Symbol('A')
    W = sp.Poly(list(reversed(w)), A, modulus=p)
    print('p', p, 'ideal degree', deg, 'elim degree', W.degree())
    fac = sp.factor_list(W.as_expr(), modulus=p)
    print('factorization degrees mod p:', [(sp.Poly(f, A, modulus=p).degree(), e) for f, e in fac[1]])
    Den = sp.Poly(list(reversed(den)), A, modulus=p)
    roots = [r for r in range(p) if W.eval(r) % p == 0] if W.degree() <= 20 and p < 2*10**6 else []
    print('F_p roots:', roots)
    for r in roots:
        dv = Den.eval(r) % p
        pt = []
        for item in params:
            coeffs = item[0][1] if isinstance(item[0], list) else item[1]
            P = sp.Poly(list(reversed(coeffs)), A, modulus=p)
            val = (-P.eval(r) * pow(int(dv), -1, p)) % p
            pt.append(val if val < p//2 else val - p)
        print('root', r, '->', dict(zip(vars_[:len(pt)], pt)))
