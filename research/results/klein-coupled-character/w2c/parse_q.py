import sys, sympy as sp
def parse_q(name):
    txt = open(name).read().strip().rstrip(':')
    data = sp.sympify(txt)
    body = data[1]
    deg = body[2]; par = body[5]
    w = par[1][0][1]
    A = sp.Symbol('A'); W = sp.Poly(list(reversed(list(w))), A)
    fac = [(sp.Poly(f,A).degree(), e) for f,e in sp.factor_list(W.as_expr())[1]]
    return int(deg), W.degree(), fac
if __name__ == '__main__':
    for name in sys.argv[1:]:
        try:
            deg, d, fac = parse_q(name)
            print(name, '| ideal degree', deg, '| #points over Qbar', d, '| factor degrees over Q', fac)
        except Exception as ex:
            print(name, 'ERR', str(ex)[:80], '| raw head:', open(name).read()[:40].replace('\n',' '))
