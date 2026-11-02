"""Direct full-system tests of explicit denominators. Usage: test_list.py <listname>"""
import sys, time
from common import *
S = add(mul(X0, X0), mul(Y0, Y0)); Pm = mul(X0, Y0)
Q1 = add(add(mul(mul(X1, X1), mul(Y1, Y1)), mul(mul(X1, X2), mul(Y1, Y2))), mul(mul(X2, X2), mul(Y2, Y2)))
Q2 = mul(mul(X1, Y1), add(mul(X2, X2), scale(mul(Y2, Y2), 2)))
GS = add(add(scale(mul(S, S), 3), scale(mul(Pm, Pm), 2)), scale(Q1, 5))
GP = add(scale(mul(Pm, S), 4), scale(Q2, Fraction(5, 7)))
D8 = sub(scale(mul(GS, GS), 4), mul(GP, GP))
lists = {
 'd9': [('D8*X1', mul(D8, X1)), ('D8*Y1', mul(D8, Y1)), ('D8*X2', mul(D8, X2)), ('D8*Y2', mul(D8, Y2)),
        ('D8*X0', mul(D8, X0)), ('D8*Y0', mul(D8, Y0)),
        ('D8*(X1+Y1)', mul(D8, add(X1, Y1))), ('D8*(X1-Y1)', mul(D8, sub(X1, Y1))), ('D8*(X2+Y2)', mul(D8, add(X2, Y2))),
        ('D8*(X0+Y0)', mul(D8, add(X0, Y0))), ('D8*(X0-Y0)', mul(D8, sub(X0, Y0))), ('D8*(X1+X2)', mul(D8, add(X1, X2))), ('D8*(X0+X2)', mul(D8, add(X0, X2)))],
}
for name, D in lists[sys.argv[1]]:
    t0 = time.time()
    r = test_D(D, 'D8*linear', note=name, verbose=True)
    print(name, '->', {k: v for k, v in r.items() if k not in ('system', 'res', 'out')}, f'({time.time()-t0:.0f}s)', flush=True)
