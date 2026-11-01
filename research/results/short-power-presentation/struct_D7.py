import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
from fractions import Fraction
F, A, q, rec = load()
D = core.parse(json.load(open(sys.argv[1]))[0]['D'])
dY = [core.pderiv(D, i) for i in range(3, 6)]
print('dY1 == 3 dY0:', core.sub(dY[1], core.scale(dY[0], 3)) == {}, ' dY2 == 5 dY0:', core.sub(dY[2], core.scale(dY[0], 5)) == {})
dX = [core.pderiv(D, i) for i in range(3)]
print('dX1 == 2 dX0:', core.sub(dX[1], core.scale(dX[0], 2)) == {}, ' dX2 == 3 dX0:', core.sub(dX[2], core.scale(dX[0], 3)) == {})
# rewrite as polynomial in X0,X1,X2,M: substitute Y1=Y2=0, Y0 -> M (works iff D depends on Y only via M)
P = core.substitute_scalar(core.substitute_scalar(D, 4, 0), 5, 0)
# P is now a poly in X0,X1,X2,Y0 with Y0 standing for M
print('P(X, M) with Y0 := M, %d terms:' % len(P))
print(core.to_str(P))
# bidegree pattern in (X, M)
from collections import Counter
print('M-degrees present:', sorted(set(e[3] for e in P)))
# gcd / factor via sympy
import sympy as sp
X0, X1, X2, Mv = sp.symbols('X0 X1 X2 M')
expr = sum(sp.Rational(c.numerator, c.denominator) * X0**e[0] * X1**e[1] * X2**e[2] * Mv**e[3] for e, c in P.items())
print('factor:', sp.factor(expr))
