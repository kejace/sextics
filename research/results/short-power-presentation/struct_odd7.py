import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from collections import Counter
F, A, q, rec = load()
J = json.load(open('kernel_d7_p1.json'))
D = core.parse(J['basis_rational'][0])
print('odd d=7 kernel element: terms', len(D), 'deg', core.degree(D), 'parity', core.parity_of(D), flush=True)
dX = [core.pderiv(D, i) for i in range(3)]
print('dX1 == 2 dX0:', core.sub(dX[1], core.scale(dX[0], 2)) == {}, ' dX2 == 3 dX0:', core.sub(dX[2], core.scale(dX[0], 3)) == {}, flush=True)
print('X-degrees present:', sorted(Counter(sum(e[:3]) for e in D).items()), flush=True)
print('Y-degrees present:', sorted(Counter(sum(e[3:]) for e in D).items()), flush=True)
print('max num bits', max(abs(c.numerator).bit_length() for c in D.values()), 'max den bits', max(c.denominator.bit_length() for c in D.values()), flush=True)
# compare with D7 even: substitute to see if D7odd(Y, L) relates to D7even(X, M) by X<->Y, L<->M with coefficient maps
D7 = core.parse(json.load(open('kernel_d7_p0.json'))['basis_rational'][0])
print('D7 even: terms', len(D7), 'max num bits', max(abs(c.numerator).bit_length() for c in D7.values()), 'max den bits', max(c.denominator.bit_length() for c in D7.values()), flush=True)
# D7 as polynomial in (X0,X1,X2,M): coefficients via Y1=Y2=0, Y0->M
P = core.substitute_scalar(core.substitute_scalar(D7, 4, 0), 5, 0)
print('D7(X, M) terms:', len(P), flush=True)
import sympy as sp
X0, X1, X2, Mv = sp.symbols('X0 X1 X2 M')
expr = sum(sp.Rational(c.numerator, c.denominator) * X0**e[0] * X1**e[1] * X2**e[2] * Mv**e[3] for e, c in P.items())
expr = sp.expand(expr / sp.Rational(list(P.values())[0].numerator, list(P.values())[0].denominator))
print('D7(X,M) (normalized by first coeff), M^4 part:', sp.factor(expr.coeff(Mv, 4)), flush=True)
print('M^2 part:', sp.factor(expr.coeff(Mv, 2)), flush=True)
print('M^0 part factor:', sp.factor(expr.coeff(Mv, 0)), flush=True)
Q = core.substitute_scalar(core.substitute_scalar(D, 1, 0), 2, 0)   # X1=X2=0, X0 -> L
Y0, Y1, Y2, Lv = sp.symbols('Y0 Y1 Y2 L')
expr2 = sum(sp.Rational(c.numerator, c.denominator) * Y0**e[3] * Y1**e[4] * Y2**e[5] * Lv**e[0] for e, c in Q.items())
expr2 = sp.expand(expr2 / sp.Rational(list(Q.values())[0].numerator, list(Q.values())[0].denominator))
print('odd element as poly in (Y, L): L-degrees', sorted(set(e[0] for e in Q)), flush=True)
print('L^4 part:', sp.factor(expr2.coeff(Lv, 4)), flush=True)
print('L^2 part:', sp.factor(expr2.coeff(Lv, 2)), flush=True)
print('L^0 part factor:', sp.factor(expr2.coeff(Lv, 0)), flush=True)
