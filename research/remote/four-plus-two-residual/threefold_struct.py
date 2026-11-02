"""Structured core-only denominators in the reduced (threefold) model."""
import sys, random
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual')
from threefold import *
random.seed(21)
p = core.PRIMES[0]
P = core.parse
f = P('X0**6 + 5*X0**5*X1/7 + 10*X0*X1**5/7 + X1**6')
g = P('Y0**6 + 5*Y0**5*Y1/7 + 15*Y0*Y1**5/7 + Y1**6')
K = P('5*X0**3*X1*Y0*Y1/7 + 10*X0*X1*Y0**3*Y1/7')
Hi = {i: pderiv(H, i) for i in CORE}
X0, X1, Y0, Y1 = P('X0'), P('X1'), P('Y0'), P('Y1')
tests = []
fam = 'threefold-structural'
tests += [(f, 'f(X0,X1)'), (g, 'g(Y0,Y1)'), (K, 'coupling K'), (add(f, g), 'f+g'), (add(f, scale(g, -1)), 'f-g'),
          (add(f, scale(g, 2)), 'f+2g'), (add(f, scale(g, Fraction(1, 2))), 'f+g/2'), (add(add(f, g), scale(K, -1)), 'f+g-K'),
          (add(add(f, g), scale(K, 2)), 'f+g+2K (H+K)'), (add(add(f, g), scale(K, Fraction(1, 2))), 'f+g+K/2'),
          (mul(f, g), 'f*g (d=12)'), (mul(X0, mul(X1, mul(Y0, Y1))), 'X0X1Y0Y1'),
          (P('X0**2*Y0**2'), 'X0^2 Y0^2'), (P('X0**3*Y0**3'), 'a itself'), (P('X0**6+Y0**6'), 'X0^6+Y0^6'),
          (P('X0**6+X1**6+Y0**6+Y1**6'), 'pure powers'), (P('X0**6+X1**6+Y0**6'), 'cone X0^6+X1^6+Y0^6'),
          (add(f, P('Y0**6')), 'cone f+Y0^6'), (add(g, P('X0**6')), 'cone g+X0^6'),
          (P('X0**4+X1**4+Y0**4'), 'cone quartic'), (P('X0**3*X1+X1**4+Y0**4'), 'cone quartic 2'),
          (P('X0**2+Y0**2'), 'X0^2+Y0^2'), (P('X0**2-Y0**2'), 'X0^2-Y0^2'), (P('X0*X1+Y0*Y1'), 'X0X1+Y0Y1'),
          (P('X0**2+2*Y0**2'), 'X0^2+2Y0^2 (quadratic factor of K)'), (P('X0*Y0'), 'X0Y0')]
# Jacobian multiples (pass E1 automatically): partials, random linear/quadratic combos, products of partials
for i in CORE: tests.append((Hi[i], f'H_{core.VARS[i]} (d=5)'))
for r in range(3):
    D = {}
    for i in CORE:
        l = {}
        for m in cmon(1): l[m] = Fraction(random.choice([-2, -1, 0, 1, 2, 3]))
        D = add(D, mul(l, Hi[i]))
    tests.append((D, f'random linear combo of partials sum l_i H_i (d=6) #{r}'))
for r in range(3):
    D = {}
    for i in CORE:
        l = {}
        for m in cmon(2): l[m] = Fraction(random.choice([-2, -1, 0, 1, 2, 3]))
        D = add(D, mul(l, Hi[i]))
    tests.append((D, f'random quadratic combo sum Q_i H_i (d=7) #{r}'))
tests.append((mul(Hi[0], Hi[3]), 'H_X0 * H_Y0 (d=10)'))
tests.append((mul(Hi[0], Hi[1]), 'H_X0 * H_X1 (d=10)'))
tests.append((add(mul(Hi[0], Hi[0]), mul(Hi[3], Hi[3])), 'H_X0^2 + H_Y0^2 (d=10)'))
# 2x2 Hessian minors (d=8) and Hessian 3x3 minors (d=12)
def hess(i, j): return pderiv(Hi[i], j)
Hm = {(i, j): hess(i, j) for i in CORE for j in CORE}
for (i, j) in [(0, 3), (0, 1), (3, 4), (1, 4)]:
    m2 = add(mul(Hm[(i, i)], Hm[(j, j)]), scale(mul(Hm[(i, j)], Hm[(i, j)]), -1))
    tests.append((m2, f'2x2 Hessian minor ({core.VARS[i]},{core.VARS[j]}) (d=8)'))
# radial derivatives
for (i, j, name) in [(0, 3, 'X0,Y0'), (1, 4, 'X1,Y1'), (0, 1, 'X0,X1'), (3, 4, 'Y0,Y1')]:
    R = add(mul(P(core.VARS[i]), Hi[i]), mul(P(core.VARS[j]), Hi[j]))
    tests.append((R, f'({name}) radial derivative of H (d=5)'))
# discriminant-like: H as a cubic in X0^2? skip. Perturbations of f+g by monomials
for mname in ['X0**3*X1*Y0*Y1', 'X0*X1*Y0**3*Y1', 'X0**2*X1**2*Y0**2', 'X0**4*Y0**2', 'X0**2*Y0**4']:
    for t in [1, -1, Fraction(5, 7)]:
        tests.append((add(add(f, g), scale(P(mname), t)), f'f+g+({t})*{mname}'))
for D, note in tests:
    if core.degree(D) > 12: continue
    run3(D, note, fam, p)
