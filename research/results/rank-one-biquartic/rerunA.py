import lib
from lib import *
import random
Fi = [core.pderiv(F, i) for i in range(6)]
items = [(RX, '6f+4tP (eps_X F)', 'structural-v2', 'ramification divisor R; engine v2 rerun'),
         (RY, '6g+2tP (eps_Y F)', 'structural-v2', 'engine v2 rerun')]
# J_6^{even} basis elements: X_j F_i (i,j<3), Y_j F_i (i,j>=3)
names = ['X0','X1','X2','Y0','Y1','Y2']
for i in range(6):
    for j in range(6):
        if (i < 3) == (j < 3):
            items.append((core.mul(P(names[j]), Fi[i]), f'{names[j]}*F_{i}', 'jac6', 'J_6 even basis element'))
random.seed(1)
basis = [core.mul(P(names[j]), Fi[i]) for i in range(6) for j in range(6) if (i < 3) == (j < 3)]
for k in range(4):
    D = {}
    cs = [random.randint(-3, 3) for _ in basis]
    for c, b in zip(cs, basis):
        D = core.add(D, core.scale(b, c))
    items.append((D, f'random J_6 combo {cs}', 'jac6-random', ''))
for D, expr, fam, note in items:
    r = run_D(D, expr, fam, note + ' [engine v2]')
    show(r)
