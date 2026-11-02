import sys; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation')
from hlib import *
from hodge.core import pderiv, mul, add, sub, scale, to_str
V = ['X0','X1','X2','Y0','Y1','Y2']
Fi = [pderiv(F, i) for i in range(6)]
tests = []
# d=5: partials (X: Y-even, Y: Y-odd), and symmetric combos
for i in range(6):
    tests.append((f'd/d{V[i]} F', Fi[i], 'J5-partial'))
tests.append(('F_X0+F_X1+F_X2', add(add(Fi[0],Fi[1]),Fi[2]), 'J5-sym'))
tests.append(('F_X0-F_X1', sub(Fi[0],Fi[1]), 'J5-std'))
tests.append(('F_Y0-F_Y1', sub(Fi[3],Fi[4]), 'J5-std-odd'))
tests.append(('F_Y0+F_Y1+F_Y2', add(add(Fi[3],Fi[4]),Fi[5]), 'J5-sym-odd'))
# d=6: X_i F_Xj, Y_i F_Yj, block Euler, pair radial, std-type
X = [core.parse(v) for v in V]
for i in range(3):
    for j in range(3):
        tests.append((f'{V[i]}*F_{V[j]}', mul(X[i], Fi[j]), 'J6-basis'))
        tests.append((f'{V[i+3]}*F_{V[j+3]}', mul(X[i+3], Fi[j+3]), 'J6-basis'))
EX = {}
for i in range(3): EX = add(EX, mul(X[i], Fi[i]))
tests.append(('E_X(F)=sum X_i F_Xi', EX, 'J6-blockEuler'))
for i in range(3):
    tests.append((f'pair radial {i}: X{i}F_X{i}+Y{i}F_Y{i}', add(mul(X[i],Fi[i]), mul(X[i+3],Fi[i+3])), 'J6-pairradial'))
gA = sub(add(mul(X[0],Fi[0]), mul(X[3],Fi[3])), add(mul(X[1],Fi[1]), mul(X[4],Fi[4])))
tests.append(('g_A(F) = (X0d0+Y0d3-X1d1-Y1d4)F', gA, 'J6-std'))
tests.append(('(X0d0-X1d1)F', sub(mul(X[0],Fi[0]), mul(X[1],Fi[1])), 'J6-std'))
tests.append(('(Y0d3-Y1d4)F', sub(mul(X[3],Fi[3]), mul(X[4],Fi[4])), 'J6-std'))
# transposition-type: (X0 d1 + X1 d0 + Y0 d4 + Y1 d3) F  (swap generator) and (X0 d1 - X1 d0 + Y0 d4 - Y1 d3)F (rotation)
tests.append(('(X0d1+X1d0+Y0d4+Y1d3)F', add(add(mul(X[0],Fi[1]), mul(X[1],Fi[0])), add(mul(X[3],Fi[4]), mul(X[4],Fi[3]))), 'J6-swap'))
tests.append(('(X0d1-X1d0+Y0d4-Y1d3)F', add(sub(mul(X[0],Fi[1]), mul(X[1],Fi[0])), sub(mul(X[3],Fi[4]), mul(X[4],Fi[3]))), 'J6-rot'))
tests.append(('(X0d1-X1d0)F', sub(mul(X[0],Fi[1]), mul(X[1],Fi[0])), 'J6-rotX'))
tests.append(('(Y0d4-Y1d3)F', sub(mul(X[3],Fi[4]), mul(X[4],Fi[3])), 'J6-rotY'))
# 3-cycle-type: (X1 d0 + X2 d1 + X0 d2 + Y1 d3 + Y2 d4 + Y0 d5)F
c3 = {}
for i in range(3):
    c3 = add(c3, mul(X[(i+1)%3], Fi[i])); c3 = add(c3, mul(X[3+(i+1)%3], Fi[3+i]))
tests.append(('cyclic shift field (X_{i+1}d_i + Y_{i+1}d_{i+3})F', c3, 'J6-cyc'))
for name, D, fam in tests:
    e = test_D(name, fam, D=D)
    print(name, '-> d', e['d'], 'e1', e['e1'], 'cons', e.get('consistent'), 'rank', e['rank'], 'null', e['nullity'], e['seconds'], 's', flush=True)
