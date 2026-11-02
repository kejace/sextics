import sys; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation')
from hlib import *
cands = []
# degree 2
for s in ['X0**2+X1**2+X2**2', 'Y0**2+Y1**2+Y2**2', 'X0**2+Y0**2', 'X0**2-Y0**2', 'X0*X1+Y0*Y1', 'X0*X1-Y0*Y1',
          '(X0-X1)**2+(Y0-Y1)**2', '(X0-X1)**2-(Y0-Y1)**2', '(X0+X1)**2+(Y0+Y1)**2', '(X0-X1)**2+(X1-X2)**2+(X0-X2)**2',
          '(Y0-Y1)**2+(Y1-Y2)**2+(Y0-Y2)**2', 'X2**2+Y2**2', 'X0**2+X1**2+Y0**2+Y1**2', 'X0**2', 'X0*X1', 'Y0*Y1', 'X0**2+X1**2',
          'X0**2+Y1**2', 'X0*X1+X1*X2+X2*X0', 'Y0*Y1+Y1*Y2+Y2*Y0', 'X0**2+X1**2+X2**2+Y0**2+Y1**2+Y2**2', 'X0*Y1+X1*Y0', 'X0*Y1-X1*Y0']:
    cands.append((s, 'quadric'))
# degree 3 (Y-even cubics)
for s in ['X0**3+X1**3+X2**3', 'X0*X1*X2', '(X0-X1)*(X1-X2)*(X2-X0)', 'X0**3+X1**3', 'X0*(Y1**2-Y2**2)', 'X0*Y0**2+X1*Y1**2+X2*Y2**2',
          'X0*Y0**2', 'X0*Y1*Y2+X1*Y0*Y2+X2*Y0*Y1', 'X0**3+X1**3+X2**3-3*X0*X1*X2', 'X0*Y0**2-X1*Y1**2', 'X0**2*X1-X1**2*X0',
          '(X0+X1+X2)*(Y0**2+Y1**2+Y2**2)', '(X0+X1+X2)**3', 'X0*(X1**2+Y1**2)', 'X2*((X0-X1)**2+(Y0-Y1)**2)', 'X0**3+Y0**2*X1']:
    cands.append((s, 'cubic'))
for s, fam in cands:
    e = test_D(s, fam)
    print(s, '->', 'e1' if e['e1'] else 'E1-FAIL', e.get('consistent'), e['rank'], e['nullity'], e['seconds'], flush=True)
