"""Norms of 'Y2^3 = g' cones: 5-var denominators f + lam*g^2 for binary cubics g (predicted active family)."""
import sys, pickle, os
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from five import *
from fractions import Fraction as Fr
os.makedirs('served', exist_ok=True)
P = core.parse
gs = {'s': 'Y0**3-Y1**3', 'S': 'Y0**3+Y1**3', 'Y0^3': 'Y0**3', 'Y1^3': 'Y1**3', 'Y0^2Y1': 'Y0**2*Y1', 'Y0Y1^2': 'Y0*Y1**2',
      'Y0^3+2Y1^3': 'Y0**3+2*Y1**3', 'Y0^3-2Y1^3': 'Y0**3-2*Y1**3', 'Y0^2Y1+Y0Y1^2': 'Y0**2*Y1+Y0*Y1**2',
      'Y0^3+Y0Y1^2': 'Y0**3+Y0*Y1**2', 'Y0^3-3Y0Y1^2': 'Y0**3-3*Y0*Y1**2', 'Y0^3+Y0^2Y1': 'Y0**3+Y0**2*Y1',
      'Y1^3+Y0^2Y1': 'Y1**3+Y0**2*Y1', 'Y0^3+Y1^3+Y0Y1^2': 'Y0**3+Y1**3+Y0*Y1**2',
      'X0^3': 'X0**3', 'X0^2Y0': 'X0**2*Y0', 'X0Y0^2': 'X0*Y0**2', 'X0^2Y0+Y1^3': 'X0**2*Y0+Y1**3', 'X0X1Y0': 'X0*X1*Y0',
      'c/X.. X0^3+X1^3': 'X0**3+X1**3', 'X0Y0Y1': 'X0*Y0*Y1', 'X1Y0^2': 'X1*Y0**2', 'Y0^3+X0^2X1': 'Y0**3+X0**2*X1'}
lams = [Fr(1), Fr(-1), Fr(1,3), Fr(-1,3), Fr(3), Fr(-3), Fr(2), Fr(-2), Fr(1,7), Fr(-1,7), Fr(1,4), Fr(4), Fr(1,2)]
if len(sys.argv) > 1: gs = {k: v for k, v in gs.items() if k in sys.argv[1].split(',')}
if len(sys.argv) > 2: lams = [Fr(x) for x in sys.argv[2].split(',')]
for gname, gexpr in gs.items():
    g = P(gexpr)
    for lam in lams:
        D = add(f, scale(mul(g, g), lam))
        if rparity_of(D) != 0: print('skip (rho-odd)', gname); break
        name = f'f+({lam})*({gname})^2'
        r = test5(name, 'cone-norm:f+lam*g^2', note=f'norm of Y2^3 = sqrt(-lam) g; g={gexpr}', Dpoly=D)
        fn = name.replace('/', '_').replace('*', 'x').replace('|', 'I') + '.pkl'
        with open('served/' + fn, 'wb') as fh:
            pickle.dump({'name': name, 'family': 'cone-norm:f+lam*g^2', 'D': {e: str(c) for e, c in D.items()}, 'p': r['p'], 'd': r['d'],
                         'served_basis': r['served_basis'], 'served_dim': r['served_dim'], 'served_mod_J': r['served_mod_J'],
                         'serves_A': r['serves_A'], 'serves_A+': r['serves_A+'], 'serves_A-': r['serves_A-'], 'sigma_sector': r['sigma_sector'],
                         'rank_N': r['rank_N'], 'nullity_N': r['nullity_N'], 'rows': r['rows'], 'cols': r['cols'], 'seconds': r['seconds']}, fh)
print('DONE')
