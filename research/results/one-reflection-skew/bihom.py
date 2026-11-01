"""Random bihomogeneous D of bidegree (a,b) (X-degree a, Y-degree b even): generic activity per bidegree."""
import sys, pickle, random, os
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from five import *
from fractions import Fraction
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
os.makedirs('served', exist_ok=True)
degs = [(7,0),(5,2),(3,4),(1,6),(8,0),(6,2),(4,4),(2,6),(0,8)]
if len(sys.argv) > 2: degs = [tuple(int(x) for x in s.split(':')) for s in sys.argv[2].split(',')]
for (a, b) in degs:
    for k in range(2):
        D = {}
        for m in mons5(a + b, 0):
            if m[0] + m[1] + m[2] == a:
                D[m] = Fraction(random.randint(-3, 3))
        D = {e: c for e, c in D.items() if c}
        name = f'random-bihom-({a},{b})-#{k}'
        r = test5(name, 'random-bihomogeneous', note=f'random bihomogeneous bidegree ({a},{b})', Dpoly=D)
        with open('served/' + name + '.pkl', 'wb') as fh:
            pickle.dump({'name': name, 'family': 'random-bihomogeneous', 'D': {e: str(c) for e, c in D.items()}, 'p': r['p'], 'd': r['d'],
                         'served_basis': r['served_basis'], 'served_dim': r['served_dim'], 'served_mod_J': r['served_mod_J'],
                         'serves_A': r['serves_A'], 'serves_A+': r['serves_A+'], 'serves_A-': r['serves_A-'],
                         'rank_N': r['rank_N'], 'nullity_N': r['nullity_N'], 'rows': r['rows'], 'cols': r['cols'], 'seconds': r['seconds']}, fh)
print('DONE')
