"""Tangency-along-a-surface family: D = P0^2 + lam*f*P1^2 (norm of P0 + Y2^3 P1 up to scaling), random P0, P1."""
import sys, pickle, random, os
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from five import *
from fractions import Fraction as Fr
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
k = int(sys.argv[2]) if len(sys.argv) > 2 else 1   # deg P1 = k, deg P0 = k+3, deg D = 2k+6
n = int(sys.argv[3]) if len(sys.argv) > 3 else 4
os.makedirs('served', exist_ok=True)
for trial in range(n):
    # sigma-pure choice: P0 and P1 sigma-pure with P0^2, P1^2 even automatically; rho: P0 rho-odd or even, P1 likewise
    sp0 = random.randint(0, 1); rp0 = random.randint(0, 1); sp1 = random.randint(0, 1); rp1 = random.randint(0, 1)
    P0 = {m: Fr(random.randint(-2, 2)) for m in mons5(k + 3, rp0) if spar(m) == sp0}
    P1 = {m: Fr(random.randint(-2, 2)) for m in mons5(k, rp1) if spar(m) == sp1}
    P0 = {e: c for e, c in P0.items() if c}; P1 = {e: c for e, c in P1.items() if c}
    lam = random.choice([Fr(1), Fr(-1), Fr(3), Fr(-3), Fr(1, 3), Fr(-1, 3)])
    D = add(mul(P0, P0), scale(mul(f, mul(P1, P1)), lam))
    name = f'tangent:P0^2+({lam})*f*P1^2 k={k} #{trial} chars({rp0}{sp0},{rp1}{sp1})'
    r = test5(name, f'tangent-norm:k{k}', note=f'P0={core.to_str(P0)[:80]}; P1={core.to_str(P1)[:60]}', Dpoly=D)
    fn = name.replace('/', '_').replace('*', 'x').replace('|', 'I').replace(' ', '_') + '.pkl'
    with open('served/' + fn, 'wb') as fh:
        pickle.dump({'name': name, 'family': f'tangent-norm:k{k}', 'D': {e: str(c) for e, c in D.items()}, 'p': r['p'], 'd': r['d'],
                     'served_basis': r['served_basis'], 'served_dim': r['served_dim'], 'served_mod_J': r['served_mod_J'],
                     'serves_A': r['serves_A'], 'serves_A+': r['serves_A+'], 'serves_A-': r['serves_A-'], 'sigma_sector': r['sigma_sector'],
                     'rank_N': r['rank_N'], 'nullity_N': r['nullity_N'], 'rows': r['rows'], 'cols': r['cols'], 'seconds': r['seconds']}, fh)
print('DONE')
