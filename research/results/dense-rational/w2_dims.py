"""E1-admissible denominator spaces per degree: ker([A] : (R/J)_d -> (R/J)_{d+6}) via hodge.jacobian.mult_kernel,
both parity sectors, d = 4..12 (mod p, two primes for d >= 10).  Logs one batch line per (d, parity) and saves the
kernel bases (standard-monomial representatives mod p) for d >= 10 to w2_kernel_d{d}_par{par}_p{p}.json.
Usage: python3 w2_dims.py dmin dmax [parities] [primes]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from w2lib import *
dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
pars = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0, 1]
primes = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else [core.PRIMES[0]]
F, A, q = load()
for d in range(dmin, dmax + 1):
    for par in pars:
        for p in primes:
            t0 = time.time()
            out = jacobian.mult_kernel(F, A, d, p=p, parity=par, verbose=True)
            rec = {'batch': 'e1-kernel-dims', 'family': 'mult_kernel', 'd': d, 'sector': 'even' if par == 0 else 'odd', 'p': p,
                   'dim_S_d': len(monomials(d, par)), 'dim_J_d': out['jacobian_dim'], 'dim_RJ_d': out['ring_dim'],
                   'kernel_mod_J': out['dim'], 'kernel_full': out['kernel_dim'], 'rows': out['rows'], 'cols': out['cols'],
                   'seconds': round(time.time() - t0, 1), 'note': 'E1-admissible D space = J_d + kernel; kernel_full = dim{D in S_d : A D in J}'}
            print(json.dumps(rec), flush=True)
            log_batch(rec)
            if d >= 10 and out['dim']:
                with open(f'{LANE}/w2_kernel_d{d}_par{par}_p{p}.json', 'w') as f:
                    json.dump({'p': p, 'd': d, 'par': par, 'dim': out['dim'],
                               'basis': [{' '.join(map(str, k)): v for k, v in b.items()} for b in out['basis']]}, f)
