#!/usr/bin/env python3
"""First-order deformation of the 12 Fermat subset-Euler certificates (D_I = sum_{i in I} Z_i^6, N0 = A E_I/6, Q0 = 0)
along F_eps = F0 + eps F1 (F1 = 9973(F - F0)): consistency of the linear order-eps system in (D1, N1, Q1), mod p.
Re-run of wave-1 deform.py (which finished 4 of 12 seeds) with a new prime, all 12 seeds, parallel.
Usage: python3 wave2_excl_deform.py prime workers"""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
import deform
from multiprocessing import Pool
p = int(sys.argv[1]); workers = int(sys.argv[2])
deform.p = p
deform.F0p = core.reduce_mod(deform.F0, p); deform.F1p = core.reduce_mod(deform.F1, p); deform.Ap = core.reduce_mod(deform.A, p)
def one(I):
    t0 = time.time()
    r = deform.run(I, allow_D1=True)
    r.update({'exp': 'fermat-first-order-deformation', 'prime': p, 'seconds': round(time.time() - t0, 1),
              'result': 'first-order system CONSISTENT mod p' if r['consistent'] else 'first-order system inconsistent mod p (rank_aug = rank + 1)'})
    log2(**r); print(json.dumps(r), flush=True)
    return r
if __name__ == '__main__':
    subsets = [I for I in itertools.combinations(range(6), 3) if (0 in I) != (3 in I)]
    with Pool(workers) as pool: res = pool.map(one, subsets)
    with open(f'{OUT}/deform_first_order_p{p}.json', 'w') as f: json.dump(res, f, indent=1)
