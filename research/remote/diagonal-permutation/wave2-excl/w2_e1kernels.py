"""JOB excl-dp-e1k: re-verify the E1 kernel picture at fresh primes and extend it to (01)-sectors and d <= 12.

Outputs: e1k_<job>.json (all results), log lines exp='e1kernel'.
Rigor: ker_p = 0  ==>  ker_Q = 0 (modular rank is a lower bound for the rational rank);
       dim R_n mod p == Koszul value  ==>  good reduction in degree n (dim J_n mod p = dim J_n over Q)."""
import sys, os, time, json
import w2lib as W
from hodge.core import monomials
import numpy as np

job = os.environ.get('W2JOB', 'e1k')
OUT = f'{W.PKT}/e1k_{job}.json'
results = json.load(open(OUT)) if os.path.exists(OUT) else {}


def save():
    json.dump(results, open(OUT, 'w'), indent=1, default=str)


def run(d, p, par, sigma):
    key = f'd{d}_p{p}_par{par}_sig{sigma}'
    if key in results:
        return results[key]
    t0 = time.time()
    k = W.e1_kernel(d, p, par=par, sigma=sigma, verbose=True)
    good_d = (sigma is None and k['ring_dim'] == W.koszul_hilbert(d, par))
    good_d6 = (sigma is None and k['target_dim'] == W.koszul_hilbert(d + 6, 1 - par))
    rec = {'d': d, 'p': p, 'par': par, 'sigma': sigma, 'ring_dim': k['ring_dim'], 'target_dim': k['target_dim'],
           'rank': k['rank'], 'kernel_dim': k['dim'], 'jac_dim_d': k['jac_dim_d'], 'jac_dim_d6': k['jac_dim_d6'],
           'koszul_d': W.koszul_hilbert(d, par), 'koszul_d6': W.koszul_hilbert(d + 6, 1 - par),
           'good_reduction': (good_d and good_d6) if sigma is None else None, 'seconds': round(time.time() - t0, 1)}
    if 0 < k['dim'] <= 40:
        rec['basis_modp'] = [W.core.to_str({e: W.Fraction(v if v <= p // 2 else v - p) for e, v in b.items()}) for b in k['basis']]
    results[key] = rec
    save()
    W.log('e1kernel', d=d, prime=p, par=par, sigma=sigma, result=f"ker={k['dim']}", rank=k['rank'],
          ring_dim=k['ring_dim'], target_dim=k['target_dim'], seconds=rec['seconds'],
          note=('good reduction' if rec['good_reduction'] else '') + (' sector' if sigma is not None else ''))
    print(key, rec, flush=True)
    return rec


mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
if mode in ('all', 'full'):
    # full parity sectors, fresh primes: P1 for d=1..9, P2 for d=1..8
    for p, dmax in ((W.P1, 9), (W.P2, 8)):
        for d in range(1, dmax + 1):
            for par in (0, 1):
                run(d, p, par, None)
if mode in ('all', 'sectors'):
    # (01)-sectors at P1 for d = 8..12 (both parities, both signs); kernel dims add up to the full ones
    for d in range(8, 13):
        for par in (0, 1):
            for sigma in (1, -1):
                run(d, W.P1, par, sigma)
if mode == 'd9p2':
    for par in (0, 1):
        run(9, W.P2, par, None)
print('DONE', flush=True)
