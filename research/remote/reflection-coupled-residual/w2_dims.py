"""E1-admissible spaces per (d, chi): two primes (+ third if non-Jacobian part), CRT lift of canonical reps.
Usage: python3 w2_dims.py <d> <a> <b>   -> writes e1space-d<d>-chi<a><b>.json and logs a batch line."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual')
import w2rcr
from hodge import core
d = int(sys.argv[1]); chi = (int(sys.argv[2]), int(sys.argv[3]))
F, A, q = w2rcr.load()
t0 = time.time()
r1 = w2rcr.e1_space_modp_v2(F, A, d, chi, w2rcr.PRIMES[0], verbose=True)
out = {k: v for k, v in r1.items() if k not in ('reps_modp', 'Dmons', 'Jpivots')}
out['primes'] = [w2rcr.PRIMES[0]]
if r1['dim_nonjac'] > 0:
    r2 = w2rcr.e1_space_modp_v2(F, A, d, chi, w2rcr.PRIMES[1], verbose=True)
    r3 = w2rcr.e1_space_modp_v2(F, A, d, chi, w2rcr.PRIMES[2], verbose=True)
    out['primes'] = w2rcr.PRIMES[:3]
    out['dims_agree'] = (r1['dim_colon'], r1['dim_J']) == (r2['dim_colon'], r2['dim_J']) == (r3['dim_colon'], r3['dim_J'])
    polys = w2rcr.crt_lift_reps(r1, r2, r3)
    if polys is None:
        out['rational_lift'] = False
        out['reps_modp'] = {'p': r1['p'], 'Dmons': r1['Dmons'], 'reps': r1['reps_modp']}
    else:
        out['rational_lift'] = True
        out['reps'] = [core.to_str(P) for P in polys]
        out['reps_nterms'] = [len(P) for P in polys]
        out['reps_maxden'] = [max(c.denominator for c in P.values()) for P in polys]
out['total_seconds'] = round(time.time() - t0, 1)
json.dump(out, open(f'e1space-d{d}-chi{chi[0]}{chi[1]}.json', 'w'))
w2rcr.log_batch({'step': 'e1-dims', 'd': d, 'chi': list(chi), 'dim_S': out['dim_S'], 'dim_colon': out['dim_colon'],
                 'dim_J': out['dim_J'], 'dim_nonjac': out['dim_nonjac'], 'primes': out['primes'],
                 'rational_lift': out.get('rational_lift'), 'seconds': out['total_seconds']})
print(json.dumps({k: v for k, v in out.items() if k not in ('reps', 'reps_modp')}))
