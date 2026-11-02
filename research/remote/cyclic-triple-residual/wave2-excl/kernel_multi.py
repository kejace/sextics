"""Exact E1-kernel elements (mod J) in every sector at degree d, with exact E1 solutions, then sealed exclusion tests.
usage: python3 kernel_multi.py d nprimes threads"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
from sect import e1_kernel
from multiprocessing import Pool
d = int(sys.argv[1]); nprimes = int(sys.argv[2]); threads = int(sys.argv[3])
S = Setup(); F, A = S.F, S.A
primes = core.PRIMES[:nprimes]
t0 = time.time()
res, gens = e1_kernel(F, A, d, p=primes[0], want_vectors=False)
sectors = [c for c, r in res.items() if r['ker_mod_J']]
print(f'd={d}: sectors with kernel: {[(c, res[c]["ker_mod_J"]) for c in sectors]} ({time.time()-t0:.1f}s)', flush=True)
elements = []; summary = {}
for c in sectors:
    results, info = kernel_exact(S, d, c, primes)
    summary[str(c)] = info
    print(f'  sector {c}: {info}', flush=True)
    if results is None: continue
    for i, r in enumerate(results):
        if r is None:
            elements.append({'sector': list(c), 'index': i, 'Kz': None, 'Pz': None, 'verified': False}); continue
        elements.append({'sector': list(c), 'index': i, 'Kz': {','.join(map(str, e)): v for e, v in r['Kz'].items()},
                         'Pz': [{','.join(map(str, e)): v for e, v in Pi.items()} for Pi in r['Pz']] if r.get('Pz') else None,
                         'verified': bool(r.get('verified')), 'L': r.get('L')})
    log({'exp': 'e1-kernel-exact-elements', 'family': f'degree {d} sector {c}', 'd': d, 'prime': ','.join(map(str, primes)),
         'result': f"dim {info.get('dim')} (mod p), reconstructed {info.get('reconstructed')}, exactly verified A*K in J: {info.get('verified')}",
         'rank': info.get('dim'), 'seconds': info.get('seconds'), 'job': f'excl-ctr-kernel{d}', 'note': 'CRT over primes + ratrecon + exact fmpq E1 solve; verified elements give a lower bound on the Q-kernel'})
json.dump({'d': d, 'primes': primes, 'summary': summary, 'elements': elements}, open(f'{PK}/kernel_exact_d{d}.json', 'w'))
print(f'saved kernel_exact_d{d}.json with {len(elements)} elements ({time.time()-t0:.1f}s)', flush=True)
# sealed exclusion of every verified element
def test_one(k):
    ent = elements[k]
    if not ent['verified']: return k, None
    K = {tuple(map(int, e.split(','))): Fraction(v) for e, v in ent['Kz'].items()}
    out = sealed_exclusion(F, A, 1, K)
    return k, out
with Pool(threads) as pool:
    for k, out in pool.imap_unordered(test_one, range(len(elements))):
        if out is None: continue
        ent = elements[k]
        log({'exp': 'e1-kernel-element-sealed-test', 'D': f"kernel element d={d} sector {ent['sector']} #{ent['index']}", 'd': d,
             'prime': ','.join(map(str, out['primes'])), 'result': ('SEALED INCONSISTENT over Q' if out['sealed_inconsistent'] else
             ('CONSISTENT mod p ***' if out['consistent_mod_p'] else 'inconsistent mod p, not sealed')),
             'rank': out['rank'][0], 'seconds': out['seconds'], 'job': f'excl-ctr-kernel{d}',
             'note': f"sector system rows={out['rows']} cols={out['cols']} nullity={out['nullity']} rank_aug={out['rank_aug']} verified_null={out['verified_null']}"})
        print(f"element {k} sector {ent['sector']}: {out['sealed_inconsistent']} rank={out['rank']} aug={out['rank_aug']} null={out['nullity']} sealed={out['verified_null']} {out['seconds']}s", flush=True)
