"""Exact (LLL-small) E1-kernel elements in every sector at degree d + sealed exclusion of each.  usage: d threads"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
from sect import e1_kernel
from multiprocessing import Pool
d = int(sys.argv[1]); threads = int(sys.argv[2])
S = Setup(); F, A = S.F, S.A
t0 = time.time()
res, gens = e1_kernel(F, A, d, p=1000003, want_vectors=False)
sectors = [c for c, r in res.items() if r['ker_mod_J']]
print(f'd={d}: sectors with kernel (mod p dims): {[(c, res[c]["ker_mod_J"]) for c in sectors]}', flush=True)
def one_sector(c):
    try:
        return c, kernel_exact_lll(S, d, c)
    except Exception as ex:
        return c, (None, {'error': str(ex)[:300]})
elements = []; summary = {}
with Pool(min(threads, len(sectors))) as pool:
    for c, (results, info) in pool.imap_unordered(one_sector, sectors):
        summary[str(c)] = info
        print(f'  sector {c}: {info}', flush=True)
        log({'exp': 'e1-kernel-exact-lll', 'family': f'degree {d} sector {c}', 'd': d, 'prime': None,
             'result': f"exact kernel dim {info.get('dim')} (mod-p dim {res[c]['ker_mod_J']}), verified A*K in J: {info.get('verified')}, max |coeff| {info.get('max_abs_coeff')}",
             'rank': info.get('dim'), 'seconds': info.get('seconds'), 'job': f'excl-ctr-kexact{d}', 'note': 'exact fmpz nullspace + LLL + exact fmpq E1 solve'})
        if results is None: continue
        for i, r in enumerate(results):
            elements.append({'sector': list(c), 'index': i, 'Kz': {','.join(map(str, e)): v for e, v in r['Kz'].items()},
                             'Pz': [{','.join(map(str, e)): v for e, v in Pi.items()} for Pi in r['Pz']] if r.get('Pz') else None,
                             'verified': bool(r.get('verified')), 'L': r.get('L')})
json.dump({'d': d, 'summary': summary, 'elements': elements}, open(f'{PK}/kernel_exact_d{d}.json', 'w'))
print(f'saved kernel_exact_d{d}.json with {len(elements)} elements ({time.time()-t0:.1f}s)', flush=True)
def test_one(k):
    ent = elements[k]
    if not ent['verified']: return k, None
    K = {tuple(map(int, e.split(','))): Fraction(v) for e, v in ent['Kz'].items()}
    return k, sealed_exclusion(F, A, 1, K)
with Pool(threads) as pool:
    for k, out in pool.imap_unordered(test_one, range(len(elements))):
        if out is None: continue
        ent = elements[k]
        res_ = ('SEALED INCONSISTENT over Q' if out['sealed_inconsistent'] else ('CONSISTENT mod p ***' if out['consistent_mod_p'] else 'inconsistent mod p, not sealed'))
        log({'exp': 'e1-kernel-element-sealed-test', 'D': f"exact kernel element d={d} sector {ent['sector']} #{ent['index']}", 'd': d,
             'prime': ','.join(map(str, out['primes'])), 'result': res_, 'rank': out['rank'][0], 'seconds': out['seconds'], 'job': f'excl-ctr-kexact{d}',
             'note': f"sector system rows={out['rows']} cols={out['cols']} nullity={out['nullity']} rank_aug={out['rank_aug']} verified_null={out['verified_null']}"})
        print(f"element {k} sector {ent['sector']}: {res_} rank={out['rank']} aug={out['rank_aug']} null={out['nullity']} sealed={out['verified_null']} {out['seconds']}s", flush=True)
