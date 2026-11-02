"""Job: E1-admissible denominator spaces per degree and G-sector (dims mod p, exact bases over Q).
usage: python3 w2_e1space.py dmin dmax [procs] [exact|dims]"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import time, json
import multiprocessing as mp
from w2core import *
dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
procs = int(sys.argv[3]) if len(sys.argv) > 3 else 8
mode = sys.argv[4] if len(sys.argv) > 4 else 'exact'
JOB = f'w2e1space-{dmin}-{dmax}'
F, A, q = load()
G = full_group(F, A)
S = Sect(F, A, G)

def work(args):
    d, c = args
    t0 = time.time()
    try:
        r = e1_space(S, F, A, q, d, c, exact=(mode == 'exact'))
        if mode == 'exact':
            save_space(r, f'{WD}/w2-e1space-d{d}-c{csec(c)}.json')
        out = {k: v for k, v in r.items() if k not in ('basis', 'Dm')}
        out['n_basis'] = len(r.get('basis', []))
        return out
    except Exception as ex:
        import traceback
        return {'d': d, 'sector': list(c), 'error': repr(ex), 'tb': traceback.format_exc()[-800:], 'seconds': round(time.time() - t0, 1)}

for d in range(dmin, dmax + 1):
    t0 = time.time()
    secs = sorted(S.sectors(d))
    with mp.get_context('fork').Pool(procs) as pool:
        results = pool.map(work, [(d, c) for c in secs])
    ok = [r for r in results if 'error' not in r]
    tot = {k: sum(r[k] for r in ok) for k in ('dim_S', 'dim_J', 'kerA_in_S', 'ker_mod_J', 'forced', 'dim_R_d')}
    tot['sporadic'] = tot['ker_mod_J'] - tot['forced']
    entry = {'batch': 'e1space', 'job': JOB, 'd': d, 'mode': mode, 'totals': tot, 'n_sectors': len(secs),
             'sectors': {csec(r['sector']): {k: r[k] for k in ('dim_S', 'dim_J', 'kerA_in_S', 'ker_mod_J', 'dim_R_d', 'dim_R_t', 'forced', 'seconds') if k in r} for r in ok},
             'errors': [r for r in results if 'error' in r], 'seconds': round(time.time() - t0, 1)}
    w2log(entry)
    print(f'd={d} {tot} ({entry["seconds"]}s) errors={len(entry["errors"])}', flush=True)
    for r in sorted(ok, key=lambda r: r['sector']):
        if r['ker_mod_J']:
            print('   sector', r['sector'], {k: r[k] for k in ('dim_S', 'dim_J', 'kerA_in_S', 'ker_mod_J', 'dim_R_d', 'dim_R_t', 'forced', 'seconds')}, flush=True)
    for r in entry['errors']:
        print('   ERROR', r['sector'], r['error'], r['tb'][-300:], flush=True)
print('done', flush=True)
