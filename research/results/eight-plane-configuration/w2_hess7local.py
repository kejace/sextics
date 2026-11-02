import sys, json, time; sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
for pD in (0, 1):
    res = hess_quadrics(7, pD, p=PRIME, verbose=True)
    write_msolve(res, f'{WD}/hess7local-quadrics-d7-p{pD}.ms')
    rows = quadric_rows(res)
    with open(f'{WD}/hess7local-quadrics-d7-p{pD}.json', 'w') as f:
        json.dump({'d': 7, 'par': pD, 'p': PRIME, 'nW': res['nW'], 'Wm': [[i, list(m)] for i, m in res['Wm']],
                   'rows': [[[list(k), v] for k, v in r.items()] for r in rows]}, f)
    jlog(f'{WD}/wave2-search.log.jsonl', {'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'wave': 2, 'family': 'hessian-relaxation-quadrics (local)', 'd': 7, 'par': pD, 'p': PRIME,
         'nW': res['nW'], 'npairs': res['npairs'], 'dimR_2d_odd': res['dimR_2d_odd'], 'rank_quadrics': res['rank_quadrics'], 't_null': res['t_null'], 't_quad': res['t_quad'],
         'note': 'spurious subspace g*Euler (dim 3) lies in the variety; emptiness of the complement not decidable by Macaulay ranks at this size'})
    print('d=7 par', pD, 'rank', res['rank_quadrics'], 'dimR', res['dimR_2d_odd'], flush=True)
