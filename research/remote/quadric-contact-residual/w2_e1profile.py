#!/usr/bin/env python3
"""Exact (over Q) E1-admissible denominator spaces per degree and G_aff character sector.

  python3 w2_e1profile.py --dmin 4 --dmax 12 [--cones]

For each d and each character c occurring in degree d: dim S_d[c], dim J_d[c], dim (J:A)_d[c], #extras,
gauge dim F*S_{d-6}[c].  With --cones also the kernels restricted to coordinate variable subsets
(3, 4 and 5 variables; cones with vertex a coordinate plane / line / point over curves / surfaces / 3-folds)
per sector.  One JSON line per degree in wave2-search.log.jsonl; full tables in e1_sectors_wave2.json.
"""
import sys, os, json, time, argparse, itertools
sys.path.insert(0, '/tmp/hodge/results/quadric-contact-residual')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import w2lib as W

ap = argparse.ArgumentParser()
ap.add_argument('--dmin', type=int, default=4)
ap.add_argument('--dmax', type=int, default=12)
ap.add_argument('--cones', action='store_true')
ap.add_argument('--out', default='e1_sectors_wave2.json')
a = ap.parse_args()

path = os.path.join(W.RES, a.out)
table = json.load(open(path)) if os.path.exists(path) else {}
for d in range(a.dmin, a.dmax + 1):
    t0 = time.time()
    key = str(d)
    table.setdefault(key, {'sectors': {}, 'cones': {}})
    tot = {'dimS': 0, 'dimJ': 0, 'kernel': 0, 'extra': 0, 'gauge': 0}
    for c in W.chars_in_degree(d):
        cs = W.chi_str(c)
        if cs in table[key]['sectors']:
            r = table[key]['sectors'][cs]
        else:
            fam = W.e1_family(d, c)
            info = fam['info']
            r = {'dimS': info['dimS'], 'dimJ': info['dimJ'], 'kernel': info['dim_kernel'], 'extra': info['n_extra'],
                 'gauge': info['dim_gauge'], 'seconds': info['seconds'],
                 'extras': [W.core.to_str({e: v for e, v in E.items()}) for E in fam['extras']][:40]}
            table[key]['sectors'][cs] = r
            json.dump(table, open(path, 'w'))
        for k in tot:
            tot[k] += r[k]
        print(f'd={d} chi={cs} dimS={r["dimS"]} J={r["dimJ"]} ker={r["kernel"]} extra={r["extra"]} gauge={r["gauge"]} ({r["seconds"]}s)', flush=True)
    rec = {'task': 'e1-profile', 'd': d, 'n_sectors': len(table[key]['sectors']), **{'tot_' + k: v for k, v in tot.items()},
           'per_sector_extra': {cs: r['extra'] for cs, r in table[key]['sectors'].items() if r['extra']},
           'seconds': round(time.time() - t0, 1)}
    W.log_line(rec)
    print(json.dumps(rec), flush=True)
    if a.cones and d >= 6:
        t1 = time.time()
        found = {}
        for nv in (3, 4, 5):
            for sub in itertools.combinations(range(6), nv):
                mask = [i in sub for i in range(6)]
                name = ''.join(W.core.VARS[i] for i in sub)
                for c in W.chars_in_degree(d):
                    cs = W.chi_str(c)
                    k2 = f'{name}:{cs}'
                    if k2 in table[key]['cones']:
                        r = table[key]['cones'][k2]
                    else:
                        mons = W.mons_sector(d, c, mask)
                        if not mons:
                            continue
                        fam = W.e1_family(d, c, vars_mask=mask)
                        info = fam['info']
                        r = {'dimS': info['dimS'], 'dimJ': info['dimJ'], 'kernel': info['dim_kernel'], 'extra': info['n_extra'],
                             'gauge': info['dim_gauge']}
                        table[key]['cones'][k2] = r
                    if r['kernel'] > r['gauge']:
                        found[k2] = r
            json.dump(table, open(path, 'w'))
        rec = {'task': 'e1-cones', 'd': d, 'n_nontrivial': len(found), 'nontrivial': found, 'seconds': round(time.time() - t1, 1)}
        W.log_line(rec)
        print(json.dumps(rec), flush=True)
print('done')
