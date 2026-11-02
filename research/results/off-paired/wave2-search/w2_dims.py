#!/usr/bin/env python3
"""E1-admissible denominator spaces per degree and character: dim (J:A)_{d,chi}, gauge F*R_{d-6}, Jacobian part
J_{d,chi}, non-Jacobian quotient; d = 4..dmax, all 36 characters (mod p = PRIMES[0])."""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); os.chdir(os.path.dirname(os.path.abspath(__file__)))
import json, time, itertools
import w2lib
from hodge import core
from flint import nmod_mat
F, A, q = w2lib.load()
p = w2lib.P0
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
out = {}
t0 = time.time()
for d in range(4, dmax + 1):
    for a in range(6):
        for b in range(6):
            c = (a, b)
            nD = len(w2lib.monos_char(d, c))
            vecs = w2lib.colon_vecs_modp(F, A, q, d, c, [{m: 1} for m in w2lib.monos_char(d, c)], p) if nD else []
            colon = nmod_mat(len(vecs), nD, [v % p for row in vecs for v in row], p).rank() if vecs else 0
            gens = w2lib.jacobian_generators_exact(F, d, c)
            mons = w2lib.monos_char(d, c); idx = {m: i for i, m in enumerate(mons)}
            if gens:
                rows = []
                for P in gens:
                    Pp = core.reduce_mod(P, p); row = [0] * nD
                    for e, v in Pp.items(): row[idx[e]] = v
                    rows.append(row)
                jdim = nmod_mat(len(rows), nD, [v for r in rows for v in r], p).rank()
            else:
                jdim = 0
            gauge = len(w2lib.monos_char(d - 6, c)) if d >= 6 else 0
            out[f'{d}:{a},{b}'] = {'d': d, 'c': [a, b], 'nD': nD, 'colon': colon, 'jac': jdim, 'nonjac': colon - jdim,
                                   'gauge': gauge, 'eff': colon - gauge, 'yparity': a % 2}
    tot = sum(v['colon'] for k, v in out.items() if v['d'] == d)
    totJ = sum(v['jac'] for k, v in out.items() if v['d'] == d)
    print(f'd={d}: total colon {tot} (J {totJ}, non-J {tot-totJ}); (0,0): {out[f"{d}:0,0"]}  [{time.time()-t0:.0f}s]', flush=True)
    json.dump(out, open('wave2-dims.json', 'w'), indent=0)
# markdown table
lines = ['| d | chi | nD | colon | J | non-J | gauge | eff |', '|---|---|---|---|---|---|---|---|']
for k, v in out.items():
    if v['colon'] > v['gauge']:
        lines.append(f"| {v['d']} | ({v['c'][0]},{v['c'][1]}) | {v['nD']} | {v['colon']} | {v['jac']} | {v['nonjac']} | {v['gauge']} | {v['eff']} |")
open('wave2-dims.md', 'w').write('\n'.join(lines) + '\n')
w2lib.append_json('wave2-search.remote.log.jsonl', {'ts': w2lib.now(), 'step': 'dims', 'dmax': dmax, 'file': 'wave2-dims.json',
                  'total_colon_by_d': {d: sum(v['colon'] for v in out.values() if v['d'] == d) for d in range(4, dmax + 1)},
                  'total_nonjac_by_d': {d: sum(v['nonjac'] for v in out.values() if v['d'] == d) for d in range(4, dmax + 1)},
                  'wall': round(time.time() - t0, 1)})
print('done', flush=True)
