import sys, json, time; sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
res = hess_quadrics(d, None, p=PRIME, verbose=True)
rows = quadric_rows(res)
nW = res['nW']; Wm = res['Wm']
ev = [1 if (m[i] == 1 and sum(m) == 1) else 0 for (i, m) in Wm]
vals = [sum(v * ev[k] * ev[l] for (k, l), v in row.items()) % PRIME for row in rows]
print('mixed d', d, 'nW', nW, 'rank', res['rank_quadrics'], 'dimR_2d(all)', res['dimR_2d_odd'], 'nonzero at Euler:', sum(1 for v in vals if v))
write_msolve(res, f'{WD}/hessmixed-quadrics-d{d}.ms')
jlog(f'{WD}/wave2-search.log.jsonl', {'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'wave': 2, 'family': 'hessian-relaxation-quadrics (mixed parity, local)', 'd': d, 'par': 'mixed', 'p': PRIME, 'nW': nW, 'npairs': res['npairs'], 'dimR_2d_both': res['dimR_2d_odd'], 'rank_quadrics': res['rank_quadrics'], 't_null': res['t_null'], 't_quad': res['t_quad']})
