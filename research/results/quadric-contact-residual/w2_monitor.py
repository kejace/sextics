#!/usr/bin/env python3
"""Compact one-line summaries of new JSON log records read from stdin (dedup via a seen-file). Used by the Monitor."""
import sys, json, hashlib
seen_path = sys.argv[1]
try:
    seen = set(open(seen_path).read().split())
except FileNotFoundError:
    seen = set()
new = []
for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    h = hashlib.md5(line.encode()).hexdigest()
    if h in seen:
        continue
    seen.add(h); new.append(h)
    try:
        r = json.loads(line)
    except Exception:
        print('LOG:', line[:220]); continue
    t = r.get('task'); d = r.get('d')
    if t == 'varpro' and r.get('tag') in ('cone', 'uvcone') and (r.get('best_rel') or 1.0) > 1e-4 and r.get('status') == 'done':
        continue   # tiny cone families with ordinary minima: not worth a notification (they are in the report)
    if t == 'varpro':
        print(f"varpro d={d} chi={r.get('chi')} vars={r.get('vars')} tag={r.get('tag')} K={r.get('K')} dimP={r.get('dimP')} rows={r.get('rows')} "
              f"t_eval={r.get('t_eval')} restarts={r.get('restarts_done')} best={r.get('best_rel')} med={r.get('median_rel')} "
              f"n<1e-6={r.get('n_rel_lt_1e-6')} n<1e-10={r.get('n_rel_lt_1e-10')} s={r.get('seconds')} status={r.get('status')} "
              f"WITNESS={r.get('WITNESS')} err={str(r.get('error', ''))[:80]}")
    elif t == 'hess-relax':
        print(f"hess d={d} chi={r.get('chi')} members={r.get('n_members')} dimP={r.get('dimP')} dimR={r.get('dimR_2d')} restarts={r.get('restarts_done')} "
              f"best_frel={r.get('best_frel')} n<1e-12={r.get('n_frel_lt_1e-12')} polished={[p.get('full_rel') for p in r.get('polished', [])]} "
              f"s={r.get('seconds')} WITNESS={r.get('WITNESS')}")
    elif t == 'hess-gb':
        print(f"hess-gb d={d} chi={r.get('chi')} unknowns={r.get('unknowns')} quadrics={r.get('n_quadrics')} rank={r.get('quadric_rank')}/{r.get('n_pair_monomials')} "
              f"dim={r.get('dim')} vdim={r.get('vdim')} status={r.get('status')} s={r.get('seconds')}")
    elif t in ('struct', 'struct-exact'):
        print(f"struct d={d} {r.get('name')} chi={r.get('chi')} e1={r.get('e1')} consistent={r.get('consistent')} status={r.get('status')} "
              f"rank={r.get('rank')} null={r.get('nullity')} s={r.get('seconds')}")
    else:
        print(f"{t} d={d} " + json.dumps({k: v for k, v in r.items() if k not in ('task', 'd', 'ts', 'host', 'engine', 'nontrivial', 'per_sector_extra')})[:260])
with open(seen_path, 'a') as f:
    f.write(''.join(h + '\n' for h in new))
