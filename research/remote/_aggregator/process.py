#!/usr/bin/env python3
"""Fold verify-*.out.jsonl results into /tmp/hodge/notes/wave2-verification.log.jsonl and print a summary.

usage: python3 process.py <results.jsonl> [<results2.jsonl> ...]
"""
import sys, json, os, time
LOG = '/tmp/hodge/notes/wave2-verification.log.jsonl'

have = {}
if os.path.exists(LOG):
    for line in open(LOG):
        try:
            o = json.loads(line); have[o['check']] = o
        except Exception:
            pass

new = []
for path in sys.argv[1:]:
    for line in open(path):
        try:
            o = json.loads(line)
        except Exception:
            continue
        name = o['name']
        if name in have:
            continue
        spec = o.get('spec', {})
        res = o.get('result', {})
        entry = {
            'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'check': name, 'lane': o['lane'], 'kind': o['kind'],
            'claim': spec.get('source', ''),
            'method': {'e1kernel': 'jacobian.mult_kernel (dense nullspace mod p)', 'fixedD': f"fixed_d {spec.get('method','dense')} solver",
                       'e1': 'fixed_d.e1_solvable', 'a_in_j': 'jacobian.in_jacobian', 'hilbert': 'jacobian.jacobian_ring_basis',
                       'diag_group': 'symmetry.diagonal_group', 'poly_eq': 'exact polynomial identity'}.get(o['kind'], o['kind']),
            'prime': o.get('p'),
            'expected': o.get('expect'),
            'got': {k: v for k, v in res.items() if k not in ('basis', 'N', 'Q', 'generators', 'trace')},
            'agree': o.get('agree'),
            'ok': o.get('ok'),
            'error': o.get('error'),
            'remote_ts': o.get('ts'),
        }
        have[name] = entry
        new.append(entry)

with open(LOG, 'a') as f:
    for e in new:
        f.write(json.dumps(e, default=str) + '\n')

# summary
allv = list(have.values())
n_ok = sum(1 for e in allv if e['ok'])
n_err = sum(1 for e in allv if not e['ok'])
n_agree = sum(1 for e in allv if e['agree'] is True)
n_dis = [e for e in allv if e['agree'] is False]
n_noexp = sum(1 for e in allv if e['ok'] and e['agree'] is None)
print(f'verification log: {len(allv)} checks ({len(new)} new): ok {n_ok}, errors {n_err}, agree {n_agree}, DISAGREE {len(n_dis)}, no-expectation {n_noexp}')
for e in n_dis:
    print('  DISAGREE:', e['check'], 'expected', e['expected'], 'got', {k: e['got'].get(k) for k in (e['expected'] or {})})
for e in allv:
    if not e['ok']:
        print('  ERROR:', e['check'], e['error'])
by_lane = {}
for e in allv:
    by_lane.setdefault(e['lane'], []).append(e)
print('per lane (ok/agree/disagree/new-info/errors):')
for ln in sorted(by_lane):
    es = by_lane[ln]
    print(f"  {ln:30s} {sum(1 for e in es if e['ok']):3d} {sum(1 for e in es if e['agree'] is True):3d} {sum(1 for e in es if e['agree'] is False):3d} {sum(1 for e in es if e['ok'] and e['agree'] is None):3d} {sum(1 for e in es if not e['ok']):3d}")
if '--full' in sys.argv or True:
    print('--- all results (compact):')
    for e in sorted(allv, key=lambda e: e['check']):
        g = e['got']
        key = ({k: g.get(k) for k in ('dim', 'consistent', 'e1', 'in_jacobian', 'equal', 'status', 'size', 'order')} if isinstance(g, dict) else {'raw': g}) if g else {}
        key = {k: v for k, v in key.items() if v is not None}
        extra = {k: g.get(k) for k in ('kernel_dim', 'rank', 'nullity', 'ring_dim_par0', 'ring_dim_par1', 'seconds') if isinstance(g, dict) and g.get(k) is not None}
        flag = 'AGREE' if e['agree'] else ('DISAGREE' if e['agree'] is False else ('ERR' if not e['ok'] else 'info'))
        print(f"  {flag:8s} {e['check']:62s} {key} {extra}")
