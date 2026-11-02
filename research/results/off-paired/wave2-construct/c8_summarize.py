#!/usr/bin/env python3
"""Summarise the remote outputs after `rjob pull --lane off-paired`: served scan + Mordell-Weil msolve runs -> markdown rows."""
import json, os, glob, re
R = '/tmp/hodge/remote/off-paired/wave2-construct'
L = '/tmp/hodge/results/off-paired/wave2-construct'
rows = []
for base in (R, L):
    f = os.path.join(base, 'wave2-construct-served.jsonl')
    if os.path.exists(f):
        for line in open(f):
            r = json.loads(line)
            nontriv = sum(v[2] for v in r.get('served', {}).values()) if r.get('served') else None
            rows.append((r['name'], r.get('d'), r.get('char'), r.get('parity'), r.get('e1'), r.get('status'), nontriv, r.get('time')))
        break
seen = set()
print('| D | deg | char | parity | E1 | status | nontrivial served (all 18 Y-odd sectors) | s |')
print('|---|---|---|---|---|---|---|---|')
for row in rows:
    if row[0] in seen: continue
    seen.add(row[0]); print('| ' + ' | '.join(str(x) for x in row) + ' |')
for base in (R, L):
    for e in ('-1', '+1'):
        out = os.path.join(base, f'mw_sections_eps{e}.out'); log = os.path.join(base, f'mw_sections_eps{e}.log')
        if os.path.exists(out):
            txt = open(out).read()[:400].replace('\n', ' ')
            print(f'MW eps={e}: out={txt!r}')
            if os.path.exists(log):
                lg = open(log).read()
                m = re.findall(r'(dimension|#solutions|overall|elapsed)[^\n]*', lg)
                print('   log:', ' | '.join(m[-6:]))
