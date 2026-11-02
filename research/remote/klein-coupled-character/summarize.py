import json, collections
rows = [json.loads(l) for l in open('/tmp/hodge/results/klein-coupled-character/log.jsonl')]
print('total log entries', len(rows))
byfam = collections.defaultdict(list)
for r in rows: byfam[r['family']].append(r)
for fam, rs in byfam.items():
    cons = sum(1 for r in rs if r.get('consistent') is True)
    inc = sum(1 for r in rs if r.get('consistent') is False)
    e1f = sum(1 for r in rs if r.get('e1') is False)
    secs = sum(r.get('seconds', 0) or 0 for r in rs)
    print(f'{fam:28s} n={len(rs):4d} consistent={cons} inconsistent={inc} e1-fail={e1f} degrees={sorted(set(r["d"] for r in rs))} time={secs:.0f}s')
print()
# detailed tables
def table(fam, cols=('D','d','k','e1','consistent','rank','nullity','rows','cols','seconds')):
    rs = byfam.get(fam, [])
    if not rs: return
    print(f'### {fam}')
    print('| ' + ' | '.join(cols) + ' |'); print('|' + '---|' * len(cols))
    for r in rs:
        print('| ' + ' | '.join(str(r.get(c, '')) for c in cols) + ' |')
    print()
for fam in ['polar', 'polar-times-coordinate', 'polar-product', 'hessian-block', 'companion', 'e1-kernel-rep', 'e1-kernel-random', 'e1-kernel-plus-J-random']:
    table(fam)
table('J_d-random', ('D','d','k','consistent','rank','nullity','rows','cols','seconds'))
for fam in ['scan-1param', 'scan-2param', 'served-targets', 'served-targets-cone']:
    rs = byfam.get(fam, [])
    if rs:
        print(f'### {fam}')
        for r in rs: print(f"- {r['D']} (d={r['d']}, k={r.get('k')}): consistent={r.get('consistent')} rank={r.get('rank')} [{r.get('seconds')}s] -- {r['note']}")
        print()
