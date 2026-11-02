"""Recover sweep results from the remote job's printed output (w2_sweep.out in the mirror) when the shared remote log
was overwritten by another agent's rsync. Adds minimal 'sweep' entries for candidates missing from the canonical log."""
import json, re, os, time
WD = '/tmp/hodge/results/cyclic-triple-residual'
out = '/tmp/hodge/remote/cyclic-triple-residual/w2_sweep.out'
log = f'{WD}/wave2-search.log.jsonl'
rows = [json.loads(l) for l in open(log) if l.strip()]
have = {r['name'] for r in rows if r.get('batch') == 'sweep'}
pat = re.compile(r'^d=\s*(\d+) e1=(\S+) cons=(\S+) rank=(\S+)/(\S+) (\S+)s (\S+)\s+(.*?)\s*(\*\*\* WITNESS \*\*\*)?\s*$')
added = 0
if os.path.exists(out):
    for line in open(out):
        m = pat.match(line.rstrip('\n'))
        if not m: continue
        d, e1, cons, rank, raug, secs, fam, name, wit = m.groups()
        name = name.strip()
        if name in have: continue
        def conv(x):
            return None if x == 'None' else (True if x == 'True' else False if x == 'False' else int(x))
        rec = {'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'engine': 2, 'lane': 'cyclic-triple-residual', 'batch': 'sweep', 'job': 'w2sweep',
               'family': fam, 'name': name, 'd': int(d), 'terms': None, 'e1': conv(e1), 'consistent': conv(cons), 'rank': conv(rank),
               'rank_aug': conv(raug), 'seconds': float(secs), 'recovered_from': 'w2_sweep.out (shared remote log overwritten by rsync)',
               'WITNESS': bool(wit)}
        open(log, 'a').write(json.dumps(rec) + '\n'); have.add(name); added += 1
print(f'recovered {added} sweep entries')
