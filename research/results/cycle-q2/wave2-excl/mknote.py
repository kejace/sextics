"""Rewrite /tmp/hodge/notes/wave2-cycle-q2-excl.md from the log + commentary.md."""
import json, time, collections, os
LOG = '/tmp/hodge/results/cycle-q2/wave2-excl.log.jsonl'
rows = [json.loads(l) for l in open(LOG)] if os.path.exists(LOG) else []
witness = any(r.get('status') == 'SOLVED-EXACT' and r.get('instance', 'cycle-q2') == 'cycle-q2' for r in rows)
status = 'WITNESS FOUND' if witness else 'no witness'
comm = open('/tmp/hodge/results/cycle-q2/wave2-excl/commentary.md').read() if os.path.exists('/tmp/hodge/results/cycle-q2/wave2-excl/commentary.md') else ''
out = [f"# Wave 2 — lane `cycle-q2` — lens EXCLUSION / RIGOR\n",
       f"STATUS: {status} (in progress; {len(rows)} logged experiments; last update {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())})\n",
       "Log: `/tmp/hodge/results/cycle-q2/wave2-excl.log.jsonl`. Packet: `/tmp/hodge/results/cycle-q2/wave2-excl/` (independent solver `xlib.py`, `receipts/`, `farkas/`, `README.md`).\n",
       comm, "\n## Experiment log (auto-generated from the jsonl)\n"]
by = collections.OrderedDict()
for r in rows:
    by.setdefault(r.get('exp', '?'), []).append(r)
for exp, rs in by.items():
    out.append(f"\n### {exp} ({len(rs)} records)\n")
    for r in rs[-400:]:
        keys = [k for k in ('instance', 'family', 'D', 'd', 'status', 'ranks', 'primes', 'farkas', 'result', 'note', 'seconds') if k in r]
        out.append('- ' + '; '.join(f"{k}={str(r[k])[:90]}" for k in keys) + '\n')
open('/tmp/hodge/notes/wave2-cycle-q2-excl.md', 'w').write(''.join(out))
