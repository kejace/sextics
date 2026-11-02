#!/usr/bin/env python3
"""Harvest finished wave-2 msolve outputs (leading ideal or reduced GB) into wave-3 receipts.
Usage: wave3_harvest.py TAG  (reads lines "<input.ms> <output>" from stdin). Appends one JSON line per pair to
/tmp/hodge/results/_closure/wave3-receipts.jsonl and <dir>/wave3/receipts.jsonl; prints a compact summary."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wave3_verdict import verdict

REC = '/tmp/hodge/results/_closure/wave3-receipts.jsonl'
tag = sys.argv[1]
for line in sys.stdin:
    parts = line.split()
    if len(parts) != 2:
        continue
    ms, out = parts
    ok = os.path.exists(out) and os.path.getsize(out) > 0
    v = verdict(ms, out, 0 if ok else 1)
    lane = ms.split('/tmp/hodge/results/')[-1].split('/')[0] if '/tmp/hodge/results/' in ms else '?'
    rec = {'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'tag': tag, 'lane': lane, 'file': ms, 'lead': out,
           'engine': 'wave-2 msolve run (harvested output)', 'output_mtime': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(os.path.getmtime(out))) if ok else None,
           'result': v}
    with open(REC, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    d = os.path.join(os.path.dirname(ms), 'wave3'); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, 'receipts.jsonl'), 'a') as f:
        f.write(json.dumps(rec) + '\n')
    print(lane, os.path.basename(ms), v.get('n_vars'), v.get('n_polys'), v.get('output_kind'), v.get('n_lead'), v.get('max_pure_power'), v['verdict'])
