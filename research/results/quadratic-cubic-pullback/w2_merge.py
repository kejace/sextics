#!/usr/bin/env python3
"""Merge remote per-job JSONL outputs (mirror /tmp/hodge/remote/<lane>/w2_*.jsonl) into the lane's wave2-search.log.jsonl,
stripping bulky vectors ('v') from search records (kept in the mirror files).  Idempotent (dedupes by content)."""
import json, glob, os, sys
LANE = 'quadratic-cubic-pullback'
MIR = f'/tmp/hodge/remote/{LANE}'
LOG = f'/tmp/hodge/results/{LANE}/wave2-search.log.jsonl'
have = set()
if os.path.exists(LOG):
    have = set(l.strip() for l in open(LOG) if l.strip())
added = 0
with open(LOG, 'a') as out:
    for path in sorted(glob.glob(f'{MIR}/w2_*.jsonl')) + sorted(glob.glob(f'/tmp/hodge/results/{LANE}/w2_*.jsonl')):
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except Exception:
                continue
            r.pop('v', None)
            r.setdefault('src', os.path.basename(path))
            s = json.dumps(r, sort_keys=True)
            if s in have:
                continue
            have.add(s); out.write(s + '\n'); added += 1
print(f'merged {added} new records into {LOG} (total {len(have)})')
