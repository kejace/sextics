"""Logging helpers for wave-2 construction lens (lane polarized-k3-remainder)."""
import json, time, os
LANE = '/tmp/hodge/results/polarized-k3-remainder'
LOG = os.path.join(LANE, 'wave2-construct.log.jsonl')
def log(**kw):
    kw.setdefault('time', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(LOG, 'a') as f:
        f.write(json.dumps(kw, default=str) + '\n')
    return kw
