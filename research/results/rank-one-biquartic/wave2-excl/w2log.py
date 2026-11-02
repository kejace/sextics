"""Incremental logging for wave-2 EXCLUSION lens (lane rank-one-biquartic)."""
import json, time, os
LANE = 'rank-one-biquartic'
WORK = f'/tmp/hodge/results/{LANE}'
LOG = f'{WORK}/wave2-excl.log.jsonl'

def log(entry, path=None):
    """Append one JSON line (one experiment) to the wave-2 log; also mirror into the job's own receipt file."""
    entry = dict(entry)
    entry.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    entry.setdefault('lens', 'excl')
    entry.setdefault('engine', 2)
    with open(path or LOG, 'a') as fh:
        fh.write(json.dumps(entry, default=str) + '\n')
    return entry
