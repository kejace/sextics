import json, time, os
LOG = '/tmp/hodge/results/polarized-k3-remainder/wave2-excl.log.jsonl'
def log(**rec):
    rec.setdefault('time', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('lens', 'excl'); rec.setdefault('lane', 'polarized-k3-remainder'); rec.setdefault('engine', 2)
    with open(LOG, 'a') as fh:
        fh.write(json.dumps(rec, default=str) + '\n')
    print('LOGGED', json.dumps(rec, default=str)[:400], flush=True)
