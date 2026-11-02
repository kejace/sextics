"""Append-as-you-go JSONL logging for the LEAN CAPACITY program."""
import json, time, sys
LOG = '/tmp/hodge/lean/capacity/log.jsonl'
def log(step, **kw):
    rec = {'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'step': step}
    rec.update(kw)
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec
if __name__ == '__main__':
    step = sys.argv[1]
    note = ' '.join(sys.argv[2:])
    print(log(step, note=note))
