"""append one JSON line per check to the wave2-verify log (import or CLI)."""
import json, sys, time
LOG = '/tmp/hodge/results/quadratic-cubic-pullback/wave2-verify.log.jsonl'
def log(claim, verdict, evidence, **kw):
    rec = {'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'claim': claim, 'verdict': verdict, 'evidence': evidence}
    rec.update(kw)
    with open(LOG, 'a') as f: f.write(json.dumps(rec) + '\n')
    print('[vlog]', verdict, '|', claim, '|', evidence, flush=True)
if __name__ == '__main__':
    log(sys.argv[1], sys.argv[2], sys.argv[3])
