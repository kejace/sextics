import json, time, sys
LOG='/tmp/hodge/results/cyclic-triple-residual/wave4/log.jsonl'
def log(**kw):
    kw.setdefault('t', time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    kw.setdefault('route','fano-lines'); kw.setdefault('wave',4)
    with open(LOG,'a') as f: f.write(json.dumps(kw, default=str)+'\n')
