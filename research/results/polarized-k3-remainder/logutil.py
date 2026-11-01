import json, time, os, sys
sys.path.insert(0, '/home/kejace/devel/sextics-work')
LOG = '/tmp/hodge/results/polarized-k3-remainder/log.jsonl'
def log(D_expr, d, family, e1, consistent, rank, nullity, seconds, note=''):
    rec = {"D": D_expr, "d": d, "family": family, "e1": e1, "consistent": consistent,
           "rank": rank, "nullity": nullity, "seconds": round(seconds, 2), "note": note, "engine": 2,
           "time": time.strftime('%Y-%m-%dT%H:%M:%S')}
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec
