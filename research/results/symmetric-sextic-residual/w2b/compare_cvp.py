"""Compare the minima found by this agent's CharFamily varpro with the other agent's 'cvp' runs on the same sector."""
import json, glob, collections, sys
MIR = '/tmp/hodge/remote/symmetric-sextic-residual'
def minima(recs, key='rel'):
    c = collections.Counter(f"{r[key]:.3g}" for r in recs if key in r and not r.get('spurious'))
    return sorted(c.items(), key=lambda kv: float(kv[0]))
for tag, mine_glob, theirs in [('triv d8 p0', f'{MIR}/w2b/vp1_triv-8-0_w*.jsonl', f'{MIR}/w2vp_S3S3_triv_d8_p0_cvp_all.jsonl'),
                               ('sgnX d8 p1', f'{MIR}/w2b/vp3_sgnX-8-1_w*.jsonl', f'{MIR}/w2vp_S3S3_sgnX_d8_p1_cvp_all.jsonl')]:
    mine = [json.loads(l) for fn in glob.glob(mine_glob) for l in open(fn) if l.strip()]
    try: other = [json.loads(l) for l in open(theirs) if l.strip()]
    except FileNotFoundError: other = []
    print(tag, '| mine:', len(mine), 'restarts, minima', minima(mine)[:6], '| cvp:', len(other), 'restarts, minima', minima(other)[:6])
