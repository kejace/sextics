"""w2s_report.py -- summarize wave2-search.log.jsonl (+ mirrored restart files) into markdown tables for the report."""
import json, os, glob, sys
import numpy as np
LANE = '/tmp/hodge/results/dense-rational'
MIR = '/tmp/hodge/remote/dense-rational'
LOG = os.path.join(LANE, 'wave2-search.log.jsonl')

def load_log():
    out = []
    for line in open(LOG):
        try: out.append(json.loads(line))
        except Exception: pass
    return out

def restart_stats(name):
    fn = None
    for base in (MIR, LANE):
        p = os.path.join(base, f'{name}.restarts.jsonl')
        if os.path.exists(p): fn = p
    if fn is None: return None
    rs = [json.loads(l) for l in open(fn) if l.strip()]
    if not rs: return None
    k = 'rel' if 'rel' in rs[0] else 'f'
    v = np.array([r[k] for r in rs]); it = np.array([r.get('iters', 0) for r in rs]); sec = np.array([r.get('seconds', 0) for r in rs])
    q = lambda t: float(np.quantile(v, t))
    return dict(n=len(rs), key=k, best=float(v.min()), q10=q(0.1), median=q(0.5), q90=q(0.9), worst=float(v.max()),
                n_lt_1e6=int((v < 1e-6).sum()), n_lt_1e10=int((v < 1e-10).sum()), iters_med=float(np.median(it)), sec_med=float(np.median(sec)))

def main():
    recs = load_log()
    print('## Batches (from log)')
    print('| event | name | d | sector | family / mode | restarts | outer dim | inner / eqs | best | median | <1e-6 | <1e-10 | seconds |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for r in recs:
        ev = r.get('event')
        if ev == 'nvp-batch':
            print(f"| E2 varpro | {r['name']} | {r['d']} | pD={r['pD']} | {r['kind']} | {r['restarts']} | {r['outer_dim']} | {r['inner_dim']} / {r['rows_modF']} | {r['best_rel']:.2e} | {r['median_rel']:.2e} | {r['n_below_1e6']} | {r['n_below_1e10']} | {r['seconds']} |")
        elif ev == 'hessrelax-batch':
            print(f"| relaxation | {r['name']} | {r['d']} | {r['sector']} | {r['mode']} | {r['restarts']} | {r['dimP']} | {r['neq']} eqs | {r['best_f']:.2e} | {r['quantiles'][2]:.2e} | {r['n_below_1e6']} | {r['n_below_1e9']} (<1e-9) | {r['seconds']} |")
    print()
    print('## Live restart statistics (mirrored files)')
    print('| name | restarts | best | q10 | median | q90 | <1e-6 | <1e-10 | iters med | s/restart |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    names = sorted({os.path.basename(p).split('.restarts')[0] for base in (MIR, LANE) for p in glob.glob(os.path.join(base, '*.restarts.jsonl'))})
    for n in names:
        s = restart_stats(n)
        if s: print(f"| {n} | {s['n']} | {s['best']:.2e} | {s['q10']:.2e} | {s['median']:.2e} | {s['q90']:.2e} | {s['n_lt_1e6']} | {s['n_lt_1e10']} | {s['iters_med']:.0f} | {s['sec_med']:.0f} |")
    print()
    print('## Other events')
    for r in recs:
        ev = r.get('event')
        if ev in ('e1-space', 'kernel-structure', 'structured', 'nvp-confirm', 'hessrelax-confirm', 'hessrelax-promote', 'WITNESS', 'e1-exact-modp', 'family-built', 'J-complement', 'K_full'):
            print('-', json.dumps({k: v for k, v in r.items() if k not in ('ts', 'host', 'lane', 'engine')})[:400])

if __name__ == '__main__':
    main()
