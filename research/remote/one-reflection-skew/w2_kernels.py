"""E1-admissible denominator spaces (5-var reduction) per degree and (rho,sigma) sector; saves float bases."""
import sys, json, pickle, time
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
import w2skew as W
ds = [int(x) for x in sys.argv[1:]] or list(range(4, 13))
for d in ds:
    for rD in (0, 1):
        for sD in (0, 1):
            k = W.e1_kernel5(d, rD, sD, want_float=True)
            info = {kk: v for kk, v in k.items() if kk not in ('mons', 'Kbasis', 'Jbasis', 'Cbasis')}
            print(json.dumps(info), flush=True)
            with open(f'/tmp/hodge/results/one-reflection-skew/w2ker5_d{d}_r{rD}s{sD}.pkl', 'wb') as fh:
                pickle.dump(k, fh)
            W.logbatch(batch='e1-kernel5', family='E1 kernel (5-var)', degree=d, sector=f'rho{rD}sig{sD}', restarts=0,
                       best_rel=None, seconds=k['seconds'],
                       note=f"nD={k['nD']} kerD={k['kerD']} dimJ_d={k['dimJ']} nonJ={k['nonJ']} matrix {k['rows']}x{k['cols']} sv_gap={k.get('sv_gap')}")
    W.rewrite_report()
print('DONE')
