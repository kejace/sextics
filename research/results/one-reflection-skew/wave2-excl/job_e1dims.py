"""Exact (over Q) E1-kernel dimensions K_d[chi] = {D in S_d[chi] : A D in J}, all 24 sectors, d in argv.
Writes wave2-excl/e1dims_d<d>.json (per-sector dims, with and without the F-gauge) and e1basis_d<d>.pkl (bases, gauge-fixed)."""
import sys, os, json, time, pickle, multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x

def work(args):
    d, chi = args
    t0 = time.time()
    r0 = w2x.e1_kernel(d, chi, gauge=False)
    r1 = w2x.e1_kernel(d, chi, gauge=True)
    basis = [P.to_dict() for P in r1['basis']]
    return {'d': d, 'chi': chi, 'nD': r0['nD'], 'dimK': r0['dimK'], 'dimJ': r0['dimJ_raw'], 'extra': r0['dimK'] - r0['dimJ_raw'],
            'dimK_gauged': r1['dimK'], 'dimJ_gauged': r1['dimJ'], 'gauge_removed': r1['gauge_removed'],
            'dimSJ_target': r0['dimSJ_target'], 'seconds': round(time.time() - t0, 1)}, basis

if __name__ == '__main__':
    ds = [int(a) for a in sys.argv[1:]]
    for d in ds:
        t0 = time.time()
        tasks = [(d, chi) for chi in w2x.ALLCHR]
        with mp.Pool(min(int(os.environ.get("NPROC", "8")), len(tasks))) as pool:
            res = pool.map(work, tasks)
        dims = [r for r, _ in res]
        bases = {tuple(r['chi']): b for r, b in res}
        tot = sum(r['dimK'] for r in dims); totJ = sum(r['dimJ'] for r in dims); totG = sum(r['dimK_gauged'] for r in dims)
        summary = {'d': d, 'dimK_total': tot, 'dimJ_total': totJ, 'extra_total': tot - totJ, 'dimK_gauged_total': totG,
                   'dimS_d': len(w2x.mons(d)), 'seconds': round(time.time() - t0, 1), 'sectors': dims, 'exact': True}
        with open(f'e1dims_d{d}.json', 'w') as fh: json.dump(summary, fh)
        with open(f'e1basis_d{d}.pkl', 'wb') as fh: pickle.dump(bases, fh)
        w2x.logexp(exp='e1-kernel-exact', family='E1 kernel {D: AD in J}, all sectors', d=d, prime=None, method='exact fmpz nullspace (python-flint), sector by sector',
                   result=f'dimK={tot} dimJ={totJ} extra={tot-totJ} gauged={totG}', rank=tot, seconds=summary['seconds'], job=f'excl-ors-e1dims',
                   note='exact over Q; no prime involved; per-sector dims in wave2-excl/e1dims_d%d.json' % d)
        print(json.dumps({k: v for k, v in summary.items() if k != 'sectors'}), flush=True)
