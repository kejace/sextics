"""Wave-2 multi-process Koszul variable-projection search driver (klein-coupled-character).

python3 w2_search.py --d 9 --par 0 --chi all --family full --restarts 300 --procs 8 --seed 1000 --tag w2a
  family: full (= (J:A)_d sector), w (= J_d sector with lifts A w), hess (w-space, no Koszul unknowns)
Per sector: build once in the parent, fork --procs workers (1 BLAS thread each), each runs its share of
restarts with distinct seeds, streams one JSON line per restart to wave2-<tag>-d..p..c..-w<k>.jsonl.
Parent aggregates (histogram, distinct minima, rank stats), polishes the best non-spurious points,
attempts rational reconstruction + exact confirmation for rel < 1e-8, appends ONE batch line to
wave2-search.log.jsonl and regenerates the report (w2_report.py).
"""
import sys, os, json, time, argparse, math
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import numpy as np
import multiprocessing as mp
import w2fast as wf
from hodge import core
LANE_DIR = wf.LANE_DIR

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--par', type=str, default='0')          # 0, 1 or both
ap.add_argument('--chi', type=str, default='all')
ap.add_argument('--family', default='full', choices=['full', 'w', 'hess'])
ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--procs', type=int, default=8)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--maxiter', type=int, default=300)
ap.add_argument('--polish', type=int, default=3)
ap.add_argument('--tag', default='w2')
ap.add_argument('--confirm-below', type=float, default=1e-8)
args = ap.parse_args()

pars = [0, 1] if args.par == 'both' else [int(args.par)]
chis = list(range(7)) if args.chi == 'all' else [int(c) for c in args.chi.split(',')]
F, A, q = wf.load()
KVP = None   # global for forked workers


def worker(wid, nrest, seed, path):
    np.random.seed(seed)
    with open(path, 'a') as fh:
        def cb(c):
            rec = {k: v for k, v in c.items() if k != 'a'}
            if 'a' in c and c['rel'] < 1e-6:
                rec['a'] = [float(x) for x in c['a']]
            fh.write(json.dumps(rec) + '\n'); fh.flush()
        KVP.search(n_restarts=nrest, seed=seed, maxiter=args.maxiter, callback=cb, stop_rel=1e-13)


def hist(rels):
    edges = [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0, 10.0]
    h = {}
    for r in rels:
        for e in edges:
            if r < e:
                h[f'<{e:g}'] = h.get(f'<{e:g}', 0) + 1
                break
        else:
            h['>=10'] = h.get('>=10', 0) + 1
    return h


for par in pars:
    for chi in chis:
        t0 = time.time()
        sec = wf.Sector(F, A, args.d, par, chi, verbose=True)
        if sec.nD == 0:
            wf.batch_log({'lens': 'search', 'family': args.family, 'd': args.d, 'parD': par, 'chi': chi, 'note': 'empty D sector'})
            continue
        if args.family == 'full':
            Bc, Bn, info = sec.e1_space()
        else:
            Bc, Bn, info = sec.w_space()
        if info['K'] == 0:
            wf.batch_log({'lens': 'search', 'family': args.family, 'd': args.d, 'parD': par, 'chi': chi, 'K': 0,
                          'note': 'no E1-admissible D in sector', 'dims': sec.dims()})
            print('empty family', par, chi, flush=True)
            continue
        sec.build(verbose=True)
        KVP = wf.KVP(sec, Bc, Bn, use_koszul=(args.family != 'hess'))
        if KVP.dimP == 0:
            wf.batch_log({'lens': 'search', 'family': args.family, 'd': args.d, 'parD': par, 'chi': chi, 'K': KVP.K,
                          'excluded_F_multiples': KVP.excluded_dim, 'note': 'family consists of F-multiples only', 'dims': sec.dims()})
            print('only F multiples', par, chi, flush=True)
            continue
        base = f'wave2-{args.tag}-d{args.d}p{par}c{chi}'
        paths = [os.path.join(LANE_DIR, f'{base}-w{w}.jsonl') for w in range(args.procs)]
        share = [args.restarts // args.procs + (1 if w < args.restarts % args.procs else 0) for w in range(args.procs)]
        t1 = time.time()
        ctx = mp.get_context('fork')
        procs = []
        for w in range(args.procs):
            if share[w] == 0:
                continue
            p = ctx.Process(target=worker, args=(w, share[w], args.seed + 7919 * w + 1000003 * chi + 31 * par, paths[w]))
            p.start(); procs.append(p)
        for p in procs:
            p.join()
        t_search = time.time() - t1
        # aggregate
        recs = []
        for pth in paths:
            if os.path.exists(pth):
                for line in open(pth):
                    try:
                        recs.append(json.loads(line))
                    except Exception:
                        pass
        rels = [r['rel'] for r in recs if not r['spurious']]
        nsp = sum(1 for r in recs if r['spurious'])
        best = sorted([r for r in recs if not r['spurious']], key=lambda r: r['rel'])[:max(args.polish, 5)]
        distinct = sorted({round(r['rel'], 4) for r in recs if not r['spurious']})[:12]
        # polish the best ones (need a vectors: re-run from stored a if present, else from their seed/restart)
        polished = []
        for r in best[:args.polish]:
            if 'a' in r:
                a0 = np.array(r['a'])
            else:
                # regenerate the start point deterministically and re-minimize
                rng = np.random.default_rng(r['seed'])
                u0 = None
                for it in range(r['restart'] + 1):
                    u0 = rng.standard_normal(KVP.dimP)
                u0 /= np.linalg.norm(u0)
                from scipy.optimize import minimize
                res = minimize(KVP.value_and_grad, u0, jac=True, method='L-BFGS-B',
                               options={'maxiter': args.maxiter, 'gtol': 1e-10, 'ftol': 1e-16, 'maxcor': 20})
                a0 = KVP.P @ res.x; a0 /= np.linalg.norm(a0)
            a1, f1, rel1, nit = KVP.polish(a0, maxiter=2000)
            f_, bb, xP, rr, rank = KVP.value(a1)
            entry = {'rel_before': r['rel'], 'rel_after': rel1, 'f_after': f1, 'polish_iters': nit, 'rankT': rank,
                     'nullityT': sec.nP - rank if KVP.use_koszul else None}
            if rel1 < args.confirm_below:
                c = KVP.Bc @ a1
                cands = wf.reconstruct_D(sec, c)
                entry['reconstructions'] = len(cands)
                import w2lib
                for D in cands[:3]:
                    st, data = w2lib.confirm_exact(F, A, D)
                    entry.setdefault('confirm', []).append({'status': st, 'rank': data.get('rank'), 'nullity': data.get('nullity'),
                                                            'D': core.to_str(D)[:300]})
                    if st == 'solved':
                        wf.save_witness_candidate(F, A, D, data['N'], data['Q'])
                        print('WITNESS FOUND', core.to_str(D), flush=True)
                        wf.batch_log({'lens': 'search', 'WITNESS': True, 'd': args.d, 'parD': par, 'chi': chi, 'D': core.to_str(D)})
                np.save(os.path.join(LANE_DIR, f'{base}-best-c.npy'), c)
            polished.append(entry)
        entry = {'lens': 'search', 'family': args.family, 'd': args.d, 'parD': par, 'chi': chi, 'tag': args.tag,
                 'K': KVP.K, 'dimP': KVP.dimP, 'excluded_F_multiples': KVP.excluded_dim, 'family_info': {k: v for k, v in info.items() if k != 'sv_gap'},
                 'dims': sec.dims(), 'restarts': len(recs), 'spurious': nsp, 'seed': args.seed, 'maxiter': args.maxiter,
                 'best_rel': min(rels) if rels else None, 'best_f': min((r['f'] for r in recs if not r['spurious']), default=None),
                 'hist_rel': hist(rels), 'distinct_minima': distinct, 'median_iters': float(np.median([r['iters'] for r in recs])) if recs else None,
                 'mean_restart_seconds': float(np.mean([r['seconds'] for r in recs])) if recs else None,
                 'rankT_at_best': best[0]['rankT'] if best else None, 'nP': sec.nP, 'nQ': sec.nQ, 'rows': sec.nrows,
                 'polished': polished, 'search_seconds': round(t_search, 1), 'total_seconds': round(time.time() - t0, 1),
                 'procs': args.procs, 'loadavg': open('/proc/loadavg').read().split()[0]}
        wf.batch_log(entry)
        print(json.dumps({k: v for k, v in entry.items() if k not in ('polished', 'family_info', 'dims')}), flush=True)
        try:
            os.system('python3 /tmp/hodge/results/klein-coupled-character/w2_report.py > /dev/null 2>&1')
        except Exception:
            pass
print('SEARCH_DONE', flush=True)
