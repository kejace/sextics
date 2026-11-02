"""Remote varpro driver.  python3 w2_varpro.py --job NAME --sectors triv:10:0,triv:10:1 [--restarts N] [--per-dim R]
   [--seed S] [--maxiter M] [--workers W] [--sparse-frac f] [--stop-after-hours H]
Per sector: build the admissible family (E1 kernel minus F-multiples), fork W workers, each runs its restarts and
appends one JSON line per restart to w2/<job>_<sector>_w<i>.jsonl; the parent summarises to w2/<job>.log.jsonl and
confirms (rational reconstruction + exact re-solve) every non-spurious candidate with rel < 1e-10."""
import os, sys, os, time, json, argparse, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2lib import *
import multiprocessing as mp

ap = argparse.ArgumentParser()
ap.add_argument('--job', required=True)
ap.add_argument('--sectors', required=True)
ap.add_argument('--restarts', type=int, default=200)
ap.add_argument('--per-dim', type=int, default=0, help='if >0: restarts = min(restarts, per_dim * dimP)')
ap.add_argument('--seed', type=int, default=1000)
ap.add_argument('--maxiter', type=int, default=400)
ap.add_argument('--workers', type=int, default=8)
ap.add_argument('--sparse-frac', type=float, default=0.3)
ap.add_argument('--stop-after-hours', type=float, default=100.0)
ap.add_argument('--struct', default=None, help='python file defining families (used by w2_struct)')
ap.add_argument('--dump-rel', type=float, default=1e-6, help='save the parameter vector a for restarts with rel below this')
args = ap.parse_args()
JOBLOG = f'{RES}/w2b/{args.job}.log.jsonl'
os.makedirs(f'{RES}/w2b', exist_ok=True)
F, A, q = load()
T0 = time.time()


def worker(fam, seeds, out_path, maxiter, sparse_frac):
    rng = np.random.default_rng(seeds[0] * 7919 + 17)
    with open(out_path, 'a') as f:
        for s in seeds:
            inits = ()
            if rng.random() < sparse_frac:
                # sparse init in family coordinates (1-3 random directions of the ambient basis), projected
                v = np.zeros(fam.K); idx = rng.choice(fam.K, size=int(rng.integers(1, 4)), replace=False)
                v[idx] = rng.choice([-1.0, 1.0], size=len(idx)); inits = (v,)
            def cb(c):
                rec = {'seed': c['seed'], 'rel': c['rel'], 'iters': c['iters'], 'nfev': c['nfev'],
                       'spurious': bool(c['spurious']), 'div_rel': c['div_rel'], 'sparse_init': bool(inits)}
                if c['rel'] < args.dump_rel:
                    rec['a'] = [float(x) for x in c['a']]
                f.write(json.dumps(rec) + '\n'); f.flush()
            fam.search(n_restarts=1, seed=s, maxiter=maxiter, inits=inits, callback=cb)


def summarize(files):
    recs = []
    for fn in files:
        if os.path.exists(fn):
            for line in open(fn):
                try: recs.append(json.loads(line))
                except Exception: pass
    if not recs:
        return {'n': 0}, []
    rels = np.array([r['rel'] for r in recs]); spur = np.array([r['spurious'] for r in recs])
    good = rels[~spur] if (~spur).any() else rels
    hist = {}
    for r in good:
        b = int(math.floor(math.log10(max(r, 1e-30))))
        hist[b] = hist.get(b, 0) + 1
    # distinct minima (rounded to 3 significant digits)
    minima = {}
    for r in good:
        key = f'{r:.3g}'
        minima[key] = minima.get(key, 0) + 1
    top = sorted(minima.items(), key=lambda kv: float(kv[0]))[:8]
    summ = {'n': len(recs), 'n_spurious': int(spur.sum()), 'best_rel': float(good.min()), 'median_rel': float(np.median(good)),
            'n_below_1e-4': int((good < 1e-4).sum()), 'n_below_1e-6': int((good < 1e-6).sum()),
            'n_below_1e-8': int((good < 1e-8).sum()), 'n_below_1e-10': int((good < 1e-10).sum()),
            'hist_log10': dict(sorted(hist.items())), 'top_minima': top,
            'mean_iters': float(np.mean([r['iters'] for r in recs]))}
    cands = [r for r in recs if (not r['spurious']) and r['rel'] < 1e-10 and 'a' in r]
    return summ, cands


for sec in args.sectors.split(','):
    if time.time() - T0 > args.stop_after_hours * 3600:
        break
    chi, d, pD = sec.split(':'); d = int(d); pD = int(pD)
    label = f'{chi}/d{d}/p{pD}'
    t0 = time.time()
    try:
        if chi == 'jac':
            fam, info = build_jac_family(F, A, q, d, pD, with_kernel=True, verbose=True)
            info.update(chi='jac', nD=info.get('K'), dimK=info.get('K'), dimJ=info.get('nW'), ker_beyond_J=info.get('n_kernel_extra'), nN=info.get('nT'))
        elif chi == 'jacJ':
            fam, info = build_jac_family(F, A, q, d, pD, with_kernel=False, verbose=True)
            info.update(chi='jacJ', nD=info.get('K'), dimK=info.get('K'), dimJ=info.get('nW'), ker_beyond_J=0, nN=info.get('nT'))
        else:
            fam, info = build_sector_family(F, A, q, chi, d, pD, verbose=True, label=label)
    except Exception as ex:
        wlog({'job': args.job, 'family': 'varpro-sector', 'chi': chi, 'd': d, 'pD': pD, 'error': repr(ex)}, JOBLOG)
        print('ERROR building', label, repr(ex), flush=True); continue
    base = {'job': args.job, 'family': 'varpro-sector', 'chi': chi, 'd': d, 'pD': pD, 'nD': info.get('nD'),
            'dimK': info.get('dimK'), 'dimJ': info.get('dimJ'), 'ker_beyond_J': info.get('ker_beyond_J'),
            'dimF': info.get('dimF'), 'dimP': info.get('dimP'), 'rows': info.get('rows'), 'nN': info.get('nN'),
            'nQ': info.get('nQ'), 'build_s': info.get('build_s'), 'maxiter': args.maxiter}
    print('family', label, info, flush=True)
    if fam is None or fam.dimP == 0:
        wlog(dict(base, restarts=0, note='admissible subspace empty (nothing to search)'), JOBLOG); continue
    n = args.restarts
    if args.per_dim > 0:
        n = min(n, args.per_dim * fam.dimP)
    n = max(n, 1)
    W = min(args.workers, n)
    seeds = list(range(args.seed, args.seed + n))
    files = [f'{RES}/w2b/{args.job}_{chi}-{d}-{pD}_w{i}.jsonl' for i in range(W)]
    ctx = mp.get_context('fork')
    procs = []
    for i in range(W):
        pr = ctx.Process(target=worker, args=(fam, seeds[i::W], files[i], args.maxiter, args.sparse_frac))
        pr.start(); procs.append(pr)
    # periodic checkpoints while running
    pollwait = 10
    while any(p.is_alive() for p in procs):
        time.sleep(pollwait); pollwait = min(120, pollwait + 10)
        summ, cands = summarize(files)
        if summ.get('n'):
            print(f'[{label}] progress n={summ["n"]} best={summ["best_rel"]:.3e} ({(time.time()-t0)/60:.1f} min)', flush=True)
    for p in procs: p.join()
    summ, cands = summarize(files)
    rec = dict(base, restarts=n, seeds=[seeds[0], seeds[-1]], seconds=round(time.time() - t0, 1), **summ)
    # Gauss-Newton polish of the best saved point (distinguishes a shrinking nonzero minimum from a true zero)
    try:
        saved = [r for fn in files if os.path.exists(fn) for r in map(json.loads, open(fn)) if 'a' in r and not r.get('spurious')]
        if saved:
            b = min(saved, key=lambda r: r['rel'])
            pa, px, prel, phist = polish(fam, np.array(b['a']), iters=10)
            rec['best_polished_rel'] = prel; rec['polish_hist'] = [f'{h:.2e}' for h in phist]
            if prel < 1e-10 and not any(c['rel'] < 1e-10 for c in cands):
                cands.append({'a': [float(v) for v in pa], 'rel': prel})
    except Exception as ex:
        rec['polish_error'] = repr(ex)[:200]
    # confirmation of candidates
    conf = []
    for c in cands[:20]:
        ok, res = try_confirm(fam, {'a': np.array(c['a']), 'rel': c['rel']}, dict(rec, event='confirm'))
        conf.append({'rel': c['rel'], 'ok': ok, 'res': res})
        if ok:
            rec['WITNESS'] = True
    rec['confirmations'] = conf
    wlog(rec, JOBLOG)
    print('SECTOR DONE', label, json.dumps({k: v for k, v in rec.items() if k not in ('confirmations',)}, default=str), flush=True)
print('ALL DONE', args.job, f'{(time.time()-T0)/3600:.2f} h', flush=True)
