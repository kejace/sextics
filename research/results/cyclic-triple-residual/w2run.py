"""Parallel search driver for wave-2 families (multiprocessing over seeds; workers single-threaded BLAS)."""
import os, sys, time, json, math
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import multiprocessing as mp
import numpy as np
from w2core import *

_FAM = None


def _init(fam):
    global _FAM
    _FAM = fam


def _task(args):
    seed, n, maxiter, relax, time_limit = args
    t0 = time.time()
    cands = _FAM.search(n_restarts=n, seed=seed, maxiter=maxiter, relax=relax, time_limit=time_limit)
    light = [{'rel': c['rel'], 'value': c['value'], 'bnorm2': c['bnorm2'], 'iters': c['iters'], 'spurious': c['spurious'],
              'div_rel': c['div_rel'], 'seed': seed, 'restart': c['restart'], 'relax': relax,
              'a': c['a'].tolist() if c['rel'] < 1e-4 or c is cands[0] else None} for c in cands]
    return {'seed': seed, 'n': len(cands), 'seconds': round(time.time() - t0, 1), 'cands': light}


def histogram(rels):
    bins = {}
    for r in rels:
        k = 'rel<1e-14' if r < 1e-14 else f'1e{int(math.floor(math.log10(max(r, 1e-300))))}'
        bins[k] = bins.get(k, 0) + 1
    return dict(sorted(bins.items()))


def search_batch(fam, label, F, A, q, job, n_restarts, procs=8, maxiter=400, relax=False, per_task=None,
                 promote_thresh=1e-8, time_limit=None, seed0=0, extra=None, best_path=None, log_path=W2LOG):
    """Run n_restarts random restarts split over `procs` workers; log one batch line; promote good candidates."""
    t0 = time.time()
    per_task = per_task or max(1, min(25, n_restarts // (4 * procs) or 1))
    tasks = []
    seed = seed0
    left = n_restarts
    while left > 0:
        n = min(per_task, left)
        tasks.append((seed, n, maxiter, relax, time_limit)); seed += 1; left -= n
    results = []
    if procs > 1:
        ctx = mp.get_context('fork')
        with ctx.Pool(procs, initializer=_init, initargs=(fam,)) as pool:
            for res in pool.imap_unordered(_task, tasks):
                results.append(res)
                if time_limit and time.time() - t0 > time_limit:
                    pool.terminate(); break
    else:
        _init(fam)
        for tk in tasks:
            results.append(_task(tk))
            if time_limit and time.time() - t0 > time_limit:
                break
    cands = [c for r in results for c in r['cands']]
    cands.sort(key=lambda c: (c['spurious'], c['rel']))
    rels = [c['rel'] for c in cands if not c['spurious']]
    nspur = sum(1 for c in cands if c['spurious'])
    best = cands[0] if cands else None
    entry = {'batch': 'search', 'job': job, 'label': label, 'relax': relax, 'd': fam.d, 'sector': list(fam.c),
             'group_order': fam.S.order, 'K': fam.K, 'dimP': fam.dimP, 'excluded_F_dim': fam.excluded_dim,
             'rows': fam.nrows, 'cols': fam.ncols, 'nC': fam.nC, 'nQ': fam.nQ, 'nN': fam.nN,
             'restarts': len(cands), 'restarts_requested': n_restarts, 'maxiter': maxiter, 'procs': procs,
             'best_rel': best['rel'] if best else None, 'best_abs': best['value'] if best else None,
             'median_rel': float(np.median(rels)) if rels else None, 'n_spurious': nspur,
             'n_rel_below_1e-6': sum(1 for r in rels if r < 1e-6), 'n_rel_below_1e-10': sum(1 for r in rels if r < 1e-10),
             'hist': histogram(rels), 'mean_iters': float(np.mean([c['iters'] for c in cands])) if cands else None,
             'seconds': round(time.time() - t0, 1), 'seed0': seed0, 'seeds_used': seed - seed0}
    if relax and fam._Pperp is not None:
        entry['rank_J2d'] = int(fam.rank_J2d); entry['dim_R2d'] = int(fam.nrows - fam.rank_J2d)
    if extra:
        entry.update(extra)
    w2log(entry, log_path)
    if best_path and best is not None and best.get('a') is not None:
        Df = fam.D_float(best['a'])
        with open(best_path, 'w') as f:
            json.dump({'label': label, 'relax': relax, 'rel': best['rel'], 'a': best['a'],
                       'D_float': {core.to_str({m: 1}): v for m, v in Df.items()},
                       'top5': [{'rel': c['rel'], 'seed': c['seed'], 'restart': c['restart'], 'a': c['a']} for c in cands[:5]]}, f)
    print(f'[batch {label}{" relax" if relax else ""}] restarts={len(cands)} best_rel={entry["best_rel"]:.3e} '
          f'median={entry["median_rel"]:.3e} spur={nspur} hist={entry["hist"]} {entry["seconds"]}s', flush=True)
    promoted = []
    if not relax:
        for c in cands[:5]:
            if c['spurious'] or c['rel'] > promote_thresh or c.get('a') is None:
                continue
            promoted.extend(promote(fam, {'a': np.array(c['a']), 'rel': c['rel']}, F, A, q, label, job))
    return entry, cands, promoted


def refine_relax_to_full(fam, cands, label, F, A, q, job, top=5, maxiter=600):
    """Take the best relaxation minimizers as initial points for the full search."""
    out = []
    for c in cands[:top]:
        if c.get('a') is None or c['spurious']:
            continue
        res = fam.search(n_restarts=1, init=np.array(c['a']), maxiter=maxiter)
        r = res[0]
        out.append({'from_relax_rel': c['rel'], 'full_rel': r['rel'], 'iters': r['iters'], 'spurious': r['spurious']})
        if r['rel'] < 1e-8 and not r['spurious']:
            promote(fam, r, F, A, q, label + '-refined', job)
    w2log({'batch': 'refine', 'job': job, 'label': label, 'd': fam.d, 'sector': list(fam.c), 'results': out})
    return out
