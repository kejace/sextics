#!/usr/bin/env python3
"""Wave-2 joint (D,N,Q) search driver: Koszul variable projection over E1-admissible families.

Usage (run inside the lane dir, engine on PYTHONPATH):
  python3 w2_varpro.py --d 8 --chi all --restarts 200 --procs 8 --seed 1000 --tag d8all
  python3 w2_varpro.py --d 7 --parity 0 --restarts 500 --procs 8 --tag d7even      # character-mixed, Y-even D
  python3 w2_varpro.py --d 7 --mixed --restarts 200 --procs 8 --tag d7mixed         # mixed Y-parity too
  python3 w2_varpro.py --d 10 --chi 0000 --vars X0,Y0,X1,Y1 --mode basis ...        # cone family (variable subset)

Per (degree, sector): builds the exact E1 family, the Koszul varpro system, runs `restarts` random
restarts spread over `procs` worker processes (distinct seeds), logs ONE JSON line per sector batch to
wave2-search.log.jsonl, saves the best parameter vectors to varpro-<tag>-d<d>-c<chi>.json, and for every
candidate with rel < --recon attempts rational reconstruction + exact re-solve (witness-candidate.json
on success; 'WITNESS FOUND' printed).
"""
import sys, os, json, time, argparse, math
sys.path.insert(0, '/tmp/hodge/results/quadric-contact-residual')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1'); os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import multiprocessing as mp
import w2lib as W

VP = None   # shared (fork) family system


def _worker(args):
    seed, n, maxiter, stop_rel = args
    out = VP.search(n_restarts=n, seed=seed, maxiter=maxiter, stop_rel=stop_rel)
    return [{'rel': c['rel'], 'f': c['f'], 'bnorm2': c['bnorm2'], 'a': c['a'].tolist(), 'spurious': c['spurious'],
             'div_rel': c['div_rel'], 'iters': c['iters'], 'seed': seed, 'restart': c['restart']} for c in out]


def run_sector(d, c, parity, mixed, mode, vars_mask, restarts, procs, seed, maxiter, tag, recon, keep, time_limit, per_dim=6, min_restarts=16, sym=False):
    global VP
    t0 = time.time()
    chi_s = 'mixed' if (c is None and parity is None) else (f'par{parity}' if c is None else W.chi_str(c))
    base = {'task': 'varpro', 'd': d, 'chi': chi_s, 'mode': mode + ('+sym' if sym else ''), 'tag': tag,
            'vars': None if vars_mask is None else ''.join(v for v, m in zip(W.core.VARS, vars_mask) if m)}
    try:
        members, info, fam = W.build_family(d, c, parity=parity, vars_mask=vars_mask, mode=mode, verbose=True, sym=sym)
    except Exception as ex:
        W.log_line({**base, 'status': 'family-error', 'error': str(ex)[:300]})
        return
    rec = {**base, 'dimS': info['dimS'], 'dim_kernel': info['dim_kernel'], 'dimJ': info['dimJ'],
           'n_extra': info['n_extra'], 'dim_gauge': info['dim_gauge'], 'n_members': info['n_members'],
           'max_lift_res': info['max_lift_res'], 't_family': info['t_family']}
    if info['n_members'] == 0:
        W.log_line({**rec, 'status': 'empty-family'})
        print('empty family', base, flush=True)
        return
    VP = W.KoszulVarPro(members, d, c, parity=parity, verbose=True)
    rec.update({'rows': VP.rows, 'nK': VP.nK, 'nQ': VP.nQ, 'K': VP.K, 'dimP': VP.dimP, 'excluded': VP.excluded,
                'mstack_nnz': int(VP.Mstack.nnz) if VP.Mstack is not None else 0, 't_build': round(VP.build_time, 1)})
    if VP.dimP == 0:
        W.log_line({**rec, 'status': 'only-gauge'})
        print('only gauge directions', base, flush=True)
        return
    # timing probe
    tp = time.time(); VP.value_and_grad(np.ones(VP.dimP)); rec['t_eval'] = round(time.time() - tp, 4)
    # adaptive number of restarts: per_dim * dimP clamped to [min_restarts, restarts]
    n_restarts = restarts if per_dim <= 0 else max(min_restarts, min(restarts, per_dim * VP.dimP))
    rec['restarts_planned'] = n_restarts
    procs = max(1, min(procs, n_restarts))
    jobs = [(seed + 7919 * i, 1, maxiter, recon * 1e-2) for i in range(n_restarts)]
    cands = []
    timed_out = False
    if procs == 1:
        for jb in jobs:
            cands.extend(_worker(jb))
            if time_limit and time.time() - t0 > time_limit:
                timed_out = True; break
    else:
        ctx = mp.get_context('fork')
        pool = ctx.Pool(processes=procs)
        try:
            for res in pool.imap_unordered(_worker, jobs):
                cands.extend(res)
                if time_limit and time.time() - t0 > time_limit:
                    print(f'time limit reached after {len(cands)} restarts; terminating pool', flush=True)
                    timed_out = True
                    break
            pool.close()
        finally:
            pool.terminate(); pool.join()
    rec['timed_out'] = timed_out
    cands.sort(key=lambda c: (c['spurious'], c['rel']))
    rels = np.array([c['rel'] for c in cands if not c['spurious']]) if cands else np.array([])
    rec.update({'restarts_done': len(cands), 'n_spurious': int(sum(1 for c in cands if c['spurious'])),
                'best_rel': float(rels.min()) if rels.size else None,
                'median_rel': float(np.median(rels)) if rels.size else None,
                'n_rel_lt_1e-6': int((rels < 1e-6).sum()), 'n_rel_lt_1e-10': int((rels < 1e-10).sum()),
                'hist_log10': {str(k): int(v) for k, v in zip(*np.unique(np.floor(np.log10(np.maximum(rels, 1e-300))).astype(int), return_counts=True))} if rels.size else {},
                'mean_iters': float(np.mean([c['iters'] for c in cands])) if cands else None,
                'seconds': round(time.time() - t0, 1), 'status': 'done'})
    # save best candidates
    outp = f'{W.RES}/varpro-{tag}-d{d}-c{chi_s}.json'
    best = cands[:keep]
    with open(outp, 'w') as f:
        json.dump({'rec': rec, 'members': [m[2] for m in members], 'best': best}, f)
    rec['out'] = os.path.basename(outp)
    # reconstruction attempts
    recon_att = []
    for cnd in cands:
        if cnd['spurious'] or cnd['rel'] > recon:
            continue
        att = W.try_reconstruct_and_confirm(VP, np.array(cnd['a']), verbose=True, tag=f'{tag}-d{d}-c{chi_s}-s{cnd["seed"]}r{cnd["restart"]}')
        recon_att.append({'rel': cnd['rel'], 'attempts': att})
        if any(a.get('witness_path') for a in att):
            rec['WITNESS'] = True
            break
    rec['recon'] = recon_att
    W.log_line(rec)
    print(json.dumps({k: v for k, v in rec.items() if k not in ('hist_log10', 'recon')}), flush=True)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--chi', default=None, help="character string like 0000, or 'all' (every sector of the degree), or comma list")
    ap.add_argument('--parity', type=int, default=None)
    ap.add_argument('--mixed', action='store_true')
    ap.add_argument('--mode', default='gens', choices=['gens', 'basis'])
    ap.add_argument('--vars', default=None, help='comma list of variables allowed in D (cone families)')
    ap.add_argument('--restarts', type=int, default=100)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--maxiter', type=int, default=300)
    ap.add_argument('--tag', default='w2')
    ap.add_argument('--recon', type=float, default=1e-10)
    ap.add_argument('--keep', type=int, default=10)
    ap.add_argument('--time-limit', type=float, default=0.0, help='seconds per sector batch (0 = none)')
    ap.add_argument('--per-dim', type=int, default=6, help='restarts = per_dim * dimP clamped to [min-restarts, restarts]; 0 = always --restarts')
    ap.add_argument('--min-restarts', type=int, default=16)
    ap.add_argument('--skip-parity', type=int, default=None, help='skip sectors of this Y-parity')
    ap.add_argument('--sym', action='store_true', help='restrict to <s0,r1>-invariant denominators (sector 0000 only is meaningful)')
    a = ap.parse_args()
    vars_mask = None
    if a.vars:
        allowed = set(a.vars.split(','))
        vars_mask = [v in allowed for v in W.core.VARS]
    if a.mixed:
        sectors = [(None, None)]
    elif a.chi is None and a.parity is not None:
        sectors = [(None, a.parity)]
    elif a.chi == 'all':
        sectors = [(c, None) for c in W.chars_in_degree(a.d) if a.skip_parity is None or c[3] != a.skip_parity]
    else:
        sectors = [(W.parse_chi(s), None) for s in a.chi.split(',')]
    print(f'{len(sectors)} sector(s) to run', flush=True)
    for c, par in sectors:
        try:
            run_sector(a.d, c, par, a.mixed, a.mode, vars_mask, a.restarts, a.procs, a.seed, a.maxiter, a.tag, a.recon,
                       a.keep, a.time_limit, a.per_dim, a.min_restarts, a.sym)
        except Exception as ex:
            import traceback; traceback.print_exc()
            W.log_line({'task': 'varpro', 'd': a.d, 'chi': 'mixed' if c is None and par is None else (W.chi_str(c) if c is not None else f'par{par}'),
                        'tag': a.tag, 'status': 'family-error', 'error': str(ex)[:300]})


if __name__ == '__main__':
    main()
