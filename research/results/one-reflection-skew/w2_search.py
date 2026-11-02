"""Wave-2 joint search driver (5-variable reduction, variable projection, multiprocessing).

Usage (run inside the lane dir, engine on PYTHONPATH):
  python3 w2_search.py --job NAME --d 7 --rho 0 --sigma 0 --family kernel --target A \
        --restarts 400 --procs 8 --method lm --maxiter 60 --seed0 1000 [--inst one-reflection-skew]
Families: kernel (E1 kernel minus f-multiples), jacobian (J_d minus f-multiples), cone:X0,X2,Y0 (kernel ∩ C[vars]),
          kernel+jac (same as kernel), full (all monomials of the sector; no E1 restriction, diagnostics only).
Outputs in the lane dir: NAME.restarts.jsonl (one line per restart), NAME.status.json (rewritten every batch),
NAME.cands.json (all restarts with rel < 1e-6, with monomial coefficient vectors); batches appended to
wave2-search.log.jsonl.  Exact confirmation is attempted for every restart with rel < 1e-7.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import sys, time, json, argparse, traceback
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
import numpy as np
import multiprocessing as mp
import w2lib as W
from w2lib import core

VP = None  # set in parent before forking


def _work(task):
    seed, init, method, maxiter = task
    try:
        out = VP.restart(seed, maxiter=maxiter, init=init, method=method)
    except Exception as ex:
        out = {'seed': seed, 'init': init, 'method': method, 'rel': float('nan'), 'error': f'{ex.__class__.__name__}: {ex}',
               'c': None, 'secs': 0, 'iters': 0, 'nfev': 0, 'phi': float('nan'), 'msg': 'error'}
    if isinstance(out, list):
        for o in out:
            o.pop('u', None)
        return out
    out.pop('u', None)
    return [out]


def build_family(prob, op, fam, d, pD, sD, T, exclude_fmult=False):
    info = {}
    if fam in ('kernel', 'kernel+jac') or fam.startswith('cone:'):
        Kb, kinfo = W.e1_kernel_exact(prob, d, pD, sD, T=T)
        info['kernel'] = {k: v for k, v in kinfo.items() if k != 'exact_basis'}
        B = Kb
        if fam.startswith('cone:'):
            vars_ = fam.split(':', 1)[1].split(',')
            keep = {core.VARS.index(v) for v in vars_}
            Bc = W.monomial_subfamily(prob, d, pD, sD, lambda m: all(m[i] == 0 for i in range(5) if i not in keep))
            B = W.intersect(Kb, Bc)
            info['cone_monomials'] = int(Bc.shape[1])
    elif fam == 'jacobian':
        B = W.jacobian_family(prob, d, pD, sD)
    elif fam == 'full':
        B = np.eye(op.nD)
    else:
        raise SystemExit(f'unknown family {fam}')
    # NOTE (wave-2 finding): do NOT remove the f-multiple subspace from the family.  A solution D* = D0 + f*G
    # generally has a nonzero component along f*S_{d-6}; restricting to the orthogonal complement loses it
    # (this is exactly what hid the needle in the first calibration runs; the engine's varpro has the same flaw).
    # In the strict 5-variable problem there are no spurious f-multiple solutions: (fG, fW) solves iff (G, W) does.
    if exclude_fmult and d >= 6 and fam != 'full':
        Bf = W.fmult_family(prob, d, pD, sD)
        if Bf.shape[1]:
            inter = W.intersect(B, Bf)
            if inter.shape[1]:
                B = W.complement_within(B, inter)
            info['fmult_removed'] = int(inter.shape[1])
    info['K'] = int(B.shape[1])
    return B, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inst', default='one-reflection-skew')
    ap.add_argument('--job', required=True)
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--rho', type=int, default=0)
    ap.add_argument('--sigma', default='0', help='0, 1 or none')
    ap.add_argument('--family', default='kernel')
    ap.add_argument('--target', default='A', help='A or Aplus or an expression')
    ap.add_argument('--restarts', type=int, default=100)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--method', default='lm', help='lm, lbfgs or both')
    ap.add_argument('--maxiter', type=int, default=0)
    ap.add_argument('--seed0', type=int, default=1000)
    ap.add_argument('--batch', type=int, default=16)
    ap.add_argument('--sparse-frac', type=float, default=0.5)
    ap.add_argument('--time-limit', type=float, default=0, help='seconds; stop submitting new restarts after this')
    ap.add_argument('--cand-tol', type=float, default=1e-6)
    ap.add_argument('--confirm-tol', type=float, default=1e-7)
    ap.add_argument('--defl-runs', type=int, default=10, help='sequential deflated LM runs per task (method defl)')
    ap.add_argument('--exclude-fmult', action='store_true', help='(not recommended) remove f*S_{d-6} from the family')
    ap.add_argument('--defl-perturb', type=float, default=0.0)
    a = ap.parse_args()
    W.DEFL_RUNS = a.defl_runs
    W.DEFL_PERTURB = a.defl_perturb
    global VP
    t_start = time.time()
    prob = W.Problem.from_instance(a.inst)
    sD = None if a.sigma == 'none' else int(a.sigma)
    if not prob.use_sigma:
        sD = None
    if a.target == 'A':
        T = prob.A
    elif a.target == 'Aplus':
        T = core.parse('X0**3*(Y0**3+Y1**3)/2')
    elif a.target == 'Aminus':
        T = core.parse('X0**3*(Y0**3-Y1**3)/2')
    else:
        T = core.parse(a.target)
    op = W.Operator5(prob, a.d, a.rho, sD, T, verbose=True)
    B, finfo = build_family(prob, op, a.family, a.d, a.rho, sD, T, exclude_fmult=a.exclude_fmult)
    K = B.shape[1]
    hdr = {'event': 'search-start', 'job': a.job, 'inst': a.inst, 'd': a.d, 'rho': a.rho, 'sigma': sD, 'family': a.family,
           'target': a.target, 'K': K, 'rows': op.m, 'cols': op.n, 'nD': op.nD, 'nnz': int(op.nnz), 'method': a.method,
           'restarts_planned': a.restarts, 'procs': a.procs, 'family_info': finfo}
    print(json.dumps(hdr), flush=True)
    W.log_event(hdr)
    if K == 0:
        W.log_event({'event': 'search-end', 'job': a.job, 'note': 'empty family', 'd': a.d, 'rho': a.rho, 'sigma': sD, 'family': a.family})
        print('empty family; nothing to search'); return
    VP = W.VarPro5(op, B)
    maxiter = a.maxiter or (150 if a.method in ('lm', 'lmk', 'defl') else 300)
    tasks = []
    rng = np.random.default_rng(a.seed0)
    for k in range(a.restarts):
        init = 'sparse' if rng.random() < a.sparse_frac else 'gauss'
        if a.method == 'both':
            meth = 'lm' if k % 2 == 0 else 'lbfgs'
            mi = 60 if meth == 'lm' else 300
        else:
            meth, mi = a.method, maxiter
        tasks.append((a.seed0 + k, init, meth, mi))
    if a.method == 'defl':
        # each task is a chain of defl-runs restarts; keep the total number of restarts as requested
        tasks = tasks[:max(1, a.restarts // a.defl_runs)]
    rest_path = f'{W.WD}/{a.job}.restarts.jsonl'
    cand_path = f'{W.WD}/{a.job}.cands.json'
    stat_path = f'{W.WD}/{a.job}.status.json'
    rels, cands, confirmed = [], [], []
    best = None
    n_done = 0
    batch_rels = []
    t_batch = time.time()

    def flush_status(final=False):
        arr = np.array([r for r in rels if np.isfinite(r)])
        hist = {}
        if arr.size:
            edges = [0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0, 10.0]
            cnt = np.histogram(arr, bins=edges)[0]
            hist = {f'<{edges[i+1]:g}': int(cnt[i]) for i in range(len(cnt))}
        st = {'job': a.job, 'd': a.d, 'rho': a.rho, 'sigma': sD, 'family': a.family, 'target': a.target, 'K': K,
              'method': a.method, 'done': n_done, 'planned': a.restarts, 'best_rel': (None if best is None else best['rel']),
              'median_rel': (float(np.median(arr)) if arr.size else None), 'hist': hist, 'n_cands': len(cands),
              'confirmed': confirmed, 'elapsed': round(time.time() - t_start, 1), 'final': final, 'ts': W.now()}
        with open(stat_path, 'w') as fh:
            json.dump(st, fh, indent=1)
        return st

    def handle(out):
        nonlocal best, n_done
        n_done += 1
        rel = out.get('rel', float('nan'))
        rels.append(rel); batch_rels.append(rel)
        rec = {k: out[k] for k in out if k != 'c'}
        rec.update({'job': a.job, 'd': a.d, 'rho': a.rho, 'sigma': sD, 'family': a.family, 'target': a.target})
        with open(rest_path, 'a') as fh:
            fh.write(json.dumps(rec, default=str) + '\n')
        if out.get('c') is not None and np.isfinite(rel) and (best is None or rel < best['rel']):
            best = out
        if out.get('c') is not None and np.isfinite(rel) and rel < a.cand_tol:
            c = np.asarray(out['c'], float); c = c / c[int(np.argmax(np.abs(c)))]
            cands.append({'seed': out['seed'], 'rel': rel, 'method': out.get('method'), 'c': c.tolist(),
                          'D_float': core.to_str({m: round(float(v), 10) for m, v in zip(op.Dm, c) if abs(v) > 1e-9})})
            with open(cand_path, 'w') as fh:
                json.dump({'Dm': [list(m) for m in op.Dm], 'cands': cands}, fh)
            if rel < a.confirm_tol and not confirmed:
                try_confirm(out)          # one exact confirmation per job is enough (further hits are only counted)

    def try_confirm(out):
        c = np.asarray(out['c'], float)
        try:
            zdim, svals = W.zero_set_dim(VP, c)
        except Exception:
            zdim, svals = None, None
        W.log_event({'event': 'candidate', 'job': a.job, 'inst': a.inst, 'd': a.d, 'rho': a.rho, 'sigma': sD, 'family': a.family,
                     'target': a.target, 'seed': out['seed'], 'rel': out['rel'], 'zero_set_dim_est': zdim, 'J_svals_tail': svals})

        def logger(rec):
            r2 = dict(rec); r2.update({'event': 'exact-confirm', 'job': a.job, 'inst': a.inst, 'd': a.d, 'rho': a.rho, 'sigma': sD,
                                        'family': a.family, 'target': a.target, 'seed': out['seed'], 'rel_num': out['rel']})
            W.log_event(r2)
            print(json.dumps(r2)[:400], flush=True)
        st, inf = W.snap_confirm(VP, prob, c, a.rho, sD, T, logger=logger)
        if st == 'solved':
            Dc = inf['D5']
            path = f'{W.WD}/witness-candidate.json' if a.inst == W.LANE else f'{W.WD}/{a.job}.witness.json'
            ok = W.save_witness_candidate(prob, inf['D6'], inf['N6'], {}, path=path,
                                          extra={'job': a.job, 'D5': core.to_str(Dc), 'N5': [core.to_str(n) for n in inf['N5']], 'target': a.target})
            confirmed.append({'D5': core.to_str(Dc), 'verified6': ok, 'seed': out['seed']})
            W.log_event({'event': 'WITNESS' if a.inst == W.LANE else 'witness-other-instance', 'job': a.job, 'inst': a.inst, 'd': a.d,
                         'D5': core.to_str(Dc), 'verified6': ok, 'target': a.target, 'path': path})
            print(('!!!!!!!! WITNESS FOUND (5-var -> 6-var verified=%s) !!!!!!!!' % ok) if a.inst == W.LANE else
                  f'(calibration instance {a.inst}: exact solution recovered, verified={ok})', flush=True)
        else:
            W.log_event({'event': 'unconfirmed-candidate', 'job': a.job, 'inst': a.inst, 'd': a.d, 'seed': out['seed'], 'rel': out['rel'],
                         'attempts': inf.get('attempts', [])[-6:]})

    ctx = mp.get_context('fork')
    with ctx.Pool(a.procs) as pool:
        it = pool.imap_unordered(_work, tasks)
        for outs in it:
            for out in outs:
                handle(out)
            if n_done % a.batch < len(outs) or n_done >= len(tasks):
                arr = np.array([r for r in batch_rels if np.isfinite(r)])
                rec = {'event': 'batch', 'job': a.job, 'inst': a.inst, 'd': a.d, 'rho': a.rho, 'sigma': sD, 'family': a.family,
                       'target': a.target, 'K': K, 'method': a.method, 'restarts': int(len(batch_rels)), 'done_total': n_done,
                       'best_rel_batch': (float(arr.min()) if arr.size else None), 'median_rel_batch': (float(np.median(arr)) if arr.size else None),
                       'best_rel_total': (None if best is None else best['rel']), 'batch_secs': round(time.time() - t_batch, 1),
                       'n_cands': len(cands)}
                W.log_event(rec)
                print(json.dumps(rec), flush=True)
                batch_rels = []; t_batch = time.time()
                flush_status()
            if a.time_limit and time.time() - t_start > a.time_limit:
                print('time limit reached; terminating pool', flush=True)
                pool.terminate()
                break
    st = flush_status(final=True)
    W.log_event({'event': 'search-end', 'job': a.job, 'inst': a.inst, 'd': a.d, 'rho': a.rho, 'sigma': sD, 'family': a.family,
                 'target': a.target, 'K': K, 'method': a.method, 'restarts': n_done, 'best_rel': st['best_rel'],
                 'median_rel': st['median_rel'], 'hist': st['hist'], 'n_cands': len(cands), 'confirmed': confirmed,
                 'elapsed': st['elapsed']})
    print(json.dumps(st), flush=True)


if __name__ == '__main__':
    main()
