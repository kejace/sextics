"""Full-E2 variable-projection search (Kaufman-LM) over the Jacobian-type w-space (or a saved E1 family) at degree d.
Usage: python3 w2_varpro.py --d 7 --sector even|odd|mixed --restarts N --seed S --workers W [--threads T]
       [--maxiter 60] [--batch 20] [--family jacobian|npz:<file>] [--init npy:<file>]
One JSON line per restart in w2_vp_d{d}_{sector}.jsonl; one batch line per --batch restarts in the lane log.
Any restart with relres < 1e-9 triggers rational reconstruction + exact re-solve through the engine; a verified
witness is saved to witness-candidate.json (and the batch line says WITNESS).
"""
import sys, os, time, json, argparse
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True); ap.add_argument('--sector', default='even')
ap.add_argument('--restarts', type=int, default=100); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--workers', type=int, default=4); ap.add_argument('--threads', type=int, default=1)
ap.add_argument('--maxiter', type=int, default=60); ap.add_argument('--batch', type=int, default=20)
ap.add_argument('--family', default='jacobian'); ap.add_argument('--init', default=None)
ap.add_argument('--tag', default='')
args = ap.parse_args()
os.environ['OMP_NUM_THREADS'] = str(args.threads); os.environ['OPENBLAS_NUM_THREADS'] = str(args.threads)
import numpy as np
import multiprocessing as mp
from w2lib import *

d = args.d; sector = args.sector
pD = {'even': 0, 'odd': 1, 'mixed': None}[sector]
famname = args.family if args.family == 'jacobian' else os.path.basename(args.family)
tag = f'vp_d{d}_{sector}_{famname}{args.tag}'
S = None


def build():
    global S
    F, A, q = load()
    sp = spaces(d, pD)
    if args.family == 'jacobian':
        fam = jacobian_family(F, A, d, pD, sp)
    else:
        z = np.load(args.family.split(':', 1)[1], allow_pickle=True)
        fam = Family(d, pD, z['Dmat'], [z[f'N{i}'] for i in range(NV)], list(z['names']), str(z['kind']))
    S = NormE2(F, A, q, fam, sp, verbose=True)
    S.Fexact, S.Aexact, S.qexact = F, A, q


def run_seed(seed):
    rng = np.random.default_rng(seed)
    u0 = rng.standard_normal(S.dimU)
    g0 = rng.standard_normal(S.ng) * 0.3 if S.ng else None
    if args.init and seed == args.seed:
        w0 = np.load(args.init.split(':', 1)[1])
        u0 = S.Pp.T @ w0; g0 = S.G.T @ w0 if S.ng else None
    res = S.lm(u0, g0, maxiter=args.maxiter)
    out = {'seed': int(seed), 'relres': res['relres'], 'f': res['f'], 'bnorm2': res['bnorm2'], 'iters': res['iters'],
           'rank': res['rank'], 'time': res['time'], 'gnorm': float(np.linalg.norm(res['g'])) if S.ng else 0.0}
    if res['relres'] < 1e-6:
        np.save(f'{LANE}/w2_{tag}_hit_{seed}.npy', res['w'])
    if res['relres'] < 1e-9:
        # rational reconstruction + exact confirmation
        Dc = res['Dc']
        tries = []
        for fr in S.reconstruct(Dc):
            D = S.D_exact(fr)
            status, data = confirm_exact(S.Fexact, S.Aexact, S.qexact, D)
            tries.append({'D': core.to_str(D)[:200], 'status': status})
            if status == 'solved':
                rec = save_candidate(d, D, data['N'], data['Q'], f'{LANE}/witness-candidate.json')
                out['WITNESS'] = True
                break
        out['exact_tries'] = tries
    return out


if __name__ == '__main__':
    build()
    chk = S.fd_check(seed=1)
    print(f'[check] gradient fd: {chk}', flush=True)
    seeds = list(range(args.seed, args.seed + args.restarts))
    t0 = time.time(); allres = []
    with mp.Pool(args.workers) as pool:
        for b0 in range(0, len(seeds), args.batch):
            chunk = seeds[b0:b0 + args.batch]
            outs = pool.map(run_seed, chunk, chunksize=1)
            for o in outs:
                append_json(f'{LANE}/w2_{tag}.jsonl', o)
            vals = np.array([o['relres'] for o in outs]); allres += list(vals)
            rec = {'batch': 'full-E2-varpro', 'family': f'{famname} (normalized E2, Kaufman-LM, u+g gauge chart)', 'd': d,
                   'sector': sector, 'nw': S.nw, 'gauge': S.ng, 'rows': S.R, 'nK': S.nK, 'nC': S.nC,
                   'restarts': len(chunk), 'seeds': [chunk[0], chunk[-1]], 'best_relres': float(vals.min()),
                   'median_relres': float(np.median(vals)), 'q10_relres': float(np.quantile(vals, 0.1)),
                   'mean_iters': float(np.mean([o['iters'] for o in outs])), 'mean_time': float(np.mean([o['time'] for o in outs])),
                   'ranks': sorted({o['rank'] for o in outs}), 'hits_below_1e-6': int((vals < 1e-6).sum()),
                   'witness': any(o.get('WITNESS') for o in outs), 'cum_restarts': len(allres), 'cum_best': float(min(allres)),
                   'seconds': round(time.time() - t0, 1)}
            if rec['witness']:
                rec['note'] = 'WITNESS'
            print(json.dumps(rec), flush=True)
            log_batch(rec)
