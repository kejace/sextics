"""Steps 2-5: sector-restricted variable-projection search jobs.
Usage:
  python3 w2s_search.py --mode {hessJ,wspace,family,cone} --d D [--keys all|k1,k2,...] [--restarts N] [--cpus C]
                        [--seed S] [--maxiter M] [--job NAME] [--min-dim 1] [--vars X0,X1,Y0,Y1] [--warm]
hessJ : Hessian relaxation (necessary): quad(w) in J_{2d+6} over the w-space (very cheap)
wspace: full E2 on the Jacobian part (N = A w/q + Koszul, Q free)
family: full E2 on the whole E1 kernel K_d^key (N, Q free in the sector)
cone  : family restricted to polynomials in the given variable subset (vertex a line/plane)
"""
import sys, time, json, argparse, os
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import numpy as np
import w2s_lib as L
from hodge import core

ap = argparse.ArgumentParser()
ap.add_argument('--mode', required=True); ap.add_argument('--d', type=int, required=True)
ap.add_argument('--keys', default='all'); ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--cpus', type=int, default=1); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--maxiter', type=int, default=300); ap.add_argument('--job', default=None)
ap.add_argument('--min-dim', type=int, default=1); ap.add_argument('--vars', default=None)
ap.add_argument('--warm', action='store_true', help='wspace: warm-start from hessJ near-zeros')
ap.add_argument('--chunk', type=int, default=None); ap.add_argument('--ls', default='auto')
ap.add_argument('--max-rows', type=int, default=10**9)
args = ap.parse_args()
job = args.job or f'search-{args.mode}-d{args.d}-s{args.seed}'
I = L.Inst()
d = args.d
keys = I.keys_of_degree(d, 0) if args.keys == 'all' else [tuple(int(x) for x in k.split('-')) for k in args.keys.split(',')]
t00 = time.time()
for key in keys:
    kd = L.load_kernel(I, d, key, parity=0)
    tag = f'{args.mode}-d{d}-{L.keystr(key)}'
    extra = {'kernel_dim': kd['kernel_dim'], 'jac_dim': kd['jac_dim'], 'nonjac_dim': kd['nonjac_dim']}
    try:
        if args.mode in ('hessJ', 'wspace'):
            wb = L.wspace_basis(I, d, key)
            if len(wb) < args.min_dim:
                L.log_line({'step': 'skip', 'job': job, 'tag': tag, 'reason': f'wdim={len(wb)}<{args.min_dim}'}); continue
            if args.mode == 'hessJ':
                svp = L.HessRelax(I, key, wb, verbose=True)
            else:
                svp = L.SectorVarPro(I, key, mode='wspace', wbasis=wb, verbose=True, ls=args.ls)
        elif args.mode in ('family', 'cone'):
            basis = kd['basis']
            if args.mode == 'cone':
                vs = [int(v[1]) + (3 if v[0] == 'Y' else 0) for v in args.vars.split(',')]
                basis = L.cone_subbasis(kd, vs)
                tag += '-' + args.vars.replace(',', '')
                extra['vars'] = args.vars
            if len(basis) < args.min_dim:
                L.log_line({'step': 'skip', 'job': job, 'tag': tag, 'reason': f'dim={len(basis)}<{args.min_dim}'}); continue
            svp = L.SectorVarPro(I, key, mode='family', basis=basis, verbose=True, ls=args.ls)
        else:
            raise SystemExit('unknown mode')
        if svp.rows > args.max_rows:
            L.log_line({'step': 'skip', 'job': job, 'tag': tag, 'reason': f'rows={svp.rows}>{args.max_rows}'}); continue
        if svp.dimP == 0:
            L.log_line({'step': 'skip', 'job': job, 'tag': tag, 'reason': 'dimP=0 (all directions spurious)', **extra}); continue
        L.run_search(svp, tag, job, args.restarts, cpus=args.cpus, seed0=args.seed, chunk=args.chunk,
                     maxiter=args.maxiter, extra=extra)
    except Exception as ex:
        import traceback; traceback.print_exc()
        L.log_line({'step': 'error', 'job': job, 'tag': tag, 'error': repr(ex)[:500]})
print('DONE', round(time.time() - t00, 1), 's', flush=True)
