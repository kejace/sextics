"""Sector-restricted E1-eliminated variable-projection search over the E1-admissible family (kernel + J_d^c).

usage: python3 w2_search.py --d 8 --sectors 04,00,... (or 'all', 'even', 'odd') --restarts N --workers W
                            [--relax] [--seed S] [--maxiter 400] [--tag name] [--time-limit sec per sector]
Needs wave2/e1space_d{d}_c{ab}.json (from w2_dims.py).  Appends one batch line per (d, sector, mode) to
wave2-search.log.jsonl and per-restart records to wave2/runs/<tag>_d{d}_c{ab}_{mode}.jsonl.
Candidates with rel < 1e-6 are rationally reconstructed (increasing denominators) and re-solved exactly.
"""
import sys, time, json, os, argparse
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2')
from w2lib import *
import numpy as np
from multiprocessing import Pool

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--sectors', default='all')
ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--workers', type=int, default=8)
ap.add_argument('--relax', action='store_true')
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--maxiter', type=int, default=400)
ap.add_argument('--tag', default='s')
ap.add_argument('--time-limit', type=float, default=None)
ap.add_argument('--recon-tol', type=float, default=1e-6)
ap.add_argument('--min-family', type=int, default=1)
args = ap.parse_args()

inst = Inst()
os.makedirs(f'{W2}/runs', exist_ok=True)
if args.sectors == 'all':
    secs = inst.chars
elif args.sectors == 'even':
    secs = [c for c in inst.chars if c[1] % 2 == 0]
elif args.sectors == 'odd':
    secs = [c for c in inst.chars if c[1] % 2 == 1]
else:
    secs = [(int(s[0]), int(s[1])) for s in args.sectors.split(',')]

V = None


def init_worker(d, c, relax):
    global V
    meta = load_space(d, c)
    basis, w_of = basis_from_meta(inst, meta)
    V = SectorVP(inst, basis, c, relax=relax, w_of_basis=w_of if relax else None, verbose=False)


def run_seed(seed):
    global V
    t0 = time.time()
    res = V.search(n_restarts=1, seed=seed, maxiter=args.maxiter, verbose=False)
    r = res[0]
    return {'seed': seed, 'rel': r['rel'], 'f': r['f'], 'bnorm': r['bnorm'], 'iters': r['iters'],
            'a': r['a'].tolist(), 'seconds': round(time.time() - t0, 2)}


for c in secs:
    d = args.d
    fn = f'{W2}/e1space_d{d}_c{c[0]}{c[1]}.json'
    if not os.path.exists(fn):
        print(f'missing {fn}', flush=True); continue
    meta = load_space(d, c)
    famdim = meta['dimJ'] + meta['ker_mod_J']
    if famdim < args.min_family:
        print(f'd={d} c={c}: family dim {famdim} < {args.min_family}, skip', flush=True); continue
    mode = 'relax' if args.relax else 'full'
    tag = f'{args.tag}_d{d}_c{c[0]}{c[1]}_{mode}'
    t0 = time.time()
    basis, w_of = basis_from_meta(inst, meta)
    Vm = SectorVP(inst, basis, c, relax=args.relax, w_of_basis=w_of if args.relax else None, verbose=True)
    if Vm.dimP == 0:
        line = {'batch': 'varpro', 'd': d, 'sector': list(c), 'mode': mode, 'family_dim': Vm.K, 'excluded_F': Vm.excluded,
                'note': 'family consists of F-multiples only (nothing to search)', 'seconds': round(time.time()-t0, 1)}
        log_batch(line); print(json.dumps(line), flush=True); continue
    sizes = {'K': Vm.K, 'family_rank': Vm.family_rank, 'n1': Vm.n1, 'n2': Vm.n2, 'nN': Vm.nN, 'nQ': Vm.nQ, 'rank1': Vm.rank1,
             'dimK': Vm.dimK, 'rQ': Vm.rQ, 'e1_resid': Vm.e1_resid, 'excluded_F': Vm.excluded, 'dimP': Vm.dimP}
    print(f'[{tag}] sizes {sizes}', flush=True)
    seeds = [args.seed * 100000 + i for i in range(args.restarts)]
    results = []
    tstart = time.time()
    with Pool(args.workers, initializer=init_worker, initargs=(d, c, args.relax)) as pool:
        for r in pool.imap_unordered(run_seed, seeds, chunksize=1):
            results.append(r)
            with open(f'{W2}/runs/{tag}.jsonl', 'a') as f:
                f.write(json.dumps(r) + '\n')
            if args.time_limit and time.time() - tstart > args.time_limit:
                pool.terminate(); break
    rels = np.array([r['rel'] for r in results]) if results else np.array([1.0])
    hist = {f'<1e-{k}': int(np.sum(rels < 10.0 ** (-k))) for k in (2, 3, 4, 6, 8, 10, 12)}
    best = sorted(results, key=lambda r: r['rel'])[:5]
    # reconstruction + exact confirmation for small residuals
    confirmed = []; attempts = []
    for r in sorted(results, key=lambda r: r['rel']):
        if r['rel'] > args.recon_tol:
            break
        a = np.array(r['a'])
        for fr in Vm.reconstruct(a):
            status, data = Vm.confirm(fr)
            attempts.append({'seed': r['seed'], 'rel': r['rel'], 'status': status,
                             'rank': data.get('rank'), 'rank_aug': data.get('rank_aug'), 'D': data.get('D', '')[:300]})
            print(f'[{tag}] exact confirm seed {r["seed"]} rel {r["rel"]:.2e} -> {status}', flush=True)
            if status == 'solved':
                path = save_witness_candidate(data['D'], data['N'], data['Q'], d, family=f'varpro {tag}')
                confirmed.append({'D': data['D'], 'path': path})
                print('*** WITNESS FOUND ***', data['D'], flush=True)
                break
        if confirmed:
            break
    line = {'batch': 'varpro', 'd': d, 'sector': list(c), 'yparity': c[1] % 2, 'mode': mode, 'family_dim': Vm.K,
            'sizes': sizes, 'restarts': len(results), 'seed_base': args.seed, 'maxiter': args.maxiter,
            'best_rel': float(rels.min()), 'median_rel': float(np.median(rels)), 'hist': hist,
            'best': [{'seed': b['seed'], 'rel': b['rel'], 'iters': b['iters']} for b in best],
            'recon_attempts': attempts, 'confirmed': confirmed, 'seconds': round(time.time() - t0, 1),
            'runs_file': f'{W2}/runs/{tag}.jsonl'}
    log_batch(line)
    print(json.dumps({k: v for k, v in line.items() if k != 'recon_attempts'}), flush=True)
    if confirmed:
        print('WITNESS FOUND; stopping', flush=True)
        break
