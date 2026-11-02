"""Multiprocess E1-reduced variable-projection search over a denominator family (one degree, one Klein sector).

python3 w2_search.py --d 7 --chi 0 0 --family jac|colon --mode ob|full --restarts 200 --procs 8 --seed 1 --tag NAME
        [--maxiter 300] [--stop-rel 1e-13] [--sub RANGE]   (--sub 'k0:k1' restricts the family to a slice of the basis)
Per-restart records: wave2-<tag>.w<id>.jsonl ; summary line appended to wave2-search.log.jsonl ;
candidates with rel < 1e-8 saved to wave2-<tag>.cands.json ; rel < 1e-10 -> rational reconstruction + exact re-solve.
"""
import sys, os, json, time, argparse, math
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual')
import numpy as np
import multiprocessing as mp
import w2rcr
from hodge import core, fixed_d, jacobian
from fractions import Fraction

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--chi', type=int, nargs=2, required=True)
ap.add_argument('--family', default='jac')
ap.add_argument('--mode', default='full')
ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--procs', type=int, default=4)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--tag', required=True)
ap.add_argument('--maxiter', type=int, default=300)
ap.add_argument('--stop-rel', type=float, default=1e-13)
ap.add_argument('--sub', default=None)
ap.add_argument('--cands-file', default=None, help='JSON list of poly strings to use as the family (family=list)')
args = ap.parse_args()
d = args.d; chi = tuple(args.chi)
F, A, q = w2rcr.load()
LDIR = w2rcr.LDIR
t0 = time.time()


def build_basis():
    polys, gens = [], []
    if args.family in ('jac', 'colon'):
        polys, gens = w2rcr.jacobian_family(F, d, chi)
        polys = list(polys)
    reps = []
    if args.family == 'colon':
        fn = f'{LDIR}/e1space-d{d}-chi{chi[0]}{chi[1]}.json'
        info = json.load(open(fn))
        if info.get('dim_nonjac', 0) > 0:
            assert info.get('rational_lift'), 'no rational lift available'
            reps = [core.parse(s) for s in info['reps']]
    if args.family == 'list':
        reps = [core.parse(s) for s in json.load(open(args.cands_file))]
    basis = reps + polys
    if args.sub:
        k0, k1 = [int(x) for x in args.sub.split(':')]
        basis = basis[k0:k1]
    return basis, len(reps), len(polys)


basis, nreps, njac = build_basis()
print(f'[search] d={d} chi={chi} family={args.family} K={len(basis)} (nonjac {nreps}, jac {njac}) mode={args.mode}', flush=True)
fam = w2rcr.ReducedFamily(F, A, basis, mode=args.mode, verbose=True)
print(f'[search] built in {fam.build_time:.1f}s; e1 relres {fam.e1_relres:.2e}; dimP={fam.dimP}', flush=True)
gc = fam.check_grad(0)
print(f'[search] gradcheck f={gc[0]:.3e} analytic={gc[1]:.6e} numeric={gc[2]:.6e}', flush=True)

per = [args.restarts // args.procs + (1 if w < args.restarts % args.procs else 0) for w in range(args.procs)]


def worker(wid):
    rng_seed = args.seed * 100003 + wid
    path = f'{LDIR}/wave2-{args.tag}.w{wid}.jsonl'
    fh = open(path, 'a')
    best = None

    def cb(c):
        nonlocal best
        rec = {'restart': c['restart'], 'rel': c['rel'], 'f': c['value'], 'iters': c['iters'], 'nfev': c['nfev'],
               'rank': c['rank'], 'spurious': c['spurious'], 'div_rel': c['div_rel'], 'seconds': c['seconds']}
        if c['rel'] < 1e-6 and not c['spurious']:
            rec['a'] = [float(x) for x in c['a']]
        fh.write(json.dumps(rec) + '\n'); fh.flush()
        if best is None or c['rel'] < best:
            best = c['rel']
    res = fam.search(n_restarts=per[wid], seed=rng_seed, maxiter=args.maxiter, stop_rel=args.stop_rel, callback=cb)
    fh.close()
    out = [{'rel': c['rel'], 'a': [float(x) for x in c['a']], 'spurious': c['spurious'], 'rank': c['rank'],
            'iters': c['iters']} for c in res[:5]]
    return out


if __name__ == '__main__':
    ctx = mp.get_context('fork')
    with ctx.Pool(args.procs) as pool:
        results = pool.map(worker, range(args.procs))
    # merge per-restart records
    rels, ranks, iters, spur = [], [], [], 0
    for wid in range(args.procs):
        for line in open(f'{LDIR}/wave2-{args.tag}.w{wid}.jsonl'):
            r = json.loads(line)
            rels.append(r['rel']); iters.append(r['iters']); spur += int(r['spurious'])
            if r.get('rank') is not None: ranks.append(r['rank'])
    rels = np.array(rels)
    lr = np.log10(np.maximum(rels, 1e-300))
    hist = {f'1e{int(b)}': int(n) for b, n in zip(*np.unique(np.floor(lr), return_counts=True))}
    cands = [c for res in results for c in res if c['rel'] < 1e-8 and not c['spurious']]
    cands.sort(key=lambda c: c['rel'])
    summary = {'step': 'varpro', 'tag': args.tag, 'd': d, 'chi': list(chi), 'family': args.family, 'mode': args.mode,
               'K': fam.K, 'dimP': fam.dimP, 'rows': fam.SB.n,
               'nK': getattr(fam, 'nK', None), 'nQ': getattr(fam, 'nQ', None), 'ob_dim': getattr(fam, 'rdim', None),
               'restarts': int(len(rels)), 'best_rel': float(rels.min()) if len(rels) else None,
               'median_rel': float(np.median(rels)) if len(rels) else None, 'hist_log10': hist,
               'spurious': spur, 'rank_min': min(ranks) if ranks else None, 'rank_max': max(ranks) if ranks else None,
               'iters_mean': float(np.mean(iters)) if iters else None, 'e1_relres': fam.e1_relres,
               'seconds': round(time.time() - t0, 1), 'procs': args.procs, 'seed': args.seed, 'engine': 3}
    if cands:
        json.dump(cands, open(f'{LDIR}/wave2-{args.tag}.cands.json', 'w'))
        summary['n_cands_1e-8'] = len(cands)
    # reconstruction attempts for rel < 1e-10
    recon = []
    for c in cands:
        if c['rel'] > 1e-10:
            continue
        a = np.array(c['a'])
        for md, err, P in fam.reconstruct_D(a, max_den=100000):
            if not P:
                continue
            ok, info = fixed_d.e1_solvable(F, A, 1, P)
            rec = {'max_den': md, 'err': err, 'e1': bool(ok), 'nterms': len(P)}
            if ok:
                status, data = fam.confirm(P, verbose=False)
                rec['status'] = status
                if status == 'solved':
                    N, Q = data['N'], data['Q']
                    okv, vinfo = core.verify_witness(F, A, 1, P, N, Q)
                    divF, _ = jacobian.divisible_by(F, P)
                    rec['verified'] = bool(okv); rec['F_divides_D'] = bool(divF)
                    if okv and not divF:
                        w = {'d': d, 'D': core.to_str(P), 'N': [core.to_str(x) for x in N], 'Q': core.to_str(Q),
                             'tag': args.tag, 'instance': w2rcr.LANE}
                        json.dump(w, open(f'{LDIR}/witness-candidate.json', 'w'))
                        print('WITNESS FOUND', flush=True)
                        summary['WITNESS'] = True
            recon.append(rec)
            if rec.get('status') == 'solved':
                break
    if recon:
        summary['reconstruction'] = recon
    w2rcr.log_batch(summary)
    print(json.dumps(summary), flush=True)
