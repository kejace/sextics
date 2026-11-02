"""Multi-seed character-twisted varpro search over an exact family (kernel + J_d^chi or sub-family).

usage: python3 w2_varpro.py --chi triv --d 10 --pD 0 --mode cvp|kvp|hess --sub all|kernel|jac --seeds 0-7
                            --restarts 25 --maxiter 400 [--blocks S3S3] [--nproc 8] [--tag T]
Per seed: random-restart L-BFGS; every candidate with rel < 1e-8 is rationally reconstructed and re-solved exactly
(reduced solver) immediately.  Output: w2vp_<tag>.jsonl (one line per restart), batch line in wave2-search.log.jsonl.
"""
import sys, os, time, json, argparse, pickle
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import numpy as np
import w2lib as W
from hodge import core
from multiprocessing import Pool

ap = argparse.ArgumentParser()
ap.add_argument('--chi', default='triv'); ap.add_argument('--d', type=int, required=True); ap.add_argument('--pD', type=int, default=0)
ap.add_argument('--blocks', default='S3S3'); ap.add_argument('--mode', default='cvp'); ap.add_argument('--sub', default='all')
ap.add_argument('--seeds', default='0-3'); ap.add_argument('--restarts', type=int, default=20); ap.add_argument('--maxiter', type=int, default=400)
ap.add_argument('--nproc', type=int, default=int(os.environ.get('W2_NPROC', '4'))); ap.add_argument('--tag', default='')
ap.add_argument('--rel_confirm', type=float, default=1e-8)
args = ap.parse_args()
s0, s1 = [int(x) for x in args.seeds.split('-')]
seeds = list(range(s0, s1 + 1))
tag = args.tag or f'{args.blocks}_{args.chi}_d{args.d}_p{args.pD}_{args.mode}_{args.sub}'
OUT = f'/tmp/hodge/results/symmetric-sextic-residual/w2vp_{tag}.jsonl'

F, A, q = W.load()
G = W.group(args.chi, args.blocks)
fam = W.family(F, A, q, args.d, args.pD, args.chi, args.blocks, verbose=True)
kernel, jac = fam['kernel'], fam['jac']
if args.sub == 'kernel': basis = kernel
elif args.sub == 'jac': basis = jac
else: basis = kernel + jac
if not basis:
    rec = {'batch': 'varpro', 'tag': tag, 'family': f'{args.blocks}/{args.chi}', 'd': args.d, 'pD': args.pD, 'mode': args.mode,
           'sub': args.sub, 'K': 0, 'note': 'empty family'}
    W.wlog(rec); print(json.dumps(rec)); sys.exit(0)
t0 = time.time()
if args.mode == 'cvp':
    vp = W.CharVarPro(F, A, q, basis, G, verbose=True)
else:
    pairs = []
    if args.sub in ('all', 'kernel'):
        for P in kernel:
            N0, relres = W.e1_lift_float(F, A, q, P, G)
            if relres > 1e-8:
                print(f'[w2_varpro] float lift residual {relres:.2e} for a kernel element; skipping it', file=sys.stderr); continue
            pairs.append((P, N0))
        print(f'[w2_varpro] {len(pairs)} kernel lifts', file=sys.stderr)
    if args.sub in ('all', 'jac'):
        pairs += W.jac_family_with_lifts(F, A, q, args.d, args.pD, G)
    vp = W.KoszulVarPro(F, A, q, pairs, G, koszul=(args.mode == 'kvp'), verbose=True)
build = time.time() - t0
info = {'K': vp.K, 'dimP': vp.dimP, 'rows': vp.rows, 'cols': vp.ncols, 'nN': vp.nN, 'nQ': vp.nQ, 'excluded': vp.excluded_dim, 'build_s': round(build, 1)}
print('[w2_varpro]', tag, info, flush=True)

def run_seed(seed):
    t0 = time.time(); out = []; hits = []
    def on_cand(c):
        rec = {'seed': seed, 'restart': c['restart'], 'rel': c['rel'], 'f': c['value'], 'spurious': bool(c['spurious']),
               'div_rel': c['div_rel'], 'iters': c['iters']}
        if c['rel'] < 1e-4 and not c['spurious']:
            rec['a'] = [float(v) for v in c['a']]
        if c['rel'] < args.rel_confirm and not c['spurious']:
            st, data, fr = W.try_reconstruct_confirm(vp, c['a'], verbose=True, log_tag=f' seed{seed}')
            rec['confirm'] = st; rec['frac'] = [str(x) for x in fr] if fr else None
            if st == 'solved':
                D = data['D']; W.save_witness_candidate(D, data['N'], data['Q'], extra={'tag': tag, 'seed': seed, 'a': [str(x) for x in fr]})
                rec['WITNESS'] = True; hits.append(rec)
                with open('/tmp/hodge/results/symmetric-sextic-residual/WITNESS_FOUND.txt', 'a') as f:
                    f.write(json.dumps(rec) + '\n')
        with open(OUT, 'a') as f:
            f.write(json.dumps({'tag': tag, **rec}) + '\n')
        out.append(rec)
    res = vp.search(n_restarts=args.restarts, seed=seed, maxiter=args.maxiter, on_candidate=on_cand, stop_rel=1e-13)
    best = res[0]
    # rank statistics at the best point
    a = best['a']; M, b = vp.assemble(a / np.linalg.norm(a))
    sv = np.linalg.svd(M, compute_uv=False); rk = int(np.sum(sv > 1e-9 * sv[0]))
    return {'seed': seed, 'best_rel': best['rel'], 'best_spurious': bool(best['spurious']), 'n': len(out),
            'rels': [r['rel'] for r in out], 'rank': rk, 'nullity': M.shape[1] - rk, 'rows': M.shape[0],
            'seconds': round(time.time() - t0, 1), 'hits': len(hits)}

if __name__ == '__main__':
    t1 = time.time()
    with Pool(min(args.nproc, len(seeds))) as pool:
        results = pool.map(run_seed, seeds)
    rels = sorted(r for res in results for r in res['rels'])
    hist = {'<1e-10': sum(r < 1e-10 for r in rels), '<1e-6': sum(r < 1e-6 for r in rels), '<1e-3': sum(r < 1e-3 for r in rels),
            '<1e-2': sum(r < 1e-2 for r in rels), '<1e-1': sum(r < 1e-1 for r in rels), 'total': len(rels)}
    rec = {'batch': 'varpro', 'tag': tag, 'family': f'{args.blocks}/{args.chi}', 'd': args.d, 'pD': args.pD, 'mode': args.mode, 'sub': args.sub,
           **info, 'seeds': args.seeds, 'restarts_per_seed': args.restarts, 'maxiter': args.maxiter, 'restarts': len(rels),
           'best_rel': rels[0] if rels else None, 'median_rel': rels[len(rels)//2] if rels else None, 'hist': hist,
           'rank': results[0]['rank'], 'nullity': results[0]['nullity'], 'hits': sum(r['hits'] for r in results),
           'seconds': round(time.time() - t1, 1), 'quartiles': [rels[int(f*(len(rels)-1))] for f in (0, .25, .5, .75, 1)] if rels else None}
    W.wlog(rec); print(json.dumps(rec), flush=True)
