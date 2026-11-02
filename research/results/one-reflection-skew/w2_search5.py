"""Multi-restart variable-projection search (5-var reduction or 6-var sector) for one-reflection-skew.

usage: python3 w2_search5.py --d 8 --sector 0,0 --family full|J|nonJ|mons|cone4:i|cone3:i,j|bihom:a,b
                             [--target A|A+|A-] [--restarts 200] [--seed 0] [--procs 8] [--maxiter 400] [--tag NAME]
       python3 w2_search5.py --six --d 8 --chi 0,0,0 --family full|J|nonJ|mons ...   (6-var, D of pure character chi)
Writes per-restart lines to w2run_<tag>.jsonl (crash-safe), a batch line to wave2-search.log.jsonl at the end,
candidates (rel < 1e-6) to w2cand_<tag>.pkl, and tries rational reconstruction + exact confirmation on each candidate.
"""
import os, sys, time, json, argparse, pickle, math
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1'); os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
import numpy as np
import multiprocessing as mp

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--sector', default='0,0')
ap.add_argument('--chi', default='0,0,0')
ap.add_argument('--six', action='store_true')
ap.add_argument('--family', default='full')
ap.add_argument('--target', default='A')
ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--procs', type=int, default=8)
ap.add_argument('--maxiter', type=int, default=400)
ap.add_argument('--tag', default=None)
ap.add_argument('--keep', type=float, default=1e-6)
ap.add_argument('--confirm', type=float, default=1e-8)
args = ap.parse_args()

import w2skew as W
from fractions import Fraction
from hodge import core

tag = args.tag or (f"{'six' if args.six else 'five'}_d{args.d}_{'c' + args.chi.replace(',', '') if args.six else 's' + args.sector.replace(',', '')}_{args.family.replace(':', '').replace(',', '')}_{args.target.replace('+','p').replace('-','m')}_seed{args.seed}")
RUN = f'{W.WD}/w2run_{tag}.jsonl'
t0 = time.time()

def build_family():
    if args.six:
        chi = tuple(int(x) for x in args.chi.split(','))
        k = W.load_kernel6(args.d, chi)
        mons = k['mons']
        if args.family == 'full': B = k['Kbasis']
        elif args.family == 'J': B = k['Jbasis']
        elif args.family == 'nonJ': B = k['Cbasis']
        elif args.family == 'mons': B = np.eye(len(mons))
        else: raise SystemExit('family')
        info = {k_: v for k_, v in k.items() if k_ in ('nD', 'kerD', 'dimJ', 'nonJ')}
        vp = W.VP6(args.d, mons, B, verbose=True)
        return vp, info, f'chi{chi[0]}{chi[1]}{chi[2]}'
    rD, sD = (int(x) for x in args.sector.split(','))
    k = W.load_kernel(args.d, rD, sD)
    mons = k['mons']
    base = {'full': k['Kbasis'], 'J': k['Jbasis'], 'nonJ': k['Cbasis'], 'mons': np.eye(len(mons))}
    fam = args.family
    if fam in base:
        B = base[fam]
    else:
        B = W.restrict_basis(k['Kbasis'], W.family_mask(mons, fam))
    tgt = {'A': W.A, 'A+': W.A_plus, 'A-': W.A_minus}[args.target]
    info = {k_: v for k_, v in k.items() if k_ in ('nD', 'kerD', 'dimJ', 'nonJ')}
    info['K'] = int(B.shape[0])
    if B.shape[0] == 0:
        W.logbatch(batch='varpro5', family=f'{fam}/{args.target}', degree=args.d, sector=f'rho{rD}sig{sD}', restarts=0,
                   best_rel=None, seconds=0, note=f'empty family (E1 kernel has no vectors in this sub-family); {info}')
        W.rewrite_report()
        raise SystemExit(0)
    vp = W.VP5(args.d, rD, sD, mons, B, target=tgt, verbose=True)
    return vp, info, f'rho{rD}sig{sD}'

VP, INFO, SECT = build_family()
print(f'[driver] tag={tag} family K={VP.K} rows={VP.rows} cols={getattr(VP, "ncols", getattr(VP, "nN", None))} build {VP.build_time:.1f}s', file=sys.stderr, flush=True)

def worker(seed):
    recs = []
    def on_restart(c):
        rec = {'seed': seed, 'restart': c['restart'], 'rel': c['rel'], 'iters': c['iters'], 'secs': c['secs']}
        if c.get('a') is not None:
            rec['a'] = [float(v) for v in c['a']]
        if 'spurious' in c:
            rec['spurious'] = bool(c['spurious'])
        with open(RUN, 'a') as fh:
            fh.write(json.dumps(rec) + '\n')
        recs.append(rec)
    VP.search(n_restarts=1, seed=seed, maxiter=args.maxiter, keep_rel=args.keep, on_restart=on_restart)
    return recs

seeds = list(range(args.seed, args.seed + args.restarts))
done = set()
if os.path.exists(RUN):
    for line in open(RUN):
        try: done.add(json.loads(line)['seed'])
        except Exception: pass
seeds = [s for s in seeds if s not in done]
print(f'[driver] {len(seeds)} restarts to run ({len(done)} already in {RUN})', file=sys.stderr, flush=True)
if seeds:
    with mp.get_context('fork').Pool(args.procs) as pool:
        for _ in pool.imap_unordered(worker, seeds, chunksize=1):
            pass
# collect
recs = [json.loads(l) for l in open(RUN)]
rels = np.array([r['rel'] for r in recs]); iters = [r['iters'] for r in recs]
best = float(rels.min()) if len(rels) else None
qs = np.quantile(rels, [0, 0.1, 0.5, 0.9, 1]).tolist() if len(rels) else []
cands = [r for r in recs if r['rel'] < args.keep and r.get('a') is not None and not r.get('spurious')]
with open(f'{W.WD}/w2cand_{tag}.pkl', 'wb') as fh:
    pickle.dump({'tag': tag, 'args': vars(args), 'info': INFO, 'cands': cands, 'mons': VP.mons, 'Basis': VP.Basis}, fh)
note = f"K={VP.K} rows={VP.rows} cols={getattr(VP,'ncols',getattr(VP,'nN',None))} {INFO} rel quantiles(0,.1,.5,.9,1)={['%.2e'%q for q in qs]} n_cand(rel<{args.keep})={len(cands)}"
# exact confirmation of candidates
confirmed = []
for c in sorted(cands, key=lambda r: r['rel'])[:20]:
    Dfl = VP.D_float(np.array(c['a']))
    for Dq in W.reconstruct_D(Dfl):
        if args.six:
            from hodge import fixed_d
            ok, info = fixed_d.e1_solvable(W.F6, W.A6, 1, Dq)
            st, data = ('inconsistent', {'note': 'e1 fails'}) if not ok else fixed_d.solve_fixed_D(W.F6, W.A6, 1, Dq)
            if st == 'solved':
                okv, inf = W.save_witness6(Dq, data['N'], data['Q'], note=f'6-var varpro {tag}')
                confirmed.append((core.to_str(Dq), okv)); print('!!!!!!!! WITNESS FOUND (6-var) !!!!!!!!', okv, file=sys.stderr, flush=True)
        else:
            tgt = {'A': W.A, 'A+': W.A_plus, 'A-': W.A_minus}[args.target]
            st, data = W.confirm_exact5(Dq, target=tgt)
            if st == 'solved':
                okv, inf = W.save_witness6(data['D6'], data['N6'], {}, note=f'5-var varpro {tag} target {args.target}')
                confirmed.append((core.to_str(Dq), okv)); print('!!!!!!!! WITNESS FOUND (5-var lift) !!!!!!!!', okv, file=sys.stderr, flush=True)
        W.logbatch(batch='confirm', family=f'{args.family}/{args.target}', degree=args.d, sector=SECT, restarts=1, best_rel=c['rel'],
                   seconds=round(time.time() - t0, 1), note=f'rational candidate {core.to_str(Dq)[:120]} -> {st} {str({k: v for k, v in data.items() if k in ("served_dim","served_mod_J","rank","nullity","note")})[:200]}')
        if st == 'solved':
            break
W.logbatch(batch='varpro6' if args.six else 'varpro5', family=f'{args.family}/{args.target}', degree=args.d, sector=SECT,
           restarts=len(recs), best_rel=best, iters=f'{min(iters) if iters else None}/{int(np.median(iters)) if iters else None}',
           seconds=round(time.time() - t0, 1), note=note + (f' CONFIRMED {confirmed}' if confirmed else ''), tag=tag)
W.rewrite_report()
print('[driver] done', note, file=sys.stderr, flush=True)
