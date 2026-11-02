#!/usr/bin/env python3
"""Multiprocess sector-restricted VarPro search (LM, random restarts) for one (d, chi, family).

  python3 w2_varpro.py --d 9 --c 0,0 --kind colon --restarts 2000 --workers 8 --seed0 1000 --job NAME
         [--maxiter 120] [--time-limit 7000] [--F fermat|offpaired] [--A expr] [--opt lm|lbfgs] [--store-rel 1e-3]

Outputs in cwd: wave2-<job>.w<k>.jsonl (per-seed records), wave2-<job>.summary.json, and one batch line appended to
wave2-search.remote.log.jsonl.  Candidates with relres < 1e-9 are rationally reconstructed and re-solved exactly;
a verified witness is written to witness-candidate-<job>.json and witness-candidate.json.
"""
import os, sys
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1'); os.environ['OMP_NUM_THREADS'] = os.environ.get('W2_THREADS', '1')
os.environ['OPENBLAS_NUM_THREADS'] = os.environ['OMP_NUM_THREADS']
import argparse, json, time, math, glob
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import w2lib
from hodge import core

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--c', default='0,0')
ap.add_argument('--kind', default='colon')
ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--workers', type=int, default=8)
ap.add_argument('--seed0', type=int, default=0)
ap.add_argument('--maxiter', type=int, default=120)
ap.add_argument('--time-limit', type=float, default=7000.0)
ap.add_argument('--F', default='offpaired')
ap.add_argument('--A', default=None)
ap.add_argument('--opt', default='lm')
ap.add_argument('--store-rel', type=float, default=1e-3)
ap.add_argument('--confirm-rel', type=float, default=1e-9)
ap.add_argument('--job', required=True)
ap.add_argument('--cond', type=float, default=1e-9)
args = ap.parse_args()

c = tuple(int(x) for x in args.c.split(','))
F, A, q = w2lib.load()
if args.F == 'fermat':
    F = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
elif args.F != 'offpaired':
    F = core.parse(args.F)
if args.A:
    A = core.parse(args.A)
T0 = time.time()
print(f'[{w2lib.now()}] job {args.job}: d={args.d} c={c} kind={args.kind} F={args.F} restarts={args.restarts} workers={args.workers}', flush=True)
basis, gauge, info = w2lib.build_family(F, A, q, args.d, c, args.kind, verbose=True)
if not basis or len(basis) <= len(gauge):
    rec = {'ts': w2lib.now(), 'job': args.job, 'step': 'varpro', 'F': args.F, 'd': args.d, 'sector': list(c), 'family': args.kind,
           'K': len(basis), 'gauge': len(gauge), 'note': 'empty family after gauge removal', 'restarts': 0}
    w2lib.append_json('wave2-search.remote.log.jsonl', rec); print(rec, flush=True); sys.exit(0)
fam = w2lib.SectorVarPro(F, A, q, basis, c, gauge=gauge, verbose=True, cond=args.cond)
print(f'[{w2lib.now()}] family built: K={fam.K} rank={fam.rankB} gauge={fam.ng} dimP={fam.dimP} rows={fam.R} N={fam.nN} Q={fam.nQ} '
      f'build {fam.build_time:.1f}s', flush=True)
try:
    g1, g2 = fam.fd_check(0)
    print(f'gradient check: analytic {g1:.6e} fd {g2:.6e}', flush=True)
except Exception as ex:
    print('gradient check failed:', ex, flush=True)


def worker(k, seeds):
    out = f'wave2-{args.job}.w{k}.jsonl'
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try: done.add(json.loads(line)['seed'])
            except Exception: pass
    for s in seeds:
        if s in done:
            continue
        if time.time() - T0 > args.time_limit:
            break
        try:
            if args.opt == 'lm':
                r = fam.lm_from(s, maxiter=args.maxiter)
            else:
                r = fam.minimize_from(s, maxiter=max(args.maxiter, 300))
        except Exception as ex:
            w2lib.append_json(out, {'seed': int(s), 'error': str(ex)}); continue
        rec = {k_: v for k_, v in r.items() if k_ not in ('u', 'Dc')}
        if r['relres'] < args.store_rel:
            rec['Dc'] = [float(x) for x in r['Dc']]
        w2lib.append_json(out, rec)


seeds = list(range(args.seed0, args.seed0 + args.restarts))
chunks = [seeds[k::args.workers] for k in range(args.workers)]
ctx = mp.get_context('fork')
procs = [ctx.Process(target=worker, args=(k, chunks[k])) for k in range(args.workers)]
for p in procs: p.start()
for p in procs: p.join()

# aggregate
recs = []
for fn in sorted(glob.glob(f'wave2-{args.job}.w*.jsonl')):
    for line in open(fn):
        try: recs.append(json.loads(line))
        except Exception: pass
rels = sorted(r['relres'] for r in recs if 'relres' in r)
hist = {}
for r in rels:
    b = 'lt1e-12' if r < 1e-12 else ('lt1e-9' if r < 1e-9 else ('lt1e-6' if r < 1e-6 else ('lt1e-4' if r < 1e-4 else ('lt1e-2' if r < 1e-2 else 'ge1e-2'))))
    hist[b] = hist.get(b, 0) + 1
# distinct local minima (rounded)
mins = {}
for r in rels:
    key = f'{r:.3e}'
    mins[key] = mins.get(key, 0) + 1
top_minima = sorted(mins.items(), key=lambda kv: float(kv[0]))[:12]
summary = {'ts': w2lib.now(), 'job': args.job, 'step': 'varpro', 'F': args.F, 'd': args.d, 'sector': list(c), 'family': args.kind,
           'K': fam.K, 'rank': fam.rankB, 'gauge': fam.ng, 'dimP': fam.dimP, 'rows': fam.R, 'nN': fam.nN, 'nQ': fam.nQ,
           'restarts': len(rels), 'best_relres': rels[0] if rels else None, 'median_relres': rels[len(rels)//2] if rels else None,
           'hist': hist, 'top_minima': top_minima, 'opt': args.opt, 'maxiter': args.maxiter,
           'mean_time_per_restart': round(float(np.mean([r['time'] for r in recs if 'time' in r])), 2) if recs else None,
           'statuses': {s: sum(1 for r in recs if r.get('status') == s) for s in set(r.get('status') for r in recs if 'status' in r)},
           'wall': round(time.time() - T0, 1)}
# confirmation attempts
witness = None
cands = sorted([r for r in recs if 'Dc' in r and r['relres'] < args.confirm_rel], key=lambda r: r['relres'])
attempts = []
for r in cands[:20]:
    att = w2lib.try_reconstruct_and_confirm(fam, np.array(r['Dc']), verbose=False)
    attempts.append({'seed': r['seed'], 'relres': r['relres'], 'attempts': att})
    if any(a.get('status') == 'solved' for a in att):
        sol = [a for a in att if a.get('status') == 'solved'][0]
        witness = {'d': args.d, 'D': sol['D'], 'N': sol['N'], 'Q': sol['Q'], 'job': args.job, 'seed': r['seed'], 'relres': r['relres'],
                   'F': args.F, 'sector': list(c), 'family': args.kind}
        json.dump(witness, open(f'witness-candidate-{args.job}.json', 'w'))
        if args.F == 'offpaired' and args.A is None:
            json.dump(witness, open('witness-candidate.json', 'w'))
        print('!!!!!!!! WITNESS FOUND (exact, verified) !!!!!!!!', flush=True)
        break
summary['confirm_attempts'] = attempts
summary['witness'] = bool(witness)
json.dump(summary, open(f'wave2-{args.job}.summary.json', 'w'), indent=1)
w2lib.append_json('wave2-search.remote.log.jsonl', {k_: v for k_, v in summary.items() if k_ != 'confirm_attempts'})
print(json.dumps({k_: v for k_, v in summary.items() if k_ not in ('confirm_attempts',)}, indent=1), flush=True)
