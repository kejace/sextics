#!/usr/bin/env python3
"""Multi-restart joint search (variable projection, LM) over an E1-admissible family for polarized-k3-remainder.

  python3 w2_search.py --d 8 --kind jac --mode full --restarts 500 --workers 8 --seed0 1000 --job NAME
        [--maxiter 60] [--time-limit 7000] [--relax-restarts 0] [--cond 1e-10] [--sketch-margin 300]

mode: full  -> complete E2 with Koszul-correction unknowns (c) and quotient (C)
      A1    -> relaxation c = 0 (N = A w only), C free                      (sufficient sub-search)
      A2    -> relaxation B(a) in L_{2d}, L = (F) + J^2 + (R_ijk)           (necessary quadric system)
Outputs (cwd = lane dir): wave2-<job>.w<k>.jsonl (per seed), wave2-<job>.summary.json, one batch line in
wave2-search.log.jsonl; candidates with relres < confirm-rel are rationally reconstructed and re-solved
exactly; a verified witness goes to witness-candidate-<job>.json and witness-candidate.json.
"""
import os, sys
os.environ.setdefault('OMP_NUM_THREADS', os.environ.get('W2_THREADS', '1'))
os.environ['OPENBLAS_NUM_THREADS'] = os.environ['OMP_NUM_THREADS']
import argparse, json, time, math, glob
import multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np
import w2k3 as W
from hodge import core, fixed_d, jacobian

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--kind', default='jac')              # jac | full
ap.add_argument('--mode', default='full')             # full | A1 | A2
ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--workers', type=int, default=8)
ap.add_argument('--seed0', type=int, default=0)
ap.add_argument('--maxiter', type=int, default=60)
ap.add_argument('--time-limit', type=float, default=7000.0)
ap.add_argument('--job', required=True)
ap.add_argument('--cond', type=float, default=1e-10)
ap.add_argument('--sketch-margin', type=int, default=300)
ap.add_argument('--store-rel', type=float, default=1e-2)
ap.add_argument('--confirm-rel', type=float, default=1e-9)
ap.add_argument('--family-note', default='')
args = ap.parse_args()

F, A, q = W.load()
T0 = time.time()
print(f'[{W.now()}] job {args.job}: d={args.d} kind={args.kind} mode={args.mode} restarts={args.restarts} '
      f'workers={args.workers} seed0={args.seed0}', flush=True)
kernel_data = None
if args.kind == 'full':
    kf = f'wave2-kernel-float-d{args.d}.json'
    if os.path.exists(kf):
        kernel_data = json.load(open(kf))
        print(f'loaded float kernel from {kf}: dim {kernel_data["kernel_dim"]}', flush=True)
fam = W.Family(F, A, args.d, kind=args.kind, verbose=True, cond=args.cond, sketch_margin=args.sketch_margin,
               kernel_data=kernel_data, with_relaxations=(args.mode == 'A2'), seed=args.seed0)
if args.kind == 'full' and kernel_data is None and fam.kernel_info is not None:
    pass
fam_desc = {'d': args.d, 'kind': args.kind, 'K': fam.K, 'Kw': len(fam.wcols), 'Kker': fam.Kker, 'rankD': fam.rankD,
            'gauge': fam.ngauge, 'removed': fam.nbad, 'dimP': fam.dimP, 'rows': fam.SR.n, 'nc': fam.nc, 'nC': fam.nC,
            'ks': fam.ks, 'build_time': round(fam.build_time, 1)}
if args.mode == 'A2':
    fam_desc.update({'rankL': fam.rankL, 'codimL': fam.codimL})
print(json.dumps(fam_desc), flush=True)
if fam.dimP == 0:
    rec = {'step': 'search', 'job': args.job, 'mode': args.mode, 'note': 'empty family after gauge removal', **fam_desc}
    W.log_batch(rec); print(rec); sys.exit(0)


def run_one(seed):
    if args.mode == 'full':
        return fam.lm(seed=seed, maxiter=args.maxiter, tol=1e-13)
    else:
        return fam.relax_minimize(seed=seed, mode=('L' if args.mode == 'A2' else 'F'), maxiter=max(args.maxiter, 300))


def worker(k, seeds):
    out = f'wave2-{args.job}.w{k}.jsonl'
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try:
                done.add(json.loads(line)['seed'])
            except Exception:
                pass
    for s in seeds:
        if s in done:
            continue
        if time.time() - T0 > args.time_limit:
            break
        try:
            r = run_one(s)
        except Exception as ex:
            W.append_json(out, {'seed': int(s), 'error': repr(ex)}); continue
        rec = {k_: v for k_, v in r.items() if k_ not in ('u', 'a', 'hist')}
        rec['seed'] = int(s)
        if r['relres'] < args.store_rel:
            rec['a'] = [float(x) for x in r['a']]
        W.append_json(out, rec)


seeds = list(range(args.seed0, args.seed0 + args.restarts))
chunks = [seeds[k::args.workers] for k in range(args.workers)]
ctx = mp.get_context('fork')
procs = [ctx.Process(target=worker, args=(k, chunks[k])) for k in range(args.workers)]
for p in procs:
    p.start()
for p in procs:
    p.join()

# ---- aggregate ------------------------------------------------------------------------------------------
recs = []
for fn in sorted(glob.glob(f'wave2-{args.job}.w*.jsonl')):
    for line in open(fn):
        try:
            recs.append(json.loads(line))
        except Exception:
            pass
rels = sorted(r['relres'] for r in recs if 'relres' in r)
bins = [(1e-12, 'lt1e-12'), (1e-9, 'lt1e-9'), (1e-6, 'lt1e-6'), (1e-4, 'lt1e-4'), (1e-2, 'lt1e-2'), (1e-1, 'lt1e-1'), (1e9, 'ge1e-1')]
hist = {}
for r in rels:
    for b, name in bins:
        if r < b:
            hist[name] = hist.get(name, 0) + 1; break
mins = {}
for r in rels:
    k_ = f'{r:.2e}'; mins[k_] = mins.get(k_, 0) + 1
top = sorted(mins.items(), key=lambda kv: float(kv[0]))[:15]
summary = {'step': 'search', 'job': args.job, 'mode': args.mode, 'family': args.kind + (' ' + args.family_note if args.family_note else ''),
           **fam_desc, 'restarts': len(rels), 'best_relres': rels[0] if rels else None,
           'median_relres': rels[len(rels) // 2] if rels else None, 'hist': hist, 'top_minima': top,
           'maxiter': args.maxiter, 'mean_time_per_restart': round(float(np.mean([r['time'] for r in recs if 'time' in r])), 1) if recs else None,
           'statuses': {str(s): sum(1 for r in recs if str(r.get('status')) == str(s)) for s in set(r.get('status') for r in recs if 'status' in r)},
           'errors': sum(1 for r in recs if 'error' in r), 'wall': round(time.time() - T0, 1)}

# ---- candidates: diagnostics + reconstruction + exact confirmation ---------------------------------------
witness = None
cands = sorted([r for r in recs if 'a' in r and r['relres'] < max(args.confirm_rel, 1e-6)], key=lambda r: r['relres'])
attempts = []
fU = core.parse("(X0+Y0)**6 + 5*(X0+Y0)**5*(X1+Y1)/7 + 5*(X0+Y0)**2*(X1+Y1)**3*(X2+Y2)/7 + 15*(X0+Y0)*(X2+Y2)**5/7 + (X1+Y1)**6 + 10*(X1+Y1)**5*(X2+Y2)/7 + (X2+Y2)**6")
fV = core.parse("(X0-Y0)**6 + 5*(X0-Y0)**5*(X1-Y1)/7 + 5*(X0-Y0)**2*(X1-Y1)**3*(X2-Y2)/7 + 15*(X0-Y0)*(X2-Y2)**5/7 + (X1-Y1)**6 + 10*(X1-Y1)**5*(X2-Y2)/7 + (X2-Y2)**6")
for r in cands[:25]:
    a = np.array(r['a'])
    att = {'seed': r['seed'], 'relres': r['relres']}
    try:
        if args.mode == 'full':
            att['full_relres'] = fam.full_residual(a)
        att['div_fU'] = fam.divisibility(a, fU, 6); att['div_fV'] = fam.divisibility(a, fV, 6)
        att['div_F'] = fam.divisibility(a, F, 6)
        att['div_Fi'] = [fam.divisibility(a, fam.Fi[i], 5) for i in range(6)]
        att['support'] = int(np.sum(np.abs(fam.D_of(a)) > 1e-8 * np.max(np.abs(fam.D_of(a)))))
    except Exception as ex:
        att['diag_error'] = repr(ex)
    # reconstruction attempts (only meaningful if relres small)
    att['recon'] = []
    if r['relres'] < args.confirm_rel:
        tried = set()
        for mode in ('w', 'D'):
            for md in (1, 2, 6, 12, 60, 120, 1000, 10000, 100000, 1000000):
                fr, err = fam.reconstruct(a, max_den=md, mode=mode)
                if err > 1e-6:
                    att['recon'].append({'mode': mode, 'max_den': md, 'err': err, 'status': 'no-fit'}); continue
                if mode == 'w' and fam.Kker == 0:
                    w = [{} for _ in range(6)]
                    for j, (i, m) in enumerate(fam.wcols):
                        if fr[j]:
                            w[i][m] = fr[j]
                    D = {}
                    for i in range(6):
                        D = core.add(D, core.mul(w[i], fam.Fi[i]))
                else:
                    D = fam.SD.poly_exact(fr) if mode == 'D' else None
                    if D is None:
                        continue
                    w = W.w_from_D_exact(F, args.d, D) if fam.Kker == 0 else None
                Ds = core.to_str(D)
                if Ds in tried or not D:
                    continue
                tried.add(Ds)
                t0 = time.time()
                if not W.coprime_to_F(F, D):
                    att['recon'].append({'mode': mode, 'max_den': md, 'status': 'spurious (F | D)'}); continue
                if w is not None:
                    st, data = W.exact_reduced_solve(F, A, args.d, w, verbose=False, max_primes=12)
                else:
                    out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=False)
                    st = 'solved' if out.get('exact') else ('inconsistent' if not out['consistent'] else 'lift-failed')
                    data = out
                rec_ = {'mode': mode, 'max_den': md, 'err': err, 'status': st, 'D': Ds[:400], 'time': round(time.time() - t0, 1),
                        'rank': data.get('rank'), 'nullity': data.get('nullity')}
                att['recon'].append(rec_)
                if st == 'solved':
                    ok, info = core.verify_witness(F, A, q, data['D'] if 'D' in data else D, data['N'], data['Q'])
                    if ok and W.coprime_to_F(F, D):
                        witness = {'d': args.d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']],
                                   'Q': core.to_str(data['Q']), 'job': args.job, 'seed': r['seed'], 'relres': r['relres']}
                        json.dump(witness, open(f'witness-candidate-{args.job}.json', 'w'))
                        json.dump(witness, open('witness-candidate.json', 'w'))
                        print('!!!!!!!! WITNESS FOUND (exact, verified) !!!!!!!!', flush=True)
                        break
            if witness:
                break
    attempts.append(att)
    if witness:
        break
summary['candidates'] = attempts
summary['witness'] = bool(witness)
json.dump(summary, open(f'wave2-{args.job}.summary.json', 'w'), indent=1)
W.log_batch({k_: v for k_, v in summary.items() if k_ != 'candidates'} | {'n_candidates_examined': len(attempts),
             'best_full_relres': min([a_.get('full_relres', 1.0) for a_ in attempts], default=None)})
print(json.dumps({k_: v for k_, v in summary.items() if k_ != 'candidates'}, indent=1), flush=True)
for a_ in attempts[:5]:
    print(json.dumps(a_)[:1500], flush=True)
